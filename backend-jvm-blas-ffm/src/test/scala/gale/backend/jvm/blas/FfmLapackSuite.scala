package gale.backend.jvm.blas

import gale.backend.{Capability, PureBackend}
import gale.linalg.*
import gale.spectral.SpectralCapability

import scala.concurrent.duration.*
import scala.concurrent.{Await, Future}

class FfmLapackSuite extends munit.FunSuite:
  private lazy val backend = FfmBlasBackend
    .load(thresholds = Some(FfmBlasThresholds(
      nativeLuMinSize = 0,
      nativeCholeskyMinSize = 0,
      nativeQrMinSize = 0
    )))
    .fold(throw _, identity)

  override def afterAll(): Unit =
    if backend != null then backend.close()

  private def requireLapack(): Unit =
    assume(backend.libraryInfo.hasLapack, s"${backend.libraryInfo.name} has no complete LAPACK symbol set")
    assert(backend.capabilities.contains(Capability.NativeLapack))
    assert(backend.spectral.exists(_.capabilities.contains(SpectralCapability.DenseSymmetricEigen)))

  private def assertMatrixClose(actual: DMat, expected: DMat, tolerance: Double = 2e-11): Unit =
    assertEquals(actual.shape, expected.shape)
    var i = 0
    while i < actual.rows do
      var j = 0
      while j < actual.cols do
        val scale = math.max(1.0, math.max(math.abs(actual(i, j)), math.abs(expected(i, j))))
        assert(math.abs(actual(i, j) - expected(i, j)) <= tolerance * scale,
          s"($i,$j): ${actual(i, j)} != ${expected(i, j)}")
        j += 1
      i += 1

  test("NativeLapack metadata exposes typed factorizations and symmetric eigen"):
    requireLapack()
    assert(backend.denseFactorizations.isDefined)
    assert(backend.spectral.isDefined)

  test("facade-routed LU solves and reports singular pivots with Gale errors"):
    requireLapack()
    val a = Matrix.dense(4, 4)(
      0.0, 2.0, 3.0, 1.0,
      4.0, 5.0, 6.0, -2.0,
      7.0, 8.0, 10.0, 3.0,
      2.0, -1.0, 1.0, 8.0
    )
    val b = Vec(3.0, -2.0, 5.0, 7.0)
    val lu = a.lu(using backend).orThrow
    val x = lu.solve(b).orThrow
    val residual = a.*(x)(using PureBackend) - b
    assert(residual.norm2 <= 2e-12 * math.max(1.0, b.norm2))

    val singular = Matrix.dense(3, 3)(1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 0.0, 1.0, 1.0)
    assert(backend.denseFactorizations.get.lu(singular).left.exists(_.isInstanceOf[LinAlgError.SingularMatrix]))

  test("native Cholesky reconstructs across scale and preserves typed non-SPD failure"):
    requireLapack()
    for scale <- Seq(1e-150, 1.0, 1e150) do
      val spd = Matrix.dense(3, 3)(
        6.0 * scale, 2.0 * scale, 1.0 * scale,
        2.0 * scale, 5.0 * scale, 2.0 * scale,
        1.0 * scale, 2.0 * scale, 4.0 * scale
      )
      val chol = spd.cholesky(using backend).orThrow
      assertMatrixClose(chol.lower.*(chol.lower.t)(using PureBackend), spd, 5e-12)

    val indefinite = Matrix.dense(2, 2)(1.0, 2.0, 2.0, 1.0)
    assert(backend.denseFactorizations.get.cholesky(indefinite).left.exists(_.isInstanceOf[LinAlgError.NotPositiveDefinite]))

  test("native QR reconstructs tall and wide matrices and supports least squares"):
    requireLapack()
    for (rows, cols) <- Seq((9, 4), (4, 9), (7, 7)) do
      val a = Matrix.tabulate(rows, cols)((i, j) =>
        math.sin((i + 1).toDouble * (j + 2).toDouble) + (if i == j then 2.0 else 0.0)
      )
      val qr = backend.denseFactorizations.get.qr(a).orThrow
      assertMatrixClose(qr.q.*(qr.r)(using PureBackend), a, 2e-10)
      assertMatrixClose(qr.q.t.*(qr.q)(using PureBackend), DMat.eye(rows), 2e-10)

    val tall = Matrix.tabulate(10, 4)((i, j) => math.cos(i * 0.4 + j * 0.7) + (if i == j then 1.0 else 0.0))
    val truth = Vec(1.0, -2.0, 0.5, 3.0)
    val observations = tall.*(truth)(using PureBackend)
    val solved = tall.qr(using backend).solveLeastSquares(observations).orThrow
    assert((solved - truth).norm2 <= 2e-10)

  test("symmetric eigen obeys residual, orthogonality, values-only, and scaling contracts"):
    requireLapack()
    val a = Matrix.dense(4, 4)(
      5.0, 2.0, -1.0, 0.5,
      2.0, 4.0, 1.5, -0.25,
      -1.0, 1.5, 3.0, 0.75,
      0.5, -0.25, 0.75, 2.0
    )
    val spectral = backend.spectral.get
    val raw = spectral.denseSymmetricEigen(a, wantVectors = true).orThrow
    assertEquals(raw.values.length, 4)
    assert(raw.values.toSeq.sliding(2).forall(pair => pair(0) <= pair(1)))
    assertMatrixClose(raw.vectors.t.*(raw.vectors)(using PureBackend), DMat.eye(4), 2e-10)
    val diagonal = Matrix.tabulate(4, 4)((i, j) => if i == j then raw.values(i) else 0.0)
    assertMatrixClose(a.*(raw.vectors)(using PureBackend), raw.vectors.*(diagonal)(using PureBackend), 3e-10)

    val valuesOnly = spectral.denseSymmetricEigen(a, wantVectors = false).orThrow
    assertEquals(valuesOnly.vectors.shape, Shape(Rows(4), Cols(0)))
    val scaledInput = Matrix.tabulate(a.rows, a.cols)((i, j) => 7.0 * a(i, j))
    val scaled = spectral.denseSymmetricEigen(scaledInput, wantVectors = false).orThrow
    var i = 0
    while i < 4 do
      assertEqualsDouble(scaled.values(i), 7.0 * raw.values(i), 2e-11, s"eigenvalue $i")
      i += 1

  test("symmetric eigen and Cholesky read the lower triangle only"):
    requireLapack()
    val clean = Matrix.dense(3, 3)(4.0, 1.0, 0.5, 1.0, 3.0, -0.25, 0.5, -0.25, 2.0)
    val garbled = Matrix.tabulate(3, 3)((i, j) => if i < j then 1000.0 + i * 10 + j else clean(i, j))
    val cleanEig = backend.spectral.get.denseSymmetricEigen(clean, wantVectors = false).orThrow.values
    val garbledEig = backend.spectral.get.denseSymmetricEigen(garbled, wantVectors = false).orThrow.values
    assertEquals(cleanEig.toSeq, garbledEig.toSeq)
    assertMatrixClose(
      backend.denseFactorizations.get.cholesky(clean).orThrow.lower,
      backend.denseFactorizations.get.cholesky(garbled).orThrow.lower
    )

  test("shared backend supports concurrent independent LAPACK calls"):
    requireLapack()
    import scala.concurrent.ExecutionContext.Implicits.global
    val jobs = (0 until 16).map: seed =>
      Future:
        val a = Matrix.tabulate(12, 12)((i, j) =>
          if i == j then 20.0 + seed else ((i * 7 + j * 3 + seed) % 11 - 5).toDouble / 20.0
        )
        val chol = backend.denseFactorizations.get.cholesky(a).orThrow
        chol.lower(0, 0)
    val results = Await.result(Future.sequence(jobs), 30.seconds)
    assert(results.forall(_.isFinite))
