package gale.spectral

import gale.linalg.*

/** Bit-exact regression guard for the dense nonsymmetric kernel
  * ([[DenseSpectralKernels.nonsymmetricEigen]] and
  * [[DenseSpectralKernels.hessenberg]]) on the JVM.
  *
  * Fixture transcendental functions use StrictMath so identical inputs and objective values
  * are used on ARM and x86 JVMs. This constrains the test corpus, not production arithmetic.
  *
  * Each case of a fixed corpus is reduced to one 64-bit digest of the raw bits
  * the kernel returns (eigenvalue parts, packed vectors, or the failure value;
  * `H` and `Q` for the Hessenberg reduction) and compared with the digests in
  * `nonsymmetric-kernel-digests.txt`. The stored digests were generated from
  * the kernel as of breeze/integration 64ba31e, so any change to its
  * arithmetic or summation order fails here; such a change must regenerate the
  * file deliberately and justify the new bits. JVM-only: the kernels are
  * deterministic per platform, not across platforms.
  */
class NonsymmetricKernelDigestSuite extends munit.FunSuite:
  import NonsymmetricKernelDigestSuite.*

  private lazy val expected: Map[String, String] =
    val stream = getClass.getResourceAsStream(DigestResource)
    assert(stream != null, s"missing resource $DigestResource")
    val source = scala.io.Source.fromInputStream(stream, "UTF-8")
    try
      source.getLines
        .filter(line => line.nonEmpty && !line.startsWith("#"))
        .map: line =>
          val tab = line.lastIndexOf('\t')
          line.substring(0, tab) -> line.substring(tab + 1)
        .toMap
    finally source.close()

  test("nonsymmetricEigen matches the stored per-case bit digests") {
    val actual = eigenDigests(c => (vec, budget) => DenseSpectralKernels.nonsymmetricEigen(c, vec, budget))
    assert(actual.size > 2000, s"corpus shrank to ${actual.size} cases")
    assertEquals(actual.keySet, expected.keySet.filter(_.startsWith("eig ")))
    val mismatches = actual.filter((key, digest) => !expected.get(key).contains(digest)).keys.toSeq.sorted
    assert(mismatches.isEmpty, s"${mismatches.size} cases changed bits, e.g. ${mismatches.take(5).mkString("; ")}")
  }

  test("hessenberg matches the stored per-case bit digests") {
    val actual = hessenbergDigests((a, wantQ) =>
      val result = DenseSpectralKernels.hessenberg(a, wantQ)
      (result.h, result.q)
    )
    assert(actual.size > 1000, s"corpus shrank to ${actual.size} cases")
    assertEquals(actual.keySet, expected.keySet.filter(_.startsWith("hessenberg ")))
    val mismatches = actual.filter((key, digest) => !expected.get(key).contains(digest)).keys.toSeq.sorted
    assert(mismatches.isEmpty, s"${mismatches.size} cases changed bits, e.g. ${mismatches.take(5).mkString("; ")}")
  }

object NonsymmetricKernelDigestSuite:
  val DigestResource = "/gale/spectral/nonsymmetric-kernel-digests.txt"

  /** The corpus: Gaussian, small-integer, triangular, companion, rotation-block,
    * Jordan, zero, `1e±300`-scaled, cyclic-permutation and graded matrices for
    * `n = 0..40, 47, 64, 65`.
    */
  def corpus: Seq[(String, DMat)] =
    val out = scala.collection.mutable.ArrayBuffer.empty[(String, DMat)]
    for n <- (0 to 40) ++ Seq(47, 64, 65) do
      for seed <- 0 until (if n <= 40 then 4 else 1) do
        val r = new scala.util.Random(1000L * n + seed)
        out += (s"gauss n=$n s=$seed" -> Matrix.tabulate(n, n)((_, _) => r.nextGaussian()))
      val r = new scala.util.Random(77L + n)
      out += (s"ints n=$n" -> Matrix.tabulate(n, n)((_, _) => (r.nextInt(5) - 2).toDouble))
      out += (s"upper n=$n" -> Matrix.tabulate(n, n)((i, j) => if j >= i then (i + 2 * j + 1).toDouble else 0.0))
      out += (s"companion n=$n" -> Matrix.tabulate(n, n)((i, j) =>
        if i == 0 then -1.0 / (j + 1) else if i == j + 1 then 1.0 else 0.0
      ))
      out += (s"rotations n=$n" -> Matrix.tabulate(n, n)((i, j) =>
        if i == j then 0.5 * (i / 2) else if i / 2 == j / 2 then (if i < j then 1.0 else -1.0) else 0.0
      ))
      out += (s"jordan n=$n" -> Matrix.tabulate(n, n)((i, j) => if i == j then 2.0 else if j == i + 1 then 1.0 else 0.0))
      out += (s"zero n=$n" -> Matrix.zeros(n, n))
      val r3 = new scala.util.Random(31L + n)
      out += (s"huge n=$n" -> Matrix.tabulate(n, n)((_, _) => r3.nextGaussian() * 1e300))
      out += (s"tiny n=$n" -> Matrix.tabulate(n, n)((_, _) => r3.nextGaussian() * 1e-300))
      out += (s"perm n=$n" -> Matrix.tabulate(n, n)((i, j) => if j == (i + 1) % math.max(n, 1) then 1.0 else 0.0))
      val r2 = new scala.util.Random(9L + n)
      out += (s"scaled n=$n" -> Matrix.tabulate(n, n)((i, j) => r2.nextGaussian() * StrictMath.pow(10.0, (i - j) * 3.0)))
    out.toSeq

  /** FNV-1a over 64-bit words, each word's bytes folded high to low. */
  final class Digest:
    private var state = 0xcbf29ce484222325L
    def long(x: Long): Digest =
      var shift = 56
      while shift >= 0 do
        state = (state ^ ((x >>> shift) & 0xffL)) * 0x100000001b3L
        shift -= 8
      this
    def double(x: Double): Digest = long(java.lang.Double.doubleToRawLongBits(x))
    def vec(v: DVec): Digest =
      long(v.length.toLong)
      var i = 0
      while i < v.length do
        double(v(i))
        i += 1
      this
    def mat(m: DMat): Digest =
      long(m.rows.toLong).long(m.cols.toLong)
      var i = 0
      while i < m.rows do
        var j = 0
        while j < m.cols do
          double(m(i, j))
          j += 1
        i += 1
      this
    def hex: String = f"$state%016x"

  def eigenDigest(
      result: Either[DenseSpectralKernels.SpectralKernelFailure, DenseSpectralKernels.NonsymmetricEigen]
  ): String =
    result match
      case Right(e) =>
        val d = new Digest().long(1L).vec(e.re).vec(e.im)
        e.vectors match
          case Some(v) => d.long(1L).mat(v).hex
          case None    => d.long(0L).hex
      case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(sweeps)) =>
        new Digest().long(2L).long(sweeps.toLong).hex

  /** `key -> digest` for every corpus case, with vectors on and off and with the
    * default and a 3-sweep budget.
    */
  def eigenDigests(
      solve: DMat => (Boolean, Int) => Either[
        DenseSpectralKernels.SpectralKernelFailure,
        DenseSpectralKernels.NonsymmetricEigen
      ]
  ): Map[String, String] =
    (for
      (label, a) <- corpus
      vec <- Seq(true, false)
      budget <- Seq(-1, 3)
    yield s"eig $label vec=$vec budget=$budget" -> eigenDigest(solve(a)(vec, budget))).toMap

  /** `key -> digest` of `H` and `Q` for every nonempty corpus case. */
  def hessenbergDigests(reduce: (DMat, Boolean) => (DMat, Option[DMat])): Map[String, String] =
    (for
      (label, a) <- corpus if a.rows > 0
      wantQ <- Seq(true, false)
    yield
      val (h, q) = reduce(a, wantQ)
      val d = new Digest().mat(h)
      s"hessenberg $label wantQ=$wantQ" -> q.fold(d.long(0L))(m => d.long(1L).mat(m)).hex
    ).toMap
