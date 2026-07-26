package gale.laws

import gale.linalg.*
import gale.solvers.Preconditioner
import gale.spectral.*
import munit.Assertions

/** Reusable capability law for iterative generalized spectral providers.
  *
  * The law exercises a small analytic SPD pencil through both the raw provider
  * method and Gale's public facade. It rejects a capability that merely inherits
  * the default `UnsupportedOperation`, malformed raw carrier metadata/shapes,
  * non-finite factors, and factors that the facade cannot validate against the
  * true operator residual and B inner product.
  */
object SpectralBackendLaws extends Assertions:

  def iterativeGeneralizedCapability(backend: SpectralBackend): Unit =
    if backend.capabilities.contains(SpectralCapability.IterativeGeneralized) then
      val n = 4
      val k = 2
      val a = diagonalOperator(IndexedSeq(1.0, 2.0, 4.0, 8.0))
      val b = diagonalOperator(IndexedSeq.fill(n)(1.0))
      val order = EigenOrder.SmallestAlgebraic
      val options =
        GeneralizedSpectralOptions(tolerance = 1e-8, maxIterations = 40)

      backend.iterativeGeneralizedEigen(
        a,
        b,
        n,
        k,
        order,
        options,
        Preconditioner.Identity
      ) match
        case Left(_: LinAlgError.UnsupportedOperation) =>
          fail(
            s"backend '${backend.name}' advertises IterativeGeneralized but inherits or returns UnsupportedOperation"
          )
        case Left(_) =>
          // A backend may decline a representative fixture for another typed
          // reason; the public facade's fallback policy handles that case.
          ()
        case Right(raw) =>
          assertCarrier(raw, backend.name, n, k)
          try
            Eigen.eigSymmetricGeneralized(
              a.assumeSymmetricOperator,
              b.assumePositiveDefiniteOperator,
              n,
              EigenSelection.Count(k, order),
              options
            )(using backend) match
              case Left(error) =>
                fail(
                  s"backend '${backend.name}' returned raw factors but the public facade returned $error"
                )
              case Right(result) =>
                assertEquals(result.diagnostics.requested, k)
                assertEquals(result.diagnostics.converged, result.size)
                var i = 0
                while i < result.size do
                  assert(
                    result.diagnostics.residuals(i) <= options.tolerance,
                    s"backend '${backend.name}' residual $i exceeds tolerance"
                  )
                  i += 1
                if result.eigenvectors.cols > 0 then
                  val block = result.eigenvectors
                  val gram = block.t * block
                  assert(
                    frobeniusIdentityError(gram) <= 1e-8,
                    s"backend '${backend.name}' returned non-B-orthonormal vectors"
                  )
          catch
            case error: LinAlgError.InvalidArgument =>
              fail(
                s"backend '${backend.name}' violates the IterativeGeneralized raw carrier contract: ${error.getMessage}"
              )

  private def assertCarrier(
      raw: RawIterativeGeneralizedEigen,
      backendName: String,
      n: Int,
      k: Int
  ): Unit =
    val convergence = raw.convergence
    assertEquals(
      convergence.requested,
      k,
      s"backend '$backendName' raw requested count"
    )
    assert(
      convergence.converged >= 0 && convergence.converged <= k,
      s"backend '$backendName' raw converged count ${convergence.converged}"
    )
    assert(
      convergence.iterations >= 0,
      s"backend '$backendName' returned negative iterations"
    )
    assertEquals(
      raw.values.length,
      convergence.converged,
      s"backend '$backendName' raw value count"
    )
    assertEquals(raw.vectors.rows, n, s"backend '$backendName' raw vector rows")
    assert(
      raw.vectors.cols >= convergence.converged,
      s"backend '$backendName' raw vector columns ${raw.vectors.cols}"
    )

    var col = 0
    while col < convergence.converged do
      assert(raw.values(col).isFinite, s"backend '$backendName' non-finite value $col")
      var row = 0
      while row < n do
        assert(
          raw.vectors(row, col).isFinite,
          s"backend '$backendName' non-finite vector entry ($row, $col)"
        )
        row += 1
      col += 1

  private def diagonalOperator(diagonal: IndexedSeq[Double]): DoubleLinearOperator =
    LinearOperator.fromFunction(diagonal.length, diagonal.length): (x, into) =>
      var i = 0
      while i < diagonal.length do
        into(i) = diagonal(i) * x(i)
        i += 1

  private def frobeniusIdentityError(gram: DMat): Double =
    var sum = 0.0
    var row = 0
    while row < gram.rows do
      var col = 0
      while col < gram.cols do
        val delta = gram(row, col) - (if row == col then 1.0 else 0.0)
        sum += delta * delta
        col += 1
      row += 1
    math.sqrt(sum)
