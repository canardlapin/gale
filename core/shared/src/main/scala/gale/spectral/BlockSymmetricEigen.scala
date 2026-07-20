package gale.spectral

import gale.linalg.DMat
import gale.linalg.DVec
import gale.linalg.DoubleLinearOperator
import gale.linalg.LinAlgError
import gale.linalg.MutableDVec
import scala.collection.mutable.ArrayBuffer

/** Portable multiplicity-safe block Krylov solver for partial symmetric
  * eigenproblems.
  *
  * A block at least as wide as the requested result is used so a repeated
  * eigenspace can contribute more than one independent direction. Each restart
  * retains the wanted Ritz space (soft locking) and expands it again through the
  * operator (a thick restart). Full reorthogonalization keeps the projected
  * problem symmetric and the returned basis orthonormal. If a Krylov component
  * becomes invariant before the requested subspace is full, a deterministic
  * orthogonal probe replenishes the block; this is important for disconnected
  * operators and exact multiplicities.
  *
  * The implementation only invokes [[DoubleLinearOperator.applyTo]] and never
  * materializes the input operator. The small Rayleigh-Ritz projection is dense,
  * as it is for every block Krylov method.
  */
private[spectral] object BlockSymmetricEigen:

  private final case class Basis(
      vectors: Array[DVec],
      images: Array[DVec]
  )

  private final case class RitzPair(
      value: Double,
      vector: DVec,
      image: DVec,
      residual: DVec,
      residualNorm: Double
  )

  def solve(
      op: DoubleLinearOperator,
      n: Int,
      k: Int,
      order: EigenOrder,
      options: SpectralOptions,
      wantVectors: Boolean,
      firstStart: DVec
  ): Either[LinAlgError, EigenDecomposition] =
    var subspaceSize =
      math.min(
        n,
        math.max(
          options.subspaceDimension.getOrElse(math.max(2 * k + 1, 20)),
          k + 1
        )
      )
    val blockWidth = math.min(n, math.max(k, 2))
    val growBy = math.max(k, 16)
    val maxRestarts = math.max(options.maxIterations, 1)
    val tolerance = options.tolerance

    var restart = 0
    var seeds = initialSeeds(firstStart, n, blockWidth)
    var done = false
    var failure: Option[LinAlgError] = None
    var result = assemble(
      n,
      k,
      Array.empty[RitzPair],
      wantVectors,
      iterations = 0,
      extremalityCertified = false
    )

    while restart < maxRestarts && !done && failure.isEmpty do
      val basis = buildBasis(op, n, subspaceSize, blockWidth, seeds, restart)
      val projected = projectedMatrix(basis)
      DenseSpectralKernels.symmetricEigen(projected, wantVectors = true) match
        case Left(DenseSpectralKernels.SpectralKernelFailure.DidNotConverge(iterations)) =>
          failure = Some(LinAlgError.DidNotConverge(iterations, 0.0))
        case Right(eigen) =>
          val projectedVectors = eigen.vectors.get
          val wanted = selectExtremeIndices(eigen.values, math.min(k, basis.vectors.length), order)
          val pairs = wanted.map: index =>
            ritzPair(basis, eigen.values(index), projectedVectors, index)
          val scale = math.max(1.0, maxAbs(eigen.values))
          val converged = pairs.filter(_.residualNorm <= tolerance * scale)
          result = assemble(
            n,
            k,
            converged,
            wantVectors,
            restart + 1,
            extremalityCertified = basis.vectors.length == n
          )

          val fullyConverged = pairs.length == k && converged.length == k
          // A full orthonormal basis makes the Rayleigh-Ritz problem equivalent
          // to the original operator. Repeating it cannot improve a deliberately
          // impossible tolerance such as zero.
          if fullyConverged || basis.vectors.length == n then done = true
          else
            val retainCount =
              math.min(
                basis.vectors.length - 1,
                math.max(k, math.min(2 * k, basis.vectors.length - 1))
              )
            val retained = selectExtremeIndices(eigen.values, retainCount, order)
            seeds = retained.map: index =>
              ritzVector(basis.vectors, projectedVectors, index)
            subspaceSize = math.min(n, subspaceSize + growBy)
      restart += 1

    failure match
      case Some(error) => Left(error)
      case None        => Right(result)

  /** Build a block Krylov basis, retaining supplied thick-restart vectors first.
    * Every stored operator image aligns with the basis vector at the same index.
    */
  private def buildBasis(
      op: DoubleLinearOperator,
      n: Int,
      target: Int,
      blockWidth: Int,
      seeds: Array[DVec],
      restart: Int
  ): Basis =
    val vectors = ArrayBuffer.empty[DVec]
    val images = ArrayBuffer.empty[Option[DVec]]

    def add(candidate: DVec): Boolean =
      orthonormalized(candidate, vectors) match
        case Some(vector) =>
          vectors += vector
          images += None
          true
        case None => false

    var seedIndex = 0
    while seedIndex < seeds.length && vectors.length < target do
      add(seeds(seedIndex))
      seedIndex += 1

    var stream = restart * math.max(blockWidth, 1) + 1
    while vectors.length < math.min(blockWidth, target) do
      if !add(deterministicVector(n, stream)) then
        addCanonicalComplement(n, vectors, images)
      stream += 1

    var frontier = 0
    while vectors.length < target do
      if frontier < vectors.length then
        val image = apply(op, vectors(frontier))
        images(frontier) = Some(image)
        add(image)
        frontier += 1
      else
        // The current Krylov component is invariant. Probe its orthogonal
        // complement instead of declaring a happy breakdown that could hide a
        // repeated eigenvector, then continue the block expansion there.
        val added = add(deterministicVector(n, stream))
        stream += 1
        if !added && !addCanonicalComplement(n, vectors, images) then
          // This is reachable only after spanning R^n.
          return finishImages(op, vectors, images)

    finishImages(op, vectors, images)

  private def finishImages(
      op: DoubleLinearOperator,
      vectors: ArrayBuffer[DVec],
      images: ArrayBuffer[Option[DVec]]
  ): Basis =
    var i = 0
    while i < vectors.length do
      if images(i).isEmpty then images(i) = Some(apply(op, vectors(i)))
      i += 1
    Basis(vectors.toArray, images.map(_.get).toArray)

  private def addCanonicalComplement(
      n: Int,
      vectors: ArrayBuffer[DVec],
      images: ArrayBuffer[Option[DVec]]
  ): Boolean =
    var coordinate = 0
    while coordinate < n do
      val unit = DVec.tabulate(n)(i => if i == coordinate then 1.0 else 0.0)
      orthonormalized(unit, vectors) match
        case Some(vector) =>
          vectors += vector
          images += None
          return true
        case None => ()
      coordinate += 1
    false

  /** Classical Gram-Schmidt twice. The second pass is essential for clustered
    * Ritz vectors and exact repeated roots.
    */
  private def orthonormalized(
      candidate: DVec,
      basis: collection.IndexedSeq[DVec]
  ): Option[DVec] =
    val work = candidate.mutableCopy
    var pass = 0
    while pass < 2 do
      var i = 0
      while i < basis.length do
        val coefficient = basis(i).dot(work.asVec)
        work.axpyInPlace(-coefficient, basis(i))
        i += 1
      pass += 1
    val norm = work.asVec.norm2
    val threshold = 1e-12 * math.max(1.0, candidate.norm2)
    if !norm.isFinite || norm <= threshold then None
    else Some((work.asVec * (1.0 / norm)).copy)

  private def projectedMatrix(basis: Basis): DMat =
    val m = basis.vectors.length
    DMat.tabulate(m, m): (i, j) =>
      // Symmetrize the two mathematically equal dot products so roundoff or a
      // slightly asymmetric operator cannot leak an asymmetric projected input
      // into the dense symmetric kernel.
      0.5 * (
        basis.vectors(i).dot(basis.images(j)) +
          basis.vectors(j).dot(basis.images(i))
      )

  private def ritzPair(
      basis: Basis,
      value: Double,
      projectedVectors: DMat,
      column: Int
  ): RitzPair =
    val vector0 = linearCombination(basis.vectors, projectedVectors, column)
    val image0 = linearCombination(basis.images, projectedVectors, column)
    val norm = vector0.norm2
    val vector = if norm == 0.0 then vector0 else vector0 * (1.0 / norm)
    val image = if norm == 0.0 then image0 else image0 * (1.0 / norm)
    val residualMutable = image.mutableCopy
    residualMutable.axpyInPlace(-value, vector)
    val residual = residualMutable.asVec.copy
    RitzPair(value, vector.copy, image.copy, residual, residual.norm2)

  private def ritzVector(
      basis: Array[DVec],
      projectedVectors: DMat,
      column: Int
  ): DVec =
    val vector = linearCombination(basis, projectedVectors, column)
    val norm = vector.norm2
    if norm == 0.0 then vector else (vector * (1.0 / norm)).copy

  private def linearCombination(
      basis: Array[DVec],
      coefficients: DMat,
      column: Int
  ): DVec =
    val out = MutableDVec.zeros(basis(0).length)
    var j = 0
    while j < basis.length do
      out.axpyInPlace(coefficients(j, column), basis(j))
      j += 1
    out.asVec.copy

  private def apply(op: DoubleLinearOperator, vector: DVec): DVec =
    val out = MutableDVec.zeros(op.rows)
    op.applyTo(vector, out)
    out.asVec.copy

  private def initialSeeds(first: DVec, n: Int, width: Int): Array[DVec] =
    val out = new Array[DVec](width)
    out(0) = first
    var i = 1
    while i < width do
      out(i) = deterministicVector(n, i)
      i += 1
    out

  /** Portable deterministic block probe. Distinct streams have distinct LCG
    * initial states; integer wraparound is identical on JVM and Scala.js.
    */
  private def deterministicVector(n: Int, stream: Int): DVec =
    var state = 123456789 ^ (stream * 0x9e3779b9)
    DVec.tabulate(n): _ =>
      state = state * 1103515245 + 12345
      ((state >>> 9) & 0x7fffff).toDouble / 0x800000.toDouble * 2.0 - 1.0

  /** Ascending indices for the requested extreme. Projected values are already
    * ascending; ties retain projected index order and therefore multiplicity.
    */
  private def selectExtremeIndices(values: DVec, k: Int, order: EigenOrder): Array[Int] =
    val n = values.length
    val selected =
      order match
        case EigenOrder.SmallestAlgebraic => (0 until k).toArray
        case EigenOrder.LargestAlgebraic  => ((n - k) until n).toArray
        case EigenOrder.SmallestMagnitude =>
          (0 until n).sortBy(i => (math.abs(values(i)), i)).take(k).toArray
        case EigenOrder.LargestMagnitude =>
          (0 until n).sortBy(i => (-math.abs(values(i)), i)).take(k).toArray
        case EigenOrder.BothEnds =>
          val high = (k + 1) / 2
          val low = k / 2
          ((0 until low) ++ ((n - high) until n)).toArray
        case EigenOrder.LargestRealPart | EigenOrder.SmallestRealPart =>
          Array.empty[Int]
    java.util.Arrays.sort(selected)
    selected

  private def assemble(
      n: Int,
      requested: Int,
      pairs: Array[RitzPair],
      wantVectors: Boolean,
      iterations: Int,
      extremalityCertified: Boolean
  ): EigenDecomposition =
    val order = pairs.indices.sortBy(i => (pairs(i).value, i)).toArray
    val values = DVec.tabulate(order.length)(i => pairs(order(i)).value)
    val residuals = DVec.tabulate(order.length)(i => pairs(order(i)).residualNorm)
    val vectors =
      if wantVectors && order.nonEmpty then
        DMat.tabulate(n, order.length)((row, col) => pairs(order(col)).vector(row))
      else DMat.zeros(n, 0)
    val diagnostics = SpectralDiagnostics(
      requested = requested,
      converged = order.length,
      residuals = residuals,
      orthogonalityError = if wantVectors then orthogonalityError(vectors) else 0.0,
      iterations = iterations,
      rank = None,
      extremalityCertified = extremalityCertified
    )
    EigenDecomposition(values, vectors, diagnostics)

  private def orthogonalityError(vectors: DMat): Double =
    val gram = vectors.t * vectors
    var sum = 0.0
    var i = 0
    while i < gram.rows do
      var j = 0
      while j < gram.cols do
        val delta = if i == j then gram(i, j) - 1.0 else gram(i, j)
        sum += delta * delta
        j += 1
      i += 1
    math.sqrt(sum)

  private def maxAbs(values: DVec): Double =
    var maximum = 0.0
    var i = 0
    while i < values.length do
      maximum = math.max(maximum, math.abs(values(i)))
      i += 1
    maximum
