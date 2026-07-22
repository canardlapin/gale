package gale.solvers

import gale.linalg.*

/** Shared configuration for the iterative solvers.
  *
  *   - `tolerance` is the convergence threshold on a residual norm. What residual,
  *     and whether the threshold is absolute or relative, is set by the solver's
  *     [[ToleranceMode]]: `cg`/`bicgstab`/`gmres` test the linear residual
  *     `||b − A x||`, while `cgnr`/`lsqr` test the least-squares (normal-equation)
  *     residual `||Aᵀ(b − A x)||`. `RelativeToRhs` scales the threshold by the
  *     norm of the relevant right-hand side (`||b||`, or `||Aᵀ b||` for the
  *     least-squares solvers).
  *   - `maxIterations` caps the total iteration count; on reaching it the solver
  *     returns [[SolverResult.NotConverged]] with the best iterate so far.
  *   - `restart` is the GMRES restart length (inner Krylov dimension); ignored by
  *     the other solvers.
  */
final case class SolverConfig(
    tolerance: Double = 1e-10,
    maxIterations: Int = 1000,
    restart: Int = 30
)

/** How [[SolverConfig.tolerance]] is compared against the residual norm.
  *
  * `Absolute` tests `||r|| <= tolerance` directly (the historical behaviour);
  * `RelativeToRhs` tests `||r|| <= tolerance * ||b||`, i.e. a relative residual.
  */
enum ToleranceMode:
  case Absolute, RelativeToRhs

sealed trait SolverResult:
  def x: DVec
  def iterations: Int
  def residual: Double
  def converged: Boolean

  def orThrow: DVec =
    if converged then x else throw LinAlgError.DidNotConverge(iterations, residual)

object SolverResult:
  final case class Converged(x: DVec, iterations: Int, residual: Double) extends SolverResult:
    val converged: Boolean = true

  final case class NotConverged(x: DVec, iterations: Int, residual: Double) extends SolverResult:
    val converged: Boolean = false

/** Caller-owned storage for repeated conjugate-gradient solves of one fixed
  * dimension.
  *
  * [[solution]] is an immutable-facing view of mutable workspace storage. It is
  * allocation-free to read, but a later [[IterativeSolvers.cgWith]] call using
  * this workspace overwrites the values visible through earlier views. Use
  * [[solutionCopy]] when a result must outlive workspace reuse.
  */
final class CgWorkspace private (val size: Int):
  private[solvers] val x = MutableDVec.zeros(size)
  private[solvers] val ax = MutableDVec.zeros(size)
  private[solvers] val r = MutableDVec.zeros(size)
  private[solvers] val z = MutableDVec.zeros(size)
  private[solvers] val p = MutableDVec.zeros(size)
  private[solvers] val ap = MutableDVec.zeros(size)

  private[solvers] val xView = x.asVec
  private[solvers] val rView = r.asVec
  private[solvers] val zView = z.asVec
  private[solvers] val pView = p.asVec
  private[solvers] val apView = ap.asVec

  private var iterationCount = 0
  private var residualValue = Double.PositiveInfinity
  private var convergedValue = false

  def solution: DVec = xView
  def solutionCopy: DVec = xView.copy
  def iterations: Int = iterationCount
  def residual: Double = residualValue
  def converged: Boolean = convergedValue

  private[solvers] def finish(iterations: Int, residual: Double, converged: Boolean): Unit =
    iterationCount = iterations
    residualValue = residual
    convergedValue = converged

object CgWorkspace:
  def apply(size: Int): CgWorkspace =
    require(size >= 0, "CG workspace size must be non-negative")
    new CgWorkspace(size)

trait Preconditioner:
  def solve(r: DVec, into: MutableVec[Double]): Unit

  def apply(r: DVec): DVec =
    val out = MutableDVec.zeros(r.length)
    solve(r, out)
    out.asVec

object Preconditioner:
  object Identity extends Preconditioner:
    def solve(r: DVec, into: MutableVec[Double]): Unit =
      if r.length != into.length then
        throw LinAlgError.VectorLengthMismatch(r.length, into.length)
      var i = 0
      while i < r.length do
        into(i) = r(i)
        i += 1

  // Storage note: `inverseDiagonal` is a plain `Array[Double]` held privately in
  // shared code. P4 keeps raw `Double` arrays off the public surface, but this one
  // never escapes — it is private, read-only after construction, and only ever
  // consumed element-by-element by `solve` — so it stays as a plain array here
  // rather than a platform `DoubleArray`.
  final class Jacobi private[solvers] (private val inverseDiagonal: Array[Double]) extends Preconditioner:
    def solve(r: DVec, into: MutableVec[Double]): Unit =
      if r.length != inverseDiagonal.length then
        throw LinAlgError.VectorLengthMismatch(inverseDiagonal.length, r.length)
      if into.length != inverseDiagonal.length then
        throw LinAlgError.VectorLengthMismatch(inverseDiagonal.length, into.length)
      var i = 0
      while i < inverseDiagonal.length do
        into(i) = inverseDiagonal(i) * r(i)
        i += 1

  object Jacobi:
    def apply(A: Matrix[Double]): Jacobi =
      if A.rows != A.cols then
        throw LinAlgError.NonSquareMatrix(A.shape)
      val inv = new Array[Double](A.rows)
      var i = 0
      while i < A.rows do
        val diag = A(i, i)
        if diag == 0.0 then
          throw LinAlgError.SingularMatrix(i)
        inv(i) = 1.0 / diag
        i += 1
      new Jacobi(inv)

  /** Block-diagonal preconditioner: `M` is the block diagonal of `A`, applied
    * by exactly solving each diagonal block. Reduces to point [[Jacobi]] when
    * `blockSize == 1` and to a direct solve of `A` when `blockSize >= n`.
    */
  final class BlockJacobi private[solvers] (
      private val blockStarts: Array[Int],
      private val blockFactors: Array[LU],
      private val dimension: Int
  ) extends Preconditioner:
    def solve(r: DVec, into: MutableVec[Double]): Unit =
      if r.length != dimension then
        throw LinAlgError.VectorLengthMismatch(dimension, r.length)
      if into.length != dimension then
        throw LinAlgError.VectorLengthMismatch(dimension, into.length)
      val rData = r.data
      val rBase = r.offset.value
      val rStep = r.stride.value
      var b = 0
      while b < blockStarts.length do
        val start = blockStarts(b)
        val end = if b + 1 < blockStarts.length then blockStarts(b + 1) else dimension
        val size = end - start
        val rhs = new Array[Double](size)
        var i = 0
        var idx = rBase + start * rStep
        while i < size do
          rhs(i) = rData(idx)
          idx += rStep
          i += 1
        val solution =
          DenseDecompositions.solve(blockFactors(b), DVec.fromArray(rhs)) match
            case Right(x)    => x.toArray
            case Left(error) => throw error
        i = 0
        while i < size do
          into(start + i) = solution(i)
          i += 1
        b += 1

  object BlockJacobi:
    def apply(A: Matrix[Double], blockSize: Int): Preconditioner =
      require(blockSize > 0, "blockSize must be positive")
      if A.rows != A.cols then
        throw LinAlgError.NonSquareMatrix(A.shape)
      val n = A.rows
      if n == 0 then
        new BlockJacobi(Array.empty[Int], Array.empty[LU], 0)
      else
        val numBlocks = (n + blockSize - 1) / blockSize
        val starts = new Array[Int](numBlocks)
        val factors = new Array[LU](numBlocks)
        var b = 0
        while b < numBlocks do
          val start = b * blockSize
          val end = math.min(start + blockSize, n)
          val size = end - start
          starts(b) = start
          // Extract the diagonal block into a row-major dense array and LU-factor
          // it once; a singular block is fatal at construction with a global index.
          val block = new Array[Double](size * size)
          var i = 0
          while i < size do
            var j = 0
            while j < size do
              block(i * size + j) = A(start + i, start + j)
              j += 1
            i += 1
          DenseDecompositions.lu(DMat.fromArrayRowMajor(size, size, block)) match
            case Right(factor) =>
              factors(b) = factor
            case Left(LinAlgError.SingularMatrix(k)) =>
              throw LinAlgError.SingularMatrix(start + k)
            case Left(error) =>
              throw error
          b += 1
        new BlockJacobi(starts, factors, n)

  final class SymmetricGaussSeidel private[solvers] (A: Matrix[Double]) extends Preconditioner:
    def solve(r: DVec, into: MutableVec[Double]): Unit =
      if r.length != A.rows || into.length != A.rows then
        throw LinAlgError.DimensionMismatch(Shape(Rows(A.rows), Cols(1)), Shape(Rows(r.length), Cols(1)))

      var i = 0
      while i < A.rows do
        var sum = r(i)
        var j = 0
        while j < i do
          sum -= A(i, j) * into(j)
          j += 1
        into(i) = sum / A(i, i)
        i += 1

      i = A.rows - 1
      while i >= 0 do
        var sum = into(i) * A(i, i)
        var j = i + 1
        while j < A.cols do
          sum -= A(i, j) * into(j)
          j += 1
        into(i) = sum / A(i, i)
        i -= 1

  object SymmetricGaussSeidel:
    def apply(A: Matrix[Double]): SymmetricGaussSeidel =
      if A.rows != A.cols then
        throw LinAlgError.NonSquareMatrix(A.shape)
      var i = 0
      while i < A.rows do
        if A(i, i) == 0.0 then
          throw LinAlgError.SingularMatrix(i)
        i += 1
      new SymmetricGaussSeidel(A)

object IterativeSolvers:
  def cg(
      A: DoubleLinearOperator,
      b: DVec,
      config: SolverConfig = SolverConfig(),
      preconditioner: Preconditioner = Preconditioner.Identity,
      initial: Option[DVec] = None,
      toleranceMode: ToleranceMode = ToleranceMode.Absolute
  ): SolverResult =
    requireSquare(A)
    if b.length != A.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(A.rows), Cols(1)), Shape(Rows(b.length), Cols(1)))
    val tol = effectiveTolerance(toleranceMode, config.tolerance, b)
    val x = initialGuess(initial, A.cols)
    val ax = MutableDVec.zeros(A.rows)
    A.applyTo(x.asVec, ax)
    val r = (b - ax.asVec).mutableCopy
    var residual = r.asVec.norm2
    if residual <= tol then
      return SolverResult.Converged(x.asVec, 0, residual)

    val z = MutableDVec.zeros(A.rows)
    preconditioner.solve(r.asVec, z)
    val p = z.asVec.mutableCopy
    var rzOld = r.asVec.dot(z.asVec)
    var iteration = 0
    while iteration < config.maxIterations do
      val ap = A * p.asVec
      val denom = p.asVec.dot(ap)
      if denom == 0.0 then
        return SolverResult.NotConverged(x.asVec, iteration, residual)
      val alpha = rzOld / denom
      x.axpyInPlace(alpha, p.asVec)
      r.axpyInPlace(-alpha, ap)
      residual = r.asVec.norm2
      iteration += 1
      if residual <= tol then
        return SolverResult.Converged(x.asVec, iteration, residual)
      preconditioner.solve(r.asVec, z)
      val rzNew = r.asVec.dot(z.asVec)
      if rzOld == 0.0 then
        return SolverResult.NotConverged(x.asVec, iteration, residual)
      val beta = rzNew / rzOld
      p *= beta
      p += z.asVec
      rzOld = rzNew
    SolverResult.NotConverged(x.asVec, iteration, residual)

  /** Allocation-controlled conjugate gradient using caller-owned workspace.
    *
    * The recurrence, tolerance modes, preconditioner semantics, and breakdown
    * behavior match [[cg]]. The returned object is the supplied workspace; its
    * [[CgWorkspace.solution]] view is overwritten by the next solve using the
    * same workspace.
    */
  def cgWith(
      A: DoubleLinearOperator,
      b: DVec,
      workspace: CgWorkspace,
      config: SolverConfig = SolverConfig(),
      preconditioner: Preconditioner = Preconditioner.Identity,
      initial: Option[DVec] = None,
      toleranceMode: ToleranceMode = ToleranceMode.Absolute
  ): CgWorkspace =
    requireSquare(A)
    if b.length != A.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(A.rows), Cols(1)), Shape(Rows(b.length), Cols(1)))
    if workspace.size != A.cols then
      throw LinAlgError.VectorLengthMismatch(A.cols, workspace.size)

    initial match
      case None => workspace.x.clear()
      case Some(guess) =>
        if guess.length != A.cols then
          throw LinAlgError.VectorLengthMismatch(A.cols, guess.length)
        workspace.x := guess

    val tolerance = effectiveTolerance(toleranceMode, config.tolerance, b)
    A.applyTo(workspace.xView, workspace.ax)
    var i = 0
    while i < A.rows do
      workspace.r(i) = b(i) - workspace.ax(i)
      i += 1
    var residual = workspace.rView.norm2
    if residual <= tolerance then
      workspace.finish(0, residual, converged = true)
      return workspace

    preconditioner.solve(workspace.rView, workspace.z)
    workspace.p := workspace.zView
    var rzOld = workspace.rView.dot(workspace.zView)
    var iteration = 0
    while iteration < config.maxIterations do
      A.applyTo(workspace.pView, workspace.ap)
      val denominator = workspace.pView.dot(workspace.apView)
      if denominator == 0.0 then
        workspace.finish(iteration, residual, converged = false)
        return workspace
      val alpha = rzOld / denominator
      workspace.x.axpyInPlace(alpha, workspace.pView)
      workspace.r.axpyInPlace(-alpha, workspace.apView)
      residual = workspace.rView.norm2
      iteration += 1
      if residual <= tolerance then
        workspace.finish(iteration, residual, converged = true)
        return workspace
      preconditioner.solve(workspace.rView, workspace.z)
      val rzNew = workspace.rView.dot(workspace.zView)
      if rzOld == 0.0 then
        workspace.finish(iteration, residual, converged = false)
        return workspace
      val beta = rzNew / rzOld
      workspace.p *= beta
      workspace.p += workspace.zView
      rzOld = rzNew

    workspace.finish(iteration, residual, converged = false)
    workspace

  def bicgstab(
      A: DoubleLinearOperator,
      b: DVec,
      config: SolverConfig = SolverConfig(),
      preconditioner: Preconditioner = Preconditioner.Identity,
      initial: Option[DVec] = None,
      toleranceMode: ToleranceMode = ToleranceMode.Absolute
  ): SolverResult =
    requireSquare(A)
    if b.length != A.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(A.rows), Cols(1)), Shape(Rows(b.length), Cols(1)))
    val tol = effectiveTolerance(toleranceMode, config.tolerance, b)
    val x = initialGuess(initial, A.cols)
    val r = (b - (A * x.asVec)).mutableCopy
    val rHat = r.toVec
    val p = MutableDVec.zeros(A.cols)
    val v = MutableDVec.zeros(A.rows)
    var rhoOld = 1.0
    var alpha = 1.0
    var omega = 1.0
    var residual = r.asVec.norm2
    if residual <= tol then
      return SolverResult.Converged(x.asVec, 0, residual)

    var iteration = 0
    while iteration < config.maxIterations do
      val rhoNew = rHat.dot(r.asVec)
      if rhoNew == 0.0 then
        return SolverResult.NotConverged(x.asVec, iteration, residual)
      val beta = (rhoNew / rhoOld) * (alpha / omega)
      // p := r + beta * (p - omega * v)
      p *= beta
      p.axpyInPlace(-beta * omega, v.asVec)
      p += r.asVec
      // Left preconditioning: apply M^{-1} to the search direction before A.
      // With Identity this is a copy, recovering the unpreconditioned recurrence.
      val pHat = preconditioner(p.asVec)
      A.applyTo(pHat, v)
      val denom = rHat.dot(v.asVec)
      if denom == 0.0 then
        return SolverResult.NotConverged(x.asVec, iteration, residual)
      alpha = rhoNew / denom
      val s = r.asVec - (v.asVec * alpha)
      val sNorm = s.norm2
      if sNorm <= tol then
        x.axpyInPlace(alpha, pHat)
        return SolverResult.Converged(x.asVec, iteration + 1, sNorm)
      val sHat = preconditioner(s)
      val t = A * sHat
      val tt = t.dot(t)
      if tt == 0.0 then
        return SolverResult.NotConverged(x.asVec, iteration, residual)
      omega = t.dot(s) / tt
      x.axpyInPlace(alpha, pHat)
      x.axpyInPlace(omega, sHat)
      // r := s - omega * t keeps r the true residual b - A x.
      r := s
      r.axpyInPlace(-omega, t)
      residual = r.asVec.norm2
      iteration += 1
      if residual <= tol then
        return SolverResult.Converged(x.asVec, iteration, residual)
      if omega == 0.0 then
        return SolverResult.NotConverged(x.asVec, iteration, residual)
      rhoOld = rhoNew
    SolverResult.NotConverged(x.asVec, iteration, residual)

  def gmres(
      A: DoubleLinearOperator,
      b: DVec,
      config: SolverConfig = SolverConfig(),
      preconditioner: Preconditioner = Preconditioner.Identity,
      initial: Option[DVec] = None,
      toleranceMode: ToleranceMode = ToleranceMode.Absolute
  ): SolverResult =
    requireSquare(A)
    if b.length != A.rows then
      throw LinAlgError.DimensionMismatch(Shape(Rows(A.rows), Cols(1)), Shape(Rows(b.length), Cols(1)))
    val tol = effectiveTolerance(toleranceMode, config.tolerance, b)
    val restart = math.max(1, config.restart)
    val x = initialGuess(initial, A.cols)
    var totalIterations = 0
    // Left-preconditioned GMRES on M^{-1} A x = M^{-1} b. The Givens recurrence
    // tracks the preconditioned residual ||M^{-1}(b - A x)||, which equals the
    // true residual when M is the identity. Each inner step costs one matvec and
    // one preconditioner apply — no per-step dense QR and no extra residual matvec.
    var residual = preconditionedResidualNorm(A, preconditioner, b, x.asVec)
    if residual <= tol then
      return SolverResult.Converged(x.asVec, 0, residual)

    while totalIterations < config.maxIterations do
      val ax = MutableDVec.zeros(A.rows)
      A.applyTo(x.asVec, ax)
      val r0 = preconditioner(b - ax.asVec)
      val beta = r0.norm2
      if beta <= tol then
        return SolverResult.Converged(x.asVec, totalIterations, beta)
      val innerLimit = math.min(restart, config.maxIterations - totalIterations)
      val basis = Array.fill(innerLimit + 1)(MutableDVec.zeros(A.cols))
      scaleInto(basis(0), r0, 1.0 / beta)
      val cs = new Array[Double](innerLimit)
      val sn = new Array[Double](innerLimit)
      val g = new Array[Double](innerLimit + 1)
      g(0) = beta
      // rMat holds the rotated upper-triangular factor: rMat(i * innerLimit + col).
      val rMat = new Array[Double](innerLimit * innerLimit)
      var j = 0
      var innerSteps = 0
      var converged = false
      var brokeDown = false
      while j < innerLimit && !converged && !brokeDown do
        val av = MutableDVec.zeros(A.rows)
        A.applyTo(basis(j).asVec, av)
        val w = preconditioner(av.asVec).mutableCopy
        val hcol = new Array[Double](j + 2)
        var i = 0
        while i <= j do
          val hij = basis(i).asVec.dot(w.asVec)
          hcol(i) = hij
          w.axpyInPlace(-hij, basis(i).asVec)
          i += 1
        val hNext = w.asVec.norm2
        hcol(j + 1) = hNext
        // Happy breakdown: the next basis vector is (numerically) zero, so the
        // Krylov space is exhausted at this step.
        val happyBreakdown = hNext <= 1e-14 * beta
        if !happyBreakdown && j + 1 < basis.length then
          scaleInto(basis(j + 1), w.asVec, 1.0 / hNext)
        // Apply the accumulated Givens rotations to the new Hessenberg column.
        i = 0
        while i < j do
          val temp = cs(i) * hcol(i) + sn(i) * hcol(i + 1)
          hcol(i + 1) = -sn(i) * hcol(i) + cs(i) * hcol(i + 1)
          hcol(i) = temp
          i += 1
        val rr = math.hypot(hcol(j), hcol(j + 1))
        if rr == 0.0 then
          // Diagonal and subdiagonal both vanished: R is singular, a hard
          // breakdown that no further step can repair.
          brokeDown = true
        else
          val c = hcol(j) / rr
          val s = hcol(j + 1) / rr
          cs(j) = c
          sn(j) = s
          i = 0
          while i < j do
            rMat(i * innerLimit + j) = hcol(i)
            i += 1
          rMat(j * innerLimit + j) = rr
          val gtemp = c * g(j)
          g(j + 1) = -s * g(j)
          g(j) = gtemp
          residual = math.abs(g(j + 1))
          innerSteps = j + 1
          if residual <= tol then converged = true
          if happyBreakdown then brokeDown = true
        j += 1

      // Back-substitute R y = g over the completed steps, then x += Σ y_i v_i.
      val k = innerSteps
      if k > 0 then
        val y = new Array[Double](k)
        var singular = false
        var ii = k - 1
        while ii >= 0 do
          var sum = g(ii)
          var jj = ii + 1
          while jj < k do
            sum -= rMat(ii * innerLimit + jj) * y(jj)
            jj += 1
          val diag = rMat(ii * innerLimit + ii)
          if diag == 0.0 then
            singular = true
            y(ii) = 0.0
          else
            y(ii) = sum / diag
          ii -= 1
        if singular then
          return SolverResult.NotConverged(x.asVec, totalIterations + k, residual)
        var i = 0
        while i < k do
          x.axpyInPlace(y(i), basis(i).asVec)
          i += 1

      totalIterations += innerSteps
      if converged then
        return SolverResult.Converged(x.asVec, totalIterations, residual)
      if brokeDown then
        return SolverResult.NotConverged(x.asVec, totalIterations, residual)
    SolverResult.NotConverged(x.asVec, totalIterations, residual)

  def cgnr(
      A: DoubleLinearOperator,
      b: DVec,
      config: SolverConfig = SolverConfig(),
      toleranceMode: ToleranceMode = ToleranceMode.Absolute
  ): SolverResult =
    val rhs = MutableDVec.zeros(A.cols)
    A.transposeApplyTo(b, rhs)
    val normal =
      LinearOperator.fromFunction(A.cols, A.cols): (x, into) =>
        val tmp = MutableDVec.zeros(A.rows)
        A.applyTo(x, tmp)
        A.transposeApplyTo(tmp.asVec, into)
    // The normal-equation right-hand side is Aᵀb, so RelativeToRhs scales the
    // tolerance by ||Aᵀb|| — the natural reference for the residual CG tracks.
    cg(normal, rhs.asVec, config, toleranceMode = toleranceMode)

  /** LSQR (Paige & Saunders, 1982): Golub–Kahan bidiagonalization driving the
    * Paige–Saunders recurrences to solve `min ||A x - b||₂` for any `m x n` `A`.
    *
    * Unlike [[cgnr]] it never forms `AᵀA`, so it is numerically well-conditioned
    * for ill-conditioned `A`. The stopping test is a simplified single-tolerance
    * variant of the paper's rules: convergence on the least-squares (normal-
    * equation) residual `||Aᵀ(b − A x)||`, which drives to zero for both
    * consistent and inconsistent systems. As in [[cgnr]], `RelativeToRhs` scales
    * the tolerance by `||Aᵀ b||`; the reported `residual` is `||Aᵀ r||`.
    */
  def lsqr(
      A: DoubleLinearOperator,
      b: DVec,
      config: SolverConfig = SolverConfig(),
      toleranceMode: ToleranceMode = ToleranceMode.Absolute
  ): SolverResult =
    val m = A.rows
    val n = A.cols
    if b.length != m then
      throw LinAlgError.DimensionMismatch(Shape(Rows(m), Cols(1)), Shape(Rows(b.length), Cols(1)))

    val x = MutableDVec.zeros(n)
    // beta_1 u_1 = b.
    val u = b.mutableCopy
    var beta = u.asVec.norm2
    if beta == 0.0 then
      // b = 0 ⇒ x = 0 solves it exactly.
      return SolverResult.Converged(x.asVec, 0, 0.0)
    u *= (1.0 / beta)

    // alpha_1 v_1 = Aᵀ u_1.
    val v = MutableDVec.zeros(n)
    A.transposeApplyTo(u.asVec, v)
    var alpha = v.asVec.norm2
    // ||Aᵀ b|| = alpha_1 * beta_1: the natural reference for the LS residual.
    val arnorm0 = alpha * beta
    if alpha == 0.0 then
      // Aᵀ b = 0 ⇒ x = 0 is already the least-squares solution.
      return SolverResult.Converged(x.asVec, 0, 0.0)
    v *= (1.0 / alpha)

    val w = v.asVec.mutableCopy
    var phibar = beta
    var rhobar = alpha
    val tol = toleranceMode match
      case ToleranceMode.Absolute      => config.tolerance
      case ToleranceMode.RelativeToRhs => config.tolerance * arnorm0

    val av = MutableDVec.zeros(m)
    val atu = MutableDVec.zeros(n)
    var arnorm = arnorm0
    var iteration = 0
    var converged = false
    while iteration < config.maxIterations && !converged do
      // Bidiagonalization: u := A v - alpha u, then normalise.
      A.applyTo(v.asVec, av)
      u *= (-alpha)
      u += av.asVec
      beta = u.asVec.norm2
      if beta != 0.0 then u *= (1.0 / beta)
      // v := Aᵀ u - beta v, then normalise.
      A.transposeApplyTo(u.asVec, atu)
      v *= (-beta)
      v += atu.asVec
      alpha = v.asVec.norm2
      if alpha != 0.0 then v *= (1.0 / alpha)

      // Givens rotation folding the new subdiagonal into the bidiagonal factor.
      val rho = math.hypot(rhobar, beta)
      val cs = rhobar / rho
      val sn = beta / rho
      val theta = sn * alpha
      rhobar = -cs * alpha
      val phi = cs * phibar
      phibar = sn * phibar

      // x := x + (phi / rho) w; w := v - (theta / rho) w.
      x.axpyInPlace(phi / rho, w.asVec)
      w *= (-(theta / rho))
      w += v.asVec

      iteration += 1
      // ||A^T r_k|| = phibar_{k+1} * alpha_{k+1} * |c_k| (Paige–Saunders, S2).
      arnorm = phibar * alpha * math.abs(cs)
      if arnorm <= tol then converged = true

    if converged then SolverResult.Converged(x.asVec, iteration, arnorm)
    else SolverResult.NotConverged(x.asVec, iteration, arnorm)

  private def requireSquare(A: DoubleLinearOperator): Unit =
    if A.rows != A.cols then
      throw LinAlgError.DimensionMismatch(Shape(Rows(A.rows), Cols(A.rows)), Shape(Rows(A.rows), Cols(A.cols)))

  /** Effective convergence threshold for the chosen [[ToleranceMode]]. */
  private def effectiveTolerance(mode: ToleranceMode, tolerance: Double, b: DVec): Double =
    mode match
      case ToleranceMode.Absolute      => tolerance
      case ToleranceMode.RelativeToRhs => tolerance * b.norm2

  /** `||M^{-1}(b - A x)||`, the preconditioned residual norm GMRES tracks. */
  private def preconditionedResidualNorm(
      A: DoubleLinearOperator,
      preconditioner: Preconditioner,
      b: DVec,
      x: DVec
  ): Double =
    val ax = MutableDVec.zeros(A.rows)
    A.applyTo(x, ax)
    preconditioner(b - ax.asVec).norm2

  /** `out := alpha * x`. */
  private def scaleInto(out: MutableDVec, x: DVec, alpha: Double): Unit =
    out := x
    out *= alpha

  /** Materialise the starting iterate: a fresh zero vector when no guess is given,
    * or a mutable copy of the supplied guess. A guess of any length other than
    * `cols` — including a length-0 vector for a nonzero-dimension system — is a
    * dimension mismatch, not a silent fallback to zero.
    */
  private def initialGuess(initial: Option[DVec], cols: Int): MutableDVec =
    initial match
      case None => MutableDVec.zeros(cols)
      case Some(guess) =>
        if guess.length == cols then guess.mutableCopy
        else throw LinAlgError.VectorLengthMismatch(cols, guess.length)

def cg(
    A: DoubleLinearOperator,
    b: DVec,
    config: SolverConfig = SolverConfig(),
    preconditioner: Preconditioner = Preconditioner.Identity,
    initial: Option[DVec] = None,
    toleranceMode: ToleranceMode = ToleranceMode.Absolute
): SolverResult =
  IterativeSolvers.cg(A, b, config, preconditioner, initial, toleranceMode)

def cgWith(
    A: DoubleLinearOperator,
    b: DVec,
    workspace: CgWorkspace,
    config: SolverConfig = SolverConfig(),
    preconditioner: Preconditioner = Preconditioner.Identity,
    initial: Option[DVec] = None,
    toleranceMode: ToleranceMode = ToleranceMode.Absolute
): CgWorkspace =
  IterativeSolvers.cgWith(A, b, workspace, config, preconditioner, initial, toleranceMode)

def bicgstab(
    A: DoubleLinearOperator,
    b: DVec,
    config: SolverConfig = SolverConfig(),
    preconditioner: Preconditioner = Preconditioner.Identity,
    initial: Option[DVec] = None,
    toleranceMode: ToleranceMode = ToleranceMode.Absolute
): SolverResult =
  IterativeSolvers.bicgstab(A, b, config, preconditioner, initial, toleranceMode)

def gmres(
    A: DoubleLinearOperator,
    b: DVec,
    config: SolverConfig = SolverConfig(),
    preconditioner: Preconditioner = Preconditioner.Identity,
    initial: Option[DVec] = None,
    toleranceMode: ToleranceMode = ToleranceMode.Absolute
): SolverResult =
  IterativeSolvers.gmres(A, b, config, preconditioner, initial, toleranceMode)

/** Conjugate gradient on the normal equations AᵀA x = Aᵀb (CGNR).
  *
  * Solves the least-squares / minimum-residual problem by running CG on the
  * symmetric positive (semi-)definite normal-equation operator. Forming AᵀA
  * squares the condition number of `A`, so this suits well-conditioned problems;
  * the `tolerance` applies to the normal-equation residual ‖Aᵀ(b − Ax)‖.
  *
  * Formerly named `lsqr`: it never implemented the Paige–Saunders LSQR
  * recurrence, so the honest name is CGNR.
  */
def cgnr(
    A: DoubleLinearOperator,
    b: DVec,
    config: SolverConfig = SolverConfig(),
    toleranceMode: ToleranceMode = ToleranceMode.Absolute
): SolverResult =
  IterativeSolvers.cgnr(A, b, config, toleranceMode)

/** LSQR (Paige–Saunders): Golub–Kahan bidiagonalization for `min ||A x − b||₂`.
  *
  * Numerically well-conditioned where [[cgnr]] struggles, because it never forms
  * `AᵀA`. Stopping is a simplified single-tolerance variant on the least-squares
  * residual `||Aᵀ(b − A x)||`; see [[IterativeSolvers.lsqr]].
  */
def lsqr(
    A: DoubleLinearOperator,
    b: DVec,
    config: SolverConfig = SolverConfig(),
    toleranceMode: ToleranceMode = ToleranceMode.Absolute
): SolverResult =
  IterativeSolvers.lsqr(A, b, config, toleranceMode)
