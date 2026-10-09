package gale.golden

/** The error units of the Breeze golden corpus, shared by the replay (`BreezeGoldenSuite`, which allows `C` units
  * between gale and Breeze) and the parity-side freshness check (which allows a larger multiple between a committed
  * and a regenerated Breeze reference).
  *
  * A forward unit has the shape of the standard forward-error bound `n · κ · ε · scale`. It is computed from the
  * case alone (inputs and Breeze references), per reference array: `n` from the operation's dimensions and `scale`
  * from that array. `κ` is Breeze's 2-norm condition number for conditioning-sensitive results and `1` for
  * backward-stable products and reductions. Exact or correctly rounded operations have unit `0`.
  */
object GoldenTolerance:

  val Eps: Double = Math.ulp(1.0)

  def maxAbs(a: Array[Double]): Double = a.foldLeft(0.0)((s, x) => math.max(s, math.abs(x)))
  def sumAbs(a: Array[Double]): Double = a.foldLeft(0.0)((s, x) => s + math.abs(x))
  def norm2(a: Array[Double]): Double = math.sqrt(a.foldLeft(0.0)((s, x) => s + x * x))

  /** `max_ij Σ_k |a_ik| |b_kj|` (or with `Aᵀ`): the scale of a rounding-error bound for `A B`. */
  def absProductMax(a: GoldenArray, b: GoldenArray, transposeA: Boolean): Double =
    val (m, k) = if transposeA then (a.cols, a.rows) else (a.rows, a.cols)
    var out = 0.0
    for i <- 0 until m; j <- 0 until b.cols do
      var s = 0.0
      for p <- 0 until k do s += math.abs(if transposeA then a(p, i) else a(i, p)) * math.abs(b(p, j))
      out = math.max(out, s)
    out

  /** Least-squares forward-error growth `κ + κ² ‖r‖ / (‖A‖ ‖x‖)`. `‖A‖_F ≥ ‖A‖₂` makes the `κ²` term no larger
    * than with the 2-norm, so this growth (and the bound) is at most the textbook one.
    */
  def lstsqGrowth(c: GoldenCase): Double =
    val k = c.scalar("kappa")
    k + k * k * c.scalar("residual") / (norm2(c.values("A")) * norm2(c.values("x")))

  /** The forward unit of reference `key` of case `c`, per element. */
  def forwardUnit(c: GoldenCase, key: String): Int => Double =
    val ref = c.values(key)
    lazy val kappa = c.scalar("kappa")
    def const(u: Double): Int => Double = _ => u
    val exact: Int => Double = _ => 0.0
    (c.family, key) match
      case ("matvec", "Ax") =>
        const(c.array("A").cols * Eps * absProductMax(c.array("A"), c.array("x"), transposeA = false))
      case ("matvec", "Atxt") =>
        const(c.array("A").rows * Eps * absProductMax(c.array("A"), c.array("xt"), transposeA = true))
      case ("matmul", "AB") =>
        const(c.array("A").cols * Eps * absProductMax(c.array("A"), c.array("B"), transposeA = false))
      case ("matmul", "AtD") =>
        const(c.array("A").rows * Eps * absProductMax(c.array("A"), c.array("D"), transposeA = true))
      case ("vector", "dot") =>
        val (x, y) = (c.values("x"), c.values("y"))
        const(x.length * Eps * x.indices.map(i => math.abs(x(i) * y(i))).sum)
      case ("vector", "axpy") =>
        val (x, y, alpha) = (c.values("x"), c.values("y"), c.scalar("alpha"))
        i => Eps * (math.abs(alpha * x(i)) + math.abs(y(i)))
      case (_, "kappa") =>
        // σ_min carries absolute error ~ n ε σ_max, i.e. relative error ~ n κ ε.
        const(c.array("A").rows.max(c.array("A").cols) * kappa * Eps * kappa)
      case ("lu" | "cholesky", "x" | "inv" | "L" | "det") =>
        const(c.array("A").rows * kappa * Eps * maxAbs(ref))
      case ("qr", "R") =>
        const(c.array("A").cols * kappa * Eps * maxAbs(ref))
      case ("qr", "x") =>
        const(c.array("A").cols * lstsqGrowth(c) * Eps * maxAbs(ref))
      case ("qr", "residual") =>
        val a = c.array("A")
        const(a.rows * Eps * (norm2(a.values) * norm2(c.values("x")) + norm2(c.values("b"))))
      case ("eigsym", "values") =>
        const(c.array("A").rows * Eps * maxAbs(ref))
      case ("svd", "sigma") =>
        const(c.array("A").rows.max(c.array("A").cols) * Eps * maxAbs(ref))
      case ("pinv", "P") =>
        const(c.array("A").rows.max(c.array("A").cols) * kappa * Eps * maxAbs(ref))
      case ("kron", "K") => exact
      case ("reduce-vector" | "reduce-matrix" | "sparse", "max" | "min" | "argmax" | "argmin") => exact
      case ("reduce-matrix", "argmaxRow" | "argmaxCol")                                          => exact
      case ("reduce-vector" | "sparse", "normInf")                                               => exact
      case ("sparse", "plus" | "minus")                                                          => exact
      case ("reduce-vector" | "reduce-matrix" | "sparse", "sum" | "norm1" | "norm2" | "frobenius" | "mean" | "normInf") =>
        val input = c.values(if c.family == "reduce-matrix" then "A" else if c.family == "sparse" then "xv" else "x")
        val n =
          if c.family == "reduce-matrix" then
            val a = c.array("A")
            key match
              case "norm1"   => a.rows
              case "normInf" => a.cols
              case _         => a.rows * a.cols
          else if c.family == "sparse" then c.values("d").length
          else input.length
        key match
          case "sum" | "norm1" => const(n * Eps * sumAbs(input))
          case "mean"          => const(Eps * sumAbs(input) + Eps * math.abs(ref(0)))
          case _               => const(n * Eps * math.abs(ref(0)))
      case ("sparse", "dot") =>
        val n = c.values("d").length
        val (x, y) = (dense(c, "xi", "xv"), dense(c, "yi", "yv"))
        const(n * Eps * x.indices.map(i => math.abs(x(i) * y(i))).sum)
      case ("sparse", "dotDense") =>
        val d = c.values("d")
        val x = dense(c, "xi", "xv")
        const(d.length * Eps * x.indices.map(i => math.abs(x(i) * d(i))).sum)
      case ("numerics", "exp" | "log" | "sigmoid") => i => Eps * math.abs(ref(i))
      case ("numerics", "logSumExp") =>
        val x = c.values("x")
        const(Eps * (math.abs(ref(0)) + maxAbs(x) + x.length))
      case ("numerics", "softmax") =>
        val x = c.values("x")
        val lse = c.scalar("logSumExp")
        i => Eps * math.abs(ref(i)) * (1 + math.abs(x(i)) + math.abs(lse) + x.length)
      case _ => throw new NoSuchElementException(s"no tolerance unit for ${c.family}/$key")

  private def dense(c: GoldenCase, indices: String, values: String): Array[Double] =
    val out = new Array[Double](c.values("d").length)
    val idx = c.values(indices)
    val v = c.values(values)
    for k <- idx.indices do out(idx(k).toInt) = v(k)
    out

  /** Partitions ascending eigenvalues into clusters a projector must treat together: consecutive values whose gap
    * is unresolvable at working precision (`≤ 64 n ε ‖A‖`) or tiny relative to their magnitude (`≤ 1e-8 |λ|`).
    * Graded spectra keep their small, well-separated (relative to themselves) eigenvalues apart.
    */
  def clustersOf(values: Array[Double]): Vector[Range] =
    val n = values.length
    val norm = maxAbs(values)
    val out = Vector.newBuilder[Range]
    var start = 0
    for i <- 1 until n do
      val gap = values(i) - values(i - 1)
      val scale = math.max(math.abs(values(i)), math.abs(values(i - 1)))
      if gap > 64 * n * Eps * norm && gap > 1e-8 * scale then
        out += (start until i)
        start = i
    out += (start until n)
    out.result()

  /** Worst `‖P_u − P_v‖_max / (n ε ‖A‖ / gap)` over eigenvalue clusters, for eigenvector matrices `u`, `v`
    * (columns aligned with the ascending `values`).
    */
  def projectorRatio(u: (Int, Int) => Double, v: (Int, Int) => Double, values: Array[Double]): Double =
    val n = values.length
    val norm = math.max(maxAbs(values), Double.MinPositiveValue)
    val clusters = clustersOf(values)
    clusters.map { cl =>
      val gap = clusters
        .filter(_ != cl)
        .map(q => math.min(math.abs(values(q.head) - values(cl.last)), math.abs(values(q.last) - values(cl.head))))
        .minOption
        .getOrElse(norm)
      var diff = 0.0
      for i <- 0 until n; j <- 0 until n do
        var pu = 0.0
        var pv = 0.0
        for k <- cl do
          pu += u(i, k) * u(j, k)
          pv += v(i, k) * v(j, k)
        diff = math.max(diff, math.abs(pu - pv))
      if diff == 0.0 then 0.0 else diff / (n * Eps * norm / gap)
    }.max
