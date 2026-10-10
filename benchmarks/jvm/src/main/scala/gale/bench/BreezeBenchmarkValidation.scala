package gale.bench

import breeze.linalg.{DenseMatrix as BDM, LU as BLU, cholesky as bChol, eig as bEig, eigSym as bEigSym, inv as bInv, pinv as bPinv, qr as bQr, svd as bSvd}
import gale.backend.{Backend, PureBackend}
import gale.bench.BreezeBenchData.*
import gale.linalg.*
import gale.spectral.*
import java.nio.charset.StandardCharsets
import java.nio.file.{Files, Paths, StandardOpenOption}
import java.security.MessageDigest
import org.openjdk.jmh.infra.BenchmarkParams
import scala.collection.mutable.ArrayBuffer
import scala.util.control.NonFatal

/** Untimed checks of the exact paired decomposition fixtures. The reference
  * arithmetic is scalar array code, independent of either library's BLAS.
  * This is a validator for finite, moderately scaled benchmark inputs, not a
  * general-purpose numerical testing API. Normalized backward/orthogonality
  * bounds use 128 * max(m,n) * eps, allowing accumulation and norm conversion.
  */
object BreezeBenchmarkValidation:
  type Data = Array[Array[Double]]
  final case class Metric(name: String, value: Double, limit: Double)
  final val SidecarProperty = "gale.bench.validationSidecar"
  private val Eps = math.ulp(1.0)

  def operation(params: BenchmarkParams): String =
    params.getBenchmark.split('.').last.stripPrefix("gale").stripPrefix("breeze")

  def column(x: Array[Double]): Data = x.map(v => Array(v))

  private def data(a: DMat): Data = Array.tabulate(a.rows, a.cols)(a(_, _))
  private def data(a: BDM[Double]): Data = Array.tabulate(a.rows, a.cols)(a(_, _))
  private def identity(n: Int): Data = Array.tabulate(n, n)((i, j) => if i == j then 1.0 else 0.0)
  private def transpose(a: Data): Data = Array.tabulate(a(0).length, a.length)((i, j) => a(j)(i))
  private def norm(a: Data): Double = a.iterator.flatMap(_.iterator).foldLeft(0.0)(math.hypot)
  private def minus(a: Data, b: Data): Data = Array.tabulate(a.length, a(0).length)((i, j) => a(i)(j) - b(i)(j))

  private def product(a: Data, b: Data): Data =
    require(a(0).length == b.length, "preflight product shape")
    val c = Array.ofDim[Double](a.length, b(0).length)
    var i = 0
    while i < a.length do
      var k = 0
      while k < b.length do
        val aik = a(i)(k)
        var j = 0
        while j < b(0).length do
          c(i)(j) += aik * b(k)(j)
          j += 1
        k += 1
      i += 1
    c

  private def ratio(error: Double, scale: Double): Double =
    require(error.isFinite && scale.isFinite && scale >= 0.0, s"non-finite/invalid residual or scale: $error / $scale")
    if scale == 0.0 && error == 0.0 then 0.0 else error / scale

  def successful[A](label: String, result: Either[LinAlgError, A]): A =
    result.fold(e => throw new IllegalStateException(s"$label failed: ${e.getMessage}"), value => value)

  private[bench] final class Audit(n: Int):
    val metrics = ArrayBuffer.empty[Metric]
    val tolerance: Double = 128.0 * math.max(1, n) * Eps

    def check(name: String, value: Double, limit: Double = tolerance): Unit =
      require(value.isFinite && value >= 0.0 && value <= limit, s"$name: error=$value limit=$limit")
      metrics += Metric(name, value, limit)

    def matrix(name: String, a: Data, rows: Int, cols: Int): Unit =
      require(a.length == rows && a.forall(_.length == cols), s"$name: expected ${rows}x$cols")
      require(a.forall(_.forall(_.isFinite)), s"$name: non-finite output")

    def orthogonal(name: String, q: Data): Unit =
      val k = q(0).length
      check(name, norm(minus(product(transpose(q), q), identity(k))))

    def solve(name: String, a: Data, x: Data, b: Data): Unit =
      matrix(name, x, a(0).length, b(0).length)
      check(name, ratio(norm(minus(product(a, x), b)), norm(a) * norm(x) + norm(b)))

    def inverse(name: String, a: Data, x: Data): Unit =
      solve(s"$name right", a, x, identity(a.length))
      solve(s"$name left", x, a, identity(a.length))

    def qr(name: String, a: Data, r: Data): Unit =
      val n = a(0).length
      matrix(name, r, math.min(a.length, n), n)
      require(r.indices.forall(i => (0 until math.min(i, n)).forall(j => r(i)(j) == 0.0)), s"$name R triangular")
      val ata = product(transpose(a), a)
      check(s"$name QR Gram", ratio(norm(minus(product(transpose(r), r), ata)), norm(ata)))

    def cholesky(name: String, a: Data, l: Data): Unit =
      val n = a.length
      matrix(name, l, n, n)
      require(l.indices.forall(i => (i + 1 until n).forall(j => l(i)(j) == 0.0)), s"$name Cholesky lower triangular")
      require(l.indices.forall(i => l(i)(i) > 0.0), s"$name Cholesky diagonal")
      check(s"$name Cholesky reconstruction", ratio(norm(minus(product(l, transpose(l)), a)), norm(a)))

    def determinant(name: String, a: Data, value: Double): Unit =
      val (sign, logAbs) = referenceLogDet(a)
      require(!value.isNaN && math.signum(value) == sign, s"$name determinant sign")
      if logAbs > math.log(Double.MaxValue) + 1.0 then
        require(value.isInfinity, s"$name determinant should overflow (reference logabs=$logAbs)")
        check(s"$name expected determinant overflow", 0.0)
      else
        require(value.isFinite && value != 0.0, s"$name determinant non-finite/zero")
        check(s"$name log determinant", math.abs(math.log(math.abs(value)) - logAbs), tolerance * math.max(1.0, math.abs(logAbs)))

    def pinv(name: String, a: Data, x: Data): Unit =
      matrix(name, x, a(0).length, a.length)
      val ax = product(a, x)
      val xa = product(x, a)
      check(s"$name AXA", ratio(norm(minus(product(ax, a), a)), norm(a)))
      check(s"$name XAX", ratio(norm(minus(product(xa, x), x)), norm(x)))
      check(s"$name AX symmetric", ratio(norm(minus(ax, transpose(ax))), norm(ax)))
      check(s"$name XA symmetric", ratio(norm(minus(xa, transpose(xa))), norm(xa)))

    def svd(name: String, a: Data, sigma: Array[Double], u: Data, vt: Data): Unit =
      val m = a.length
      val n = a(0).length
      val p = math.min(m, n)
      require(sigma.length == p && sigma.forall(v => v.isFinite && v >= 0.0), s"$name spectrum")
      require(sigma.sliding(2).forall(s => s.length < 2 || s(0) >= s(1)), s"$name spectrum ordering")
      matrix(s"$name U", u, m, p)
      matrix(s"$name VT", vt, p, n)
      val us = Array.tabulate(m, p)((i, j) => u(i)(j) * sigma(j))
      check(s"$name reconstruction", ratio(norm(minus(product(us, vt), a)), norm(a)))
      orthogonal(s"$name U orthogonality", u)
      orthogonal(s"$name V orthogonality", transpose(vt))

    def eigenpair(name: String, a: Data, re: Double, im: Double, vr: Array[Double], vi: Array[Double]): Unit =
      require(re.isFinite && im.isFinite, s"$name eigenvalue")
      require(vr.length == a.length && vi.length == a.length, s"$name eigenvector shape")
      require(vr.forall(_.isFinite) && vi.forall(_.isFinite), s"$name eigenvector non-finite")
      var vn = 0.0
      var rn = 0.0
      var i = 0
      while i < a.length do
        var ar = 0.0
        var ai = 0.0
        var j = 0
        while j < a.length do
          ar += a(i)(j) * vr(j)
          ai += a(i)(j) * vi(j)
          j += 1
        rn = math.hypot(rn, ar - re * vr(i) + im * vi(i))
        rn = math.hypot(rn, ai - re * vi(i) - im * vr(i))
        vn = math.hypot(vn, math.hypot(vr(i), vi(i)))
        i += 1
      require(vn > 0.0 && vn.isFinite, s"$name zero/invalid eigenvector")
      check(name, ratio(rn, (norm(a) + math.hypot(re, im)) * vn))

    /** Per-pair residuals alone allow n copies of one eigenpair. Reconstruct A
      * through the complete real modal basis, with 2x2 conjugate blocks. These
      * benchmark fixtures must have an adequately conditioned eigenbasis; a
      * singular basis or inaccurate reconstruction fails their qualification.
      */
    def eigenbasis(name: String, a: Data, re: Array[Double], im: Array[Double], v: Data): Unit =
      val n = a.length
      matrix(name, v, n, n)
      require(re.length == n && im.length == n && re.forall(_.isFinite) && im.forall(_.isFinite), s"$name spectrum")
      for j <- 0 until n do
        if im(j) > 0.0 then require(j + 1 < n && im(j + 1) == -im(j) && re(j + 1) == re(j), s"$name conjugate packing")
        if im(j) < 0.0 then require(j > 0 && im(j - 1) == -im(j) && re(j - 1) == re(j), s"$name conjugate packing")
      val vd = Array.tabulate(n, n) { (i, j) =>
        if im(j) > 0.0 then re(j) * v(i)(j) - im(j) * v(i)(j + 1)
        else if im(j) < 0.0 then re(j) * v(i)(j) - im(j) * v(i)(j - 1)
        else re(j) * v(i)(j)
      }
      val reconstructed = transpose(referenceSolve(transpose(v), transpose(vd)))
      matrix(s"$name reconstructed", reconstructed, n, n)
      check(s"$name complete reconstruction", ratio(norm(minus(reconstructed, a)), norm(a)))

  /** Both libraries are run, even in a single-side JMH fork. Only the selected
    * operation is validated, so e.g. a QR trial does not run nonsymmetric eig.
    */
  def validate(op: String, a: Data, rhs: Data = Array.empty, backend: Backend = PureBackend): Vector[Metric] =
    require(a.nonEmpty && a(0).nonEmpty && a.forall(_.length == a(0).length), "nonempty rectangular fixture required")
    require(a.forall(_.forall(_.isFinite)) && rhs.forall(_.forall(_.isFinite)), "finite fixture required")
    given Backend = backend
    val m = a.length
    val n = a(0).length
    val audit = new Audit(math.max(m, n))
    val ga = galeMatrix(a)
    val ba = breezeMatrix(a)
    op match
      case "Inv" =>
        audit.inverse("gale inverse", a, data(successful("gale inverse", ga.inverse)))
        audit.inverse("breeze inverse", a, data(bInv(ba)))
      case "Solve" | "LuSolve" | "CholSolve" | "Lstsq" =>
        require(rhs.length == m && rhs.nonEmpty && rhs.forall(_.length == rhs(0).length), "RHS shape")
        val gx = op match
          case "CholSolve" => data(successful("gale Cholesky solve", ga.cholesky.flatMap(_.solve(galeMatrix(rhs)))))
          case "Lstsq" => column(successful("gale least squares", ga.leastSquares(galeVector(rhs.map(_(0))))).toArray)
          case "Solve" => column(successful("gale solve", ga.solve(galeVector(rhs.map(_(0))))).toArray)
          case _ => data(successful("gale multi-RHS solve", ga.solve(galeMatrix(rhs))))
        val bx =
          if op == "CholSolve" then breezeCholeskySolve(ba, breezeMatrix(rhs))
          else if op == "Solve" || op == "Lstsq" then column((ba \ breezeVector(rhs.map(_(0)))).toArray)
          else data(ba \ breezeMatrix(rhs))
        for (label, x) <- Seq("gale" -> gx, "breeze" -> bx) do
          if op == "Lstsq" then
            audit.matrix(label, x, n, 1)
            val r = minus(product(a, x), rhs)
            audit.check(s"$label least-squares stationarity", ratio(norm(product(transpose(a), r)), norm(a) * (norm(a) * norm(x) + norm(rhs))))
          else audit.solve(s"$label $op", a, x, rhs)
      case "Lu" =>
        val g = successful("gale LU", ga.lu)
        require(g.diagnostics.isSuccess, "gale LU diagnostics")
        val gp = data(g.packed)
        audit.matrix("gale LU", gp, n, n)
        val (bp, swaps) = BLU.primitive(ba)
        val p = (0 until n).toArray
        for i <- swaps.indices do
          val j = swaps(i) - 1 // LAPACK's one-based sequential interchanges.
          require(j >= i && j < n, s"breeze LU pivot $i -> $j")
          val t = p(i); p(i) = p(j); p(j) = t
        for (label, packed, permutation) <- Seq(("gale", gp, g.pivots.toArray), ("breeze", data(bp), p)) do
          audit.matrix(label, packed, n, n)
          require(permutation.sorted.sameElements((0 until n).toArray), s"$label LU permutation")
          val l = Array.tabulate(n, n)((i, j) => if i == j then 1.0 else if i > j then packed(i)(j) else 0.0)
          val u = Array.tabulate(n, n)((i, j) => if i <= j then packed(i)(j) else 0.0)
          val pa = Array.tabulate(n, n)((i, j) => a(permutation(i))(j))
          audit.check(s"$label LU reconstruction", ratio(norm(minus(product(l, u), pa)), norm(a)))
      case "Chol" =>
        val g = successful("gale Cholesky", ga.cholesky)
        require(g.diagnostics.isSuccess, "gale Cholesky diagnostics")
        for (label, l) <- Seq("gale" -> data(g.lower), "breeze" -> data(bChol(ba))) do
          audit.cholesky(label, a, l)
      case "Qr" =>
        val g = ga.qr
        require(g.diagnostics.isSuccess, "gale QR diagnostics")
        val gr = data(g.r)
        val br = data(bQr.justR(ba))
        for (label, r) <- Seq("gale" -> gr, "breeze" -> br) do
          audit.qr(label, a, r)
        val q = data(g.q)
        audit.matrix("gale Q", q, m, m)
        audit.orthogonal("gale Q orthogonality", q)
        val paddedR = Array.tabulate(m, n)((i, j) => if i < gr.length then gr(i)(j) else 0.0)
        audit.check("gale QR reconstruction", ratio(norm(minus(product(q, paddedR), a)), norm(a)))
      case "Svd" =>
        val g = successful("gale SVD", ga.svd)
        require(g.diagnostics.allConverged, "gale SVD not converged")
        val b = bSvd(ba)
        val p = math.min(m, n)
        audit.svd("gale", a, g.singularValues.toArray, data(g.u), data(g.vt))
        audit.svd("breeze", a, b.singularValues.toArray,
          Array.tabulate(m, p)(b.leftVectors(_, _)), Array.tabulate(p, n)(b.rightVectors(_, _)))
      case "Pinv" =>
        audit.pinv("gale pinv", a, data(successful("gale pinv", ga.pinv)))
        audit.pinv("breeze pinv", a, data(bPinv(ba)))
      case "EigSym" =>
        val g = successful("gale eigSym", Eigen.eigSymmetric(ga, EigenSelection.All, EigenVectors.Right))
        require(g.size == n && g.diagnostics.allConverged, "gale eigSym count/convergence")
        val b = bEigSym(ba)
        for (label, values, v) <- Seq(("gale", g.eigenvalues.toArray, data(g.eigenvectors)), ("breeze", b.eigenvalues.toArray, data(b.eigenvectors))) do
          require(values.length == n && values.forall(_.isFinite), s"$label eigSym spectrum")
          audit.matrix(label, v, n, n)
          audit.orthogonal(s"$label eigSym orthogonality", v)
          for j <- 0 until n do audit.eigenpair(s"$label eigSym pair $j", a, values(j), 0.0, v.map(_(j)), new Array[Double](n))
      case "Eig" =>
        val g = successful("gale eig", Eigen.eigNonsymmetric(ga, EigenSelection.All, EigenVectors.Right))
        require(g.size == n && g.diagnostics.allConverged, "gale eig count/convergence")
        val gv = Array.ofDim[Double](n, n)
        for j <- 0 until n do
          val lambda = g.eigenvalue(j)
          val (vr, vi) = g.eigenvector(j)
          audit.eigenpair(s"gale eig pair $j", a, lambda.re, lambda.im, vr.toArray, vi.toArray)
          if lambda.im >= 0.0 then
            for i <- 0 until n do
              gv(i)(j) = vr(i)
              if lambda.im > 0.0 then gv(i)(j + 1) = vi(i)
        audit.eigenbasis("gale eig", a, Array.tabulate(n)(g.eigenvalue(_).re), Array.tabulate(n)(g.eigenvalue(_).im), gv)
        val b = bEig(ba)
        require(b.eigenvalues.length == n && b.eigenvaluesComplex.length == n, "breeze eig count")
        audit.matrix("breeze eigenvectors", data(b.eigenvectors), n, n)
        for j <- 0 until n do
          val im = b.eigenvaluesComplex(j)
          val rc = if im < 0.0 then j - 1 else j
          val ic = if im > 0.0 then j + 1 else j
          require(rc >= 0 && ic < n, "breeze eig packing")
          val vr = Array.tabulate(n)(b.eigenvectors(_, rc))
          val vi = Array.tabulate(n)(i => if im == 0.0 then 0.0 else b.eigenvectors(i, ic) * math.signum(im))
          audit.eigenpair(s"breeze eig pair $j", a, b.eigenvalues(j), im, vr, vi)
        audit.eigenbasis("breeze eig", a, b.eigenvalues.toArray, b.eigenvaluesComplex.toArray, data(b.eigenvectors))
      case "Det" =>
        for (label, value) <- Seq("gale" -> successful("gale determinant", ga.det), "breeze" -> breeze.linalg.det(ba)) do
          audit.determinant(label, a, value)
      case "Gemm" =>
        val expected = product(a, rhs)
        for (label, actual) <- Seq("gale" -> data(ga * galeMatrix(rhs)), "breeze" -> data(ba * breezeMatrix(rhs))) do
          audit.matrix(label, actual, m, rhs(0).length)
          audit.check(s"$label product", ratio(norm(minus(actual, expected)), norm(a) * norm(rhs)))
      case other => throw new IllegalArgumentException(s"unvalidated operation $other")
    audit.metrics.toVector

  /** Same netlib calls and copies as MultiRhsBreezeJmh.breezeCholSolve. */
  private def breezeCholeskySolve(a: BDM[Double], b: BDM[Double]): Data =
    val l = a.data.clone()
    val x = b.data.clone()
    val info = new org.netlib.util.intW(0)
    val lapack = dev.ludovic.netlib.lapack.LAPACK.getInstance()
    lapack.dpotrf("L", a.rows, l, a.rows, info)
    require(info.`val` == 0, s"breeze dpotrf info=${info.`val`}")
    lapack.dpotrs("L", a.rows, b.cols, l, a.rows, x, a.rows, info)
    require(info.`val` == 0, s"breeze dpotrs info=${info.`val`}")
    Array.tabulate(b.rows, b.cols)((i, j) => x(i + j * b.rows))

  /** Scalar partial-pivot elimination, independent of both libraries' solvers. */
  private def referenceSolve(a: Data, b: Data): Data =
    val h = a.map(_.clone())
    val x = b.map(_.clone())
    val n = a.length
    for k <- 0 until n do
      var pivot = k
      for i <- k + 1 until n do
        if math.abs(h(i)(k)) > math.abs(h(pivot)(k)) then pivot = i
      if pivot != k then
        val row = h(k); h(k) = h(pivot); h(pivot) = row
        val rhs = x(k); x(k) = x(pivot); x(pivot) = rhs
      val d = h(k)(k)
      require(d.isFinite && d != 0.0, "singular/non-finite modal basis")
      for i <- k + 1 until n do
        val factor = h(i)(k) / d
        for j <- k + 1 until n do h(i)(j) -= factor * h(k)(j)
        for j <- x(i).indices do x(i)(j) -= factor * x(k)(j)
    for i <- (0 until n).reverse; j <- x(i).indices do
      var sum = x(i)(j)
      for k <- i + 1 until n do sum -= h(i)(k) * x(k)(j)
      x(i)(j) = sum / h(i)(i)
    x

  /** Independent textbook elimination, including partial pivoting. Used instead
    * of comparing infinity with infinity for the n=256 determinant fixture.
    */
  private def referenceLogDet(a: Data): (Double, Double) =
    val x = a.map(_.clone())
    var sign = 1.0
    var logAbs = 0.0
    var k = 0
    while k < x.length do
      var pivot = k
      for i <- k + 1 until x.length do
        if math.abs(x(i)(k)) > math.abs(x(pivot)(k)) then pivot = i
      if pivot != k then
        val row = x(k); x(k) = x(pivot); x(pivot) = row
        sign = -sign
      val d = x(k)(k)
      require(d.isFinite && d != 0.0, "reference determinant requires nonsingular fixture")
      sign *= math.signum(d)
      logAbs += math.log(math.abs(d))
      for i <- k + 1 until x.length do
        val factor = x(i)(k) / d
        for j <- k + 1 until x.length do x(i)(j) -= factor * x(k)(j)
      k += 1
    (sign, logAbs)

  private def quote(s: String): String =
    "\"" + s.flatMap {
      case '\\' => "\\\\"
      case '"' => "\\\""
      case c if c < ' ' => f"\\u${c.toInt}%04x"
      case c => c.toString
    } + "\""

  private def fingerprint(a: Data, rhs: Data): String =
    val digest = MessageDigest.getInstance("SHA-256")
    for matrix <- Seq(a, rhs) do
      digest.update(java.nio.ByteBuffer.allocate(4).putInt(matrix.length).array())
      for row <- matrix do
        digest.update(java.nio.ByteBuffer.allocate(4).putInt(row.length).array())
        for v <- row do digest.update(java.nio.ByteBuffer.allocate(8).putLong(java.lang.Double.doubleToRawLongBits(v)).array())
    digest.digest().map(b => f"${b & 0xff}%02x").mkString

  /** Invoked only from Level.Trial setup; failures prevent all timed iterations. */
  def preflight(params: BenchmarkParams, a: Data, fixture: String, rhs: Data = Array.empty): Unit =
    val keys = params.getParamsKeys.toArray.map(_.toString)
    val selected = if keys.contains("backend") then params.getParam("backend") else "pure"
    val state = new GaleBackendState
    state.backend = selected
    state.selectBackend()
    val runtime = s""""forkId":${quote(forkId)},"lane":${quote(System.getProperty("gale.bench.lane", "unset"))},"jdk":${quote(System.getProperty("java.version"))},"vectorModule":${ModuleLayer.boot().findModule("jdk.incubator.vector").isPresent},"blas":${quote(dev.ludovic.netlib.blas.BLAS.getInstance().getClass.getName)},"lapack":${quote(dev.ludovic.netlib.lapack.LAPACK.getInstance().getClass.getName)}"""
    val fields = s"""$runtime,"benchmark":${quote(params.getBenchmark)},"fixture":${quote(fixture)},"fixtureSha256":${quote(fingerprint(a, rhs))},"backend":${quote(selected)},"params":${quote(keys.sorted.map(k => s"$k=${params.getParam(k)}").mkString(","))},"sourceSha":${quote(System.getProperty("gale.bench.sourceSha", "unrecorded"))}"""
    def record(status: String, metrics: Vector[Metric], message: String): Unit =
      val encoded = metrics.map(m => s"""{"name":${quote(m.name)},"value":${m.value},"limit":${m.limit}}""").mkString("[", ",", "]")
      val line = s"""{$fields,"status":${quote(status)},"metrics":$encoded,"message":${quote(message)}}\n"""
      val default = System.getProperty(NetlibSidecarProperty, DefaultNetlibSidecar) + ".validation.jsonl"
      val path = Paths.get(System.getProperty(SidecarProperty, default))
      Option(path.toAbsolutePath.getParent).foreach(Files.createDirectories(_))
      Files.write(path, line.getBytes(StandardCharsets.UTF_8), StandardOpenOption.CREATE, StandardOpenOption.APPEND)
      System.err.println(s"[breeze-validation] ${params.getBenchmark} $selected $fixture $status $message")
    try
      val metrics = validate(operation(params), a, rhs, state.selected)
      record("passed", metrics, s"${metrics.size} numerical checks")
    catch
      case NonFatal(e) =>
        record("failed", Vector.empty, Option(e.getMessage).getOrElse(e.getClass.getName))
        throw new IllegalStateException(s"Benchmark preflight failed: ${params.getBenchmark} fixture=$fixture backend=$selected", e)
