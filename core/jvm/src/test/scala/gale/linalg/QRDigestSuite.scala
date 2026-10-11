package gale.linalg

/** Bit-exact regression guard for the portable dense QR and the default least-squares facade on the JVM.
  *
  * Each case of a fixed corpus is reduced to one 64-bit digest of the raw bits the factorization returns (reflectors,
  * `tau`, `R`, rank and rank tolerance) and of `A.leastSquares(b)` (the coefficients or the error), unpivoted and
  * column-pivoted, and compared with
  * the digests in `qr-digests.txt`. The stored digests were generated from the kernel as of breeze/integration
  * 5d3950e, so any change to its arithmetic or summation order fails here; such a change must regenerate the file
  * deliberately and justify the new bits. JVM-only: the kernels are deterministic per platform, not across platforms.
  *
  * NaNs are digested canonically (`doubleToLongBits`): a NaN's sign and payload depend on which fma operand carries
  * the negation and which NaN operand the hardware propagates, neither of which is part of any numerical contract.
  * Every finite and infinite value, and whether each entry is NaN, is pinned bit for bit.
  *
  * Setting `GALE_QR_DIGEST_OUT` to a file path writes the current digests there instead of comparing.
  */
class QRDigestSuite extends munit.FunSuite:
  import QRDigestSuite.*

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

  test("QR and least squares match the stored per-case bit digests") {
    val actual = digests
    assert(actual.size > 9000, s"corpus shrank to ${actual.size} cases")
    sys.env.get("GALE_QR_DIGEST_OUT") match
      case Some(path) =>
        val lines = actual.toSeq.sortBy(_._1).map((key, digest) => s"$key\t$digest")
        java.nio.file.Files.writeString(java.nio.file.Path.of(path), lines.mkString("", "\n", "\n"))
      case None =>
        val missing = expected.keySet.diff(actual.keySet).toSeq.sorted
        assert(missing.isEmpty, s"${missing.size} stored cases vanished, e.g. ${missing.take(5).mkString("; ")}")
        val mismatches = actual.filter((key, digest) => !expected.get(key).contains(digest)).keys.toSeq.sorted
        assert(mismatches.isEmpty, s"${mismatches.size} cases changed bits, e.g. ${mismatches.take(5).mkString("; ")}")
  }

object QRDigestSuite:
  val DigestResource = "/gale/linalg/qr-digests.txt"

  private val Pivoted = QROptions(pivoting = QRPivoting.Column)

  /** Shapes on both sides of the blocked threshold, including square, wide and single-row or -column cases. */
  private val shapes: Seq[(Int, Int)] =
    val small = for
      m <- 1 to 20
      n <- 1 to 20
      if m >= n || (m + n) % 3 == 0
    yield (m, n)
    small ++ Seq((33, 7), (40, 40), (64, 16), (65, 17), (95, 95), (100, 96), (128, 32), (256, 64), (200, 120))

  /** Gaussian, small-integer (often rank-deficient), zero-column, duplicated-column, `1e±300`-scaled and graded
    * matrices, each paired with a Gaussian right-hand side.
    */
  def corpus: Seq[(String, DMat, DVec)] =
    val out = scala.collection.mutable.ArrayBuffer.empty[(String, DMat, DVec)]
    for (m, n) <- shapes do
      val big = m * n > 4000
      val r = new scala.util.Random(1000L * m + n)
      val b = Vec.tabulate(m)(_ => r.nextGaussian())
      out += ((s"gauss m=$m n=$n", Matrix.tabulate(m, n)((_, _) => r.nextGaussian()), b))
      if !big then
        out += ((s"ints m=$m n=$n", Matrix.tabulate(m, n)((_, _) => (r.nextInt(3) - 1).toDouble), b))
        out += ((s"zerocol m=$m n=$n", Matrix.tabulate(m, n)((i, j) => if j == n / 2 then 0.0 else r.nextGaussian()), b))
        out += ((s"dupcol m=$m n=$n", Matrix.tabulate(m, n)((i, j) => (i + 1).toDouble * (if j % 2 == 0 then 1.0 else 2.0)), b))
        out += ((s"huge m=$m n=$n", Matrix.tabulate(m, n)((_, _) => r.nextGaussian() * 1e300), b))
        out += ((s"tiny m=$m n=$n", Matrix.tabulate(m, n)((_, _) => r.nextGaussian() * 1e-300), b))
        out += ((s"graded m=$m n=$n", Matrix.tabulate(m, n)((i, j) => r.nextGaussian() * math.pow(10.0, (i - j) * 2.0)), b))
      if !DenseWorkspace.usesBlockedQR(m, n) then
        // Non-unit layouts reach the [A | b] gather: a transposed view, a column
        // window of a wider matrix, and a right-hand side read with stride 3 from
        // offset 4 of a longer column.
        val strided = Matrix.tabulate(m + 1, 3)((_, _) => r.nextGaussian()).col(1).slice(1, m + 1)
        out += ((s"transposed m=$m n=$n", Matrix.tabulate(n, m)((_, _) => r.nextGaussian()).t, strided))
        out += ((s"window m=$m n=$n", Matrix.tabulate(m + 2, n + 3)((_, _) => r.nextGaussian()).slice(1, m + 1, 2, n + 2), strided))
        // Non-finite entries off the leading column reach the fused norm's scaled
        // fallback; a non-finite right-hand side reaches only the transformed b.
        if m > 2 && n > 1 then
          val nan = Matrix.tabulate(m, n)((i, j) => if i == 2 && j == 1 then Double.NaN else r.nextGaussian())
          val inf = Matrix.tabulate(m, n)((i, j) => if i == 2 && j == 1 then Double.PositiveInfinity else r.nextGaussian())
          out += ((s"nancol m=$m n=$n", nan, b))
          out += ((s"infcol m=$m n=$n", inf, b))
          out += ((s"nanb m=$m n=$n", Matrix.tabulate(m, n)((_, _) => r.nextGaussian()), Vec.tabulate(m)(i => if i == 2 then Double.NaN else b(i))))
    // A right-hand side of the wrong length, once per family, on a tall and a wide
    // shape: pins the error and its precedence over the underdetermined check.
    for (m, n) <- Seq((7, 4), (3, 6)) do
      val wrong = Vec.tabulate(m + 1)(_.toDouble)
      for (key, a, _) <- out.toSeq.filter(_._1.endsWith(s" m=$m n=$n")) do
        out += ((s"${key.takeWhile(_ != ' ')} badb m=$m n=$n", a, wrong))
    out.toSeq

  /** FNV-1a over 64-bit words, each word's bytes folded high to low; every NaN folds as the canonical NaN. */
  final class Digest:
    private var state = 0xcbf29ce484222325L
    def long(x: Long): Digest =
      var shift = 56
      while shift >= 0 do
        state = (state ^ ((x >>> shift) & 0xffL)) * 0x100000001b3L
        shift -= 8
      this
    def double(x: Double): Digest = long(java.lang.Double.doubleToLongBits(x))
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

  def qrDigest(qr: QR): String =
    val d = new Digest().mat(qr.reflectors).mat(qr.r)
    var k = 0
    while k < qr.tau.length do
      d.double(qr.tau(k))
      k += 1
    d.long(qr.diagnostics.rank.fold(-1L)(_.toLong)).double(qr.diagnostics.rankTolerance.getOrElse(Double.NaN)).hex

  def solveDigest(result: Either[LinAlgError, DVec]): String =
    result match
      case Right(x)    => new Digest().long(1L).vec(x).hex
      case Left(error) => new Digest().long(2L).long(error.toString.hashCode.toLong).hex

  def digests: Map[String, String] =
    corpus.flatMap { (key, a, b) =>
      Seq(
        s"qr $key" -> qrDigest(a.qr),
        s"lstsq $key" -> solveDigest(a.leastSquares(b)),
        s"qrp $key" -> qrDigest(a.qr(Pivoted)),
        s"lstsqp $key" -> solveDigest(a.leastSquares(b, Pivoted))
      )
    }.toMap
