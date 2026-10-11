package gale.linalg

/** Bit-exact regression guard for [[Numerics.sigmoid]] on the JVM.
  *
  * A fixed corpus (random normals at several scales, the neighbourhoods of the `-36` form switch, of `-37.5` where
  * `1 + exp(x)` rounds to 1, of the `exp` overflow and underflow edges near `±709.78` and `-745.13`, subnormal and
  * extreme inputs, and the special values) is pushed through every layout the facade dispatches on: a contiguous
  * vector, an offset slice, a strided column, a contiguous matrix and a transposed (strided-row) matrix.
  *
  * Two checks pin the output bits to the hybrid kernel of breeze/integration 96c95e5:
  *
  *   - every element must equal [[SigmoidDigestSuite.reference96c95e5]], a frozen copy of that kernel evaluated with
  *     the same `Math.exp`. This holds on every platform, because both sides share the platform's `exp`.
  *   - each family's 64-bit digest must equal the one stored for this `os.arch` in `sigmoid-digests.txt`, generated
  *     from 96c95e5 itself. `Math.exp` is deterministic per platform but not across platforms, so the digests are
  *     keyed by architecture; an architecture with no stored digests skips this half.
  *
  * NaN sign and payload are unspecified (docs/user/advanced/numerical-contract.md), so digests and comparisons
  * canonicalize NaN (`doubleToLongBits`); every NaN position is still pinned.
  *
  * Setting `GALE_SIGMOID_DIGEST_OUT` to a file path writes this architecture's current digests there instead of
  * comparing.
  */
class SigmoidDigestSuite extends munit.FunSuite:
  import SigmoidDigestSuite.*

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

  test("sigmoid is bit-identical to the frozen 96c95e5 hybrid in every layout") {
    for (family, xs) <- corpus; (layout, got) <- layouts(xs) do
      var i = 0
      while i < xs.length do
        val want = reference96c95e5(xs(i))
        assert(
          bits(got(i)) == bits(want),
          s"$family/$layout: sigmoid(${xs(i)}) = ${got(i)}, the 96c95e5 hybrid gives $want"
        )
        i += 1
  }

  test("sigmoid stays bit-identical to the 96c95e5 hybrid once the kernel is JIT-compiled") {
    // A single pass over the corpus mostly runs interpreted code; repeat the mixed family until the kernel's hot
    // loops are compiled, and check the last pass of each layout.
    val xs = corpus.last._2
    val want = xs.map(x => bits(reference96c95e5(x)))
    var last = Seq.empty[(String, Array[Double])]
    for _ <- 1 to 40 do last = layouts(xs)
    for (layout, got) <- last do
      val bad = got.indices.filter(i => bits(got(i)) != want(i))
      assert(bad.isEmpty, s"$layout: ${bad.size} elements differ, e.g. sigmoid(${xs(bad.head)}) = ${got(bad.head)}")
  }

  test("sigmoid matches the stored per-family digests for this architecture") {
    val actual = digests
    assert(actual.size >= 60, s"corpus shrank to ${actual.size} cases")
    sys.env.get("GALE_SIGMOID_DIGEST_OUT") match
      case Some(path) =>
        val lines = actual.toSeq.sortBy(_._1).map((key, digest) => s"$key\t$digest")
        java.nio.file.Files.writeString(java.nio.file.Path.of(path), lines.mkString("", "\n", "\n"))
      case None =>
        val stored = expected.filter(_._1.startsWith(s"$arch "))
        assume(stored.nonEmpty, s"no stored sigmoid digests for os.arch=$arch")
        val missing = stored.keySet.diff(actual.keySet).toSeq.sorted
        assert(missing.isEmpty, s"${missing.size} stored cases vanished, e.g. ${missing.take(5).mkString("; ")}")
        val mismatches = actual.filter((key, digest) => !stored.get(key).contains(digest)).keys.toSeq.sorted
        assert(mismatches.isEmpty, s"${mismatches.size} cases changed bits, e.g. ${mismatches.take(5).mkString("; ")}")
  }

object SigmoidDigestSuite:
  val DigestResource = "/gale/linalg/sigmoid-digests.txt"

  private val arch: String = System.getProperty("os.arch")

  /** The sigmoid kernel as of breeze/integration 96c95e5, verbatim. */
  def reference96c95e5(v: Double): Double =
    if v < -36.0 then
      val t = math.exp(v)
      t / (1.0 + t)
    else 1.0 / (1.0 + math.exp(-v))

  private def bits(x: Double): Long = java.lang.Double.doubleToLongBits(x)

  /** `x` and its `k` nearest neighbours on each side. */
  private def around(x: Double, k: Int): Seq[Double] =
    val below = Iterator.iterate(x)(math.nextDown).slice(1, k + 1).toSeq
    val above = Iterator.iterate(x)(math.nextUp).slice(1, k + 1).toSeq
    below.reverse ++ (x +: above)

  private def grid(from: Double, to: Double, count: Int): Seq[Double] =
    (0 until count).map(i => from + (to - from) * i / (count - 1))

  def corpus: Seq[(String, Array[Double])] =
    val normals =
      for scale <- Seq(1e-6, 1e-2, 1.0, 4.0, 16.0, 64.0, 256.0, 1024.0)
      yield
        val r = new scala.util.Random(java.lang.Double.doubleToLongBits(scale))
        s"normal scale=$scale" -> Array.fill(4096)(r.nextGaussian() * scale)
    val uniform =
      val r = new scala.util.Random(5300L)
      Seq("uniform [-1,1)" -> Array.fill(4096)(r.nextDouble() * 2.0 - 1.0),
        "uniform [-800,800)" -> Array.fill(4096)(r.nextDouble() * 1600.0 - 800.0))
    val logMax = math.log(Double.MaxValue)
    val edges = Seq(
      "switch -36" -> (around(-36.0, 256) ++ grid(-37.0, -35.0, 513)),
      "switch +36" -> (around(36.0, 64) ++ grid(35.0, 37.0, 129)),
      "tail -37.5" -> (around(-37.5, 64) ++ around(-36.7368005696771, 64) ++ grid(-40.0, -33.0, 701)),
      "overflow ±709.78" -> (around(logMax, 64) ++ around(-logMax, 64) ++ around(-708.3964185322641, 64) ++
        grid(-712.0, -706.0, 601) ++ grid(706.0, 712.0, 61)),
      "underflow -745" -> (around(-745.1332191019411, 64) ++ around(-744.4400719213812, 64) ++
        around(745.1332191019411, 16) ++ grid(-750.0, -740.0, 1001)),
      "subnormal inputs" -> {
        val tiny = Seq(java.lang.Double.MIN_VALUE, 2 * java.lang.Double.MIN_VALUE, 1e-310, 1e-308,
          java.lang.Double.MIN_NORMAL, math.nextDown(java.lang.Double.MIN_NORMAL), 1e-300, 1e-17, 1.1102230246251565e-16)
        tiny ++ tiny.map(-_)
      },
      "extremes" -> Seq(Double.MaxValue, -Double.MaxValue, 1e300, -1e300, 1e10, -1e10, 1e5, -1e5, 800.0, -800.0,
        1000.0, -1000.0, 1e6, -1e6),
      "specials" -> Seq(0.0, -0.0, Double.PositiveInfinity, Double.NegativeInfinity, Double.NaN,
        java.lang.Double.longBitsToDouble(0xfff8000000000000L), java.lang.Double.longBitsToDouble(0x7ff0000000000001L),
        java.lang.Double.longBitsToDouble(0xfff0000000000abcL), 1.0, -1.0, 0.5, -0.5)
    ).map((name, xs) => name -> xs.toArray)
    // Every family once more with the specials interleaved, so non-finite values meet finite neighbours mid-block.
    val mixed =
      val r = new scala.util.Random(36L)
      val specials = edges.last._2
      "mixed" -> (normals ++ uniform ++ edges).flatMap(_._2).map(x => if r.nextInt(97) == 0 then specials(r.nextInt(specials.length)) else x).toArray
    normals ++ uniform ++ edges :+ mixed

  /** `sigmoid(xs)` through each layout the facade dispatches on. */
  def layouts(xs: Array[Double]): Seq[(String, Array[Double])] =
    val n = xs.length
    val contiguous = Numerics.sigmoid(Vec(xs*)).toSeq.toArray
    val padded = Vec.tabulate(n + 7)(i => if i >= 3 && i < n + 3 then xs(i - 3) else -1e9)
    val slice = Numerics.sigmoid(padded.slice(3, n + 3)).toSeq.toArray
    // Column 1 of an n × 3 matrix: stride 3, offset 1.
    val column = Numerics.sigmoid(Matrix.tabulate(n, 3)((i, j) => if j == 1 then xs(i) else Double.NaN).col(1)).toSeq.toArray
    val matrix = Numerics.sigmoid(Matrix.tabulate(1, n)((_, j) => xs(j))).valuesRowMajor.toArray
    // The transpose of an n × 2 row-major matrix: a 2 × n view whose rows have stride 2; row 0 holds xs.
    val transposed = Numerics.sigmoid(Matrix.tabulate(n, 2)((i, j) => if j == 0 then xs(i) else 0.0).t).valuesRowMajor
    val transposedRow0 = transposed.take(n).toArray
    Seq("vec" -> contiguous, "slice" -> slice, "column" -> column, "matrix" -> matrix, "transposed" -> transposedRow0)

  /** FNV-1a over 64-bit words, each word's bytes folded high to low; every NaN folds as the canonical NaN. */
  private def digest(xs: Array[Double]): String =
    var state = 0xcbf29ce484222325L
    def long(x: Long): Unit =
      var shift = 56
      while shift >= 0 do
        state = (state ^ ((x >>> shift) & 0xffL)) * 0x100000001b3L
        shift -= 8
    long(xs.length.toLong)
    xs.foreach(x => long(bits(x)))
    f"$state%016x"

  def digests: Map[String, String] =
    corpus.flatMap { (family, xs) =>
      layouts(xs).map((layout, ys) => s"$arch $family $layout" -> digest(ys)) :+ (s"$arch $family input" -> digest(xs))
    }.toMap
