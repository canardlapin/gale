package gale.linalg

import gale.kernel.DoubleKernels
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.PlatformMath.fma

/** Equivalence court for the wide column-pivoted QR screen.
  *
  * Widths past the compact regime discard trailing-column scans that a bound on the downdated norm proves cannot win.
  * That is a performance decision only: the pivot order, the tie choice, and the resulting factor must stay exactly
  * what recomputing every column would have produced, because the permutation is a published contract and the rank
  * decision reads the factor's diagonal.
  *
  * The oracle below is a deliberate transcription of the exhaustive algorithm rather than a call into it, so the suite
  * does not verify the selection logic against itself. Its Householder construction and reflector application are the
  * production kernels verbatim; only the pivot search differs, and only by scanning every trailing column
  * unconditionally.
  */
class PivotedQRScreenSuite extends munit.FunSuite:

  private case class Factored(permutation: Array[Int], r: Array[Double])

  private def exhaustive(m: Int, n: Int, source: Array[Double]): Factored =
    val limit = math.min(m, n)
    val r = DoubleArray.fromArray(source.clone())
    val reflectors = DoubleArray.alloc(m * limit)
    val tau = DoubleArray.alloc(limit)
    val scratch = DoubleArray.alloc(math.max(math.max(m, n), 1))
    val permutation = Array.tabulate(n)(i => i)

    var k = 0
    while k < limit do
      var pivot = k
      var bestNorm = -1.0
      var col = k
      while col < n do
        val norm = DoubleKernels.dnrm2(m - k, r, k * n + col, n)
        if norm > bestNorm then
          bestNorm = norm
          pivot = col
        col += 1
      val selectedNorm =
        if bestNorm >= 0.0 then bestNorm
        else DoubleKernels.dnrm2(m - k, r, k * n + pivot, n)

      if pivot != k then
        var row = 0
        while row < m do
          val left = row * n + k
          val right = row * n + pivot
          val moved = r(left)
          r(left) = r(right)
          r(right) = moved
          row += 1
        val original = permutation(k)
        permutation(k) = permutation(pivot)
        permutation(pivot) = original

      // Production `factorHouseholder`, verbatim.
      val diagIndex = k * n + k
      val x0 = r(diagIndex)
      reflectors(k * limit + k) = 1.0
      if selectedNorm > 0.0 then
        val beta = if x0 >= 0.0 then -selectedNorm else selectedNorm
        val x0OverBeta = x0 / beta
        val denominatorOverBeta = x0OverBeta - 1.0
        tau(k) = 1.0 - x0OverBeta
        r(diagIndex) = beta
        var i = k + 1
        while i < m do
          val idx = i * n + k
          reflectors(i * limit + k) = (r(idx) / beta) / denominatorOverBeta
          r(idx) = 0.0
          i += 1
      else
        tau(k) = 0.0
        var i = k + 1
        while i < m do
          r(i * n + k) = 0.0
          i += 1

      // Production `applyReflectorToColumns`, verbatim.
      val tauK = tau(k)
      val width = n - (k + 1)
      if tauK != 0.0 && width > 0 then
        var j = 0
        while j < width do
          scratch(j) = 0.0
          j += 1
        var i = k
        while i < m do
          val vi = reflectors(i * limit + k)
          val rRow = i * n + (k + 1)
          j = 0
          while j < width do
            scratch(j) = fma(vi, r(rRow + j), scratch(j))
            j += 1
          i += 1
        j = 0
        while j < width do
          scratch(j) = tauK * scratch(j)
          j += 1
        i = k
        while i < m do
          val vi = reflectors(i * limit + k)
          val rRow = i * n + (k + 1)
          j = 0
          while j < width do
            r(rRow + j) = fma(-vi, scratch(j), r(rRow + j))
            j += 1
          i += 1
      k += 1

    Factored(permutation, Array.tabulate(m * n)(i => r(i)))

  private val pivoted = QROptions(pivoting = QRPivoting.Column)

  private def checkAgreement(label: String, m: Int, n: Int, values: Array[Double]): Unit =
    val expected = exhaustive(m, n, values)
    val produced = Matrix.dense(m, n, values.toIndexedSeq).qr(pivoted)

    assertEquals(
      produced.columnPermutation.toIndexSeq.toVector,
      expected.permutation.toVector,
      s"$label: permutation diverged from exhaustive recomputation"
    )
    var index = 0
    while index < m * n do
      val row = index / n
      val col = index % n
      val actual = produced.r(row, col)
      val reference = expected.r(index)
      // Bit equality, with NaN treated as a single value: the screen must not
      // perturb the factor at all, not merely leave it close.
      assert(
        java.lang.Double.compare(actual, reference) == 0,
        s"$label: R($row, $col) was $actual, exhaustive recomputation gave $reference"
      )
      index += 1

  /** Deterministic bit-mixing source; no dependence on a platform RNG. */
  private final class Source(private var state: Long):
    def next(): Double =
      state = state * 6364136223846793005L + 1442695040888963407L
      val bits = (state >>> 11).toDouble / 9007199254740992.0
      bits * 2.0 - 1.0

  test("wide pivoted QR matches exhaustive recomputation on general designs") {
    val shapes = Vector((40, 9), (64, 12), (33, 16), (128, 24), (17, 20), (9, 9), (12, 40))
    shapes.zipWithIndex.foreach { case ((m, n), shapeIndex) =>
      val source = new Source(0x5deece66dL + shapeIndex)
      val values = Array.fill(m * n)(source.next())
      checkAgreement(s"general ${m}x$n", m, n, values)
    }
  }

  test("wide pivoted QR matches exhaustive recomputation on near-tied columns") {
    // Columns separated by less than the screen's own bound width are exactly
    // the case a downdated estimate can misjudge, and the case that made a
    // vendor `dgeqp3` disagree with this library's permutation.
    val m = 96
    val n = 14
    (0 until 24).foreach { trial =>
      val source = new Source(0x9e3779b97f4a7c15L + trial)
      val base = Array.fill(m)(source.next())
      val values = Array.ofDim[Double](m * n)
      var col = 0
      while col < n do
        val perturbation = math.pow(10.0, -(6.0 + col.toDouble))
        var row = 0
        while row < m do
          values(row * n + col) = base(row) + perturbation * source.next()
          row += 1
        col += 1
      checkAgreement(s"near-tied trial $trial", m, n, values)
    }
  }

  test("wide pivoted QR matches exhaustive recomputation on exact duplicates") {
    // Exactly equal norms must resolve to the same leftmost column both ways.
    val m = 48
    val n = 12
    val source = new Source(0x2545f4914f6cdd1dL)
    val column = Array.fill(m)(source.next())
    val values = Array.ofDim[Double](m * n)
    var row = 0
    while row < m do
      var col = 0
      while col < n do
        values(row * n + col) = if col % 3 == 0 then column(row) else source.next()
        col += 1
      row += 1
    checkAgreement("exact duplicates", m, n, values)
  }

  test("wide pivoted QR matches exhaustive recomputation on rank-deficient designs") {
    // Trailing norms collapse toward zero here, which widens the screen's
    // bounds until it stops discarding anything. The answer must not move.
    val m = 80
    val n = 16
    val rank = 5
    val source = new Source(0xda3e39cb94b95bdbL)
    val basis = Array.fill(m * rank)(source.next())
    val values = Array.ofDim[Double](m * n)
    var col = 0
    while col < n do
      val weights = Array.fill(rank)(source.next())
      var row = 0
      while row < m do
        var acc = 0.0
        var t = 0
        while t < rank do
          acc += basis(row * rank + t) * weights(t)
          t += 1
        values(row * n + col) = acc
        row += 1
      col += 1
    checkAgreement("rank deficient", m, n, values)
  }

  test("wide pivoted QR matches exhaustive recomputation across extreme scales") {
    // Columns spanning the range where a squared norm overflows or underflows,
    // so the screen must fall back rather than trust its estimate.
    val m = 40
    val n = 10
    val exponents = Vector(0, -170, 160, -320, 300, 20, -60, 120, -200, 0)
    val source = new Source(0x14057b7ef767814fL)
    val values = Array.ofDim[Double](m * n)
    var col = 0
    while col < n do
      val scale = math.pow(10.0, exponents(col).toDouble)
      var row = 0
      while row < m do
        values(row * n + col) = source.next() * scale
        row += 1
      col += 1
    checkAgreement("extreme scales", m, n, values)
  }

  test("wide pivoted QR matches exhaustive recomputation with non-finite columns") {
    val m = 24
    val n = 10
    val source = new Source(0xbf58476d1ce4e5b9L)
    val values = Array.fill(m * n)(source.next())
    values(3 * n + 4) = Double.NaN
    values(7 * n + 6) = Double.PositiveInfinity
    values(11 * n + 6) = Double.NegativeInfinity
    checkAgreement("non-finite", m, n, values)
  }

  test("wide pivoted QR reserves the downdated norm and bound regions") {
    val requirement = DenseWorkspace.qrRequirement(4096, 24, pivoted)
    assertEquals(requirement.map(_.doubleElements), Right(4096))

    val narrow = DenseWorkspace.qrRequirement(30, 24, pivoted)
    assertEquals(narrow.map(_.doubleElements), Right(72))

    val compact = DenseWorkspace.qrRequirement(4096, 6, pivoted)
    assertEquals(compact.map(_.doubleElements), Right(4096))
  }
