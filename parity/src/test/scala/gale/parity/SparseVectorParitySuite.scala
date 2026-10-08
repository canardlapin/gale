package gale.parity

import breeze.linalg.CSCMatrix as BCSC
import breeze.linalg.DenseVector as BDV
import breeze.linalg.SparseVector as BSV
import breeze.linalg.VectorBuilder
import breeze.linalg.max as bMax
import breeze.linalg.min as bMin
import breeze.linalg.norm as bNorm
import breeze.linalg.sum as bSum
import gale.linalg.*
import gale.parity.ParitySupport.*
import gale.sparse.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAllNoShrink

/** Breeze differential tests for [[gale.sparse.SparseVector]], `CSR.rowSparse`,
  * `CSC.colSparse`, and `CSR * SparseVector`, against Breeze's `SparseVector`,
  * `VectorBuilder`, and `CSCMatrix * SparseVector`.
  *
  * Inputs are unsorted `(index, value)` lists with repeated indices and explicit
  * zeros, handed unchanged to Gale's `fromEntries` and Breeze's `VectorBuilder`.
  * Lengths are `1..64`, plus a fixed length-64K case. Replay with [[ParitySeed]].
  *
  * Tolerances (`ε = 2^-52`): stored patterns, `+`, `-`, scaling, `max`/`min`,
  * and `normInf` are exact (one rounding per entry, or a selection; derived
  * operations start from bit-identical operands). A value built from `k`
  * repeated indices is within `2 k ε Σ|v|`: Breeze's builders do not sum
  * duplicates strictly in input order, so three or more can differ in the last
  * bit. `sum`, `dot`, and `A * x` use `|Δ| ≤ 2 k ε Σ|terms|`
  * for `k` terms; `norm1` and `norm2` are relative `4 k ε`.
  *
  * Agreements worth stating, since they are easy to get wrong in a port: Breeze
  * also sums duplicates, stores explicit zeros, keeps `+`/`-` cancellations and
  * `x * 0.0` entries active, includes the implicit zero in `max`/`min`, and
  * skips a NaN facing an implicit zero in `dot` and `CSCMatrix * SparseVector`.
  *
  * Divergences from Breeze 2.1.0 (asserted as Gale's documented semantics with
  * Breeze's observed behaviour pinned; listed in the Breeze migration guide):
  *   1. `max`/`min` of a length-0 vector: Gale throws `LinAlgError.EmptyInput`;
  *      Breeze's `max` returns `-Inf` and `min` throws `IllegalArgumentException`.
  *   1. Duplicates: Gale offers `DuplicatePolicy.Sum` (the default, as Breeze),
  *      `Last` and `Error`; Breeze always sums.
  *   1. `compact` returns a new vector; Breeze's `compact()` mutates in place.
  *   1. `norm2` is scaled: `1e300` entries give a finite norm; Breeze overflows.
  *   1. `CSR * SparseVector` returns a dense `DVec`; Breeze's `CSCMatrix *
  *      SparseVector` returns a `SparseVector`.
  *   1. An out-of-range `apply` throws `LinAlgError.IndexOutOfBounds`; Breeze
  *      throws `IndexOutOfBoundsException`.
  *   1. `Sparse.coo(...).toCSR()` keeps explicit zeros, so `rowSparse` and
  *      `colSparse` report them as active; Breeze's `CSCMatrix.Builder` drops
  *      them (its `VectorBuilder` keeps them).
  */
class SparseVectorParitySuite extends ScalaCheckSuite:
  override def scalaCheckInitialSeed =
    ParitySeed.initial("P2Kv_CDT0S26ns7yNc_OL6jFN88kBW8sa5Sl0pvMr24=")

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(60).withWorkers(1)

  private val Eps = math.ulp(1.0)
  private val Inf = Double.PositiveInfinity
  private val NaN = Double.NaN

  private val lengthGen = Gen.choose(1, 64)
  private val seedGen = Gen.choose(1L, 9_000_000L)
  private val caseGen: Gen[(Int, Long)] = Gen.zip(lengthGen, seedGen)
  private val scalarGen: Gen[Double] = Gen.oneOf(Gen.const(0.0), Gen.const(-1.0), Gen.chooseNum(-10.0, 10.0))

  /** Up to `2 · length` unsorted entries; roughly 15% are explicit zeros and
    * repeated indices are common for short vectors.
    */
  private def entryData(length: Int, seed: Long): Seq[(Int, Double)] =
    val rng = new scala.util.Random(seed)
    val count = rng.nextInt(2 * length + 1)
    Seq.fill(count) {
      val index = rng.nextInt(length)
      val value = if rng.nextDouble() < 0.15 then 0.0 else rng.nextDouble() * 2.0 - 1.0
      (index, value)
    }

  private def gale(length: Int, entries: Seq[(Int, Double)]): SparseVector =
    SparseVector.fromEntries(length, entries)

  private def breeze(length: Int, entries: Seq[(Int, Double)]): BSV[Double] =
    val builder = new VectorBuilder[Double](length)
    entries.foreach((i, v) => builder.add(i, v))
    builder.toSparseVector

  private def breezeActive(b: BSV[Double]): (Seq[Int], Seq[Double]) =
    (b.index.take(b.activeSize).toSeq, b.data.take(b.activeSize).toSeq)

  /** Same length, active pattern and active values. `indexTol(i)` bounds the
    * value difference at index `i`; it is zero except where duplicates were
    * summed (see [[duplicateTol]]).
    */
  private def assertSameStorage(
      g: SparseVector,
      b: BSV[Double],
      clue: => String,
      indexTol: Int => Double = _ => 0.0
  ): Unit =
    val (bIdx, bVal) = breezeActive(b)
    assertEquals(g.length, b.length, s"length $clue")
    assertEquals(g.activeIndices, bIdx.toIndexedSeq, s"active indices $clue (gale=$g breeze=$b)")
    for (index, p) <- g.activeIndices.zipWithIndex do
      val gv = g.activeValues(p)
      val bv = bVal(p)
      if !(gv == bv || (gv.isNaN && bv.isNaN) || math.abs(gv - bv) <= indexTol(index)) then
        fail(s"active value at $index $clue: gale=$gv breeze=$bv tol=${indexTol(index)}")

  /** Breeze's `VectorBuilder` and `CSCMatrix.Builder` do not sum repeated
    * indices strictly in input order, so a sum of three or more duplicates can
    * differ in the last bits: bound it by `2 k ε Σ|v|` over the `k` duplicates.
    */
  private def duplicateTol[K](entries: Seq[(K, Double)]): K => Double =
    val grouped = entries.groupBy(_._1).view.mapValues(vs => 2.0 * vs.length * Eps * vs.map(e => math.abs(e._2)).sum).toMap
    key => grouped.getOrElse(key, 0.0)

  /** A Breeze vector holding exactly Gale's stored entries, so that derived
    * operations start from bit-identical operands.
    */
  private def toBreeze(g: SparseVector): BSV[Double] =
    new BSV(g.activeIndices.toArray, g.activeValues.toArray, g.length)

  private def assertDenseEqual(g: DVec, b: BDV[Double], clue: => String): Unit =
    assertVecClose(g, b, 0.0, clue)

  private def assertSumClose(g: Double, b: Double, k: Int, sumAbs: Double, clue: => String): Unit =
    val bound = 2.0 * math.max(k, 1) * Eps * sumAbs
    if !(math.abs(g - b) <= bound) then
      fail(s"$clue: gale=$g breeze=$b |Δ|=${math.abs(g - b)} bound=2kεΣ|x|=$bound")

  private def assertRelClose(g: Double, b: Double, rel: Double, clue: => String): Unit =
    if !(g == b || math.abs(g - b) <= rel * math.max(math.abs(g), math.abs(b))) then
      fail(s"$clue: gale=$g breeze=$b tol=$rel")

  private def assertExact(g: Double, b: Double, clue: => String): Unit =
    if !(g == b || (g.isNaN && b.isNaN)) then fail(s"$clue: gale=$g breeze=$b (exact)")

  // ---------------------------------------------------------------------------
  // Construction and access
  // ---------------------------------------------------------------------------

  property("fromEntries matches VectorBuilder: sorted pattern, summed duplicates, explicit zeros") {
    forAllNoShrink(caseGen) { (sample: (Int, Long)) =>
      val (length, seed) = sample
      val entries = entryData(length, seed)
      val clue = s"length=$length seed=$seed entries=${entries.length}"
      val g = gale(length, entries)
      val b = breeze(length, entries)
      assertSameStorage(g, b, clue, duplicateTol(entries))
      val tol = duplicateTol(entries)
      val dense = g.toDense
      val bDense = b.toDenseVector
      for i <- 0 until length do
        assertScalarClose(g(i), b(i), tol(i), s"apply($i) $clue")
        assertScalarClose(dense(i), bDense(i), tol(i), s"toDense($i) $clue")
    }
  }

  property("compact drops exactly the zeros Breeze's compact() drops") {
    forAllNoShrink(caseGen) { (sample: (Int, Long)) =>
      val (length, seed) = sample
      val entries = entryData(length, seed)
      val g = gale(length, entries)
      val b = toBreeze(g)
      val compacted = b.copy
      compacted.compact()
      assertSameStorage(g.compact, compacted, s"compact length=$length seed=$seed")
      // Gale's compact is persistent: the original keeps its explicit zeros.
      assertEquals(g.activeSize, b.activeSize)
    }
  }

  // ---------------------------------------------------------------------------
  // Arithmetic, products and reductions
  // ---------------------------------------------------------------------------

  property("dot, +, -, scaling, sum, norms, max and min match Breeze") {
    forAllNoShrink(caseGen, scalarGen) { (sample: (Int, Long), alpha: Double) =>
      val (length, seed) = sample
      val leftEntries = entryData(length, seed)
      val rightEntries = entryData(length, seed + 1L)
      val clue = s"length=$length seed=$seed alpha=$alpha"
      val gx = gale(length, leftEntries)
      val gy = gale(length, rightEntries)
      val bx = toBreeze(gx)
      val by = toBreeze(gy)
      val denseData = vectorData(length, seed + 2L)
      val dense = galeVector(denseData)
      val bDense = breezeVector(denseData)
      val xd = bx.toDenseVector
      val yd = by.toDenseVector

      val productAbs = (0 until length).map(i => math.abs(xd(i) * yd(i))).sum
      assertSumClose(gx.dot(gy), bx.dot(by), gx.activeSize, productAbs, s"sparse·sparse $clue")
      val denseAbs = (0 until length).map(i => math.abs(xd(i) * denseData(i))).sum
      assertSumClose(gx.dot(dense), bx.dot(bDense), gx.activeSize, denseAbs, s"sparse·dense $clue")

      assertSameStorage(gx + gy, bx + by, s"+ $clue")
      assertSameStorage(gx - gy, bx - by, s"- $clue")
      assertSameStorage(gx * alpha, bx * alpha, s"* $clue")

      val absX = (0 until length).map(i => math.abs(xd(i))).sum
      assertSumClose(gx.sum, bSum(bx), gx.activeSize, absX, s"sum $clue")
      assertRelClose(gx.norm1, bNorm(bx, 1.0), 4.0 * math.max(gx.activeSize, 1) * Eps, s"norm1 $clue")
      assertRelClose(gx.norm2, bNorm(bx), 4.0 * math.max(gx.activeSize, 1) * Eps, s"norm2 $clue")
      assertExact(gx.normInf, bNorm(bx, Inf), s"normInf $clue")
      // Both include the implicit zero whenever activeSize < length.
      assertExact(gx.max, bMax(bx), s"max $clue")
      assertExact(gx.min, bMin(bx), s"min $clue")
    }
  }

  test("length-64K vectors at about 1% density") {
    val length = 65536
    for seed <- Seq(11L, 12L) do
      val rng = new scala.util.Random(seed)
      val left = Seq.fill(700)((rng.nextInt(length), rng.nextDouble() * 2.0 - 1.0))
      val right = Seq.fill(700)((rng.nextInt(length), rng.nextDouble() * 2.0 - 1.0))
      val gx = gale(length, left)
      val gy = gale(length, right)
      assertSameStorage(gx, breeze(length, left), s"construction seed=$seed", duplicateTol(left))
      val bx = toBreeze(gx)
      val by = toBreeze(gy)
      assertSameStorage(gx + gy, bx + by, s"+ seed=$seed")
      val absX = gx.norm1
      assertSumClose(gx.dot(gy), bx.dot(by), gx.activeSize, absX * gy.normInf, s"dot seed=$seed")
      assertSumClose(gx.sum, bSum(bx), gx.activeSize, absX, s"sum seed=$seed")
      assertRelClose(gx.norm2, bNorm(bx), 4.0 * gx.activeSize * Eps, s"norm2 seed=$seed")
      assertExact(gx.max, bMax(bx), s"max seed=$seed")
      assertExact(gx.min, bMin(bx), s"min seed=$seed")
  }

  // ---------------------------------------------------------------------------
  // Compressed rows/columns and CSR * SparseVector
  // ---------------------------------------------------------------------------

  private val matrixCaseGen: Gen[(Int, Int, Long)] =
    Gen.zip(Gen.choose(1, 32), Gen.choose(1, 32), seedGen)

  /** Unsorted COO entries with duplicates and explicit zeros. */
  private def matrixEntries(rows: Int, cols: Int, seed: Long): Seq[(Int, Int, Double)] =
    val rng = new scala.util.Random(seed)
    Seq.fill(rng.nextInt(rows * cols + 1)) {
      val value = if rng.nextDouble() < 0.15 then 0.0 else rng.nextDouble() * 2.0 - 1.0
      (rng.nextInt(rows), rng.nextInt(cols), value)
    }

  private def galeCsr(rows: Int, cols: Int, entries: Seq[(Int, Int, Double)]): CSR =
    entries.foldLeft(Sparse.coo(rows, cols))((b, e) => b.add(e._1, e._2, e._3)).toCSR()

  private def breezeCsc(rows: Int, cols: Int, entries: Seq[(Int, Int, Double)]): BCSC[Double] =
    val builder = new BCSC.Builder[Double](rows, cols)
    entries.foreach((i, j, v) => builder.add(i, j, v))
    builder.result

  private def compacted(b: BSV[Double]): BSV[Double] =
    val out = b.copy
    out.compact()
    out

  /** Distinct coordinates of `line` in the COO input, sorted: the pattern a
    * Gale row or column keeps, explicit zeros included.
    */
  private def cellsIn(entries: Seq[(Int, Int, Double)], line: Int, byColumn: Boolean): IndexedSeq[Int] =
    entries
      .collect { case (i, j, _) if (if byColumn then j else i) == line => if byColumn then i else j }
      .distinct
      .sorted
      .toIndexedSeq

  /** Column `j` of a Breeze CSC exactly as stored (explicit zeros included). */
  private def breezeStoredColumn(a: BCSC[Double], j: Int): BSV[Double] =
    val from = a.colPtrs(j)
    val until = a.colPtrs(j + 1)
    new BSV(a.rowIndices.slice(from, until), a.data.slice(from, until), until - from, a.rows)

  property("rowSparse/colSparse keep the stored pattern and CSR * SparseVector matches CSCMatrix * SparseVector") {
    forAllNoShrink(matrixCaseGen) { (sample: (Int, Int, Long)) =>
      val (rows, cols, seed) = sample
      val clue = s"${rows}x$cols seed=$seed"
      val entries = matrixEntries(rows, cols, seed)
      val csr = galeCsr(rows, cols, entries)
      val csc = csr.toCSC
      val bA = breezeCsc(rows, cols, entries)
      val bAt = breezeCsc(cols, rows, entries.map((i, j, v) => (j, i, v)))
      val cellTol = duplicateTol(entries.map((i, j, v) => ((i, j), v)))
      // Breeze's CSCMatrix.Builder drops explicit zeros (divergence 7), so the
      // stored patterns are compared after compaction; Gale's explicit zeros
      // are checked separately against the COO input.
      for j <- 0 until cols do
        val g = csc.colSparse(j)
        assertSameStorage(g.compact, compacted(breezeStoredColumn(bA, j)), s"colSparse($j) $clue", i => cellTol((i, j)))
        assertEquals(g.activeIndices, cellsIn(entries, j, byColumn = true), s"colSparse($j) pattern $clue")
      for i <- 0 until rows do
        val g = csr.rowSparse(i)
        assertSameStorage(g.compact, compacted(breezeStoredColumn(bAt, i)), s"rowSparse($i) $clue", j => cellTol((i, j)))
        assertEquals(g.activeIndices, cellsIn(entries, i, byColumn = false), s"rowSparse($i) pattern $clue")

      // The product starts from Gale's stored entries, bit for bit.
      val stored = (0 until rows).flatMap { i =>
        val r = csr.rowSparse(i)
        r.activeIndices.zip(r.activeValues).map((j, v) => (i, j, v))
      }
      val bStored = breezeCsc(rows, cols, stored)
      val xEntries = entryData(cols, seed + 3L)
      val gx = gale(cols, xEntries)
      val bx = toBreeze(gx)
      val gy = csr * gx
      val by = (bStored * bx).toDenseVector
      val ad = bStored.toDense
      val xd = bx.toDenseVector
      for i <- 0 until rows do
        val rowAbs = (0 until cols).map(j => math.abs(ad(i, j) * xd(j))).sum
        assertSumClose(gy(i), by(i), cols, rowAbs, s"(A * x)($i) $clue")
        assertSumClose(gy(i), csr.rowSparse(i).dot(gx), cols, rowAbs, s"(A * x)($i) vs rowSparse dot $clue")
    }
  }

  // ---------------------------------------------------------------------------
  // Non-finite values and implicit zeros (agreements)
  // ---------------------------------------------------------------------------

  test("a NaN facing an implicit zero is skipped by dot and A * x in both libraries") {
    val gNaN = SparseVector.fromEntries(3, Seq(0 -> NaN))
    val gOne = SparseVector.fromEntries(3, Seq(1 -> 1.0))
    val bNaN = breeze(3, Seq(0 -> NaN))
    val bOne = breeze(3, Seq(1 -> 1.0))
    assertExact(gNaN.dot(gOne), 0.0, "gale sparse·sparse")
    assertExact(bNaN.dot(bOne), 0.0, "breeze sparse·sparse")
    val denseWithNaN = Array(NaN, 1.0, 1.0)
    assertExact(gOne.dot(galeVector(denseWithNaN)), 1.0, "gale sparse·dense")
    assertExact(bOne.dot(breezeVector(denseWithNaN)), 1.0, "breeze sparse·dense")
    // An active entry, even an explicit zero, does meet the NaN.
    assert(SparseVector.fromEntries(3, Seq(0 -> 0.0)).dot(galeVector(denseWithNaN)).isNaN)
    assert(breeze(3, Seq(0 -> 0.0)).dot(breezeVector(denseWithNaN)).isNaN)

    val entries = Seq((0, 0, 1.0), (0, 2, NaN), (1, 1, 2.0))
    val x = Seq(0 -> 1.0, 1 -> 1.0)
    val gy = galeCsr(2, 3, entries) * gale(3, x)
    val by = (breezeCsc(2, 3, entries) * breeze(3, x)).toDenseVector
    assertDenseEqual(gy, by, "A * x with NaN facing an implicit zero")
    assertEquals(gy.toSeq, Seq(1.0, 2.0))
  }

  test("max and min include the implicit zero; NaN propagates") {
    val negative = Seq(0 -> -1.0, 2 -> -3.0)
    val positive = Seq(0 -> 1.0, 2 -> 3.0)
    assertExact(gale(4, negative).max, bMax(breeze(4, negative)), "max of negatives")
    assertExact(gale(4, negative).max, 0.0, "max of negatives is the implicit zero")
    assertExact(gale(4, positive).min, bMin(breeze(4, positive)), "min of positives")
    assertExact(gale(2, negative.take(1) :+ (1 -> -2.0)).max, -1.0, "fully stored: no implicit zero")
    for at <- Seq(0, 2) do
      val withNaN = Seq(0 -> 1.0, 2 -> 3.0).map((i, v) => if i == at then (i, NaN) else (i, v))
      assert(gale(4, withNaN).max.isNaN && bMax(breeze(4, withNaN)).isNaN, s"max NaN at $at")
      assert(gale(4, withNaN).min.isNaN && bMin(breeze(4, withNaN)).isNaN, s"min NaN at $at")
  }

  test("cancellations and zero scaling stay stored in both libraries") {
    val a = Seq(0 -> 1.0, 2 -> 2.0)
    val b = Seq(0 -> -1.0)
    assertSameStorage(gale(3, a) + gale(3, b), breeze(3, a) + breeze(3, b), "cancellation")
    assertEquals((gale(3, a) + gale(3, b)).activeSize, 2)
    assertSameStorage(gale(3, a) * 0.0, breeze(3, a) * 0.0, "scale by zero")
    // Summed duplicates that cancel remain an explicit zero.
    val cancelling = Seq(1 -> 1e16, 1 -> 1.0, 1 -> -1e16, 0 -> -0.0)
    assertSameStorage(gale(5, cancelling), breeze(5, cancelling), "cancelling duplicates")
  }

  test("empty and all-implicit vectors: sums and norms agree") {
    val g = SparseVector.zeros(0)
    val b = BSV.zeros[Double](0)
    assertExact(g.sum, bSum(b), "sum")
    assertExact(g.norm1, bNorm(b, 1.0), "norm1")
    assertExact(g.norm2, bNorm(b), "norm2")
    assertExact(g.normInf, bNorm(b, Inf), "normInf")
    assertExact(SparseVector.zeros(5).max, bMax(BSV.zeros[Double](5)), "max of all-implicit")
    assertExact(SparseVector.zeros(5).dot(SparseVector.zeros(5)), 0.0, "dot of all-implicit")
  }

  // ---------------------------------------------------------------------------
  // Divergences
  // ---------------------------------------------------------------------------

  test("divergence: max/min of a length-0 vector throw EmptyInput; Breeze returns -Inf or throws IAE") {
    intercept[LinAlgError.EmptyInput](SparseVector.zeros(0).max)
    intercept[LinAlgError.EmptyInput](SparseVector.zeros(0).min)
    assertEquals(bMax(BSV.zeros[Double](0)), Double.NegativeInfinity)
    intercept[IllegalArgumentException](bMin(BSV.zeros[Double](0)))
  }

  test("divergence: DuplicatePolicy.Last and Error have no Breeze counterpart (Breeze sums)") {
    val entries = Seq(1 -> 2.0, 3 -> 0.0, 1 -> 3.0)
    assertEquals(SparseVector.fromEntries(5, entries, DuplicatePolicy.Last)(1), 3.0)
    assert(SparseVector.tryFromEntries(5, entries, DuplicatePolicy.Error).isLeft)
    assertEquals(SparseVector.fromEntries(5, entries)(1), 5.0)
    assertEquals(breeze(5, entries)(1), 5.0)
    val viaApply = BSV(5)((1, 2.0), (3, 0.0), (1, 3.0))
    assertEquals(viaApply(1), 5.0)
  }

  test("divergence: compact returns a new vector; Breeze compact() mutates in place") {
    val entries = Seq(0 -> 0.0, 1 -> 1.0, 2 -> -0.0)
    val g = gale(3, entries)
    val compacted = g.compact
    assertEquals((g.activeSize, compacted.activeSize), (3, 1))
    val b = breeze(3, entries)
    b.compact()
    assertEquals(b.activeSize, 1)
  }

  test("divergence: norm2 is scaled; Breeze overflows") {
    val entries = Seq(0 -> 1e300, 2 -> 1e300)
    assertRelClose(gale(4, entries).norm2, math.sqrt(2.0) * 1e300, 4 * Eps, "gale norm2")
    assertEquals(bNorm(breeze(4, entries)), Inf)
  }

  test("divergence: CSR * SparseVector is dense; Breeze CSCMatrix * SparseVector is sparse") {
    val entries = Seq((0, 0, 1.0), (1, 1, 2.0))
    val gy: DVec = galeCsr(3, 2, entries) * gale(2, Seq(0 -> 1.0))
    val by: BSV[Double] = breezeCsc(3, 2, entries) * breeze(2, Seq(0 -> 1.0))
    assertEquals(gy.length, 3)
    assertDenseEqual(gy, by.toDenseVector, "values agree")
  }

  test("divergence: CSR/CSC keep explicit zeros from the builder; Breeze's CSCMatrix.Builder drops them") {
    val entries = Seq((0, 0, 1.0), (1, 0, 0.0), (1, 1, 2.0))
    val csr = galeCsr(2, 2, entries)
    assertEquals(csr.rowSparse(1).activeIndices, IndexedSeq(0, 1))
    assertEquals(csr.toCSC.colSparse(0).activeSize, 2)
    val bA = breezeCsc(2, 2, entries)
    assertEquals(bA.activeSize, 2)
    assertEquals(breezeStoredColumn(bA, 0).activeSize, 1)
    // Breeze's VectorBuilder, in contrast, keeps the explicit zero, as Gale does.
    assertEquals(breeze(2, Seq(0 -> 0.0)).activeSize, 1)
  }

  test("divergence: out-of-range apply throws LinAlgError.IndexOutOfBounds") {
    intercept[LinAlgError.IndexOutOfBounds](SparseVector.zeros(4)(10))
    val b = BSV.zeros[Double](4)
    intercept[IndexOutOfBoundsException](b(10))
  }
