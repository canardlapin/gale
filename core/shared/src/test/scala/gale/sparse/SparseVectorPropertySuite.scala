package gale.sparse

import gale.linalg.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

/** Every SparseVector operation against a dense reference. Values are small
  * multiples of 1/2, so sums and products are exact in any order and the
  * comparisons below are exact; one property covers general doubles with a
  * tolerance.
  */
class SparseVectorPropertySuite extends ScalaCheckSuite:
  override def scalaCheckInitialSeed =
    "6xYHq0dE2pTQyQd3pLqZr1yW8Hk7m1vB9cG2nJ5sT0D="

  override def scalaCheckTestParameters =
    super.scalaCheckTestParameters.withMinSuccessfulTests(200).withWorkers(1)

  // Zero appears often so explicit zeros and cancellations are exercised.
  private val exactValue: Gen[Double] =
    Gen.frequency(3 -> Gen.const(0.0), 10 -> Gen.choose(-16, 16).map(_ * 0.5))

  private val lengthGen: Gen[Int] =
    Gen.frequency(1 -> Gen.const(0), 2 -> Gen.const(1), 10 -> Gen.choose(2, 40))

  private def entriesGen(length: Int, value: Gen[Double]): Gen[List[(Int, Double)]] =
    if length == 0 then Gen.const(Nil)
    else
      Gen.choose(0, 2 * length).flatMap(n => Gen.listOfN(n, Gen.zip(Gen.choose(0, length - 1), value)))

  private def sparseGen(length: Int, value: Gen[Double] = exactValue): Gen[SparseVector] =
    entriesGen(length, value).map(entries => SparseVector.fromEntries(length, entries))

  private val sparseAny: Gen[SparseVector] =
    lengthGen.flatMap(n => sparseGen(n))

  private val sparsePair: Gen[(SparseVector, SparseVector)] =
    for
      n <- lengthGen
      a <- sparseGen(n)
      kind <- Gen.choose(0, 2)
      b <- kind match
        case 0 => sparseGen(n) // independent patterns: disjoint or overlapping
        case 1 => Gen.const(a.mapActive(v => v * 2.0 - 1.0)) // identical pattern
        case _ =>
          // pattern disjoint from `a` by construction
          val free = (0 until n).filterNot(a.activeIndices.toSet)
          Gen.listOf(Gen.zip(Gen.oneOf(free.toSeq :+ -1), exactValue)).map { entries =>
            SparseVector.fromEntries(n, entries.filter(_._1 >= 0))
          }
    yield (a, b)

  private def denseReference(
      length: Int,
      entries: Seq[(Int, Double)],
      duplicates: DuplicatePolicy
  ): Seq[Double] =
    val out = Array.fill(length)(0.0)
    entries.foreach { case (i, v) =>
      duplicates match
        case DuplicatePolicy.Sum  => out(i) += v
        case DuplicatePolicy.Last => out(i) = v
        case DuplicatePolicy.Error => ()
    }
    out.toSeq

  private def assertCanonical(v: SparseVector): Unit =
    val idx = v.activeIndices
    assertEquals(idx.length, v.activeSize)
    assertEquals(v.activeValues.length, v.activeSize)
    assert(idx.indices.drop(1).forall(k => idx(k - 1) < idx(k)), s"indices not strictly increasing: $idx")
    assert(idx.forall(i => i >= 0 && i < v.length), s"index out of range: $idx")

  property("fromEntries matches the dense reference under Sum and Last and keeps explicit zeros") {
    forAll(lengthGen.flatMap(n => entriesGen(n, exactValue).map(n -> _))) { case (n, entries) =>
      for policy <- Seq(DuplicatePolicy.Sum, DuplicatePolicy.Last) do
        val v = SparseVector.fromEntries(n, entries, duplicates = policy)
        assertCanonical(v)
        assertEquals(v.toDense.toSeq, denseReference(n, entries, policy))
        assertEquals(v.activeIndices, entries.map(_._1).distinct.sorted.toIndexedSeq)
      val distinct = entries.map(_._1).distinct.length == entries.length
      val strict = SparseVector.tryFromEntries(n, entries, duplicates = DuplicatePolicy.Error)
      assertEquals(strict.isRight, distinct)
    }
  }

  property("apply agrees with toDense at every position") {
    forAll(sparseAny) { v =>
      val dense = v.toDense
      (0 until v.length).foreach(i => assertEquals(v(i), dense(i)))
    }
  }

  property("fromDense round-trips and stores no zeros; compact preserves the dense value") {
    forAll(sparseAny) { v =>
      val dense = v.toDense
      val back = SparseVector.fromDense(dense)
      assertCanonical(back)
      assertEquals(back.toDense.toSeq, dense.toSeq)
      assert(back.activeValues.forall(_ != 0.0))
      assertEquals(back.activeSize, dense.toSeq.count(_ != 0.0))
      val compacted = v.compact
      assertEquals(compacted.toDense.toSeq, dense.toSeq)
      assertEquals(compacted.activeIndices, back.activeIndices)
    }
  }

  property("dot, +, - and axpyInto agree with dense arithmetic for every pattern relation") {
    forAll(sparsePair, exactValue) { case ((a, b), alpha) =>
      val da = a.toDense
      val db = b.toDense
      assertEquals(a.dot(b), da.dot(db))
      assertEquals(a.dot(db), da.dot(db))
      val sum = a + b
      val diff = a - b
      assertCanonical(sum)
      assertCanonical(diff)
      assertEquals(sum.toDense.toSeq, (da + db).toSeq)
      assertEquals(diff.toDense.toSeq, (da - db).toSeq)
      val union = (a.activeIndices ++ b.activeIndices).distinct.sorted
      assertEquals(sum.activeIndices, union)
      assertEquals(diff.activeIndices, union)
      val y = db.mutableCopy
      a.axpyInto(alpha, y)
      val expected = db.mutableCopy
      expected.axpyInPlace(alpha, da)
      assertEquals(y.toVec.toSeq, expected.toVec.toSeq)
    }
  }

  property("scaling and mapActive keep the pattern and touch active entries only") {
    forAll(sparseAny, exactValue) { (v, alpha) =>
      val scaled = v * alpha
      assertEquals(scaled.activeIndices, v.activeIndices)
      assertEquals(scaled.toDense.toSeq, (v.toDense * alpha).toSeq)
      val shifted = v.mapActive(_ + 1.0)
      assertEquals(shifted.activeIndices, v.activeIndices)
      val active = v.activeIndices.toSet
      (0 until v.length).foreach { i =>
        assertEquals(shifted(i), if active(i) then v(i) + 1.0 else 0.0)
      }
    }
  }

  property("reductions and norms equal the dense reductions, implicit zeros included") {
    forAll(sparseAny) { v =>
      val dense = v.toDense.toSeq
      assertEquals(v.sum, dense.sum)
      assertEquals(v.norm1, dense.map(math.abs).sum)
      assertEquals(v.normInf, if dense.isEmpty then 0.0 else dense.map(math.abs).max)
      assertEqualsDouble(v.norm2, math.sqrt(dense.map(x => x * x).sum), 1e-12)
      if v.length > 0 then
        assertEquals(v.max, dense.max)
        assertEquals(v.min, dense.min)
    }
  }

  property("general doubles: dot and norm2 agree with dense within rounding") {
    val gen =
      for
        n <- Gen.choose(1, 64)
        a <- sparseGen(n, Gen.choose(-1e3, 1e3))
        b <- sparseGen(n, Gen.choose(-1e3, 1e3))
      yield (a, b)
    forAll(gen) { case (a, b) =>
      val scale = a.norm2 * b.norm2 + 1.0
      assertEqualsDouble(a.dot(b), a.toDense.dot(b.toDense), 1e-12 * scale)
      assertEqualsDouble(a.norm2, a.toDense.norm2, 1e-12 * (a.norm2 + 1.0))
    }
  }

  property("CSR rows, CSC columns and CSR * SparseVector agree with dense") {
    val gen =
      for
        rows <- Gen.choose(1, 8)
        cols <- Gen.choose(1, 8)
        entries <- Gen.listOf(Gen.zip(Gen.choose(0, rows - 1), Gen.choose(0, cols - 1), exactValue))
        x <- sparseGen(cols)
      yield (rows, cols, entries, x)
    forAll(gen) { case (rows, cols, entries, x) =>
      val builder = Sparse.coo(rows, cols)
      entries.foreach { case (r, c, v) => builder.add(r, c, v) }
      val csr = builder.toCSR()
      val csc = builder.toCSC()
      (0 until rows).foreach { r =>
        val row = csr.rowSparse(r)
        assertCanonical(row)
        assertEquals(row.length, cols)
        assertEquals(row.toDense.toSeq, csr.row(r).toSeq)
      }
      (0 until cols).foreach { c =>
        val col = csc.colSparse(c)
        assertCanonical(col)
        assertEquals(col.length, rows)
        assertEquals(col.toDense.toSeq, csc.col(c).toSeq)
      }
      val y = MutableDVec.zeros(rows)
      csr.mulInto(x.toDense, y)
      assertEquals((csr * x).toSeq, y.toVec.toSeq)
    }
  }
