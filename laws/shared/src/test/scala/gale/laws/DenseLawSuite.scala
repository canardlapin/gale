package gale.laws

import gale.linalg.*
import munit.ScalaCheckSuite
import org.scalacheck.Gen
import org.scalacheck.Prop.forAll

/** Property-based algebraic laws for the dense types, driven through the reusable
  * [[VecLaws]] / [[MatrixLaws]] bundles. Lives in the laws module (not core)
  * because it depends only on the public API and on the bundles it exercises.
  */
class DenseLawSuite extends ScalaCheckSuite:
  private val scalarGen: Gen[Double] =
    Gen.chooseNum(-10.0, 10.0)

  private val vec3Gen: Gen[DVec] =
    Gen.listOfN(3, scalarGen).map(DVec.fromSeq)

  private val mat2x3Gen: Gen[DMat] =
    Gen.listOfN(6, scalarGen).map(values => Matrix.dense(2, 3, values))

  private val bigScalarGen: Gen[Double] =
    Gen.chooseNum(-1e3, 1e3)

  private val dimGen: Gen[Int] =
    Gen.choose(0, 12)

  private val smallDimGen: Gen[Int] =
    Gen.choose(0, 4)

  private def dvecGen(n: Int): Gen[DVec] =
    Gen.listOfN(n, bigScalarGen).map(DVec.fromSeq)

  private def dmatGen(rows: Int, cols: Int): Gen[DMat] =
    Gen.listOfN(rows * cols, bigScalarGen).map(values => Matrix.dense(rows, cols, values))

  private val shapedMatGen: Gen[DMat] =
    for
      m <- dimGen
      n <- dimGen
      a <- dmatGen(m, n)
    yield a

  private val matVecGen: Gen[(DMat, DVec)] =
    for
      m <- dimGen
      n <- dimGen
      a <- dmatGen(m, n)
      x <- dvecGen(n)
    yield (a, x)

  private val matVecTransposeGen: Gen[(DMat, DVec)] =
    for
      m <- dimGen
      n <- dimGen
      a <- dmatGen(m, n)
      y <- dvecGen(m)
    yield (a, y)

  private val scaleAddGen: Gen[(Double, DVec, DVec)] =
    for
      n <- dimGen
      a <- bigScalarGen
      x <- dvecGen(n)
      y <- dvecGen(n)
    yield (a, x, y)

  private val tripleGen: Gen[(DMat, DMat, DMat)] =
    for
      p <- smallDimGen
      q <- smallDimGen
      r <- smallDimGen
      s <- smallDimGen
      a <- dmatGen(p, q)
      b <- dmatGen(q, r)
      c <- dmatGen(r, s)
    yield (a, b, c)

  property("DVec addition is commutative within floating tolerance") {
    forAll(vec3Gen, vec3Gen) { (x: DVec, y: DVec) =>
      VecLaws.additionCommutes(x, y)
    }
  }

  property("DVec scalar multiplication distributes over vector addition") {
    forAll(vec3Gen, vec3Gen, scalarGen) { (x: DVec, y: DVec, alpha: Double) =>
      VecLaws.scalarDistributesOverAddition(alpha, x, y)
    }
  }

  property("DMat matrix-vector multiplication is linear in the vector") {
    forAll(mat2x3Gen, vec3Gen, vec3Gen) { (A: DMat, x: DVec, y: DVec) =>
      MatrixLaws.matVecIsLinear(A, x, y)
    }
  }

  property("DMat transpose is an involution and views match element access") {
    forAll(mat2x3Gen) { (A: DMat) =>
      MatrixLaws.transposeIsInvolution(A)
      VecLaws.assertClose(A.row(1), Vec(A(1, 0), A(1, 1), A(1, 2)))
      VecLaws.assertClose(A.col(2), Vec(A(0, 2), A(1, 2)))
    }
  }

  property("(A.t).t equals A elementwise for random shapes (0..12)") {
    forAll(shapedMatGen) { (a: DMat) =>
      MatrixLaws.transposeIsInvolution(a)
    }
  }

  property("matvec through a strided source vector matches the contiguous result") {
    forAll(matVecGen) { (pair: (DMat, DVec)) =>
      val (a, x) = pair
      VecLaws.assertCloseRel(a * strided(x), a * x)
    }
  }

  property("matvec through a transposed view matches an explicit transpose") {
    forAll(matVecTransposeGen) { (pair: (DMat, DVec)) =>
      val (a, y) = pair
      VecLaws.assertCloseRel(a.t * y, explicitTranspose(a) * y)
    }
  }

  property("scalar multiplication distributes over vector addition for large values") {
    forAll(scaleAddGen) { (triple: (Double, DVec, DVec)) =>
      val (a, x, y) = triple
      VecLaws.scalarDistributesOverAddition(a, x, y)
    }
  }

  property("matrix multiplication is associative for small shapes") {
    forAll(tripleGen) { (triple: (DMat, DMat, DMat)) =>
      val (a, b, c) = triple
      MatrixLaws.multiplicationAssociates(a, b, c)
    }
  }

  /** A stride-2 view of `x`'s values over a fresh backing buffer, built through
    * the public API: column 0 of an `n x 2` row-major matrix has row-stride 2.
    */
  private def strided(x: DVec): DVec =
    val values = new Array[Double](x.length * 2)
    var i = 0
    while i < x.length do
      values(i * 2) = x(i)
      i += 1
    Matrix.dense(x.length, 2, values.toIndexedSeq).col(0)

  private def explicitTranspose(a: DMat): DMat =
    Matrix.tabulate(a.cols, a.rows)((i, j) => a(j, i))
