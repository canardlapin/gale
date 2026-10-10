package gale.bench

import gale.bench.BreezeBenchData.*
import gale.bench.BreezeBenchmarkValidation.*
import gale.linalg.LinAlgError

class BreezeBenchmarkValidationSuite extends munit.FunSuite:
  test("a failed operation cannot become a successful timing fixture") {
    val e = intercept[IllegalStateException](successful("solve", Left(LinAlgError.SingularMatrix(1))))
    assert(e.getMessage.contains("solve failed"))
  }

  test("reject a zero eigenvector and the wrong conjugate sign") {
    val a = Array(Array(1.0, -2.0), Array(2.0, 1.0))
    val audit = new Audit(2)
    audit.eigenpair("positive", a, 1.0, 2.0, Array(1.0, 0.0), Array(0.0, -1.0))
    audit.eigenpair("negative", a, 1.0, -2.0, Array(1.0, 0.0), Array(0.0, 1.0))
    intercept[IllegalArgumentException](audit.eigenpair("zero", a, 1.0, 2.0, Array(0.0, 0.0), Array(0.0, 0.0)))
    intercept[IllegalArgumentException](audit.eigenpair("wrong conjugate", a, 1.0, -2.0, Array(1.0, 0.0), Array(0.0, -1.0)))
  }

  test("reject a corrupted R even when it remains triangular and finite") {
    val a = Array(Array(3.0, 0.0), Array(0.0, 4.0))
    val audit = new Audit(2)
    audit.qr("valid", a, a)
    intercept[IllegalArgumentException](audit.qr("corrupt", a, Array(Array(3.0, 1.0), Array(0.0, 4.0))))
  }

  test("reject a non-triangular square root and an incomplete eigenbasis") {
    val eye = Array(Array(1.0, 0.0), Array(0.0, 1.0))
    val audit = new Audit(3)
    audit.cholesky("valid", eye, eye)
    // This rotation has a positive diagonal and L L^T = I, but is not Cholesky.
    intercept[IllegalArgumentException](audit.cholesky("rotation", eye, Array(Array(0.6, -0.8), Array(0.8, 0.6))))
    val a = Array(Array(1.0, 0.0, 0.0), Array(0.0, 2.0, 0.0), Array(0.0, 0.0, 3.0))
    val basis = Array(Array(1.0, 0.0, 0.0), Array(0.0, 1.0, 0.0), Array(0.0, 0.0, 1.0))
    audit.eigenbasis("valid", a, Array(1.0, 2.0, 3.0), new Array[Double](3), basis)
    // Every individual pair (1,e1) is valid; three copies do not span A.
    val duplicated = Array(Array(1.0, 1.0, 1.0), Array(0.0, 0.0, 0.0), Array(0.0, 0.0, 0.0))
    intercept[IllegalArgumentException](audit.eigenbasis("duplicated", a, Array(1.0, 1.0, 1.0), new Array[Double](3), duplicated))
  }

  test("reject a corrupted SVD and an inaccurate or non-finite solution") {
    val a = Array(Array(4.0, 0.0), Array(0.0, 3.0))
    val eye = Array(Array(1.0, 0.0), Array(0.0, 1.0))
    val audit = new Audit(2)
    audit.svd("valid", a, Array(4.0, 3.0), eye, eye)
    intercept[IllegalArgumentException](audit.svd("corrupt", a, Array(5.0, 3.0), eye, eye))
    audit.solve("valid", a, column(Array(1.0, 2.0)), column(Array(4.0, 6.0)))
    intercept[IllegalArgumentException](audit.solve("wrong", a, column(Array(1.0, 1.0)), column(Array(4.0, 6.0))))
    intercept[IllegalArgumentException](audit.solve("NaN", a, column(Array(Double.NaN, 2.0)), column(Array(4.0, 6.0))))
  }

  test("determinant overflow needs an independent sign and magnitude witness") {
    val audit = new Audit(2)
    val a = Array(Array(1e200, 0.0), Array(0.0, 1e200))
    audit.determinant("overflow", a, Double.PositiveInfinity)
    intercept[IllegalArgumentException](audit.determinant("wrong sign", a, Double.NegativeInfinity))
    intercept[IllegalArgumentException](audit.determinant("NaN", a, Double.NaN))
    intercept[IllegalArgumentException](audit.determinant("spurious overflow", Array(Array(2.0)), Double.PositiveInfinity))
  }

  test("paired small fixtures validate through the actual public operations") {
    val a = diagonallyDominant(8, 2100L)
    val s = spd(8, 200L)
    val rhs = column(vectorData(8, 300L))
    for op <- Seq("Inv", "Lu", "Qr", "Svd", "Pinv", "Det", "Solve", "LuSolve") do
      assert(validate(op, a, rhs).nonEmpty, op)
    for op <- Seq("Chol", "CholSolve", "EigSym") do assert(validate(op, s, rhs).nonEmpty, op)
    assert(validate("Lstsq", matrixData(24, 8, 400L), column(vectorData(24, 500L))).nonEmpty)
    assert(validate("Gemm", a, matrixData(8, 3, 2400L)).nonEmpty)
    assert(validate("Eig", Array(Array(1.0, -2.0), Array(2.0, 1.0))).nonEmpty)
  }

  test("corrected full-rank pinv fixture passes at n=64 where old general fixture was deficient") {
    val old = breeze.linalg.svd(breezeMatrix(matrixData(64, 64, 2200L))).singularValues
    assert(old(63) < 64.0 * math.ulp(1.0) * old(0), "old fixture must not be described as numerically full-rank")
    assert(validate("Pinv", diagonallyDominant(64, 2100L)).nonEmpty)
  }

  test("selected Vector backend validates without changing the timed API") {
    val vector = gale.backend.jvm.vector.VectorBackend
    assert(validate("Inv", diagonallyDominant(64, 2100L), backend = vector).nonEmpty)
    assert(validate("LuSolve", diagonallyDominant(64, 1700L), matrixData(64, 16, 1900L), vector).nonEmpty)
  }
