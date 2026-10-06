package gale.spectral

import gale.linalg.{DMat, Matrix}
import munit.FunSuite

class ValidatedSvdSuite extends FunSuite:
  private def vector(values: Double*): gale.linalg.DVec = gale.linalg.DVec.tabulate(values.size)(values(_))
  private def candidate(a: DMat): SVD =
    Svds.svd(a, SingularSelection.All).fold(error => fail(error.toString), identity)

  test("enclosure contains known rectangular singular values") {
    val a = Matrix.tabulate(3, 2)((r, c) => if r == c then if r == 0 then 5.0 else 2.0 else 0.0)
    val result = ValidatedSvd.enclosure(MatrixEnclosure.exact(a).toOption.get, candidate(a)).toOption.get
    assert(result.lower(0) <= 5.0 && result.upper(0) >= 5.0)
    assert(result.lower(1) <= 2.0 && result.upper(1) >= 2.0)
    assert(result.factorError >= 0.0 && result.leftDefect < 1.0 && result.rightDefect < 1.0)
  }

  test("entry interval uncertainty expands a valid singular enclosure") {
    val a = Matrix.dense(1, 1)(3.0)
    val lower = Matrix.dense(1, 1)(2.9)
    val upper = Matrix.dense(1, 1)(3.1)
    val result = ValidatedSvd.enclosure(MatrixEnclosure.checked(lower, upper).toOption.get, candidate(a)).toOption.get
    assert(result.lower(0) <= 2.9 && result.upper(0) >= 3.1)
  }

  test("rotated tied and zero spectra remain enclosed") {
    // Columns are orthogonal with norm five; the last rectangular root is zero.
    val a = Matrix.dense(3, 3)(3.0, 4.0, 0.0, 4.0, -3.0, 0.0, 0.0, 0.0, 0.0)
    val result = ValidatedSvd.enclosure(MatrixEnclosure.exact(a).toOption.get, candidate(a)).toOption.get
    assert(result.lower(0) <= 5.0 && result.upper(0) >= 5.0)
    assert(result.lower(1) <= 5.0 && result.upper(1) >= 5.0)
    assert(result.lower(2) == 0.0 && result.upper(2) >= 0.0)
  }

  test("bad factors and nonfinite enclosure endpoints refuse certification") {
    val a = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0)
    val s = candidate(a)
    val malformed = s.copy(u = DMat.zeros(2, 2))
    assert(ValidatedSvd.enclosure(MatrixEnclosure.exact(a).toOption.get, malformed).isLeft)
    assert(MatrixEnclosure.checked(Matrix.dense(1, 1)(Double.NaN), Matrix.dense(1, 1)(1.0)).isLeft)
  }

  test("ill-conditioned and scale-extreme finite candidates either enclose or refuse") {
    val a = Matrix.dense(2, 2)(1e150, 0.0, 0.0, 1e-150)
    ValidatedSvd.enclosure(MatrixEnclosure.exact(a).toOption.get, candidate(a)) match
      case Right(result) => assert(result.lower(0) <= 1e150 && result.upper(0) >= 1e150)
      case Left(_) => () // refusal is the defined conservative outcome on an unsafe scale.
  }

  test("subnormal-scale and interval-overflow courts refuse rather than narrow") {
    val tiny = Matrix.dense(1, 1)(java.lang.Double.MIN_VALUE)
    ValidatedSvd.enclosure(MatrixEnclosure.exact(tiny).toOption.get, candidate(tiny)) match
      case Right(result) => assert(result.lower(0) <= java.lang.Double.MIN_VALUE && result.upper(0) >= java.lang.Double.MIN_VALUE)
      case Left(_) => ()
    val huge = Matrix.dense(1, 1)(Double.MaxValue)
    val invalid = MatrixEnclosure.checked(huge, huge).toOption.get
    assert(ValidatedSvd.enclosure(invalid, candidate(huge)).isLeft)
  }

  test("candidate factors are checked independently of solver metadata and estimates") {
    val a = Matrix.dense(2, 2)(5.0, 0.0, 0.0, 2.0)
    val base = candidate(a)
    // These nonorthogonal factors reconstruct a exactly: diag(1.25,1) * diag(4,2).
    val scaled = base.copy(
      u = Matrix.dense(2, 2)(1.25, 0.0, 0.0, 1.0),
      vt = Matrix.dense(2, 2)(1.0, 0.0, 0.0, 1.0),
      singularValues = vector(4.0, 2.0),
      rank = 0
    )
    Vector(scaled, scaled.copy(singularValues = vector(3.5, 2.25))).foreach { factors =>
      val result = ValidatedSvd.enclosure(MatrixEnclosure.exact(a).toOption.get, factors)
        .fold(error => fail(error.toString), identity)
      assert(result.leftDefect > 0.5 && result.leftDefect < 1.0)
      Vector(5.0, 2.0).zipWithIndex.foreach { (root, i) =>
        assert(result.lower(i) <= root && root <= result.upper(i))
      }
    }
  }

  test("wide matrices and nondegenerate matrix boxes enclose known spectra") {
    val wide = Matrix.dense(2, 3)(5.0, 0.0, 0.0, 0.0, 2.0, 0.0)
    val wideBounds = ValidatedSvd.enclosure(MatrixEnclosure.exact(wide).toOption.get, candidate(wide))
      .fold(error => fail(error.toString), identity)
    Vector(5.0, 2.0).zipWithIndex.foreach { (root, i) =>
      assert(wideBounds.lower(i) <= root && root <= wideBounds.upper(i))
    }
    val center = Matrix.dense(2, 2)(3.0, 0.0, 0.0, 1.0)
    val box = MatrixEnclosure.checked(
      Matrix.dense(2, 2)(2.9, 0.0, 0.0, 0.9),
      Matrix.dense(2, 2)(3.1, 0.0, 0.0, 1.1)
    ).toOption.get
    val bounds = ValidatedSvd.enclosure(box, candidate(center)).fold(error => fail(error.toString), identity)
    for first <- Vector(2.9, 3.0, 3.1); second <- Vector(0.9, 1.0, 1.1) do
      assert(bounds.lower(0) <= first && first <= bounds.upper(0))
      assert(bounds.lower(1) <= second && second <= bounds.upper(1))
  }

  test("malformed enclosures and invalid candidate spectra or factors refuse") {
    val a = Matrix.dense(2, 2)(3.0, 0.0, 0.0, 1.0)
    val input = MatrixEnclosure.exact(a).toOption.get
    val base = candidate(a)
    assert(MatrixEnclosure.checked(a, DMat.zeros(1, 1)).isLeft)
    assert(MatrixEnclosure.checked(a, DMat.zeros(2, 2)).isLeft)
    assert(MatrixEnclosure.checked(null, a).isLeft)
    assert(ValidatedSvd.enclosure(null, base).isLeft)
    assert(ValidatedSvd.enclosure(input, null).isLeft)
    val invalid = Vector(
      base.copy(u = DMat.zeros(1, 2)),
      base.copy(vt = DMat.zeros(2, 1)),
      base.copy(singularValues = vector(-1.0, -2.0)),
      base.copy(singularValues = vector(1.0, 3.0)),
      base.copy(singularValues = vector(Double.NaN, 1.0)),
      base.copy(u = Matrix.dense(2, 2)(Double.NaN, 0.0, 0.0, 1.0)),
      base.copy(vt = Matrix.dense(2, 2)(1.0, 0.0, 0.0, Double.PositiveInfinity)),
      base.copy(u = Matrix.dense(2, 2)(1.0, 1.0, 0.0, 0.0))
    )
    invalid.foreach(factors => assert(ValidatedSvd.enclosure(input, factors).isLeft))
  }
