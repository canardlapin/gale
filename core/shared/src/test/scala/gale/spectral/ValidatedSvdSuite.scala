package gale.spectral

import gale.linalg.{DMat, Matrix}
import munit.FunSuite

class ValidatedSvdSuite extends FunSuite:
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
