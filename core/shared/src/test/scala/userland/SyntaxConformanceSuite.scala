package userland

import gale.linalg.*

/** Conformance for the opt-in ergonomic syntax against the PRD's worked examples
  * (PRD §"Elementwise operations must be explicit" and §"Symbol aliases … opt-in
  * syntax modules", lines 614–629): the documented forms compile and compute, the
  * elementwise product is the Hadamard product (not the core matrix product), the
  * Unicode aliases are exact synonyms of the core ops, and a shape mismatch throws
  * `DimensionMismatch` (matching the core arithmetic primitives, not the
  * Either-returning solve tier).
  *
  * Namespace note: the PRD sketches `import gale.linalg.syntax.unicode.*`; gale
  * exposes these modules at the sibling top-level `gale.syntax.{all, unicode}`
  * (a deliberate, equivalent choice). The semantics below are the contract.
  */
class SyntaxConformanceSuite extends munit.FunSuite:

  import gale.syntax.all.*

  private val a = Matrix.dense(2, 2)(1.0, 2.0, 3.0, 4.0)
  private val b = Matrix.dense(2, 2)(5.0, 6.0, 7.0, 8.0)

  test("PRD elementwise examples compile, compute, and are NOT the matrix product") {
    // A.pointwise * B — Hadamard product.
    val had = a.pointwise * b
    assertEqualsDouble(had(0, 0), 5.0, 1e-12) // 1·5
    assertEqualsDouble(had(0, 1), 12.0, 1e-12) // 2·6
    assertEqualsDouble(had(1, 0), 21.0, 1e-12) // 3·7
    assertEqualsDouble(had(1, 1), 32.0, 1e-12) // 4·8
    // Distinct from the core matrix product A * B (whose (0,0) is 1·5 + 2·7 = 19).
    assertEqualsDouble((a * b)(0, 0), 19.0, 1e-12)
    assert(had(0, 0) != (a * b)(0, 0), "pointwise must not be the matrix product")

    // A.pointwise / B
    val quot = a.pointwise / b
    assertEqualsDouble(quot(0, 0), 1.0 / 5.0, 1e-15)
    assertEqualsDouble(quot(1, 1), 4.0 / 8.0, 1e-15)

    // A.pointwise.map(math.exp)
    val mapped = a.pointwise.map(math.exp)
    assertEqualsDouble(mapped(1, 1), math.exp(4.0), 1e-12)

    // A.zipMap(B)(_ + _) — agrees with the core elementwise +.
    val summed = a.zipMap(b)(_ + _)
    val coreSum = a + b
    assertEqualsDouble(summed(0, 0), coreSum(0, 0), 1e-12)
    assertEqualsDouble(summed(1, 1), coreSum(1, 1), 1e-12)
  }

  test("elementwise shape mismatch throws DimensionMismatch (arithmetic-primitive convention)") {
    val wide = Matrix.dense(2, 3)(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    intercept[LinAlgError.DimensionMismatch](a.pointwise * wide)
    intercept[LinAlgError.DimensionMismatch](a.pointwise / wide)
    intercept[LinAlgError.DimensionMismatch](a.zipMap(wide)(_ + _))
  }

  test("PRD Unicode aliases (× and ⋅) are exact synonyms of the core ops") {
    import gale.syntax.unicode.*
    val x = Vec(1.0, 1.0)
    val y = Vec(2.0, 3.0)

    // val y = A × x  (matrix-vector)
    val mv = a × x
    val mvRef = a * x
    assertEquals(mv.length, mvRef.length)
    var i = 0
    while i < mv.length do
      assertEqualsDouble(mv(i), mvRef(i), 1e-15)
      i += 1

    // A × B  (matrix-matrix)
    val mm = a × b
    assertEqualsDouble(mm(0, 0), (a * b)(0, 0), 1e-12)
    assertEqualsDouble(mm(1, 1), (a * b)(1, 1), 1e-12)

    // val d = x ⋅ y  (inner product)
    assertEqualsDouble(x ⋅ y, x.dot(y), 1e-15)
    assertEqualsDouble(x ⋅ y, 5.0, 1e-15) // 1·2 + 1·3
  }

  test("the wildcard `all` import stays free of the Unicode operators") {
    // `all` deliberately excludes ×/⋅ (they are a separate opt-in), so referencing
    // them without importing `unicode` must not compile.
    assert(compileErrors("a × b").nonEmpty, "× must require the unicode import")
    assert(compileErrors("Vec(1.0) ⋅ Vec(2.0)").nonEmpty, "⋅ must require the unicode import")
  }
