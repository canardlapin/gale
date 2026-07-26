# Breeze parity tests

The parity module runs Gale and Breeze 2.1.0 on the same numerical inputs. It
compares results that have the same mathematical meaning. For factorizations,
where signs or bases may differ, the tests compare reconstructions, residuals,
or subspace projectors instead of comparing raw factors.

Run the suite from the repository root:

```sh
sbt parityTest
```

`EverydayOpsParitySuite` uses ScalaCheck to vary matrix shapes and data seeds.
The factorization and spectral suites use fixed adversarial and
well-conditioned fixtures. A parity test should state the shared mathematical
contract and use a tolerance that accounts for the algorithms being compared.

## Migration pain points

This table records cases found while writing parity tests where the Breeze
expression is materially easier to write or read. It is an API-design input,
not a claim that Gale should copy Breeze.

| Operation | Breeze | Gale | Why the Gale expression is harder | Status |
| --- | --- | --- | --- | --- |
| Matrix literal | `DenseMatrix((1.0, 2.0), (3.0, 4.0))` | `Matrix(2, 2)(1.0, 2.0, 3.0, 4.0)` | Gale requires the dimensions and a flattened row-major value list. The compiler cannot check that visual rows have equal lengths because the rows are not present in the expression. | Consider a row-based constructor that keeps `Matrix(rows, cols)` for generated or flattened input. |
| Contiguous matrix slice | `a(1 until 4, 2 until 5)` | `a.slice(1, 4, 2, 5)` | Two range values show which endpoints belong together. Gale's four adjacent integers are easier to transpose or misread. | Consider an overload that accepts row and column ranges. |
| Exception-style solve | `a \ b` | `a.solve(b).orThrow` | Gale returns `Either[LinAlgError, A]`, which preserves the failure type but adds ceremony in programs that deliberately use exceptions. | Keep the typed default. Evaluate a clearly named throwing convenience only if migrations repeatedly add local wrappers. |

Add a row only when the shorter Breeze form also makes the operation clearer.
Do not list deliberate differences merely because their spelling differs.
