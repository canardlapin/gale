# Gale v1 acceptance audit

Audit date: 2026-07-17. Audited build: `1.0.0-SNAPSHOT` on the current local
tree. This is an engineering acceptance record, not evidence of a public
release or remote CI run.

## Acceptance result

| Criterion | Evidence | Result |
| --- | --- | --- |
| Dense, sparse, solver, spectral, kernel, and laws code cross-builds | JVM/JS tests plus full Scala.js links; packages consolidated in `gale-core` per artifact policy | pass |
| Public worked examples compile and run | `gale.examples.WorkedExamplesSuite`: 5 tests on JVM, JS, and executed Wasm suite | pass |
| Core has no forbidden numerical/effect/native/plot dependencies | generated core POM contains Scala library plus test-scoped MUnit only; source scan finds no Breeze, Spire, Cats Effect, ZIO, FFM, AWT/Swing, or plotting imports | pass |
| No public platform-storage contract | unsafe JVM/JS array-adoption helpers removed; `fromArrayCopy` / `fromFloat64ArrayCopy` and exporters are copy-only; internal storage remains `private[gale]` | pass |
| No Atlas/Zephyr/placeholder package residue | scan finds only the PRD acceptance sentence and this audit description; no declaration, import, example, or live reference | pass |
| Primitive `Vec[Double]` path | dense/type/law suites verify `Vec(...)` yields `DVec`; generic fallback is documented as correctness-oriented | pass |
| Dense solve and factorization invariants | core, laws, Breeze differential, Vector conformance, and FFM LAPACK reconstruction/typed-error suites | pass |
| Sparse v1 behavior and diagnostics | sparse/core laws, Matrix Market laws, iterative solver suites, and explicit support statement | pass |
| Optional acceleration is explicit and thresholded | no-import pure default; Vector/FFM import points; measured dashboard and declined-route control | pass |
| JS storage is typed-array-backed and private | Scala.js implementation uses private `Float64Array` abstraction; copy-only public interop | pass |
| Tiny kernels and allocation doctrine | tiny-kernel/performance-doctrine suites plus JMH/GC receipts | pass |

## Exact local gates

Run with OpenJDK 22 unless noted:

```bash
sbt testAllFull parityTest interopBreezeTest vectorBackendTest benchCompile
sbt nativeBackendTest blasFfmBackendTest benchFfmCompile
GALE_WASM=1 sbt coreJS/test benchSmokeJSFull
```

Observed test totals, all with zero failures/errors:

| Suite | Tests |
| --- | ---: |
| core JVM | 376 |
| core Scala.js | 373 |
| laws JVM / Scala.js | 30 / 30 |
| Breeze parity | 31 |
| Breeze interop/migration | 24 |
| Vector backend | 27 |
| native storage | 3 |
| FFM BLAS/LAPACK | 19 |
| Wasm-linked core | 373 |

Both JVM JMH and Scala.js benchmark projects compiled; the JDK 22 FFM JMH
project compiled separately. Core and laws test bundles also completed
`fullLinkJS`.

## Static audit commands

The audit excluded generated `target/`, the unreferenced `vendor/breeze/`
reference checkout, and tracker operation files:

```bash
rg -n -i 'atlas\.linalg|zephyr\.linalg|package (atlas|zephyr)' \
  --glob '!vendor/**' --glob '!**/target/**' --glob '!.mote/**' .
rg -n -i 'breeze|spire|cats.effect|zio|java.lang.foreign|java.awt|javax.swing|plot' \
  core/{shared,jvm,js}/src/main
rg -n 'from(Array|Float64Array)Unsafe' \
  --glob '!vendor/**' --glob '!**/target/**' --glob '!.mote/**' .
```

The only `Array[Double]`/`Float64Array` types on the public core surface are
copy-in/copy-out interop values; none aliases Gale-owned storage. Unsafe adopt
functions and raw backing access remain `private[gale]` implementation details.

## Scope and unresolved publication items

Local JDK 21 was unavailable. The checked-in workflow has required Vector legs
for JDK 21 and 22, but this audit does not present unexecuted remote CI as green.
The local JDK 22 Vector gate passed.

The code acceptance criteria pass. Public publication remains blocked by the
owner decisions listed in [the release policy](release-policy.md): project
license, SCM/homepage/developer POM metadata, publishing destination, and an
actual remote CI run on the release commit.
