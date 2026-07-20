# Gale 1.0.0-SNAPSHOT release evidence

Evidence date: 2026-07-17. Candidate identity: the commit containing this record
on `main`; its parent is `393cc60`. All commands below ran against the candidate
working tree before that evidence commit was created.

This record certifies local engineering gates. It does not claim that a public
`1.0.0` was published or that remote CI ran on the evidence commit. Apache-2.0,
the canonical GitHub repository, and POM provenance metadata were subsequently
owner-selected on 2026-07-19.

## Environment

- macOS 14.3, Apple ARM64
- sbt 1.11.7, Scala 3.3.8
- OpenJDK 22 for JVM, JS, Vector, FFM, packaging, and Scaladoc gates
- Node.js 24.1.0 for the executed Wasm profile
- runtime-discovered macOS Accelerate for local FFM tests/measurements

Local JDK 21 was unavailable. The workflow defines required Vector jobs on JDK
21 and 22 and an OpenBLAS FFM job on JDK 22, but checked-in workflow text is not
a substitute for a remote run on this candidate.

## Correctness and build gates

```bash
JAVA_HOME=<jdk22> sbt -java-home <jdk22> \
  testAllFull parityTest interopBreezeTest vectorBackendTest benchCompile

JAVA_HOME=<jdk22> sbt -java-home <jdk22> \
  nativeBackendTest blasFfmBackendTest benchFfmCompile

GALE_WASM=1 sbt coreJS/test benchSmokeJSFull
```

All commands exited zero.

| Gate | Result |
| --- | ---: |
| core JVM | 376 / 376 |
| core Scala.js JavaScript | 373 / 373 |
| laws JVM | 30 / 30 |
| laws Scala.js | 30 / 30 |
| Breeze differential parity | 31 / 31 |
| Breeze conversions/migration | 24 / 24 |
| Vector backend + conformance | 27 / 27 |
| native storage | 3 / 3 |
| FFM BLAS/LAPACK + conformance | 19 / 19 |
| Wasm-linked core | 373 / 373 |
| core/laws Scala.js full links | pass |
| JVM/JS benchmark compilation | pass |
| JDK 22 FFM JMH compilation | pass |

The [acceptance audit](v1-acceptance-audit.md) records the static dependency,
package-residue, example, and public-storage checks. It found and removed the
pre-v1 unsafe raw-array adoption API; copy-only JVM and Scala.js interop passed
all gates above.

## Artifact evidence

Binary JARs and POMs were generated successfully for the intended v1 artifact
set at version `1.0.0-SNAPSHOT`:

- `gale-core_3` and `gale-core_sjs1_3`
- `gale-laws_3` and `gale-laws_sjs1_3`
- `gale-interop-breeze_3`
- `gale-backend-jvm-vector_3`
- `gale-backend-jvm-native_3`
- `gale-backend-jvm-blas-ffm_3`

Source and Scaladoc artifacts were built for all eight. A first documentation
pass exposed unresolved Scaladoc links; those links were repaired, and focused
core/laws/interop documentation rebuilds completed without warnings. Package
inspection found no benchmark, parity, or vendor checkout content in published
core/backend/laws artifacts. The interop artifact contains only its intentional
`gale.interop.breeze` classes and declared Breeze dependency.

The generated core POM has no forbidden compile dependency. The build now adds
Apache-2.0, homepage, SCM, and developer metadata to every published module.

## Performance evidence

The canonical summary is the [backend dashboard](../benchmarks/dashboard.md).
The candidate includes full receipts for:

- Vector GEMM/GEMV versus pure Gale and Breeze VectorBLAS;
- FFM heap-copy GEMM, GEMV, LU/Cholesky/QR/eigen crossovers;
- solve/least-squares scenarios, allocation rate, and declined-route overhead;
- optimized Scala.js JavaScript versus executed Wasm.

The important shipped decisions are conservative:

- Vector ARM64/JDK 22: contiguous GEMM at `128^3`, GEMV at `128^2`;
- Accelerate/JDK 22: square GEMM at `512^3`, LU at n=128;
- FFM GEMV, QR, and Cholesky remain default-off;
- OpenBLAS/MKL automatic thresholds remain off until family-specific sweeps;
- Wasm remains experimental/default-off after two runs showed optimized
  JavaScript 20.6–42.9x faster.

## Release disposition

Engineering acceptance: **pass** for this local candidate.

Public release disposition: **not releasable yet**. Remaining actions are:

1. configure the binary publishing destination, credentials, and signing policy;
2. run required remote CI on the exact release commit, including JDK 21 Vector
   and Linux/OpenBLAS JDK 22 FFM jobs;
3. replace `1.0.0-SNAPSHOT` with `1.0.0`, build/sign/publish the same artifact
   set in the documented JDK 21/JDK 22 passes, and tag that exact commit.

No release tag, signing, or binary publication was performed by this audit.
