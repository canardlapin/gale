# Gale 0.1 compatibility and release policy

Gale is preparing its first immutable ecosystem checkpoint on the `0.1` line.
Untagged source derives a unique `0.1.0+...-SNAPSHOT` version. A snapshot names
a development state; it is not a published release or a compatibility promise.
Cross-repository validation must identify snapshots by their full Git commit.

The intended progression is:

1. `0.1.0+...-SNAPSHOT` for exact-commit development and consumer convergence;
2. `v0.1.0-M1` for the first immutable ecosystem checkpoint;
3. later `v0.1.0-Mn` or `v0.1.0-RCn` tags when another candidate is needed; and
4. `v0.1.0` only after the milestone passes its downstream soak gate.

Only an exact `v0.1.x`, `v0.1.x-Mn`, or `v0.1.x-RCn` tag may produce a
publishable version. Any source, build, dependency, or executable-documentation
change creates a new candidate and invalidates evidence tied to the old commit.

## Artifact admission

The organization is `io.github.canardlapin`; Scala 3 artifacts use the `_3`
binary suffix. The first milestone admits four coordinates:

| Artifact | Platform | Compatibility status |
| --- | --- | --- |
| `gale-core_3` | JVM | stable from `v0.1.0-M1` |
| `gale-core_sjs1_3` | Scala.js | stable from `v0.1.0-M1` |
| `gale-laws_3` | JVM | stable from `v0.1.0-M1` |
| `gale-laws_sjs1_3` | Scala.js | stable from `v0.1.0-M1` |

`gale-core` contains the public dense, sparse, solver, optimization, spectral,
sized, kernel, error, diagnostic, and backend-contract packages. Gale does not
publish empty package-shaped modules. `gale-laws` is a normal module because
downstream libraries extend its MUnit and ScalaCheck conformance suites.

The following modules are tested but provisional and are excluded from the M1
bundle: Breeze interop, the JDK Vector backend, native matrix storage, and the
FFM BLAS/LAPACK backend. Exclusion means that their tests remain part of the
candidate court while their coordinates carry no 0.1 publication or
compatibility promise. Admission requires a later explicit manifest change.

`gale-interop-ravel` is excluded because it currently depends on
`ravel-core:1.0.0-SNAPSHOT`. Parity projects, benchmarks, documentation, demos,
consumer probes, and the root aggregator are also not published. The complete
coordinate and exclusion table is the [release manifest](release-manifest.md).

## Compatibility promise

`v0.1.0-M1` establishes the compatibility baseline for the four admitted
artifacts. After that tag:

- `0.1.x` patch releases preserve the public binary and source contract and may
  contain compatible additions or fixes;
- a removal, incompatible signature change, or documented semantic-contract
  break requires `0.2.0`;
- private, `private[gale]`, internal, and explicitly experimental definitions
  remain outside the promise; and
- provisional modules acquire no promise until an admission record names their
  first baseline.

Public compatibility includes typed failure cases, result ownership,
convergence diagnostics, deterministic ordering, work-accounting fields, and
backend capability behavior where the documentation makes those observable.
It does not promise bit-identical floating-point results across legal backends.
The package-by-package inventory and semantic review boundary are recorded in
[API stability boundary](api-stability.md).

There is no compatibility baseline before M1 because no immutable artifact
exists. The build already wires MiMa/version-policy and TASTy-MiMa across the
four admitted artifacts. Before M1, `compatibilityCheck` exercises that wiring
with an empty previous-version set; it does not prove binary or TASTy
compatibility. The M1 work must therefore retain a manually reviewed API export
receipt and must not describe the candidate as binary-verified.

Immediately after publication, 0.1-line development must run
`compatibilityCheck` with
`-Dgale.compatibility.baseline=0.1.0-M1`. CI must make the compatibility job
blocking before accepting another change. A later immutable baseline may
replace M1 only through an explicit release-policy change; a moving snapshot is
never a valid baseline.

Scala 3 TASTy compatibility also constrains consumers. Gale publishes from
Scala 3.7.4 for the Scala 3 binary line. CI compiles a Scala 3.8.4 consumer
against locally published 3.7.4 artifacts as advisory evidence; that probe does
not replace compatibility checks between Gale releases.

## Required candidate evidence

One clean, pushed commit must pass:

```text
formatCheck
releaseDependencyCheck
testAllFull
compileAll
parityTest
interopBreezeTest
benchCompile
docsCheck
```

The same commit must also run `compatibilityCheck`. Before M1 this verifies the
task wiring only; after M1 it compares against the exact immutable baseline.
`coverageJVM` is advisory evidence: it runs the complete core JVM suite under
instrumentation and retains HTML, XML, and Cobertura reports. Gale does not use
a repository-wide percentage threshold as a release gate. Maintainers inspect
the report, close meaningful contract gaps, record the resulting rate, and run
`coreJVM/clean` before subsequent uninstrumented compilation or tests.

CI must also pass the maintained Vector JDK 21 and 22 lanes, the JDK 22
FFM/OpenBLAS lane, optimized Scala.js linking, and published-consumer probes.
The exact run URL and every job result belong in the release receipt. Passing a
subset of the court is not release proof.

Branch protection must require these job names:

- `Formatting and fatal-warning policy`;
- `Tests (Scala 3.7.4, JVM)` and `Tests (Scala 3.7.4, JS)`;
- `Scaladoc and executable guides`;
- `Scala.js optimized links`;
- `Breeze parity and published interop`;
- `Release dependency manifest`;
- `Benchmarks compile`;
- `Vector backend (JDK 21)` and `Vector backend (JDK 22)`; and
- `FFM BLAS/LAPACK backend (JDK 22 / OpenBLAS)`.

`JVM coverage report (advisory)`, `MiMa and TASTy-MiMa baseline (pre-M1
advisory)`, `Scala 3.8.4 consumer probe`, `Scala Next source experiment
(advisory)`, and `WebAssembly tests and profile (experimental)` remain advisory
before M1. The compatibility job loses its advisory qualifier and becomes
required immediately after the baseline is published. A red advisory lane must
be recorded and qualified; it must not be silently presented as supported or
omitted from the receipt.

The Breeze consumer probe publishes provisional Gale modules only to an
isolated local repository and compiles a separate coordinate-based consumer.
It verifies POM and transitive-dependency behavior without admitting Breeze
interop to the Central bundle.

## Ecosystem convergence

Before M1, every direct consumer must run its relevant gates against one exact
Gale commit. The receipt must record the Gale commit, consumer commit, commands,
results, source migrations, and generated dependency metadata. A consumer may
use an explicit source or local-artifact override during convergence; its
permanent dependency moves only after its gate passes.

The required direct consumers are Alder, grakern, graph4s, image4s, multivar,
reframe4s, regress4s, ScalaFIM, and signal4s. Provider publication precedes any
consumer release that needs the new coordinate.

## Supported runtimes

| Module or route | Required runtime | Status |
| --- | --- | --- |
| `gale-core` JVM and `gale-laws` JVM | JDK 21 | admitted and required |
| `gale-core` and laws for Scala.js | Node 22 in CI | admitted and required |
| Vector backend | JDK 21 or 22 plus `jdk.incubator.vector` | provisional, required test lane |
| Native and FFM BLAS/LAPACK backends | JDK 22 plus explicit native access | provisional, required test lane |
| Scala.js JavaScript in browsers | current evergreen browsers intended | not browser-certified until browser CI exists |
| Scala.js Wasm | Node 22+ experimental profile | advisory and excluded |

Release notes must distinguish the Node-tested Scala.js runtime from intended
browser support. They must not claim browser certification without a browser
gate.

## Publication and closure

The release workflow performs two different operations:

- A manual dry run validates a synthetic `0.1.x-Mn` or `0.1.x-RCn` version,
  builds the four admitted artifacts without credentials, and verifies an
  unsigned bundle without uploading it.
- A tag run derives its version from the exact Git tag, imports release secrets
  only inside publication jobs, validates the complete signed bundle, and
  uploads one `USER_MANAGED` Central Portal deployment. Publication remains a
  deliberate Portal action after validation.

Both publication aliases begin with `releaseM1Preflight`. It rejects a
snapshot-classified candidate, a non-0.1 release version, prohibited snapshot
dependencies, or any admitted project whose `publishTo` destination is not
sbt's local `target/sona-staging` repository. Bundle verification also checks
the exact Maven group and artifact set plus the MD5 and SHA-1 content digests;
the signed tag path additionally requires every primary artifact signature.

Before closing the milestone, the owner must verify the Central namespace,
project signing key, Portal credentials, and one non-publishing validation
deployment. A clean external consumer must then resolve and exercise every
admitted coordinate without local repository or sibling-checkout assistance.

The release receipt must record the exact commit and tag, toolchains, local
commands and counts, remote CI URL, direct-consumer commits and results,
generated-POM inspection, bundle manifest and signatures, external-consumer
proof, benchmark receipts, advisory failures, and owner-controlled external
state. `v0.1.0` remains a later decision after the M1 soak; this policy does not
authorize publishing it as part of the M1 epic.
