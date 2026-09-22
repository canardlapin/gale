# Gale release-readiness audit — 22 September 2026

## Verdict

**Do not publish the current candidate yet: Gale has substantial, credible numerical assurance and working artifact packaging, but this audit reproduced false convergence in preconditioned GMRES and found unresolved release-process gaps.** A public development snapshot is a reasonable next target after the blocking fixes; an immutable M1 requires the additional compatibility and downstream qualification already promised by the project.

This is an assessment, not publication approval. No production source, CI settings, credentials, tags, or public artifacts were changed. The report and its evidence are the only additions to the user's checkout.

## Candidate and scope

| Item | Audited state |
| --- | --- |
| Main candidate | `1b018d4dd082b786123b65e02cd5b81e61968e21`, including merged PR #8 |
| User's checkout | Clean `release/0.1-stabilization` at `4485cc775ae8233789b019d24a920f86391e9523` before this report; preserved |
| Main development version | `0.1.0+106-1b018d4d-SNAPSHOT` |
| Scala / build | Scala 3.7.4, sbt 1.11.7, Scala.js 1.22.0 |
| Admitted publication slice | `io.github.canardlapin:gale-core_3`, `gale-core_sjs1_3`, `gale-laws_3`, `gale-laws_sjs1_3` |
| Supported/tested core routes | JVM, Scala.js on Node; browser execution was not qualified |
| Provisional routes | Breeze interop, Vector, native storage, FFM BLAS/LAPACK; tested but excluded from the M1 bundle |
| Excluded integration | Ravel, whose current dependency is a snapshot |
| Compatibility promise | M1 establishes the four artifacts' binary, source, and documented semantic baseline |

The audit inspected build/publication definitions, live GitHub configuration and CI, laws, selected numerical and ownership implementations/tests, independent reference fixtures, benchmark receipts, documentation, and downstream pins. It also built Maven artifacts, consumed them from a separate project, challenged the bundle verifier, and ran an independent numerical counterexample. This is a risk-directed audit, not a proof of every algorithm or a security certification.

Source links below are pinned to the audited main commit because the user's release branch differs materially. [Artifact manifest][manifest]; [compatibility policy][policy].

## What has actually passed

- **All 14 jobs passed on the exact main commit** in [CI run 35728678680](https://github.com/canardlapin/gale/actions/runs/35728678680). This includes JVM/JS tests, docs, optimized links, Breeze parity and the interop consumer, dependency checks, benchmark compilation, Vector on JDK 21/22, and native/FFM tests on JDK 22 with OpenBLAS. Scala 3.8.4 and Wasm jobs also passed, but policy classifies them as advisory.
- **1,574 existing tests passed locally** on the identical source tree during PR integration: 681 core JVM, 671 core JS, 54 laws JVM, 54 laws JS, 85 parity, and 29 Breeze interop. The command also compiled the modules/benchmarks, linked optimized JS, built executable docs, and checked the published interop consumer. These are test executions across platforms, not 1,574 distinct mathematical guarantees.
- **All four admitted artifacts staged successfully** with synthetic version `0.1.0-M1`: binaries, sources, Scaladoc, POMs, and checksums. The unsigned bundle passed the existing verifier. This was local staging, with no tag, signing, upload, or public release.
- **A separate consumer resolved the staged Maven artifacts without Gale source dependencies or Ivy-local.** Depending on `gale-laws` pulled in the matching `gale-core`; published conformance laws and analytic banded/sparse Cholesky smoke tests passed on JVM, fast-linked JS, and fully optimized JS: 10 registered tests per execution. Native-only law bodies are conditional and are not exercised by this pure-backend consumer.
- POM inspection found matching core/laws versions and no snapshot dependencies. Test libraries have test scope in core and compile scope in the deliberately reusable laws module. All inspected class files had major version 52; this is packaging information, **not** a claim that the supported runtime can be lowered below JDK 21.

Local checks used JDK 22 and Node 26.7.0. Hosted CI supplies the separate JDK 21 / Node 22 evidence. [Machine-readable CI receipt](release-audit-2026-09-22/main-ci.json), [test summary](release-audit-2026-09-22/full-test-summary.txt), [artifact metadata](release-audit-2026-09-22/artifact-metadata.json), [consumer log](release-audit-2026-09-22/consumer-corrected.log).

## Prioritized findings

### F1 — Block both snapshot and M1: GMRES can report false convergence

**Confirmed on both JVM and Scala.js, using the staged public artifacts.** This uses finite inputs and a perfectly conditioned one-dimensional system; it is not an ill-conditioning edge case.

```scala
import gale.linalg.*
import gale.solvers.*

val a = Matrix.dense(1, 1)(1e12)
val b = Vec(1.0)
val result = gmres(
  a, b, SolverConfig(tolerance = 1e-10),
  Preconditioner.Jacobi(a)
)
val actualResidual = (b - a * result.x).norm2
```

| Observation | JVM and JS result |
| --- | --- |
| `result.converged` | `true` |
| `result.iterations` | `0` |
| `result.x(0)` | `0.0` |
| Reported residual | `1e-12` |
| Actual `\|\|b - Ax\|\|` | `1.0` |
| Required absolute tolerance | `1e-10` |
| Exact solution | `1e-12` |

`RelativeToRhs` fails identically because `||b|| = 1`. Four assertions failed: the two tolerance modes on each platform.

The public configuration documents the original-system residual, but GMRES measures `||M⁻¹(b-Ax)||` and compares it to a threshold based on the unpreconditioned right-hand side. The initial return already permits the false success; the restart and inner-iteration returns need the same review. [Documented contract][solver-contract]; [GMRES implementation][gmres].

The existing reusable `SolverLaws.solvesSystem` independently recomputes the residual, but its concrete law suite exercises CG on one small SPD system. GMRES's preconditioning test uses one family and checks solution accuracy; it does not exercise scalar changes in preconditioner magnitude. That explains how extensive green tests missed this case. [Solver laws][solver-laws]; [solver law suite][solver-law-suite]; [preconditioning tests][preconditioning-tests].

**Required closure:** verify the true residual before every successful GMRES return; keep absolute/relative normalization consistent with the public contract; add shared JVM/JS regression and scaling tests. Preserve bounded termination and honest work counts if extra operator applications are needed. A law that merely checks `result.residual` would not close this defect.

Evidence: [JVM failure](release-audit-2026-09-22/gmres-jvm.log), [JS failure](release-audit-2026-09-22/gmres-js.log), [reproduction harness](release-audit-2026-09-22/consumer/README.md). The deliberately failing probes are outside the normal test source tree.

### F2 — Resolve the main/release-branch split before choosing a release candidate

`main` and `release/0.1-stabilization` differ across **244 paths**. Main contains the newer numerical policy, portable sparse Cholesky, packed banded Cholesky, and PR #8 work. The release branch contains substantial QA work absent from main:

- `-Wunused:all`, `-Wvalue-discard`, and `-Wnonunit-statement`, with a documented mdoc exception;
- formatting across maintained sources and a `scalafmtSbtCheck` gate;
- scoverage reporting;
- version-policy and TASTy-MiMa wiring for a future real baseline;
- a documented `is-terminal` dependency fix for fresh-cache formatting on JDK 22;
- a bounded release performance court and associated receipts.

Main currently enables fatal deprecation/feature/unchecked warnings, but its formatter includes only four named Scala files plus an `.sbt` glob. Green formatting therefore does not mean repository-wide formatting. The earlier local JDK 22 aggregate formatter stalled; hosted JDK 21 formatting passed. The release-branch terminal fix is a credible lead for that local problem, not a reproduced diagnosis of the stalled process. [Main compiler settings][build-version]; [main formatting scope][format]; [release-branch QA configuration][release-build]; [release-branch plugin fix][release-plugins].

**Required decision:** converge one reviewed candidate from current main, selectively recovering the useful release work while preserving newer algorithms and tests. Do not publish the old release branch merely because its name sounds authoritative. Historic coverage/performance results from that branch do not certify the integrated candidate.

### F3 — Public snapshot publication is not operational yet

The build already selects the correct Sonatype snapshot destination for its normal Git-derived snapshot version. However, the checked-in release workflow handles semantic tags or unsigned milestone/RC rehearsals; it does not provide a guarded four-artifact snapshot workflow. The M1 aliases intentionally reject snapshots. Publication policy and README currently describe source-only development and tag-only publication. [Version settings][build-version]; [publication aliases][aliases]; [release workflow][release-workflow].

Live checks using the repository-owner account found:

- no GitHub releases or tags;
- no discoverable Maven metadata for any of the four coordinates in Maven Central **or** the Central snapshot repository: all eight requests returned HTTP 404;
- zero repository Actions secrets; the workflow's named PGP and Sonatype credentials are therefore not configured there;
- namespace ownership, snapshot enablement, and owner-held credentials outside GitHub remain **unverified**.

**Required closure for a public snapshot:** admit the same four coordinates explicitly, document the snapshot policy, configure owner-controlled Portal credentials and namespace snapshot access, and publish only a tested clean commit. Avoid an aggregate root `publish`: provisional modules are publishable for local probes and are excluded by the M1 aliases, not universally disabled. After upload, an anonymous external consumer must resolve and execute all four coordinates from the public resolver with local Gale resolution disabled. Record the full Git SHA, exact version, artifact hashes, and consumer receipt.

[Live configuration receipt](release-audit-2026-09-22/remote-state.json); [anonymous resolver checks](release-audit-2026-09-22/public-resolvers.json).

### F4 — Release gates are documented but not enforced at the release boundary

The owner-authenticated API reports **“Branch not protected”**, and the effective main-branch rules list is empty. This conflicts with the policy's explicit required-check list. [Required checks][policy-ci].

The tag workflow reruns `testAllFull docsCheck`, but does not itself require all the documented candidate jobs on that exact SHA. For example, parity/interop, native/Vector conformance, and formatting are separate CI jobs. The workflow also lacks an explicit Node 22 setup step despite running JS tests. Its upload step records a user-managed Portal deployment ID but does not wait for Portal validation. Upload success is not validation success or public availability. [Release preflight][release-workflow]; [upload stage][release-upload].

**Required closure:** enforce the policy's non-advisory checks, require them on the selected SHA before publication, pin the release test runtime, retain the manifest/check receipts, and verify remote validation and anonymous consumption separately. Keep the manual final-publication control for immutable releases.

### F5 — Bundle verification has narrower guarantees than its name suggests

The actual generated bundle was internally consistent and consumed successfully. The verifier checks expected directory names, version directories, required files, MD5/SHA-1 values, and optionally the existence of `.asc` files. It does not parse POM coordinates/dependencies or authenticate signatures. [Verifier][bundle-verifier].

Two controlled negative probes both returned exit code 0 and “Central bundle verified”:

1. Changing a laws POM's own and internal dependency versions to `9.9.9-SNAPSHOT`, keeping its `0.1.0-M1` directory/name, and recomputing its checksums.
2. Supplying text saying `not a PGP signature` in every required `.asc` file with `--require-signatures`.

These probes demonstrate a local verification gap, not a bypass of Sonatype's release validation. **Before M1**, parse and compare POM group/artifact/version and internal dependency versions against the manifest, reject prohibited dependencies, and verify signatures against the intended signing identity. Snapshot publication also needs explicit metadata checks because its service does not perform release-style validation. [Probe results](release-audit-2026-09-22/bundle-negative-probes.json).

### F6 — A source archive is incorrectly classified as a non-snapshot

Confirmed by loading an archive of the same commit without `.git`:

```text
version    = 0.1.0-SNAPSHOT-20260922-0856
isSnapshot = false
publishTo  = local-staging
```

The fallback puts the timestamp after `SNAPSHOT`, while classification requires the version to end in `-SNAPSHOT`. [Fallback and classification][build-version]. The guarded M1 version check still rejects this version, so this is not evidence of accidental remote release. It is a real development/publication classification defect.

**Fix before enabling snapshot automation:** keep the fallback suffix consistent and test clean Git, dirty Git, tagged Git, and no-Git archives. [Probe log](release-audit-2026-09-22/archive-version.log).

### F7 — Close the documented M1 consumer and compatibility obligations

The current workspace registry names ten direct Gale consumers, including the newer mixeff4s dependency. The release plan still enumerates nine. Current inspected builds differ from the registry and from each other: ScalaFIM pins `18d24db…`, mixeff4s pins release-branch `f869613…`, and multivar still derives a local artifact version from `83cac90…`. These are development arrangements, not evidence that they consume this candidate. [Registry and inspected build receipt](release-audit-2026-09-22/downstream-state.json).

The separate artifact consumer built during this audit is useful packaging proof. It does not replace the promised same-candidate downstream court. **Before M1**, update the actual consumer inventory, compile/test the required consumers against the selected candidate, review public API changes against their last shared pins, and retain that evidence. Set M1 as the automated binary/TASTy baseline once its immutable artifacts exist. There is no previous public baseline against which a meaningful MiMa pass can be claimed today. [API inventory and migration receipt][api-stability]; [compatibility promise][policy-compat].

A public development snapshot can precede the full M1 court if its scope and lack of a compatibility promise are explicit. It should still pass the artifact-consumer and numerical gates above.

### F8 — Documentation and assurance coverage need reconciliation

Before handing an install line to external users:

- Correct README's claim that sparse direct factorization is unimplemented; the explicit pure sparse Cholesky provider now ships in core. Release policy also still lists full dense SVD and sparse-direct factorization as deferred. [README boundary][readme-boundary]; [stale deferred list][policy-deferred].
- Reconcile numerical-contract wording that simultaneously says no sparse Cholesky provider exists and describes the pure provider. Resolve the five retained mdoc broken-link warnings and the two unresolved Scaladoc member links. Executable examples pass, but the documentation gate permits these warnings.
- Provide a verified install version/resolver and reachable documentation URL. GitHub Pages configuration exists, but `https://canardlapin.github.io/gale/` returned 404 during this audit.
- Execute the full optimized Scala.js test bundles in CI, not only `fullLinkJS`. This audit executed a small artifact consumer under full optimization; that does not certify the complete library's optimized tests. [Current alias][aliases].
- Recover coverage reporting and inspect the selected candidate's uncovered numerical/error branches. No current-main coverage or mutation result was obtained. Prior release-branch percentages must not be transplanted. Targeted mutation of convergence predicates would be valuable; a repository-wide mutation threshold is unnecessary.

Browser execution and broad hardware/performance qualification remain outside today's evidence. Either add the evidence needed for stronger claims or retain the existing Node-tested, backend-specific qualifications.

## Assurance scorecard

Ratings concern the audited main commit and actual evidence, not what another branch plans to enable.

| Dimension | Rating | Evidence | Gap or rationale |
| --- | --- | --- | --- |
| ScalaCheck use and generator quality | Present but incomplete | Seeded bounded dense/spectral suites, single workers, zero dimensions, strides, and invariant-preserving adversarial shrinkers. [Generators][generators] | Strong foundations; missing preconditioner-scaling domain exposed F1. Small spectral case budgets are not exhaustive. |
| Reusable law-test module | Strong | JVM/JS public laws and capability conformance; separate staged-artifact consumer successfully extends the suite. [Published suite][backend-laws] | Demonstrably reusable; numerical completeness is assessed separately. |
| Test framework and Discipline integration | Strong | MUnit, ScalaCheck bridge, executable custom numerical laws on both platforms. [Law tests][dense-laws] | Discipline is not required for this API; no adoption merely for branding. |
| Typeclass lawfulness and coherence | Not applicable | Core exposes concrete numerical operations rather than Cats-style Eq/Order/Monad instances. | Runtime `given Backend` selection and composition are assessed under provider conformance. |
| Backend/provider conformance | Present but incomplete | Independent decimal dense-kernel oracles, stride/alpha/beta cases, typed factor errors, capability validation, native CI. [Backend laws][backend-laws] | Coverage is especially strong for dense backends; no equivalent externally reusable conformance kit was established for every sparse/spectral provider capability. |
| Cross-platform/version CI | Present but incomplete | 14 exact-commit jobs passed; JDK 21/22, JS/Node, advisory next Scala/Wasm. | Required-check enforcement absent; full optimized test execution and browsers not comprehensively qualified. |
| Numerical/computational assurance | Present but incomplete | Analytic plants, reconstruction, residuals, clustered eigenspaces, scaling/congruence, banded storage limits. [Adversarial spectral tests][spectral-adversarial] | Confirmed false success in GMRES prevents a Strong rating. |
| Differential/independent oracles | Strong | Breeze differential court, NumPy/SciPy stored references with generator, BigDecimal kernels, analytic Laplacian. [Reference generation][numpy]; [banded tests][banded-tests] | Independence is substantial, but no claim covers every operation or scale. |
| Failure/convergence/resource contracts | Present but incomplete | Typed failures, partial spectral results, close/reuse tests, builder transfer, immutable snapshots. | F1 violates the documented success contract; non-finite solver/configuration behavior deserves a focused follow-up. |
| Work/allocation accounting | Strong | Exact A/B/preconditioner counters; matrix-free no-basis-sweep check; 100,000-row packed banded solve; retained allocation measurements. [Work tests][spectral-work]; [banded receipt][banded-receipt] | Strong for named claims; storage counts are not peak RSS, and benchmark receipts are host-specific. |
| Compiler discipline | Present but incomplete | Fatal warnings compile in exact-commit CI. | Stronger unused/discard/non-unit settings exist only on the release branch. |
| Formatting/semantic rewrites | Present but incomplete | Deterministic Scalafmt gate passes remotely. | Main has a narrow include list; no need to add Scalafix without a named rewrite problem. |
| Binary/source compatibility | Present but incomplete | Written M1 contract, API inventory, next-Scala and artifact consumer probes. | No immutable baseline yet; recover enforcement wiring and finish the export/downstream review before M1. |
| Coverage/mutation signal | Missing | No mainline coverage plugin/reporting gate or mutation result established. | Existing release-branch coverage work can be recovered; coverage supplements independent correctness checks. |
| Benchmark/performance evidence | Present but incomplete | JMH compiles, raw benchmark receipts retained, declined acceleration routes documented. [Dashboard][dashboard]; [banded receipt][banded-receipt] | No fresh integrated-candidate performance court ran here; older measurements do not certify universal speed or all completed-work equivalence. |
| Documentation/release evidence | Present but incomplete | Executable guides, API docs, manifest, metadata, local staging and artifact consumption. | Contradictions, broken links, inaccessible hosted site, absent public resolution and incomplete release enforcement. |

## Strongest numerical evidence worth preserving

The quality of the existing evidence is more important than the test total:

- Dense backend laws use independently accumulated `BigDecimal` references, rather than only comparing two implementations that could share a defect. They cover padded/transposed storage, non-unit strides, and nontrivial alpha/beta.
- Generalized eigen tests compare metric projectors for repeated/clustered roots, check ill-conditioned SPD metrics, distinguish returned residual-passing subsets from requested counts, reject non-finite geometry, and count A/B/preconditioner applications independently. These are unusually useful defenses against plausible-looking wrong answers. [Spectral adversarial tests][spectral-adversarial].
- Banded Cholesky has a closed-form discrete-Laplacian factor, determinant, inverse and solution; scaling/conditioning tests; mutation/ownership failures; and a 100,000-row narrow-band case that would make accidental dense expansion prohibitive. [Banded tests][banded-tests].
- Sparse Cholesky checks ordering, symbolic reuse, factor independence from analysis lifetime, transpose/multiple RHS, typed refusal, fill guards, and close behavior. It is explicit opt-in; the default provider remains unavailable. [Sparse provider tests][sparse-tests].
- Performance records distinguish work/storage claims from measured allocations and timings, retain slow/default-off accelerator cases, and qualify results by host/runtime. Keep those limits intact when preparing release prose.

## Smallest credible route to public availability

### First: a public development snapshot

1. Fix F1 with independent residual/scaling tests on JVM and JS. Fix F6 and reconcile the selected main/release-branch QA work into one candidate.
2. Run the full candidate court on that SHA, including parity, interop, backend checks and optimized execution. Enforce the required checks and rerun the four-artifact consumer. Retain failures as failures.
3. Make a narrow snapshot publication path for **only core and laws, JVM and JS**. Verify namespace ownership and snapshot enablement, configure credentials, and validate POMs/dependencies before upload. Reuse the existing supported publishing plugin and its verified destination.
4. Publish the exact version derived from the tested commit. Verify anonymous resolution and execution from the public endpoint with no source override or Ivy-local fallback. Only then replace the README placeholder with the actual version, resolver, and tested install command.

Sonatype's snapshot service uses `https://central.sonatype.com/repository/maven-snapshots/` and requires snapshot enablement for the namespace. Snapshots are overwritable, currently cleaned up after 90 days, and do not receive the final-release validation process. Consumers need an explicit resolver. Consequently, a public snapshot is useful for evaluation but is not a durable scientific-reproducibility checkpoint. [Sonatype snapshot documentation](https://central.sonatype.org/publish/publish-portal-snapshots/).

### Then: the immutable `0.1.0-M1` checkpoint

Complete the same-candidate downstream/API review, signing and full bundle verification, current benchmark qualification required by the release plan, a user-managed Portal validation rehearsal, and public artifact-consumer proof. Capture the baseline and enable compatibility checks against it. M1 is not just another download name: current policy makes it the start of the compatibility promise.

This sequence avoids demanding an unnecessary ecosystem rewrite. Keep MUnit/ScalaCheck and the custom laws. Do not add Cats or Cats Effect to a pure numerical core for appearances. Discipline is optional if a concrete downstream RuleSet need arises. Existing `sbt-ci-release`/site tooling is sufficient; broader `sbt-typelevel` adoption is not a publication prerequisite. Recover scoverage and compatibility tools for their specific purposes, and use narrowly targeted mutation experiments around convergence/error predicates.

## Verification boundary and evidence

The reproducible numerical failure changes the release verdict even though the existing court is green. It has **not** been fixed during this assessment.

| Check | Actual result / qualification |
| --- | --- |
| `testAllFull compileAll parityTest interopBreezeTest benchCompile docsCheck publishedInteropProbe` | Passed earlier in this session on the tree identical to audited main; 1,574 tests |
| Exact-main hosted CI | 14 jobs passed; experimental jobs remain advisory |
| Synthetic `releaseM1Unsigned` | Passed locally on JDK 22; no signing/upload |
| Unsigned bundle verification | Passed for the real staged bundle |
| Separate Maven consumer | 10 tests passed on each of JVM, fast JS, full-opt JS |
| GMRES independent contract probes | **2 failures on JVM + 2 on JS**, as described in F1 |
| Deliberately inconsistent POM / fake signatures | Both incorrectly accepted by current local verifier |
| No-Git archive version classification | Incorrect `isSnapshot = false` reproduced |
| Public metadata / documentation endpoint | Eight artifact metadata 404s; Pages URL 404 |
| Local aggregate formatter | Stalled earlier; stopped owned processes; not counted as a local pass |
| Current-candidate signing / Portal validation / publication | Not attempted |
| Full consumer ecosystem, fresh coverage/mutation/performance, browser tests | Not executed |

[Evidence directory](release-audit-2026-09-22/) contains compact logs, API receipts, artifact metadata, and the isolated consumer/reproducer. Full integration logs remain at `/private/tmp/gale-pr8-fix-tests.log`; the retained summary records their SHA-256 and the run metadata. Temporary audit work is under `/private/tmp/gale-release-audit-20260922`; staging used `/private/tmp/gale-pr8-review-20260922`. No active audit build remains.

[manifest]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/release-manifest.md#L9-L45
[policy]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/release-policy.md
[policy-compat]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/release-policy.md#L56-L80
[policy-ci]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/release-policy.md#L86-L106
[policy-deferred]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/release-policy.md#L143-L145
[api-stability]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/api-stability.md#L50-L76
[build-version]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/build.sbt#L35-L82
[aliases]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/build.sbt#L504-L547
[format]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/.scalafmt.conf#L1-L11
[release-build]: https://github.com/canardlapin/gale/blob/4485cc775ae8233789b019d24a920f86391e9523/build.sbt#L79-L121
[release-plugins]: https://github.com/canardlapin/gale/blob/4485cc775ae8233789b019d24a920f86391e9523/project/plugins.sbt#L7-L14
[release-workflow]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/.github/workflows/release.yml#L1-L95
[release-upload]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/.github/workflows/release.yml#L97-L136
[bundle-verifier]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/tools/verify-central-bundle.sh#L94-L137
[solver-contract]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/main/scala/gale/solvers/Solvers.scala#L5-L23
[gmres]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/main/scala/gale/solvers/Solvers.scala#L458-L582
[solver-laws]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/laws/shared/src/main/scala/gale/laws/SolverLaws.scala#L10-L31
[solver-law-suite]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/laws/shared/src/test/scala/gale/laws/SolverLawSuite.scala#L9-L20
[preconditioning-tests]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/test/scala/gale/solvers/PreconditionedSolverSuite.scala#L23-L52
[generators]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/laws/shared/src/test/scala/gale/laws/AdversarialGeneratorSuite.scala#L12-L140
[dense-laws]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/laws/shared/src/test/scala/gale/laws/DenseLawSuite.scala#L12-L159
[backend-laws]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/laws/shared/src/main/scala/gale/laws/BackendConformanceSuite.scala#L8-L356
[spectral-adversarial]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/test/scala/gale/spectral/GeneralizedOperatorAdversarialSuite.scala#L94-L252
[spectral-work]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/test/scala/gale/spectral/GeneralizedOperatorAdversarialSuite.scala#L254-L350
[banded-tests]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/test/scala/gale/linalg/BandedCholeskySuite.scala#L17-L177
[banded-receipt]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/docs/banded-cholesky-evidence.md#L39-L132
[sparse-tests]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/core/shared/src/test/scala/gale/sparse/direct/pure/PureSparseCholeskySuite.scala#L42-L218
[numpy]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/parity/scripts/generate_numpy_references.py#L1-L70
[dashboard]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/benchmarks/dashboard.md#L1-L67
[readme-boundary]: https://github.com/canardlapin/gale/blob/1b018d4dd082b786123b65e02cd5b81e61968e21/README.md#L94-L101
