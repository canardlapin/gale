# Gale release readiness — 26 October 2026

**Verdict: do not release this commit. Gale has a credible foundation for a first
0.1 milestone, but a false-convergence defect, failing required CI, and incomplete
release qualification block publication. No broad rewrite or feature expansion is needed.**

Assessed source: `f747d2aabd5ff21851982355f2093ac8b8f91ec5`, confirmed as both local
HEAD and GitHub `main`. Existing untracked `.agent-work/` and Python cache files
were preserved. This assessment adds this report and scratch evidence; it does not
repair product code, change GitHub configuration, tag, or publish.

The applicable target is the repository's `0.1.0-M1`: core and laws, JVM and
Scala.js, four Maven coordinates under `io.github.canardlapin`. A 1.0 compatibility
freeze is premature. Provisional JVM backends and Breeze interop need not be
admitted to make a useful first release.

## Findings in release order

### R1 — GMRES can report convergence without solving the system (blocker)

`SolverConfig` documents the linear residual `||b-Ax||` for GMRES
([Solvers.scala](../core/shared/src/main/scala/gale/solvers/Solvers.scala), lines 7–12).
GMRES instead compares `||M^-1(b-Ax)||` against that tolerance at initialization,
restart, and inner convergence (lines 472–488, 547, 578–579). Relative mode scales
the threshold by the original `||b||`, so it does not repair the inconsistency.

A decisive public-API example is `A=[1e12]`, `b=[1]`, Jacobi preconditioning,
`tolerance=1e-10`, and the default zero initial point. The preconditioned residual
is `1e-12`, while the true residual is `1`. The solver accepts `x=0` at iteration
zero. The exact solution is `1e-12`. This is a correctness defect, not a marginal
roundoff disagreement or an overly strict test.

The scratch `ReleaseContractProbe.scala` independently recomputes the residual
for both Absolute and RelativeToRhs modes. Its assertions intentionally witness
the defect; a passing witness is evidence of failure, not solver qualification.
The witness reproduced this result on both JVM and Scala.js in both modes.
See the execution receipt below. An earlier September audit had identified this
same case; current source inspection confirms that the faulty path remains.

**Acceptance:** certify the original-system residual at every successful exit,
report it consistently, and retain the preconditioned residual only as an internal
iteration estimate or separately named diagnostic. Add regressions for initialization,
restart, terminal checks, scaled/degenerate preconditioners, and both tolerance modes
on JVM and JS. Account for the extra operator applications honestly.

### R2 — Required hosted tests are red on the exact candidate (blocker)

[CI run 39128478041](https://github.com/canardlapin/gale/actions/runs/39128478041)
was freshly inspected; its source SHA matches this checkout. The JVM job executed
1,278 tests: 1,274 passed and four failed:

- QR/least-squares digest: 175 mismatching cases.
- L-BFGS family digest: 93 mismatching cases.
- Nonsymmetric eigensolver digest: 47 mismatching cases.
- Hessenberg digest: 49 mismatching cases.

The corresponding tests are in `core/jvm/src/test/scala/gale/{linalg,optim,spectral}`:
`QRDigestSuite.scala:45`, `LbfgsDigestSuite.scala:49`, and
`NonsymmetricKernelDigestSuite.scala:39,50`. The optimized-link job also fails
because its `testAllFull` alias begins with the failing JVM suite. This is not
proof of a Scala.js linker defect.

The fresh local full suite passes on macOS ARM64 with sbt-selected Homebrew
Java 25.0.1 and Node 26.7.0. Hosted CI uses Ubuntu/JDK 21 and Node 22. That
contrast narrows the investigation; it does not prove the digest differences are
harmless. Some corpora use `pow` and transcendental objectives, and their comments
acknowledge platform-sensitive bits. Root cause and numerical impact remain open.

**Acceptance:** explain the first differing inputs/outputs on the failing runtime,
check residuals/reconstructions/statuses with independent invariants, and then fix
the implementation or make the digest's platform contract explicit. Do not merely
regenerate fixtures or suppress failures. Require all mandatory jobs green at the
final candidate SHA.

### R3 — Release policy is not enforced by GitHub configuration (release gate)

The authenticated branch-protection endpoint reports `Branch not protected`;
the effective branch-rules endpoint returns `[]`. This contradicts
[release-policy.md](release-policy.md), lines 82 onward, which requires named jobs.
The release workflow runs dependency, core/laws, linking and docs checks, but does
not itself verify successful completion of all required parity, backend, formatting,
and benchmark jobs at the tagged SHA. CI's ordinary push trigger is branch-based.

**Acceptance:** configure the intended required checks and ensure the tag upload
path verifies the exact SHA's complete required check set (or runs that set through
a shared workflow). Branch protection alone is not a guarantee that an arbitrary
tag points at a qualified commit.

### R4 — A signed-bundle success message does not establish valid artifacts (release assurance)

[verify-central-bundle.sh](../tools/verify-central-bundle.sh), lines 94–134,
checks filenames, directory versions, checksums, and signature-file existence.
It does not parse POM contents, open JARs, or authenticate signatures.

The audit constructed four correctly named coordinates with matching checksums,
invalid POM/JAR bytes, and text saying `not a signature` in each `.asc` file.
Running the verifier with `--require-signatures` exited zero and printed:

```text
Central bundle verified: 4 admitted artifacts, version 0.1.0-M1
```

This does not imply that Sonatype would publish the invalid bundle. Portal
validation is a separate defense. It does mean the local release claim is too strong.
[Central's requirements](https://central.sonatype.org/publish/requirements/) cover
actual POM metadata, sources/documentation, and PGP signatures.

**Acceptance:** validate POM coordinates, dependency versions and required metadata;
check archive integrity; verify signatures against the intended release key; and
add rejection fixtures. Alternatively name a deliberately narrower structural check
honestly, but retain real signature and Portal validation before release acceptance.

### R5 — The promised API boundary is stale (settle before M1)

[api-stability.md](api-stability.md), lines 10–24, lists eight compatible core
packages but omits public `gale.numeric` and `gale.complex`. Both ship inside the
same core JAR and are advertised by the README. Numeric now includes exact sums,
centered moments, and distribution functions; complex includes owned Hermitian
factor/solve APIs. Their stability cannot be inferred from artifact exclusion.

The policy also calls full dense SVD and sparse-direct factorization deferred
(`release-policy.md:143–145`), while `DMat.svd` implements full economy SVD
(`Matrix.scala:805–818`) and `SparseDirectProvider.pure` exposes portable sparse
Cholesky (`SparseDirect.scala:195–198`). Sparse LU/QR availability must not be
inferred from that Cholesky implementation. The README's blanket sparse-direct
absence claim is therefore also stale.

**Acceptance:** review the actual public exports, explicitly admit or mark the new
surfaces experimental, and reconcile the README, stability inventory, and release
notes. M1 currently promises source/binary and semantic compatibility across 0.1.x;
that is a substantial commitment even for a milestone. Retain it deliberately or
revise the policy before publishing, not afterward.

### R6 — The current publication path has not been qualified (operational gate)

GitHub returned no tags, no release-workflow runs, and no repository-level secret
names. The workflow references `PGP_SECRET`, `PGP_PASSPHRASE`,
`SONATYPE_USERNAME`, and `SONATYPE_PASSWORD`. This audit did not inspect secret
values or establish the owner's Central namespace/key/Portal configuration.
Organization-level secret availability was not established; absence of repository
secrets alone is not proof that every possible credential source is absent.

[release-evidence.md](release-evidence.md) explicitly records an older eight-artifact
1.0 rehearsal at another SHA. It is not proof of the current four-artifact M1 path.
The current core-local consumer probes are useful, but do not cover anonymous
resolution and execution of all four final Maven artifacts.

**Acceptance:** run the current unsigned staging rehearsal, then configure owner
credentials and validate the real signed four-artifact bundle in USER_MANAGED
Portal mode. Exercise separate JVM and JS consumers against the staged coordinates,
including reuse of `gale-laws` and an optimized JS execution. After authorized
publication, prove anonymous Central resolution, update install lines, and retain
the version/SHA/CI/bundle/Portal receipts. A hosted guide/API URL is desirable for
M1; the current source docs can suffice if linked accurately.

## Assurance scorecard

Ratings assess the mechanism and inspected evidence, not universal correctness.

| Dimension | Rating | Evidence and remaining boundary |
|---|---|---|
| ScalaCheck/generators | Strong | `DenseLawSuite` spans zero/random shapes and strides; `AdversarialGeneratorSuite` has fixed seeds and invariant-preserving shrinkers. Both platforms passed. |
| Reusable law module | Strong | Eight public bundles target public APIs; JVM/JS law suites execute. Final staged law-artifact consumer remains a release check. |
| Framework/Discipline integration | Strong | Executed MUnit/ScalaCheck integration. Custom numerical assertions are appropriate; no downstream Discipline requirement was found. |
| Typeclass lawfulness/coherence | Not applicable | No Cats-style algebra/effect typeclass promise under review; backend `given` selection is assessed as provider conformance. |
| Backend/provider conformance | Strong | Shared strided GEMM/GEMV/SYRK oracles, native reconstruction and typed-error checks; hosted Vector and FFM gates passed. |
| Cross-platform/version CI | Present but incomplete | Hosted JVM digest failures, unprotected branch, and no browser execution gate. JS tests pass. Scala Next is advisory. |
| Numerical/computational assurance | Present but incomplete | Strong reconstruction, adversarial spectral and exact-arithmetic checks; GMRES false success is a decisive hole. |
| Independent oracles | Strong | Executed hosted Breeze parity/goldens; NumPy Hermitian fixtures; exact-sum rational oracles. Some moment checks share rounding helpers and are not wholly independent. |
| Failure/convergence/resource contracts | Present but incomplete | Typed errors, owned results, refusal/cancellation/work tests exist; GMRES violates its convergence contract. |
| Work/allocation accounting | Strong | Explicit spectral work tests and workspace/resource contracts; benchmark tooling distinguishes counts from timing. No fresh performance claim made. |
| Compiler discipline | Present but incomplete | Fatal deprecation/feature/unchecked warnings; unused and value-discard policies are not enabled. Not an M1 blocker. |
| Formatting/rewrites | Present but incomplete | Hosted formatting passed, but `.scalafmt.conf:5–25` intentionally limits file coverage. `scalafmtCheckAll` is not a whole-source formatting guarantee. |
| Binary/source compatibility | Present but incomplete | Explicit early-semver policy; no immutable baseline yet, so missing MiMa is reasonable. Export boundary needs R5 and post-M1 enforcement. |
| Coverage/mutation signal | Missing | No coverage/mutation gate found in build/plugins/CI. Diagnostic improvement, not a release blocker by itself. |
| Benchmark/performance evidence | Strong | JMH, allocation guidance, executable numerical preflight and receipt checks; dashboard withholds incomparable speedups. Hosted benchmark gate passed; timings not rerun. |
| Documentation/release evidence | Present but incomplete | Hosted executable guides/Scaladoc pass; stale admission claims, invalid-bundle false acceptance, and missing current publication receipt remain. |

## Verification and evidence

Fresh local `sbt testAllFull` exited zero: core JVM **1278**, core JS **1252**,
laws JVM **54**, laws JS **54**, plus both full-optimized test links. The alias
links optimized tests but does not execute them under FullOptStage. The launcher
used **Homebrew Java 25.0.1**, despite the shell's `java -version` reporting JDK 22;
Node was **26.7.0**. These are not the supported-baseline JDK 21/Node 22 receipts.
An initial sandboxed attempt could not acquire the sbt boot lock; the authorized
rerun completed. This was an environment failure, not a product failure.

The exact-SHA hosted run provides passing evidence for formatting, docs, Breeze
parity/goldens/published interop, JS tests, Ravel copies, dependency manifest,
benchmark gates, Vector JDK 21/22, and FFM/OpenBLAS. JVM and optimized-link jobs
failed as described above. The advisory Scala 3.8.4 source experiment failed;
the 3.8.4 artifact consumer and experimental Wasm jobs passed. These receipts
were retrieved now, not rerun remotely.

The isolated JVM/JS contract probes exited zero while confirming the GMRES defect.
The synthetic-version `releaseM1Preflight` also passed; it checks dependency,
version and local-staging configuration, not artifact contents or Portal acceptance.

Raw local evidence is retained under
[.agent-work/release-audit-20261026](../.agent-work/release-audit-20261026/):
`local-tests-authorized.log` and its metadata, `hosted-jvm.log`,
`hosted-jobs.json`, `contract-probes.log`, `ReleaseContractProbe.scala`,
`invalid-signed-bundle.zip`, `bundle-probe.log`, and remote configuration snapshots.
Scratch evidence is untracked; preserve it with this report if committing the audit.

No new benchmark timing, browser qualification, signed staging, Portal validation,
or downstream scientific/application acceptance is claimed. The assessment sampled
high-risk implementations and tests rather than proving every numerical routine.

## Smallest route to a defensible release

1. Fix R1 with independent true-residual regressions; diagnose and resolve R2.
2. Settle R5's exact public compatibility boundary and correct capability claims.
3. Enforce exact-candidate gates (R3), strengthen artifact validation (R4), and
   complete the four-artifact rehearsal and consumer checks (R6).
4. Cut `0.1.0-M1` only after green evidence and owner-authorized Central validation;
   use downstream published-artifact soak before promoting to `0.1.0`.

Keep the current numerical architecture, reusable laws, executable guides, and
qualified benchmark approach. Adding Cats, Cats Effect, Discipline, full
sbt-typelevel adoption, a broad lint migration, exhaustive coverage targets, or
more numerical features would not address the blockers. Keep the present
sbt-ci-release mechanism; add compatibility checking once M1 supplies a real
baseline. Scala 3.7.4 is a documented consumer floor: newer Scala compilers can
consume older TASTy, but older compiler compatibility must not be assumed
([Scala compatibility reference](https://docs.scala-lang.org/scala3/reference/language-versions/binary-compatibility.html)).
