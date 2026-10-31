# Gale development, compatibility, and release policy

The published checkpoint is `v0.1.0-M1`. `main` develops the next release, 0.2,
and derives `0.2.0+...-SNAPSHOT` versions. A snapshot names a development state;
it is not a published release or a compatibility promise. Cross-repository
validation must identify snapshots by their full Git commit.

## Branches and merging

- Use short-lived feature branches and reviewed PRs into `main`. Merge coherent
  increments once their API decisions, tests, and executable documentation are
  complete and the required hosted checks pass on the final candidate. Deferred
  features must be stated explicitly; unfinished semantics stay on the branch.
- Keep `main` usable with passing required checks. Require an up-to-date PR and
  review public API and numerical-contract changes against the published baseline.
- `release/0.1` starts at the exact `v0.1.0-M1` commit. Use it for later 0.1
  milestones, release candidates, final stabilization, and maintenance. Selectively
  backport applicable fixes from `main`; forward-port urgent maintenance fixes.
  Do not merge 0.2 feature development wholesale into that branch.
- Release tags and published artifacts are immutable. Merging or pushing a branch
  does not publish a release. Tagging, artifact admission, and Central publication
  are separate, explicitly authorized release actions.
- Keep the README's default dependency on the published version. Default released
  documentation should come from its release tag; main-branch previews must say
  they are unreleased and identify the required source version. Track changes in
  the root [changelog](../CHANGELOG.md) under `Unreleased — 0.2` until release.

`1.0.0-SNAPSHOT` was a premature development label. The numerical core is past
the PRD 0.1–0.5 feature milestones, but a `1.0` tag would freeze public APIs,
artifact admission, and SemVer for the `1.x` line before downstream soak,
publication, and a reviewed export baseline exist. This policy replaces that
label with a realistic progression.

## Version progression

1. `v0.1.0-M1` is the first published, immutable ecosystem checkpoint.
2. `0.2.0+...-SNAPSHOT` identifies current main-branch development. The 0.1
   maintenance branch retains its own 0.1 version configuration.
3. Later `v0.1.0-Mn` or `v0.1.0-RCn` tags may stabilize the maintenance line;
   `v0.1.0` requires the milestone's downstream soak gate.
4. `v0.2.0-Mn`, `v0.2.0-RCn`, and `v0.2.0` require an explicit update to release
   admission and a qualified candidate. This development transition does not
   authorize those versions for publication.
5. `v1.0.0` comes later, as the first compatibility-locked GA, after Gale has
   been depended on in published form and the remaining freeze work in this
   document is closed.

Only an exact `v0.1.x`, `v0.1.x-Mn`, or `v0.1.x-RCn` tag may produce a
publishable version today. `sbt versionPolicyCheck` checks clean admitted tags,
development and dirty states, rejected tag forms, and the no-Git fallback's
snapshot classification. CI runs it in the required quality job.
Any source, build, dependency, or executable documentation change creates a new
candidate and invalidates evidence tied to the old commit.

## Artifact admission

The organization is `io.github.canardlapin`; Scala 3 artifacts use the `_3`
binary suffix. The first milestone admits four coordinates:

| Artifact | Platforms | Compatibility status |
| --- | --- | --- |
| `gale-core` | JVM, Scala.js | stable from `v0.1.0-M1` |
| `gale-laws` | JVM, Scala.js | stable from `v0.1.0-M1` |

`gale-core` contains the public dense, sparse, solver, optimization, spectral,
sized, kernel, error, diagnostic, and backend-contract packages. Gale does not
publish empty package-shaped modules. `gale-laws` is a normal module because
downstream libraries extend its MUnit and ScalaCheck conformance suites.

The following modules are tested but provisional and are excluded from the M1
bundle: Breeze interop, the JDK Vector backend, native matrix storage, and the
FFM BLAS/LAPACK backend. Exclusion means their tests remain part of the
candidate court while their coordinates carry no 0.1 publication or
compatibility promise. Admission requires a later explicit manifest change —
likely when preparing `1.0.0`, not automatically at M1.

`gale-interop-ravel` is excluded because it currently depends on
an unpublished Ravel source revision rather than a released Central coordinate. Parity projects, benchmarks, documentation, demos,
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
The package-by-package inventory is [API stability](api-stability.md).

M1 is now an immutable published compatibility baseline. Do not describe a
candidate as binary-verified without an executed compatibility check.
`sbt -Dgale.compatibility.baseline=0.1.0-M1 compatibilityCheck` runs MiMa
(binary) and TASTy-MiMa (TASTy) for all four admitted artifacts against the
named published version; without a baseline the alias fails rather than
comparing against nothing. CI runs it against `0.1.0-M1` in the required
quality job. Intentional breaks on the 0.2 line require explicit review and
migration notes: change the CI baseline or add a reviewed MiMa filter rather
than silently resetting it. Once `0.2.0` is published it becomes the baseline
for `0.2.x`.

`1.0.0` is a later freeze, not the next tag. It may admit additional modules
and begins the `1.x` SemVer promise. It still requires a reviewed export
baseline; there will be no MiMa history for `1.0.0` itself until `1.1.0`.

## Remote CI classification

The branch-protection configuration must require these exact job names on a
release candidate commit:

- `Formatting and fatal-warning policy`
- `Tests (Scala 3.7.4, JVM)` and `Tests (Scala 3.7.4, JS)`
- `Scaladoc and executable guides`
- `Scala.js optimized links`
- `Breeze parity and published interop`
- `Release dependency manifest`
- `Benchmarks compile`
- `Vector backend (JDK 21)` and `Vector backend (JDK 22)`
- `FFM BLAS/LAPACK backend (JDK 22 / OpenBLAS)`

`Scala 3.8.4 consumer probe`, `Scala Next source experiment (advisory)`, and
`WebAssembly tests and profile (experimental)` are intentionally advisory.
They remain visible for compatibility and research evidence but cannot turn a
release candidate green or red. The required names are defined in
`.github/workflows/ci.yml`; a release record must capture the exact run URL and
classify every other check as advisory, manual, or unverified.

The Breeze job's consumer probe publishes the admitted `gale-core` and the
provisional `gale-interop-breeze` coordinates to an isolated local repository
and then compiles a separate project from those coordinates. It is a package/POM
and transitive-dependency check, not a sibling-source compilation or a claim
that Central publication has already occurred.

## Supported runtimes

| Module/route | Minimum supported runtime | Notes |
| --- | --- | --- |
| `gale-core` JVM, laws | JDK 21 | primary required CI line |
| Vector backend | JDK 21 | tested; incubating API; provisional until admitted |
| Native and BLAS/LAPACK FFM | JDK 22 | tested; provisional until admitted |
| Scala.js JavaScript | Node 22 and current evergreen browsers | optimized JS is the supported browser performance route |
| Scala.js Wasm | Node 25+ experimental profile | explicit, allow-failure, and not covered by 0.1 compatibility or performance promises |

“Current evergreen browsers” means the latest two stable major releases of
Chrome, Firefox, Safari, and Edge at the time a Gale release is cut. Browser CI
is not yet a release gate, so 0.1 release notes must state that Node is the
tested Scala.js runtime and must not overstate browser certification.

## What 1.0 still needs

Do not retarget snapshots or tags to `1.0` until the following are closed:

- an immutable 0.1 baseline that downstream consumers have actually resolved
  from a public repository;
- a decided artifact boundary for backends and Breeze interop;
- a manual API export review, including remaining source migrations such as
  `writeLinear`;
- spectral constructor tightening (`private[gale]` → `private[spectral]`) or an
  explicit decision to defer it as non-public;
- hosted Scaladoc and a real install line;
- owner Central namespace, signing, and Portal validation on the exact tag.

Full economy dense SVD and opt-in portable sparse Cholesky are implemented
and included in core. Deferred from both 0.1 and 1.0, and not release blockers:
public Schur/QZ, sparse-direct LU/QR, MINRES, Ravel publication, and
Wasm as a supported route.

## Release blockers and provenance

A public milestone requires all of the following in addition to green tests:

- a repository license chosen by the owner and represented in both a root
  `LICENSE` file and published POM metadata;
- project homepage, SCM, and developer metadata in the POM;
- a configured publishing destination and credentials outside the repository;
- clean source, binary, and documentation artifacts for every admitted module;
- the acceptance and release-evidence records for the exact release commit.

The owner selected Apache-2.0 and `https://github.com/canardlapin/gale` as the
canonical SCM repository on 2026-07-19. The root `LICENSE` and generated POMs
carry that provenance. Binary publication still requires a configured
destination, credentials, signing policy, and remote CI on the exact release
commit. `.github/workflows/release.yml` calls the complete reusable CI workflow on its
exact candidate before staging and runs independent staged JVM/JS consumers
(including optimized JS execution) before upload. It is tag-only (or an explicitly selected
tag dry-run), imports the ephemeral PGP key only in release jobs, validates the
four-artifact signed bundle against the repository variable
`RELEASE_PGP_FINGERPRINT`, and uploads it to Central Portal as
`USER_MANAGED`. The final Central publish action remains a manual Portal
operation after validation; `publishLocal` is useful packaging evidence but is
not a public release.
