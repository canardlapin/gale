# Gale 0.1 compatibility and versioning plan

Gale is on the `0.1` line. Untagged source derives a unique
`0.1.0+...-SNAPSHOT` version. A snapshot names a development state; it is not a
published release or a compatibility promise. Cross-repository validation must
identify snapshots by their full Git commit.

`1.0.0-SNAPSHOT` was a premature development label. The numerical core is past
the PRD 0.1–0.5 feature milestones, but a `1.0` tag would freeze public APIs,
artifact admission, and SemVer for the `1.x` line before downstream soak,
publication, and a reviewed export baseline exist. This policy replaces that
label with a realistic progression.

## Version progression

1. `0.1.0+...-SNAPSHOT` for exact-commit development and consumer convergence.
2. `v0.1.0-M1` for the first immutable ecosystem checkpoint.
3. Later `v0.1.0-Mn` or `v0.1.0-RCn` tags when another candidate is needed.
4. `v0.1.0` only after the milestone passes its downstream soak gate.
5. `v1.0.0` later, as the first compatibility-locked GA, after the 0.1 line has
   been depended on in published form and the remaining freeze work in this
   document is closed.

Only an exact `v0.1.x`, `v0.1.x-Mn`, or `v0.1.x-RCn` tag may produce a
publishable version today. Any source, build, dependency, or executable
documentation change creates a new candidate and invalidates evidence tied to
the old commit.

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
The package-by-package inventory is [API stability](api-stability.md).

There is no compatibility baseline before M1 because no immutable artifact
exists. Until M1 is published, maintainers must review exported API diffs
manually and must not describe a candidate as binary-verified. Immediately
after M1, 0.1-line development should add an automated binary- and
TASTy-compatibility check against `0.1.0-M1` and make that job blocking.

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

Deferred from both 0.1 and 1.0, and not release blockers: full dense SVD,
public Schur/QZ, sparse-direct factorization, MINRES, Ravel publication, and
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
commit. `.github/workflows/release.yml` is tag-only (or an explicitly selected
tag dry-run), imports the ephemeral PGP key only in release jobs, validates the
four-artifact signed bundle, and uploads it to Central Portal as
`USER_MANAGED`. The final Central publish action remains a manual Portal
operation after validation; `publishLocal` is useful packaging evidence but is
not a public release.
