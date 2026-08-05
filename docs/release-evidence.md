# Gale 1.0.0-RC1 release-foundation evidence

> **Historical record:** This rehearsal predates the decision to establish the
> first public compatibility line at `0.1`. It preserves engineering evidence
> for the old eight-artifact RC design, but it does not define the active
> release target or artifact manifest. See `docs/release-policy.md`.

Evidence date: 2026-08-03. Candidate commit:
`6b48b78` (`6b48b782d557ab120de4c14e32f8dc6eaeb96442`), on
`perf/pivoted-qr-norm-downdate`. The exact pushed SHA is recorded in Mote and
in the remote CI receipt below. No release tag was created in the real
checkout, and no Central deployment was uploaded.

This record covers the four pre-RC foundation workstreams: documentation,
artifact/dependency boundary, required CI, and signed Central Portal bundle
mechanics. It does not claim that the owner's Central namespace, project PGP
key, or Portal token has been configured.

## Exact candidate and environment

- Repository: `io.github.canardlapin/gale`
- Branch: `perf/pivoted-qr-norm-downdate`
- Exact pushed SHA: `6b48b782d557ab120de4c14e32f8dc6eaeb96442` (checked by
  `git rev-parse` and `git ls-remote` before certification)
- Local macOS Apple ARM64 court: sbt 1.11.7, Scala 3.7.4, Homebrew OpenJDK
  25.0.1, Node supplied by the local Scala.js toolchain
- All local sbt commands used isolated `/tmp/gale-*` Coursier, sbt, and Ivy
  caches; no sibling checkout supplied a dependency
- Remote CI uses the checked-in JDK 21/JDK 22 and Node 22 matrix definitions
- Mote epic: `bd-01KZ447W504M5WMBHZF5CKWYKD`

## Local correctness, compatibility, and documentation gates

The final candidate passed:

```text
git diff --check
sbt scalafmtCheckAll
sbt releaseDependencyCheck
sbt testAllFull docsCheck
sbt parityTest interopBreezeTest
```

The principal counts from the final `testAllFull docsCheck` run are:

| Gate | Result |
| --- | ---: |
| Core JVM tests | 622 / 622 |
| Core Scala.js tests | 612 / 612 |
| Laws JVM tests | 41 / 41 |
| Laws Scala.js tests | 41 / 41 |
| Scala.js optimized test links | pass |
| JVM/Scala.js Scaladoc | pass |
| mdoc pages | 26 compiled |
| Laika pages | 22 rendered |
| Breeze parity | 45 / 45 |
| Breeze conversions/migration | 24 / 24 |

`PivotedQRScreenSuite` contributes seven focused equivalence and stability
checks inside the JVM/JS totals. The required remote run below independently
executes the same public gates on Ubuntu and the JDK matrix.

## Frozen artifact and dependency boundary

The authoritative [release manifest](release-manifest.md) admits exactly these
eight coordinates:

- `gale-core_3`
- `gale-core_sjs1_3`
- `gale-laws_3`
- `gale-laws_sjs1_3`
- `gale-interop-breeze_3`
- `gale-backend-jvm-vector_3`
- `gale-backend-jvm-native_3`
- `gale-backend-jvm-blas-ffm_3`

`gale-interop-ravel` remains excluded because its `ravel-core_3:1.0.0-SNAPSHOT`
dependency cannot enter a stable or RC Gale POM. `releaseDependencyCheck`
rejects prohibited compile/runtime snapshots for every admitted project.

## Tag-derived signed bundle rehearsal

In a temporary clone of the exact candidate, a local-only `v1.0.0-RC1` tag
derived version `1.0.0-RC1`, passed `releaseVersionCheck`, and selected the
`target/sona-staging` repository. With an isolated ephemeral GPG home (not the
owner's release key), the following aliases completed:

```text
sbt releaseVersionCheck releaseJdk21Signed releaseJdk22Signed
```

The two disjoint JDK slices merged into one Maven-layout tree containing 192
files: all eight POM/JAR/sources/Scaladoc sets, checksums, and 32 detached
signatures. The exact RC bundle passed:

```text
tools/verify-central-bundle.sh --require-signatures central-bundle.zip
Central bundle verified: 8 admitted artifacts, version 1.0.0-RC1
```

An unsigned-bundle rehearsal separately failed the signature-required verifier,
confirming the gate fails closed. The release workflow remains tag-only,
imports signing secrets only in tag jobs, merges JDK 21/JDK 22 slices, and sends
one `USER_MANAGED` Central Portal deployment without automatically publishing
it.

## Remote CI receipt

GitHub run [30837452017](https://github.com/canardlapin/gale/actions/runs/30837452017)
ran on the exact candidate `6b48b782d557ab120de4c14e32f8dc6eaeb96442`. The
completed run has every required job and both required JVM/JS and Vector
matrices green:

- formatting and fatal-warning policy
- JVM and Scala.js tests
- Scaladoc and executable guides
- Scala.js optimized links
- Breeze parity and published interop
- release dependency manifest
- benchmark compilation
- Vector JDK 21 and JDK 22
- FFM BLAS/LAPACK on JDK 22/OpenBLAS

The Scala 3.8.4 consumer/source lanes are advisory and passed. The Wasm lane is
also advisory and failed; `continue-on-error: true` keeps it from masking the
required release gates. This receipt is the source-level record to carry into
the later RC tag and Central validation operation.

## Remaining RC-only gates

The engineering foundation is locally and remotely exercised, but the final
publication child remains open until the owner supplies and validates:

1. the `io.github.canardlapin` Central namespace;
2. the real project PGP key and passphrase in GitHub Actions secrets;
3. Central Portal username/token secrets and one non-publishing validation
   upload; and
4. the owner's Central validation receipt for the exact candidate.

No immutable `v1.0.0` tag, RC tag, Central upload, or public binary release was
created by this foundation work.
