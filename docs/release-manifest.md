# Gale 0.1 release manifest

This file is the authoritative artifact boundary for `v0.1.0-M1`. The first
milestone publishes the portable core and its reusable laws. JVM integrations
and acceleration backends remain part of the candidate test court, but they do
not enter the milestone bundle until a later admission decision records their
compatibility and runtime support.

## Admitted milestone artifacts

| Build project | Maven coordinate | Platform | Build pass |
| --- | --- | --- | --- |
| `coreJVM` | `io.github.canardlapin:gale-core_3` | JVM | JDK 21 |
| `coreJS` | `io.github.canardlapin:gale-core_sjs1_3` | Scala.js | JDK 21 |
| `lawsJVM` | `io.github.canardlapin:gale-laws_3` | JVM | JDK 21 |
| `lawsJS` | `io.github.canardlapin:gale-laws_sjs1_3` | Scala.js | JDK 21 |

Every admitted coordinate is released at one tag-derived version and includes
a binary JAR, sources JAR, Scaladoc JAR, POM, checksums, and signatures. The
bundle verifier rejects missing coordinates, additional coordinates, mixed
versions, snapshots, coordinates outside the admitted Maven group, missing
files or signatures, absent or incorrect MD5/SHA-1 checksums, invalid POM coordinates/metadata or
dependencies, invalid/empty JAR contents, and unexpected files. Signed validation
cryptographically verifies each payload against the configured primary key fingerprint.

## Tested but provisional modules

| Build project | Intended coordinate | Why it is provisional |
| --- | --- | --- |
| `interopBreeze` | `gale-interop-breeze_3` | JVM-only migration helper; source and POM probes remain required, but it has no 0.1 compatibility baseline |
| `vectorBackend` | `gale-backend-jvm-vector_3` | Depends on the incubating JDK Vector API and measured runtime dispatch |
| `nativeBackend` | `gale-backend-jvm-native_3` | Requires JDK 22 native access and explicit lifetime ownership |
| `blasFfmBackend` | `gale-backend-jvm-blas-ffm_3` | Requires JDK 22, a discovered native library, and platform-qualified thresholds |

These projects remain publishable to an isolated local repository so CI can
inspect their POMs and compile clean consumer probes. The `releaseM1Unsigned`
and `releaseM1Signed` aliases exclude them, and the Central bundle verifier
rejects them if they appear in the M1 bundle.

## Excluded projects

The following projects are development, integration, or evidence surfaces and
are not part of the 0.1 publication bundle:

- `interopRavelJVM` and `interopRavelJS`: the current dependency is
  an unpublished pinned Ravel source build, so this integration remains opt-in until Ravel
  has an immutable published coordinate;
- `parity`, benchmark projects, `demo`, and the documentation site;
- `scalaNextConsumer` and `publishedInteropConsumer`; and
- the root aggregator.

Excluded projects may have dedicated compile or test commands. Their presence
in the build does not admit them to publication.

## Executable checks

Run the dependency gate while converging a candidate:

```sh
sbt releaseDependencyCheck
```

The publication aliases are self-guarding. On an exact tag, run:

```sh
sbt releaseM1Preflight
```

For the credential-free manual court, supply the same synthetic version used by
the workflow:

```sh
sbt 'set ThisBuild / version := "0.1.0-M1"' releaseM1Preflight
```

The preflight repeats the dependency and version checks for all four admitted
projects and requires each `publishTo` destination to be sbt's local
`target/sona-staging` repository. `releaseM1Unsigned` and `releaseM1Signed`
invoke this preflight themselves before writing any artifact.

The gate checks compile and runtime resolution for the stable and provisional
candidate modules. It fails on external dependency revisions containing
`SNAPSHOT`; it does not treat a sibling checkout or a locally published
snapshot as release proof.

Untagged source derives a unique `0.1.0+...-SNAPSHOT` version. That version is a
development identifier, not an immutable dependency. Publication accepts only
an exact `v0.1.x`, `v0.1.x-Mn`, or `v0.1.x-RCn` tag. The first intended tag is
`v0.1.0-M1`.

The release workflow stages the four admitted modules, validates their Maven
layout, contents and signatures with `tools/verify-central-bundle.sh`, and uploads one
user-managed Central Portal deployment. A manual dry run uses a validated
synthetic milestone or RC version, creates an unsigned bundle, and never reads
release credentials or uploads it.

## Signed validation and consumer rehearsal

For owner setup, follow [the secure setup guide](release-setup.md).

Configure `PGP_SECRET`, `PGP_PASSPHRASE`, `SONATYPE_USERNAME`, and
`SONATYPE_PASSWORD` as repository secrets, and the public primary signing-key
fingerprint as repository variable `RELEASE_PGP_FINGERPRINT`. The signing job
exports only that public key for verification in the separate upload job.

```sh
tools/verify-central-bundle.sh --require-signatures \
  --public-key release-public-key.asc \
  --expected-fingerprint "$RELEASE_PGP_FINGERPRINT" central-bundle.zip
```

The release workflow runs the same required CI court as branch CI before staging.
It then verifies the bundle and runs `compat/release-consumer` from the staged Maven
repository on JVM, fast-linked JS, and fully optimized JS. That build has no sibling
source dependencies and does not use `publishLocal` to satisfy Gale coordinates.
Unsigned dry runs do not authenticate signatures and are labeled accordingly.
Only a real signed Portal validation establishes publication readiness.
