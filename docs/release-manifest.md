# Gale 1.0 release manifest

This is the authoritative artifact boundary for the first Gale release
candidate. A project listed as **published** contributes every platform variant
shown below at one version. A project listed as **excluded** remains available
for development or evidence work but must not appear in the release bundle or
in a published dependency graph.

## Published artifacts

| Build project | Maven coordinate(s) | Platform | Build pass |
| --- | --- | --- | --- |
| `coreJVM`, `coreJS` | `io.github.canardlapin:gale-core_3`, `io.github.canardlapin:gale-core_sjs1_3` | JVM, Scala.js | JDK 21 |
| `lawsJVM`, `lawsJS` | `io.github.canardlapin:gale-laws_3`, `io.github.canardlapin:gale-laws_sjs1_3` | JVM, Scala.js | JDK 21 |
| `interopBreeze` | `io.github.canardlapin:gale-interop-breeze_3` | JVM | JDK 21 |
| `vectorBackend` | `io.github.canardlapin:gale-backend-jvm-vector_3` | JVM | JDK 21 |
| `nativeBackend` | `io.github.canardlapin:gale-backend-jvm-native_3` | JVM | JDK 22 |
| `blasFfmBackend` | `io.github.canardlapin:gale-backend-jvm-blas-ffm_3` | JVM | JDK 22 |

Every published coordinate is released at one tag-derived version and must
have a binary JAR, sources JAR, Scaladoc JAR, POM, checksum, and signature.
The JDK 21 and JDK 22 passes are disjoint build jobs but one publication bundle;
neither pass may publish a partial version.

## Explicitly excluded projects

The following projects are build, integration, or evidence surfaces and are
not part of the 1.0 publication bundle:

- `interopRavelJVM` and `interopRavelJS`: the current dependency is
  `ravel-core_3:1.0.0-SNAPSHOT`, so this integration remains opt-in until Ravel
  has a stable published coordinate;
- `parity`, benchmark projects, `demo`, and the documentation site;
- `scalaNextConsumer`, the Scala-next compatibility probe; and
- the root aggregator.

The excluded projects may be compiled or tested by dedicated commands when
their external dependencies are available. They must not be pulled into
`compileAll`, `testAll`, `testAllFull`, or any release publication task.

## Executable checks

Run the release dependency gate before packaging:

```sh
sbt releaseDependencyCheck
```

The gate inspects both compile and runtime resolution for every admitted
project and fails on any dependency revision containing `SNAPSHOT`. It does not
silently treat a sibling checkout or a local snapshot as release evidence.

The current development build derives a unique `1.0.0+...-SNAPSHOT` version;
that version identifies the source tree only. A release candidate must derive
one non-snapshot `v1.0.0-RC1` (or later release) version from its tag and run
the same gate before assembling artifacts. `.github/workflows/release.yml`
stages JDK 21 and JDK 22 slices,
merges them, checks signatures and metadata with
`tools/verify-central-bundle.sh`, and uploads one user-managed Central Portal
deployment.
