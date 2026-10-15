# JVM dependency release candidate

This branch prepares `1.0.0-RC1` from source `d55fe2f97196a76ab7879e1a12f1e92403aeba06` for StoryModel's
binary dependency chain. It is an explicitly limited release candidate, not the full
stable release of every module and platform. Application sources are unchanged.

The publication set is:

- `io.github.canardlapin:gale-core_3:1.0.0-RC1`

Run the required source checks before staging:

```sh
sbt coreJVM/test parity/test
sbt stageDependencyRelease
```

`stageDependencyRelease` writes the selected artifacts into sbt's local staging
repository. It does not upload to Central. The completed bundle must contain the
POM, binary JAR, sources JAR, Javadoc JAR, detached signatures, and checksums for
each coordinate. Verify the actual bundle and fresh binary consumer, then upload
using the verified `io.github.canardlapin` Central account. Keep tokens and signing
keys outside this repository. Retain the exact commit, bundle digest, validation
result, and public Central lookup before marking publication complete.

Source and binary dependency changes on this branch require another check. The
broader module/platform stable-release policy remains applicable to a future
stable release.
