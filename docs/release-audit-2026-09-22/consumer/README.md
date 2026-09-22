# Isolated artifact consumer and GMRES regression probe

This project belongs to the release audit, not Gale's normal source/test build.
It depends only on the four staged Maven artifacts through `gale-laws` and its
transitive `gale-core` dependency. Ivy-local is excluded from `externalResolvers`.
Third-party dependencies resolve from Maven Central.

The source under `PublishedSuite.scala` checks public law-suite usability plus
analytic banded and sparse Cholesky solves. `GmresContractProbe.scala` contains
intentionally failing assertions against main `1b018d4d`, with two tolerance
modes on each platform. Do not treat a default aggregate `test` failure here as
an artifact-packaging failure or weaken these assertions to obtain green tests.

To reproduce, use an isolated clone of the audited SHA and stage without uploading:

```sh
sbt 'set ThisBuild / version := "0.1.0-M1"' releaseM1Unsigned
```

Copy this consumer directory to a temporary directory. From that copy, supply an
absolute file URL for the clone's `target/sona-staging` directory:

```sh
sbt -Dgale.audit.maven=file:///absolute/path/to/target/sona-staging/ \
  'jvm/testOnly externalconsumer.Published*' \
  'js/testOnly externalconsumer.Published*' \
  'set js / scalaJSStage := FullOptStage' \
  'js/testOnly externalconsumer.Published*'

sbt -Dgale.audit.maven=file:///absolute/path/to/target/sona-staging/ \
  'jvm/testOnly externalconsumer.GmresContractProbe'

sbt -Dgale.audit.maven=file:///absolute/path/to/target/sona-staging/ \
  'js/testOnly externalconsumer.GmresContractProbe'
```

Use Node and a compatible JDK; the audit used JDK 22 / Node 26.7.0. Normal Gale
CI separately passed on its supported JDK 21 / Node 22 lane. The synthetic M1
version above is a local rehearsal label, not an existing published release.
The retained logs show the original absolute staging path; the saved build
accepts it through `gale.audit.maven` for portability.
