# Finite symmetric solve bounds receipt

This bounded Gale primitive consumes one canonical, complete, row-major upper
triangle of finite `OutwardInterval` entries and the supplied finite candidate
once. It stores only O(n) per-row absolute sums and residual intervals. It
certifies a positive Gershgorin margin `g`, an outward residual L1 upper bound
`rho`, and `||xhat - x||_2 <= rho / g`; a nonpositive margin is an insufficient
certificate, not a conclusion that the represented matrix is indefinite.

The result applies only to the consumed rounded matrix enclosure, RHS enclosure,
and candidate. It does not enclose upstream assembly, interpolation,
coordinate conversion, whitening, identity binding, or ScalaFIM scientific
admission.

## Continuation verification, 2026-10-01

The original JVM attempt failed only because its rational oracle compared raw
case-class fields `45/324` and `5/36`. The assertion now compares the normalized
numerator/denominator pair, retaining the exact `5/36` condition. The original
failed raw log and metadata remain intact. A smallest-subnormal diagonal control
also confirms honest refusal when outward subtraction leaves no positive margin.

Using the preserved primary Eclipse Temurin JDK 21.0.12.1, the focused suite passed
5/5 on JVM and Scala.js. The broader gate passed `compileAll`, core JVM 679/679,
core JS 669/669, laws 54/54 on each platform, Breeze differential parity 84/84,
and Breeze interop 29/29. Repository `scalafmtCheckAll` passed. Because its normal
include-path policy omits these new files, both were explicitly formatted and
checked using the same formatter settings with a two-file include list.

The compile alias also publishes the candidate to the local Ivy cache for its
Scala consumer check; no remote publication occurred. Existing Scaladoc warnings
at `Matrix.scala:615` and `Sized.scala:158` and a repeated classpath flag remain,
matching the previously preserved outward-interval compile receipt. Runtime
native-provider warnings in the optional Breeze lanes are retained. This is not
a warning-free whole-build claim.

The read-only independent numerical review accepted the exact implementation
SHA256 `4b27971646311dce631ecb775613ac19b1d503a02ab03b5c6216c3e04407b8c3`
and test SHA256 `fc4d02205d5df0e033216fe96ec2672881f79b67fd3e662d8a7649dbaa48a35d`.
It verified the Gershgorin/residual inequality and O(n) auxiliary retained state.
Runtime and cumulative allocation remain O(n²). A refusal does not establish
singularity or indefiniteness; the Gershgorin certificate is sufficient only.

Full commands, raw logs, metadata and review hashes are bound by
`/private/tmp/scalafim-execution-20260929/gale-finite-resume-receipt.json`.
Raw gates are `logs/gale-finite-resume-focused-v1.log` and
`logs/gale-finite-required-v1.log` beneath that execution root; each has a
`.meta.json` sidecar. The independent review is
`gale-finite-resume-independent-review.md`. This candidate has not changed any
ScalaFIM pin, admission policy, performance claim or scientific certificate.
