# Optimization planning review evidence

Source revision: `cf1cbe0dcfbc5208e0599de7a2adaef61b890508` (initially clean).
Date: 2026-10-07.
Runtime: sbt 1.11.7, Scala 3.7.4, Homebrew Java 25.0.1, Node v26.7.0, macOS arm64.
These are local JVM and fast-linked Scala.js checks, not the JDK 21/Node 22 CI
matrix, optimized-link qualification, or performance measurements.

The [plan](../../optimization-plan.md) is the deliverable. Production code and
the normal test suites were not changed by this review.

`OptimReviewProbeSuite.scala` is a characterization fixture: its assertions
confirm current defects or limitations. **A passing probe is not an acceptance
test for a correct optimizer.** When fixing a defect, replace its characterization
with a permanent regression that asserts the corrected behavior in the normal
test source tree. Future fixes are expected to invalidate these probes.

Observed results:

- All 11 existing optimization tests passed on JVM and Scala.js.
- The initial five probes passed on both platforms in the same invocation.
- After adding the asymmetric-Rayleigh probe, all six probes passed on both
  platforms in a focused second invocation.
- Both invocations exited with status 0. Complete logs are retained as
  `existing-and-five-probes.log` and `six-probes.log`.

The six observations cover summary collisions, preceding-iterate residuals,
rounded-step false convergence, step underflow escaping `Either`, null-space
containment without completeness, and asymmetric Rayleigh false stationarity.
The containment observation is a limitation of the evidence, not a claim that
checking containment itself is incorrect.

Reproduce from the repository root at the reviewed source revision:

```sh
sbt \
  'set coreJVM / Test / unmanagedSourceDirectories += file("docs/verification/optim-review-20261007")' \
  'set coreJS / Test / unmanagedSourceDirectories += file("docs/verification/optim-review-20261007")' \
  'coreJVM/testOnly gale.optim.*' \
  'coreJS/testOnly gale.optim.*'
```

The source directory settings last only for that sbt invocation. The original
runs loaded the same fixture from `/private/tmp`; that path is visible in the
logs. The first run preceded the addition of the sixth test. The command above
runs all 17 cases per platform with the final fixture. That combined final
command is a reproduction recipe, not a claim that a third run was performed.

Source inspection, separate from executed probes, identified redundant
operator/adjoint calls, full-matrix temporaries, missing optimization benchmarks,
and assumptions requiring a mathematical audit. No speedup, exhaustive
convergence theorem, or downstream-build claim was established.
