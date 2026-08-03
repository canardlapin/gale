# Portable QR downstream admission court

Date: 2026-08-03

Decision: admit the complete portable QR stack. Against the source-pinned
pre-stack baseline, unchanged regress4s consumers improve on M1, M2, and M5,
retain their correctness receipts, satisfy the signed allocation budgets, and
retain the required batch advantages. The focused protected-control repeat
does not reproduce the short sweep's apparent M3 regression.

## Revisions and protocol

The immutable Gale baseline is commit
`83cac90a678d1b8a31c590e0c1b8fc8bf3427161`. The final candidate is commit
`bc335cf9212f9c681fd1cbd7fac3ceabb8201635`, containing the intervening
tall-skinny pivoted-QR, workspace-solve, consuming-builder, and row-scaled QR
slices. Regress4s is a non-Git source checkout; its non-target source manifest
was unchanged between runs and has SHA-256
`d05c39f8f3a3b8a1ac21c8c0f6c945b48f3d10c6c7d83f1d9ea037bc374af58c`.

The court ran on a MacBook Pro with an Apple M3 Max (14 cores, 36 GB), macOS
14.3 arm64, JDK 25.0.1, sbt 1.11.7, Scala 3.7.4, and JMH 1.37. Both broad
receipts use one thread, two forks, three 300 ms warmups, five 300 ms
measurements, average time in milliseconds, and the GC profiler. The same
regress4s tree was compiled through `-Dregress4s.gale.build=...` against the
baseline clone and final Gale checkout. The manifest hashes every non-target,
non-Mote file below the regress4s checkout in sorted path order.

```text
sbt -Dregress4s.gale.build=/absolute/path/to/gale \
  'benchmarksJVM/Jmh/run -f 2 -wi 3 -i 5 -w 300ms -r 300ms -t 1 -prof gc -rf json -rff /absolute/path/result.json regress4s.gale.*PerformanceCourtJmh.*[Ee]ndToEnd.*'
```

Raw receipts: [baseline](./2026-08-03-regress4s-portable-baseline.json),
[candidate](./2026-08-03-regress4s-portable-candidate.json),
[focused M3 baseline](./2026-08-03-regress4s-m3-baseline-exact.json), and
[focused M3 candidate](./2026-08-03-regress4s-m3-candidate-exact.json).

## Downstream results

Lower time and allocation are better. Percent changes compare the final stack
with the source-pinned baseline.

| public regress4s route | baseline ms/op | candidate ms/op | time change | baseline B/op | candidate B/op | allocation change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| M1 complete | 0.769 | 0.731 | -4.91% | 1,599,148 | 1,599,000 | -0.01% |
| M1 weighted | 1.077 | 1.041 | -3.31% | 3,024,799 | 3,024,769 | -0.00% |
| M2 batch | 0.703 | 0.687 | -2.25% | 1,514,368 | 1,514,368 | -0.00% |
| M2 sixteen scalar fits | 5.786 | 5.533 | -4.39% | 9,750,121 | 9,749,013 | -0.01% |
| M5 batch IRLS | 7.816 | 5.481 | -29.88% | 6,955,466 | 6,955,402 | -0.00% |
| M5 eight scalar fits | 18.680 | 18.645 | -0.19% | 84,765,534 | 84,765,440 | -0.00% |

The final M2 batch route is 8.05x faster than sixteen scalar fits, above the
signed 5x minimum. The final M5 batch route is 3.40x faster than eight scalar
fits, above the signed 3x minimum. M1 complete uses 1,599,000 B/op against its
2,616,722 B/op ceiling; M2 uses 1,514,368 B/op against 8,342,284 B/op; and M5
uses 6,955,402 B/op against 8,155,005 B/op. The candidate therefore satisfies
all three signed allocation budgets.

M5 is the clearest compounded downstream effect: each unchanged IRLS consumer
continues to materialize its weighted inputs, but the repeated portable
tall-skinny QR work is faster. This court does not attribute the improvement to
the new consuming or scaled APIs, because regress4s has not adopted them.

## Protected controls

The broad M3/M4/M6 controls measured time changes of +9.09%, +4.26%, and
+0.36%, respectively. All differences were inside the corresponding combined
99.9% JMH error widths, and normalized allocation changed by only +0.91%,
+0.85%, and +1.26%.

Because the short M3 point was nominally beyond the five-percent protection
rule, it was repeated with two forks, five 500 ms warmups, and ten 500 ms
measurements. The source-pinned baseline measured
`0.203125 +/- 0.002387 ms/op`; the final candidate measured
`0.202449 +/- 0.003171 ms/op`, a 0.33% improvement. Allocation changed from
615,168.86 to 615,809.65 B/op (+0.10%). The better-warmed control therefore
rejects the short sweep's apparent M3 regression.

## Correctness and scope

Each implementation slice passed Gale's full `testAllFull compileAll
parityTest interopBreezeTest benchCompile` gate. On the final stack, core passed
610 JVM and 600 Scala.js tests, laws 41 per platform, Ravel 8 per platform,
backend parity 45, Breeze interop 24, full Scala.js optimization, Scala-next
compilation, and JMH annotation processing.

Unchanged regress4s consumers were compiled and tested against the final local
Gale override. Core passed 10/10 per platform, DSL 8/8 per platform, Gale 64/64
per platform, laws 13/13 per platform, and frame4s adapters 6/6 JVM plus 5/5
Scala.js. Both JVM and Scala.js performance-court smoke runs returned
`status=ok` with unchanged M1, M2, and M5 fixture IDs, deterministic checksums,
M5 coefficients, iteration count, and convergence receipt. The nested source
override emitted only its known non-fatal Git-root, relative-glob, and existing
unused-local warnings.

This is a downstream admission receipt for unchanged consumers. The separate
microbenchmark reports retain the causal evidence for reflector traversal,
tall-skinny exact norms, workspace solves, and consuming/scaled construction.
