# Portable tall-skinny pivoted QR court

Date: 2026-08-03

Decision: select scalar-local row-first exact norm accumulation for column-pivoted
QR with `p <= 8`, and reuse the selected exact norm in Householder construction
at every width. For `p > 8`, retain independent exact column scans and reuse the
winning scan. Reject the scratch-array row-first candidate for wider designs.
No approximate norm downdate, backend dispatch, native dependency, pivot tie
rule, rank policy, ownership rule, or unpivoted QR path changes.

## Court and environment

`TallPivotedQrJmh` prepares a deterministic immutable design outside the timed
boundary. Every invocation includes Gale's mandatory row-major input copy,
exact pivot selection, compact reflector and R construction, automatic rank
decision, and owned factor results. Explicit `PureBackend` pins the shared
portable algorithm. The full court is `n={512,1024,2048,4096,10000}` and
`p={3,5,6,8,16,24}`.

The receipts were collected on a MacBook Pro with an Apple M3 Max (14 cores,
36 GB), macOS 14.3 arm64, JDK 25.0.1, sbt 1.11.7, Scala 3.7.4, and JMH 1.37.
The full court used one thread, two forks, three 300 ms warmups, five 300 ms
measurements, average time in microseconds, and the GC profiler. The focused
admission court used two forks, five 500 ms warmups, and ten 500 ms
measurements. The immutable baseline is commit
`83cac90a678d1b8a31c590e0c1b8fc8bf3427161`; candidate receipts use that commit
plus this bead's issue-scoped working-tree changes.

Representative commands:

```text
sbt 'benchmarksJVM/Jmh/run -f 2 -wi 3 -i 5 -w 300ms -r 300ms -t 1 -prof gc -rf json -rff /absolute/path/2026-08-03-tall-pivoted-qr-candidate.json gale.bench.TallPivotedQrJmh.factorPivotedQr'
sbt 'benchmarksJVM/Jmh/run -f 2 -wi 5 -i 10 -w 500ms -r 500ms -t 1 -prof gc -rf json -rff /absolute/path/2026-08-03-tall-pivoted-qr-candidate-exact.json -p n=4096 -p p=6 gale.bench.TallPivotedQrJmh.factorPivotedQr'
```

Raw receipts: [baseline sweep](./2026-08-03-tall-pivoted-qr-baseline.json),
[selected sweep](./2026-08-03-tall-pivoted-qr-candidate.json),
[focused baseline](./2026-08-03-tall-pivoted-qr-baseline-exact.json),
[focused candidate](./2026-08-03-tall-pivoted-qr-candidate-exact.json), and the
[rejected scratch candidate](./2026-08-03-tall-pivoted-qr-scratch-rejected.json).

## Regress4s shapes and admission

The final full sweep measured:

| shape | baseline us/op | candidate us/op | improvement | baseline B/op | candidate B/op |
| --- | ---: | ---: | ---: | ---: | ---: |
| `1024 x 5` | 62.077 | 52.165 | 15.97% | 90,569.43 | 90,569.20 |
| `2048 x 6` | 167.809 | 145.869 | 13.07% | 213,476.50 | 213,469.29 |
| `4096 x 6` | 335.318 | 295.489 | 11.88% | 426,562.83 | 426,558.21 |

The longer admission-critical `4096 x 6` repeat measured
`327.655 +/- 1.258 us/op` for the source-pinned baseline and
`290.528 +/- 6.903 us/op` for the candidate. This is an 11.33% improvement;
the 37.128 us/op difference is 5.29 combined JMH error widths. Normalized
allocation decreased from 426,452.56 to 426,452.03 B/op. The exact M1/M5 gate
therefore clears both the ten-percent and three-error-width requirements.

Across the full sweep, every `p=3`, `p=5`, and `p=6` point improves. Their
improvement ranges are 11.95-19.35%, 11.36-15.97%, and 11.88-14.38%,
respectively. The selected `p=8` path improves by 5.48-13.89%. Wider `p=16`
uses selected-norm reuse and improves by 2.15-6.95%. The only negative point is
`n=10000,p=24` at -2.17%, within the five-percent protection limit and well
inside its confidence interval. No normalized-allocation change is material:
the worst increase is 8.10 B/op on a 139,762 B/op result.

## Rejected wider row-first candidate

The first candidate stored every scale/sumsq pair in scratch during the
row-first matrix pass. It regressed compact designs because accumulator
load/store traffic outweighed locality: for example, `4096 x 6` measured
377.797 us/op versus the 335.318 us/op baseline. For `p=16/24`, its best full
sweep improvement was 7.70%, while protected smaller points regressed as much
as 6.69%. It therefore failed both the target-shape admission and protected
regression rules.

The selected kernel keeps up to eight scale/sumsq pairs in scalar locals, which
removes that traffic and directly covers the statistical design regime. Wider
designs keep exact column-wise scans and only reuse the selected norm. This is a
measured dispatch boundary, not an approximate numerical policy.

## Correctness and scope

The test-only reference deliberately preserves the old algorithm: each
candidate norm is recomputed in its own strided column scan and the selected
column is rescanned for Householder construction. Candidate and reference agree
exactly on permutation, rank, and R for deterministic random, tied,
rank-deficient, and overflow-sensitive fixtures. Additional tests cover safe
uniform scaling, zero columns, near ties, NaN/infinity pivot behavior,
deterministic repeats, options-aware scratch sizing and reuse, and unchanged
input/result ownership. The focused QR/workspace court passes 30/30 on both JVM
and Scala.js.

The implementation changes only portable pivoted QR and its checked workspace
requirement. Unpivoted scalar and compact-WY QR, backend routing, factorization
capabilities, and public ownership contracts are unchanged.

## Verification

The final candidate passed Gale's full `testAllFull compileAll parityTest
interopBreezeTest benchCompile` gate: core tests passed 600 JVM and 590
Scala.js cases, laws passed 41 per platform, Ravel passed 8 per platform,
backend parity passed 45, Breeze interop passed 24, full Scala.js optimization
completed, and the Scala-next and benchmark compilation gates completed.

regress4s was then compiled and tested against this checkout via
`-Dregress4s.gale.build=/Users/bbuchsbaum/code/scala/gale`. Core passed 10/10
per platform, DSL 8/8 per platform, Gale 64/64 per platform, laws 13/13 per
platform, and frame4s adapters 6/6 JVM plus 5/5 Scala.js. Both JVM and Scala.js
performance-court smoke runs returned `status=ok` with unchanged M1, M2, and M5
fixture identifiers and checksums. The nested build printed non-fatal Git
metadata warnings because it evaluates the source override outside Gale's Git
root; compilation, tests, and smoke courts all exited successfully.
