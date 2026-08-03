# QR workspace solve admission

This receipt evaluates `QR.solveLeastSquaresWith` for vector and matrix
right-hand sides. The factor and immutable inputs are built outside the timed
boundary. Each solve includes the logical RHS copy, implicit reflector
application, triangular solve, pivot unpermutation, typed result construction,
and owned coefficient storage. Workspace construction is excluded because the
API contract is explicit caller-owned reuse.

The source-pinned baseline is commit
`8e563b2ea46ecaa648b71bf30a1f9467eeb14301`; its only issue-scoped change was
the vector benchmark control. The candidate uses that commit plus this bead's
workspace implementation and benchmark methods. Raw receipts are the
[baseline](./2026-08-03-qr-workspace-solves-baseline.json) and
[candidate](./2026-08-03-qr-workspace-solves-candidate.json).

## Environment and protocol

The court ran on a MacBook Pro with an Apple M3 Max (14 cores, 36 GB), macOS
14.3 arm64, JDK 25.0.1, sbt 1.11.7, Scala 3.7.4, and JMH 1.37. It used one
thread, two forks, five 500 ms warmups, ten 500 ms measurements, average time
in microseconds, and the GC profiler. The exact regress4s M2 shape is
`n=2048,p=6,q=16`.

```text
sbt 'benchmarksJVM/Jmh/run -f 2 -wi 5 -i 10 -w 500ms -r 500ms -t 1 -prof gc -rf json -rff /absolute/path/candidate.json -p n=2048 -p p=6 -p q=16 gale.bench.QrMultiRhsJmh.solveLeastSquares.*'
```

## Results

| route | time us/op | error us/op | allocation B/op |
| --- | ---: | ---: | ---: |
| baseline matrix allocating | 147.449 | 2.250 | 263,058.044 |
| candidate matrix allocating | 146.895 | 0.800 | 263,058.037 |
| candidate matrix workspace | 147.769 | 1.008 | 938.040 |
| baseline vector allocating | 27.562 | 0.090 | 16,608.382 |
| candidate vector allocating | 28.225 | 1.408 | 16,608.392 |
| candidate vector workspace | 28.345 | 1.075 | 184.392 |

The matrix workspace route removes 262,120.004 B/op, exceeding the 200 KB/op
M2 requirement. Its time is 0.22% above the source-pinned allocating baseline
and 0.60% above the same-build allocating control. The existing matrix API is
0.38% faster than baseline. The vector workspace route removes 16,423.990 B/op
and is 2.84% above the source-pinned allocating baseline; the candidate's
allocating vector control is 2.41% above baseline. All timing comparisons remain
inside the five-percent protection rule.

The residual 938 B/op matrix and 184 B/op vector costs include the unavoidable
owned coefficient storage, result wrappers, and `Either`; the transformed
`2048x16` matrix and length-2048 vector are retained in grow-only workspace
storage across calls.

## Correctness and ownership scope

The public allocating methods are unchanged. Workspace overloads copy logical
values directly from contiguous or strided immutable inputs into checked
scratch, apply the same reflectors and triangular kernel there, and unpermute
directly into one final owned result. Tests require exact equality with the
allocating route for vector and matrix inputs, identity and nonidentity
permutations, transposed and sliced storage, `q={0,1,17}`, repeated calls, and
extreme scales. Typed dimension, underdetermined, and rank-deficient failures
are preserved and occur before scratch acquisition where their inputs make that
possible. Results never retain workspace storage, and earlier results remain
stable after reuse.

## Verification

The focused QR/workspace court passed 35/35 tests on both JVM and Scala.js.
Gale's full `testAllFull compileAll parityTest interopBreezeTest benchCompile`
gate then passed: core 605 JVM and 595 Scala.js tests, laws 41 per platform,
Ravel 8 per platform, backend parity 45, Breeze interop 24, full Scala.js
optimization, Scala-next compilation, and JMH annotation processing.
