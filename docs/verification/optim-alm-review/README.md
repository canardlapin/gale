# ALM review remediation and performance qualification

Local qualification on 2026-10-08, Apple M3 Max, 36 GiB RAM, macOS 14.3 arm64.
This report supersedes the performance figures in the preserved
[original ALM report](../optim-alm/README.md). The source is an uncommitted
working tree on the recorded Git head; `source-manifest.json` identifies the
measured files by SHA-256.

## Review items addressed

1. **Verifier consistency.** Gale `solver_success` must be a Boolean exactly
   equivalent to `status == Converged`. The independent false-convergence check
   uses the status directly. Eight mutation probes reject corrupt evidence,
   including a converged status with a false flag and nonstationary multipliers.
2. **Inner line-search overhead.** ALM inner solves interpolate cached values
   and derivatives. Oriented interior guards reject high-endpoint proposals;
   insufficient bracket contraction and nonfinite endpoints trigger bisection.
   Wolfe acceptance, evaluation caps, and final KKT tolerances remain unchanged.
3. **Repeated inner work.** Tolerances adapt to the scaled PHR update residual,
   stay nonincreasing, and retain the stationarity floor. Feasible nonstationary
   and zero-step cases force tightening. A finite positive scalar inverse
   curvature carries between subproblems, resetting when the penalty changes;
   correction pairs never carry. Bounded `innerTrace` records attempted inner
   solves, including interruptions, actual callbacks, augmented requests,
   tolerances, penalties, and the carried scalar.
4. **Allocation.** Constraint finiteness checks index the vector directly.
   All-zero gradient terms bypass `ExactSum`; nonzero severe cancellation still
   uses exact accumulation of the rounded products.

Interpolation is **specific to ALM inner solves**. Applying it to standalone
L-BFGS-B caused the rotated-32 fixture to end with `LineSearchFailed`, although
its endpoint passed the looser common accuracy gate. Public L-BFGS/B therefore
explicitly retain bisection. Final standalone qualification restores all **38
case/runtime statuses, callback counts, and iteration counts** exactly, including
rotated-32: `Converged`, 57 callbacks, 45 iterations, residual 9.078e-9.
This status comparison supplements the independent endpoint gates.

## Final measured comparison

All **960 timed endpoints** passed independent accuracy checks: 16 fixtures,
15 repeats, and four solver/runtime combinations. All 480 Gale endpoints report
`Converged` and also pass independent checks using their returned multipliers.

| Runtime | Reference | Lower median time | Geometric mean reference/Gale |
|---|---|---:|---:|
| JVM | SciPy SLSQP | 14/16 | 2.58x |
| Node | SciPy SLSQP | 13/16 | 1.91x |
| JVM | SciPy trust-constr | 16/16 | 30.62x |
| Node | SciPy trust-constr | 16/16 | 22.65x |

[Per-case timings](comparison.md) retain the losses. Both runtimes lose the tiny
mixed-box cases to SLSQP; Node also loses one ball-64 case, and its other ball-64
win is nearly tied. These synthetic, inexpensive-callback results are not a
general speed guarantee or statistical significance claim. The nonisolated host
shows timing variation; these are the final-source measurements, not the best
intermediate timings.

Summed objective calls over one solve of each fixture fell from **2,417 to
1,187 (50.9%)** on both runtimes. Constraint calls fall by the same amount.
The fixture files are byte-identical to the original comparison. Gale still
uses more objective calls than SLSQP on these cases; expensive user callbacks
can change the wall-time ranking. No expensive-callback timing claim is made.

The common independent gate remains: feasibility <=1e-8, stationarity <=1e-6,
complementarity <=1e-8, point error <=1e-5, scaled objective error <=1e-7,
exact box feasibility, and nonnegative inequality multipliers. Gale's own
convergence requires stationarity <=1e-7 and raw/scaled feasibility <=1e-8.
Scalar Python `math.fsum` recomputes all diagnostics outside timing; analytic
reference multipliers provide the common gate for all solvers, with additional
Gale checks against its own returned multipliers. No threshold was relaxed.

The timing protocol uses 50 warmups per case and 15 timed repeats, all-family
Gale warmup before measurements, sequential runtimes, prepared callbacks, and
full solve costs including final diagnostics. SciPy settings are calibrated
outside the timed sample using only accurate endpoints. JDK 21.0.12.1+1,
Node 22.18.0 with fully optimized Scala.js, Scala 3.7.4, Python 3.14.7,
NumPy 2.4.3, and SciPy 1.17.1 were used. Thread-limit environment variables
were set to 1, but `threadpoolctl` enumerated no BLAS pools, so the actual BLAS
thread count was not independently confirmed.

## Ablations and limitations

The callback-only JVM ablations enable tracing; their timings are discarded.
All rows, including failures, remain in [ablations.json](ablations.json).

| Inner policy, with ALM interpolation | Converged / 16 | Objective calls over all attempts |
|---|---:|---:|
| Fixed tightening, no scalar reuse | 14/16 | 1,800 |
| Adaptive tightening only | 16/16 | 1,274 |
| Scalar reuse only | 16/16 | 1,752 |
| Adaptive + scalar reuse, default | 16/16 | 1,187 |

The first row has real robustness regressions relative to the original
bisection/fixed baseline, which converged on 16/16. Dense-32-s0 and ball-64-s1
end with `LineSearchFailed`. The latter also fails complementarity, not just
stationarity. Failed attempts are **not qualified speed results**; their work
counts are retained for diagnosis. Adaptive-only, reuse-only, and the combined
default recover all fixtures. Relative to the successful adaptive-only policy,
scalar reuse reduces aggregate calls by another 6.8%, but is not a per-case win.

A separate nearly dependent equality fixture, with rows `(1, +/-0.001)`, is
sensitive to inner policy. The default reaches the analytic projection with
feasibility 9.203e-9 and stationarity 3.433e-11; other policies can fail honestly.
The original fixed/bisection implementation also fails that harder fixture.
The evidence demonstrates policy/path sensitivity, not a proven precision
limit or general robustness on industrial constraint collections.

The fresh reviewer reproduced a changing-activity failure in an intermediate
interpolator. The high-endpoint safeguard fixes it: all 16 combinations of
constraint slopes 1e4 through 1e7, two starts, and reuse on/off converge. The
regression tests keep this case alongside the smooth quadratic probes.

## Focused work and allocation probes

The identical one-outer-solve quadratic penalty probe gives:

| Penalty | Original objective calls | Final objective calls |
|---:|---:|---:|
| 10 | 6 | 3 |
| 1e3 | 13 | 3 |
| 1e6 | 23 | 3 |
| 1e9 | 32 | 3 |

The analytic subproblem minimizer is `1/(1+penalty)`. The first three runs end
at the intentional outer iteration limit; only the last satisfies the original
constraint tolerance. These are subproblem efficiency probes, not four fully
converged constrained problems.

A JVM thread-allocation probe on an already-solved 128-variable, zero-gradient
problem measured **169,660 -> 20,195 bytes per solve (88.1% lower)** in the final
run. Both make two user callbacks. It warms 3,000 solves and measures 1,000 via
`ThreadMXBean`; this is a focused allocation measurement, not a general heap or
throughput benchmark. Sources and baseline logs are under [probes](probes/),
and final probe output is in [checks/comparisons.log](checks/comparisons.log).

## Verification and reproduction

- Formatting, `compileAll`, `testAll`, `parityTest`, `interopBreezeTest`,
  `docsCheck`, and `benchCompile` pass.
- 926 JVM and 913 Node core tests, 54 laws tests per runtime, 85 parity tests,
  and 29 Breeze interop tests pass.
- 47 affected optimizer tests pass under fully optimized Node, covering both
  search policies, budgets/cancellation, nonfinite trials, trace accounting and
  rejected-endpoint retention, near dependence, active-set changes, scale reset,
  and solve isolation.
- Independent nine-case first-order and ten-case extensions comparisons pass
  on both runtimes, with all 38 standalone work/status records restored.
- All eight verifier mutation probes reject corrupt evidence.
- Fresh independent source/evidence review found no remaining blocker after
  restricting interpolation to ALM. Optional-policy failures remain disclosed.

Logs and exit metadata are under [checks](checks/). Existing unrelated
Scaladoc/link warnings remain; the documentation commands complete successfully.
Optional FFM/native and experimental Wasm lanes were not part of this portable
qualification. No commit or push was made for this remediation.

With JDK 21, Node 22, sbt, and the pinned Python environment on PATH:

```sh
python tools/optim-alm/run_comparison.py --output /tmp/gale-alm-review --warmups 50 --repeats 15
python tools/optim-alm/check_verifier.py /tmp/gale-alm-review
```

For each variant `fixed`, `adaptive`, `scale`, and `combined`, capture:

```sh
sbt 'benchmarksJVM/runMain gale.bench.OptimizationALMBench 5 1 docs/verification/optim-alm-review/fixtures.tsv fixed fixed trace'
```

Replace both `fixed` arguments for each variant, combine the four raw logs, then:

```sh
python tools/optim-alm/check_ablations.py --fixtures docs/verification/optim-alm-review/fixtures.json --log /tmp/ablations.log --output /tmp/ablations.json
```

The preserved original report, final comparison, source hashes, raw timing rows,
calibration data, callback traces, and artifact hashes allow the claims to be
checked without relying on this narrative.
