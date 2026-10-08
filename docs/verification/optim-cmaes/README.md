# CMA-ES qualification, 2026-10-07

Mote: `bd-01M4C84B17HGAKJ68AXEKKXSSM`. This is local qualification of an
uncommitted change on base `831009acb033e795cb116d5143071da4350157dc`.
The [executable guide](../../user/guides/cma-es.md) defines the public contracts.

## Implementation and stopping contracts

Portable `CMAES` implements positive-weight rank-one/rank-mu covariance
adaptation and cumulative step-size adaptation. It provides seeded ask/tell and
serial minimize, capped whole-vector rejection for bounds, elimination of fixed
coordinates, explicit evaluation/generation budgets, partial receipts, cancellation,
bounded progress traces, and optional increasing-population restarts. It owns
candidate/state storage and permits only one outstanding batch. Invalid receipts
are rejected atomically. Partial generations preserve observations without
adapting the distribution. Control cancellation or failure vetoes restarts.

Results distinguish target attainment, work limits, stagnation, small steps,
condition limits, unrepresentable proposals, sampling exhaustion, and callback
or numerical failures. None is a stationarity or global-optimality certificate.
The adapted covariance can lead its cached sampling factor; bounds condition
accepted samples. There is no silent eigenvalue repair or covariance reset.

## Independent recurrence and numerical evidence

The numerical oracle is `cmaes==0.12.0`, with two explicit adapters: negative
recombination weights are zeroed, and its initial generation counter is set to
-1. That package increments the counter before using `(g + 1)` in the h-sigma
transient denominator; the adapter makes the first exponent 2, as in the
canonical recurrence. No update code is replaced. An independent source review
confirmed that this counter shift otherwise affects only diagnostic history
when learning-rate adaptation is disabled.

The oracle uses even population 6 and refreshes the factor every generation.
The verifier asserts matching learning constants, an inactive external learning-rate
cap, strictly positive eigenvalues, and an inactive external sigma cap. Twelve
successive states on each of JVM and optimized Node match independent mean,
covariance, and step-size updates, including a false h-sigma gate. A separate
1D analytic test distinguishes the canonical gate from an off-by-one transient.
Package source identity and explicit adapters are recorded in [comparison.json](comparison.json).
The largest absolute covariance discrepancy in the 24 state checks is below
`1e-13`; the acceptance tolerances are `rtol=2e-11`, `atol=2e-12`.

Sixteen shared tests cover reference updates, analytic optima, rotation and
ill-conditioning, nonsmooth objectives, rank invariance, positive definiteness,
ownership, seeded repetition, invalid/stale/foreign receipts, incomplete budgets,
callback and control failures, fixed/bounded coordinates, sampling exhaustion,
unrepresentable steps, IPOP budgets, and progress vetoes at restart boundaries.
Six verifier mutation probes reject corrupted endpoints, work counts, statuses,
duplicate rows, and an incorrect step-size update.

## Multi-seed comparison with pycma

[Detailed comparison](comparison.md) records eight analytic fixture families in
4 or 8 dimensions, 20 seeds per case, and three timed repetitions per seed on
each platform: 480 timed solves each for JVM, optimized Node, and Python. Five
untimed solves per case use distinct warmup seeds. Each solve receives 15000
objective evaluations and at most 10000 generations. The common independent
clean-objective error target is `1e-8`; feasibility is exact. The noisy sphere
uses multiplicative lognormal noise and requests an observed `1e-9` target; its
success is assessed against the noiseless function at `1e-8`.

All finite endpoint, reported-objective, feasibility, receipt-budget, seeded
repetition, and recurrence checks pass. Target attainment is a separate result:

| Family | Gale JVM | Gale Node | Python |
| --- | ---: | ---: | ---: |
| Sphere, rotated ellipsoid, nonsmooth, noisy sphere, interior box, boundary box | 20/20 each | 20/20 each | 20/20 each |
| Rosenbrock | 18/20 | 18/20 | 19/20 |
| Rastrigin | 2/20 | 2/20 | 0/20 |

No restarts are used in this benchmark. Multimodal failures remain visible;
these seed counts do not establish comparative global-search superiority.
The API's restart policy has separate contract tests.

On the six families where every seed succeeds in all implementations, the
Python-to-Gale median solve-time ratios range from about **9.5 to 54** on JVM
and **9.9 to 24** on Node. Evaluation counts are similar, and are slightly higher
for Gale's median JVM run on those families. These are local full-solve timings,
including solver construction, sampling, callbacks, rejected bound proposals,
and final diagnostics/results. Data, callback, bounds and configuration preparation
are outside timing. The Python callback uses NumPy; scalar independent verification
is outside timing. Thus the comparison includes differences in callback and
runtime overhead as well as the solver implementation. No isolated-kernel speedup
or general evaluation-efficiency advantage is claimed.

Python uses `cma==4.4.4`, with active negative updates, mirrors, diagonal decoding,
and automatic conditioning alleviation disabled. Positive weights, population,
starts, initial scale, budgets and targets match Gale. External whole-vector
rejection matches the bounded policy. Covariance refresh and condition checks
are performed after each complete update in both implementations. Python retains
its own CSA constants, expected-normal-norm approximation, h-sigma threshold,
step-size safeguards, and initial covariance perturbation. This is a comparison
of closely configured standard implementations, not identical trajectories or
a benchmark of pycma's default active variant.

RNGs and noise streams are independently seeded; equal seed labels do not mean
equal draws across libraries. Repetition of a seed within each runtime is checked
exactly. Tables condition median time/evaluations on success; the Rosenbrock and
Rastrigin success sets differ and their timing ratios are not paired speedups.
The harness checks reported evaluation limits; independent callback-count tests
cover the minimize accounting, including failures and cancellation.

Hardware is Apple M3 Max / 36 GiB, using the same local machine as the earlier
optimization portfolio measurements. Runtimes are JDK 21.0.12.1, Node 22.18.0,
Scala 3.7.4, Python 3.14.7, NumPy 2.4.3, cma 4.4.4, and cmaes 0.12.0.
One-thread environment settings were requested. `threadpoolctl` did not enumerate
a BLAS pool on this NumPy build, so the actual BLAS thread count was not independently
confirmed. Runs are sequential and not OS-isolated. The raw repetitions retain timing variation; these results
are scoped to these fixtures and this host.

## Portable gates and reproduction

All required gates passed: formatting, `compileAll`, `testAllFull` (899 core JVM,
886 core JS, 54 laws per platform), 85 Breeze parity tests, 29 Breeze interop
tests, `docsCheck`, `benchCompile`, and all 16 CMA-ES tests under `FullOptStage`.
Existing unrelated documentation warnings remain. Complete commands and exit
metadata are retained in [logs](logs/final-gates-08.log).

With JDK 21, Node 22, sbt, and an isolated Python environment:

```sh
python3 -m venv /tmp/gale-cma-python
/tmp/gale-cma-python/bin/pip install -r tools/optim-cmaes/requirements.txt
python3 tools/optim-cmaes/run_comparison.py \
  --python /tmp/gale-cma-python/bin/python --output /tmp/gale-cma-rerun \
  --warmups 5 --seeds 20 --repeats 3
/tmp/gale-cma-python/bin/python tools/optim-cmaes/check_verifier.py /tmp/gale-cma-rerun
```

The runner writes to the supplied output directory and preserves checked-in
evidence. [Source and artifact hashes](source-manifest.json) identify the tested
uncommitted snapshot. Earlier exploratory logs retain an incorrect seed-test type
assertion and a Scala.js runner argument error, both fixed before qualification.
No commit/push, hosted CI, downstream migration, native backend, or Wasm
qualification is included in this result.
