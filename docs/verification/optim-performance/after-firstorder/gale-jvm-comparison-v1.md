# Gale versus Python optimization comparison

9/9 Gale harness rows meet independently recomputed accuracy and feasibility checks.

Timing ratio is Gale/Python; below 1 favors Gale. Near ties need more evidence.

| case | method | status | objective error | stationarity/KKT | valid | Gale ms | Python ms | ratio |
|---|---|---|---:|---:|:---:|---:|---:|---:|
| diag-quadratic-n16-k1e3 | LBFGS | Converged | 1.295e-15 | 8.901e-08 | True | 0.594 | 4.294 | 0.138 |
| diag-quadratic-n512-k1e3 | LBFGS | Converged | 6.920e-15 | 7.376e-08 | True | 12.358 | 20.788 | 0.594 |
| box-diag-quadratic-n16-k1e3 | AcceleratedProximalGradient | Converged | 0.000e+00 | 9.143e-08 | True | 2.449 | 0.304 | 8.044 |
| rotated-quadratic-n32-k1e3 | LBFGS | Converged | 1.242e-15 | 4.898e-08 | True | 0.777 | 6.030 | 0.129 |
| rosenbrock-n2 | LBFGS | Converged | 9.333e-18 | 6.645e-08 | True | 0.030 | 1.249 | 0.024 |
| rosenbrock-n32 | LBFGS | Converged | 2.359e-17 | 5.982e-08 | True | 0.493 | 6.944 | 0.071 |
| logistic-l2-m1024-n16 | LBFGS | Converged | 1.110e-15 | 1.005e-08 | True | 0.244 | 0.372 | 0.655 |
| lasso-m256-n32 | AcceleratedProximalGradient | Converged | 1.454e-14 | 9.169e-08 | True | 0.478 | 0.311 | 1.535 |
| lasso-m1024-n128 | AcceleratedProximalGradient | Converged | 3.686e-14 | 7.684e-08 | True | 3.592 | 1.452 | 2.473 |

## Environments

```json
{
  "fixtures": "/Users/bbuchsbaum/code/scala/gale/docs/verification/optim-python/fixtures-v1.json",
  "harness": "/Users/bbuchsbaum/code/scala/gale/.agent-work/optim-performance/qualified-firstorder/gale-jvm-v1.tsv",
  "python_baseline": {
    "numpy": "2.4.3",
    "platform": "macOS-14.3-arm64-arm-64bit-Mach-O",
    "python": "3.14.7 (main, Aug  5 2026, 10:29:49) [Clang 16.0.0 (clang-1600.0.26.6)]",
    "scikit_learn": "1.9.1",
    "scipy": "1.17.1",
    "threadpools": [
      {
        "filepath": "/private/tmp/gale-optim-python-venv/lib/python3.14/site-packages/sklearn/.dylibs/libomp.dylib",
        "internal_api": "openmp",
        "num_threads": 1,
        "prefix": "libomp",
        "user_api": "openmp",
        "version": null
      }
    ],
    "threads": {
      "OMP_NUM_THREADS": "1",
      "OPENBLAS_NUM_THREADS": "1",
      "VECLIB_MAXIMUM_THREADS": "1"
    }
  },
  "references": "/Users/bbuchsbaum/code/scala/gale/docs/verification/optim-python/reference-results-v1.json"
}
```
