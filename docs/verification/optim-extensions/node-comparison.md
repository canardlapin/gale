# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.181 | 0.110 | 1.637 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.594 | 0.181 | 3.286 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.058 | 0.185 | 0.311 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.063 | 0.059 | 1.077 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.068 | 0.051 | 1.318 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.148 | 0.204 | 0.723 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 153.688 | 9.737 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.668 | 0.736 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.238 | 0.601 | 0.396 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.281 | 0.278 | 1.008 |
