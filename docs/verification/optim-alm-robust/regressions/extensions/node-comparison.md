# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.196 | 0.154 | 1.275 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.460 | 0.277 | 1.661 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.077 | 0.255 | 0.302 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.072 | 0.083 | 0.866 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.058 | 0.072 | 0.797 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.215 | 0.340 | 0.634 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 140.380 | 12.414 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.863 | 1.036 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.333 | 0.713 | 0.467 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.376 | 0.358 | 1.051 |
