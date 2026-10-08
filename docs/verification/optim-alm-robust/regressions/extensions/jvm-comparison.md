# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.202 | 0.154 | 1.312 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.912 | 0.277 | 3.292 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.117 | 0.255 | 0.460 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.323 | 0.083 | 3.884 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.091 | 0.072 | 1.253 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.276 | 0.340 | 0.813 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 85.995 | 12.414 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 0.555 | 1.036 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.133 | 0.713 | 0.187 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.260 | 0.358 | 0.727 |
