# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.146 | 0.109 | 1.336 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.330 | 0.222 | 1.487 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.062 | 0.189 | 0.328 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.054 | 0.057 | 0.933 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.044 | 0.051 | 0.848 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.136 | 0.219 | 0.621 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 99.672 | 9.723 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.310 | 0.745 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.246 | 0.602 | 0.408 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.281 | 0.280 | 1.003 |
