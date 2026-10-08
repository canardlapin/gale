# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.235 | 0.109 | 2.149 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.240 | 0.222 | 1.083 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.071 | 0.189 | 0.376 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.133 | 0.057 | 2.320 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.037 | 0.051 | 0.727 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.104 | 0.219 | 0.473 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 62.067 | 9.723 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.023 | 0.745 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.177 | 0.602 | 0.295 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.189 | 0.280 | 0.674 |
