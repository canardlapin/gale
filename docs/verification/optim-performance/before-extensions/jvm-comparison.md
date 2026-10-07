# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.245 | 0.104 | 2.354 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 1.166 | 0.173 | 6.749 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.118 | 0.176 | 0.672 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.130 | 0.055 | 2.383 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.107 | 0.048 | 2.209 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.303 | 0.197 | 1.542 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 134.860 | 9.230 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 0.950 | 0.691 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.209 | 0.573 | 0.365 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.202 | 0.266 | 0.760 |
