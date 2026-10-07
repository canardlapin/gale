# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.186 | 0.104 | 1.790 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.572 | 0.173 | 3.308 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.055 | 0.176 | 0.314 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.064 | 0.055 | 1.171 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.064 | 0.048 | 1.324 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.158 | 0.197 | 0.806 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 148.995 | 9.230 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.533 | 0.691 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.214 | 0.573 | 0.373 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.266 | 0.266 | 1.000 |
