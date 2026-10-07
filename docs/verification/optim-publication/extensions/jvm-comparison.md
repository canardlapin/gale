# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.223 | 0.153 | 1.458 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.677 | 0.253 | 2.673 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.140 | 0.250 | 0.560 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.088 | 0.081 | 1.087 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.074 | 0.071 | 1.044 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.195 | 0.269 | 0.724 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 90.430 | 11.673 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 0.611 | 0.993 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.139 | 0.810 | 0.172 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.877 | 0.350 | 2.508 |
