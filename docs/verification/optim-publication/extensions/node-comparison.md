# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.191 | 0.153 | 1.250 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.448 | 0.253 | 1.771 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.075 | 0.250 | 0.302 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.068 | 0.081 | 0.846 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.060 | 0.071 | 0.846 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.187 | 0.269 | 0.696 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 138.041 | 11.673 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.701 | 0.993 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.317 | 0.810 | 0.391 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.384 | 0.350 | 1.098 |
