# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.340 | 0.152 | 2.240 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.574 | 0.264 | 2.174 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.120 | 0.252 | 0.478 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.189 | 0.080 | 2.349 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.063 | 0.071 | 0.890 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.173 | 0.272 | 0.635 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 87.265 | 12.238 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 0.795 | 0.935 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.201 | 0.811 | 0.248 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.860 | 0.360 | 2.388 |
