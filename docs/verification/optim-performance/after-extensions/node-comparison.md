# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.199 | 0.152 | 1.312 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.450 | 0.264 | 1.706 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.074 | 0.252 | 0.295 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.075 | 0.080 | 0.929 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.059 | 0.071 | 0.829 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.199 | 0.272 | 0.730 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 136.874 | 12.238 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 1.769 | 0.935 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.335 | 0.811 | 0.413 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.380 | 0.360 | 1.054 |
