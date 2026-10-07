# SciPy comparison

Independent scalar-loop endpoint checks; time ratio is Gale/Python.

| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |
|---|:---:|:---:|---|---:|---:|---:|---:|
| lm-exp-64 | True | True | Converged | 9.270e-12 | 0.162 | 0.110 | 1.472 |
| lm-exp-offset-256 | True | True | Converged | 3.825e-10 | 0.349 | 0.181 | 1.931 |
| lm-rosenbrock | True | True | Converged | 5.676e-10 | 0.076 | 0.185 | 0.412 |
| lm-scaled-linear-8 | True | True | Converged | 4.041e-11 | 0.135 | 0.059 | 2.307 |
| lm-rank-deficient-3 | True | True | Converged | 9.481e-10 | 0.111 | 0.051 | 2.147 |
| box-diagonal-16 | True | True | Converged | 1.493e-11 | 0.130 | 0.204 | 0.637 |
| box-diagonal-512 | True | False | Converged | 9.664e-09 | 136.256 | 9.737 | unqualified |
| box-rotated-32 | True | False | Converged | 9.078e-09 | 0.896 | 0.736 | unqualified |
| box-rosenbrock-2 | True | True | Converged | 0.000e+00 | 0.202 | 0.601 | 0.336 |
| box-logistic-1024-16 | True | True | Converged | 4.489e-09 | 0.191 | 0.278 | 0.687 |
