# Augmented Lagrangian comparison

Median full-solve milliseconds. **FAIL** means at least one endpoint failed the common accuracy gate; its time is not a qualified speed result.

| Case | JVM | Node | SciPy SLSQP | SciPy trust-constr |
|---|---:|---:|---:|---:|
| equality-diagonal-8-s0 | 0.114 | 0.218 | 0.559 | 9.091 |
| equality-diagonal-8-s1 | 0.112 | 0.207 | 0.544 | 11.009 |
| equality-dense-32-s0 | 1.025 | 2.403 | 2.337 | 45.326 |
| equality-dense-32-s1 | 1.076 | 2.344 | 2.540 | 43.710 |
| equality-dense-128-s0 | 1.696 | 3.672 | 11.070 | 19.700 |
| equality-dense-128-s1 | 1.450 | 3.026 | 18.124 | 18.790 |
| mixed-box-3-s0 | 0.326 | 0.464 | 0.145 | 4.829 |
| mixed-box-3-s1 | 0.331 | 0.368 | 0.145 | 16.952 |
| ball-8-s0 | 0.053 | 0.062 | 0.273 | 3.814 |
| ball-8-s1 | 0.058 | 0.085 | 0.250 | 5.361 |
| ball-64-s0 | 0.293 | 0.434 | 0.422 | 3.978 |
| ball-64-s1 | 0.293 | 0.563 | 0.404 | 6.492 |
| parabola-2-s0 | 0.040 | 0.057 | 0.386 | 3.150 |
| parabola-2-s1 | 0.053 | 0.079 | 0.264 | 2.491 |
| parabola-16-s0 | 0.105 | 0.204 | 0.456 | 2.810 |
| parabola-16-s1 | 0.161 | 0.327 | 0.494 | 4.830 |
