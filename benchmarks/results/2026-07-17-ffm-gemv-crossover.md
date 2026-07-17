# FFM GEMV crossover — Accelerate / JDK 22 (2026-07-17)

## Environment and protocol

- Apple ARM64, macOS 14.3
- OpenJDK 22 and runtime-discovered macOS Accelerate
- public `DMat * DVec` facade; both routes allocate a result and the FFM route
  includes all heap/native copies
- JMH 1.37, average time, one thread, 2 forks, 5 x 500 ms warmups, and
  8 x 500 ms measurements

## Result

Times are microseconds per operation. A speedup below 1 means FFM is slower.

| n | Pure Gale | FFM Accelerate | Native speedup |
| ---: | ---: | ---: | ---: |
| 64 | 0.853 | 1.708 | 0.50x |
| 128 | 3.540 | 8.330 | 0.42x |
| 256 | 14.605 | 22.798 | 0.64x |
| 512 | 64.479 | 157.623 | 0.41x |
| 1024 | 275.982 | 923.966 | 0.30x |
| 2048 | 1,148.785 | 1,354.544 | 0.85x |

There is no measured heap-copy-inclusive crossover through n=2048. Accelerate
GEMV therefore keeps `nativeGemvMinWork = Long.MaxValue`; OpenBLAS, MKL, and
reference-library defaults also remain disabled until family-specific evidence
exists. The Vector backend is copy-free and has its own separately measured
`128^2` GEMV threshold.

`dot` is deliberately not a backend-routed public seam: it is frequently called
inside algorithms where dispatch overhead would recur in the inner loop. It has
no native threshold. This is a contract decision, not a missing measurement.
