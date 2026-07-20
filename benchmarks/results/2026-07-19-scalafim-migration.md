# Scalafim migration workload receipt

Date: 2026-07-19

Machine: Apple ARM64

Runtime: OpenJDK 25.0.1, pure Gale backend

Harness: JMH 1.37, one fork, one 250 ms warmup, two 250 ms measurements

This is a focused development receipt, not a release-grade performance claim.
It exercises the new batched factorization and operator APIs at an fMRI-like
shape: a `240 x 24` design and 128 simultaneous response columns.

```bash
sbt 'benchmarksJVM/Jmh/run -i 2 -wi 1 -f 1 -r 250ms -w 250ms -prof gc gale.bench.ScalafimMigrationJmh.*'
```

| Operation | Time | Allocation |
| --- | ---: | ---: |
| Cholesky solve, 128 RHS | 51.852 us/op | 24,673 B/op |
| `X.t` operator applied to 128 columns | 466.986 us/op | 32,909 B/op |
| Pivoted QR solve, 128 RHS | 1,832.103 us/op | 270,530 B/op |
| QR `Q.t` transform, 128 RHS | 1,968.616 us/op | 245,910 B/op |

The Cholesky and QR-transform measurements are effectively the size of their
owned matrix results (`24 * 128 * 8 = 24,576` bytes and
`240 * 128 * 8 = 245,760` bytes). The operator route adds only lightweight
column-view objects to its `24 x 128` result. The pivoted QR solve intentionally
owns both the transformed `240 x 128` work matrix and the returned `24 x 128`
coefficient matrix, matching the observed allocation order.
