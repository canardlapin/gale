# Gale 0.1 API stability boundary

This document identifies the public API that becomes compatible at
`v0.1.0-M1`. It is an inventory and review boundary, not a substitute for
Scaladoc or executable examples. The binary and TASTy baselines must be
generated from the published M1 artifacts after that immutable version exists.

## Admitted `gale-core` packages

| Package | Compatible contract from M1 |
| --- | --- |
| `gale.linalg` | dense values and views, builders, shapes, operators, factorizations, solve and least-squares entry points, typed errors, destinations, and workspaces |
| `gale.sparse` | COO/CSR/CSC and structured matrices, canonicalization, Matrix Market IO, compressed patterns, and symbolic replay plans |
| `gale.solvers` | iterative solver options, results, diagnostics, preconditioners, convergence semantics, and reusable workspaces |
| `gale.spectral` | dense and partial decompositions, typed selection, result ordering, convergence and extremality diagnostics, generalized operators, and explicit metric-solve contracts |
| `gale.optim` | first-order optimization problems, options, typed failures, diagnostics, and constrained-Rayleigh helpers |
| `gale.sized` | optional compile-time sized wrappers and their checked conversion boundary |
| `gale.backend` | caller-visible capabilities, backend selection, fallback behavior, configuration, thresholds, and factorization-provider contracts |
| `gale.syntax` | the `all` and `unicode` opt-in extension modules, including `zipMapExact`, pointwise operations, matrix product, and dot aliases |

The `gale.kernel` implementation and platform-storage operations are
`private[gale]` and are not public extension points. JVM Vector, native, FFM,
and Breeze modules are separate provisional artifacts; their APIs acquire no
0.1 compatibility promise until an admission record names a baseline.

## Observable semantic contracts

Compatibility covers more than method descriptors:

- `DVec` and `DMat` results own their returned mutable storage unless a method
  is explicitly documented as a view.
- Builders, destinations, and workspaces are single-owner mutable resources.
  Their methods do not expose Gale backing arrays.
- Total factorization and solve entry points return `Either[LinAlgError, A]`.
  Throwing conveniences remain explicitly named.
- Iterative and spectral results distinguish a returned approximation from
  residual convergence and from certification of a requested global extreme.
- Sparse canonicalization, duplicate handling, stored-zero behavior, and
  structure-preserving replay retain their documented meanings.
- Selection chooses spectral membership; result types retain their documented
  canonical ordering independently of the selected backend.
- Backends may change floating-point association but must satisfy the same
  shape, validation, ownership, residual, and conformance contracts.
- Work and allocation counters retain their units and inclusion rules. A field
  must not silently change from measured work to an estimate or proxy.

An incompatible change to one of these contracts after M1 requires `0.2.0`,
even when a binary compatibility checker cannot detect it.

## Published `gale-laws` API

The M1 law module admits these reusable entry points:

- `VecLaws` and `MatrixLaws`;
- `SparseLaws` and `SolverLaws`;
- `SpectralLaws` and `GeneralizedOperatorLaws`;
- `SpectralBackendLaws`; and
- `BackendConformanceSuite`.

Their public method signatures and assertion meanings become part of the 0.1
contract. Test fixtures and private helpers do not.

## Pre-M1 API-diff receipt

Seven direct consumers still pin Gale commit
`d55fe2f97196a76ab7879e1a12f1e92403aeba06`; regress4s pins
`9c27cb5150951f97142789123594244400d971fc`. The initial 0.1 audit compared the
older common pin with `origin/main` commit
`180aa722be5c58de3dbef82a2cf0f2ec9b61187e`:

```sh
git diff --stat \
  d55fe2f97196a76ab7879e1a12f1e92403aeba06..180aa722be5c58de3dbef82a2cf0f2ec9b61187e \
  -- core/shared/src/main/scala laws/shared/src/main/scala
```

The comparison changes five core files. It adds allocation-controlled and
row-scaled QR construction and solve APIs. It also contains two deliberate
source migrations:

1. `DMatBuilder.updateRowMajor` became `writeLinear`, matching the ecosystem
   convention that sequential storage access is named explicitly.
2. `gale.syntax.all.zipMap` became `zipMapExact`, making its equal-shape
   requirement visible in the name.

The first migration still has callers in multivar and ScalaFIM. Those callers
must move to `writeLinear`; Gale will not restore a permanent alias. No direct
consumer was found using the removed Gale `zipMap` name.

Before M1, replace the comparison target with the exact candidate commit,
review every exported addition and removal, and retain the complete diff in the
release evidence. The build's pre-M1 `compatibilityCheck` has no previous
version and therefore checks task wiring rather than compatibility. After M1
publication, every admitted coordinate must run both automated checks against
the exact milestone before another 0.1-line change is accepted:

```sh
sbt -Dgale.compatibility.baseline=0.1.0-M1 compatibilityCheck
```

The corresponding CI job must become blocking at the same time. Binary and
TASTy checks complement, but do not replace, review of the semantic contracts
listed above.
