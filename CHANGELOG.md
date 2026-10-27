# Changelog

## Unreleased — 0.2

### Added

- Explicit row/column margins with `Margin.ByRow` / `Margin.ByCol`, `byRow` /
  `byCol` reductions, and zero-copy `DVec.along(axis)` tagging. Plain vectors
  retain their existing matrix-vector meaning.
- Matrix-margin arithmetic, scalar pointwise operators, and one-pass `affine`
  and subtract-then-divide `standardize` transforms with checked builder
  destinations. Standard-deviation estimation is deferred; scales are supplied.

### Changed

- Development on `main` uses 0.2 snapshots. The `release/0.1` branch starts at
  `v0.1.0-M1` for maintenance and stabilization of the published line.
- No-Git development versions retain the terminal `-SNAPSHOT` marker so they
  are classified correctly by the build.

The current published version remains **0.1.0-M1**. A main-branch merge does not
publish a release; 0.2 publication requires a separate release-admission change.

## 0.1.0-M1

First published milestone: `gale-core` and `gale-laws` for Scala 3 on the JVM
and Scala.js, under `io.github.canardlapin`.

See the [release notes](https://github.com/canardlapin/gale/releases/tag/v0.1.0-M1)
and [release manifest](docs/release-manifest.md) for the admitted scope.
