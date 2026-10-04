# Caller-selected cutoff and geometry receipt

Local implementation for Gale tickets:

- `bd-01M42AE4HXZ7122GX40PM4HCQY`: policy-selected pseudo-inverse and
  minimum-norm least squares.
- `bd-01M42AE5WV04ZJP013DEZYPM7E`: RREF and complete sparse/orthonormal null bases.
- `bd-01M42AE7ACJBQ5DJBPWH0H82H4`: basis-backed row/column spaces.

## Contract

`SvdCutoff.Default`, `Relative(tau)` and `Absolute(threshold)` keep exactly
singular values above the resolved threshold. Explicit cutoffs are finite and
non-negative; cutoff equality is discarded. The default preserves the existing
`pinv` threshold. Selected rank does not reuse `SVD.rank`'s separate policy.

`TruncatedSvd` stores retained factors and singular metadata for reuse across
vector/matrix solves, pseudo-inversion and both subspaces. Solves return
coefficients, projected residuals, rank and cutoff. Deliberate truncation solves
the minimum-norm problem for the truncated matrix. Subspaces expose basis
columns and projection, residual and vector-distance operations without square
projectors. Statistical membership/uncertainty decisions remain caller-owned.

RREF has separate contracts for an absolute elimination pivot tolerance and a
spectral cutoff. Spectral rank is selected before elimination. Null bases include
all `columns - rank` directions, including missing wide-economy-SVD directions.
Sparse output uses canonical row-reduced null equations, matching the
null-projector RREF convention in exact arithmetic; it does not promise globally
minimal support. Orthonormal output uses twice-reorthogonalized Gram-Schmidt.
Conversion failures and non-finite inputs are typed errors. No null projector is
constructed; a full null basis may itself require quadratic storage.

## Gale evidence, 2026-10-04

Base: `0c9967f1e35c46608fe3fee3500c0dd80c5ace05`. The candidate is a local working
tree change, not a published revision.

Using Eclipse Temurin JDK 21.0.12.1 and Node 22.18.0:

- `scalafmtCheckAll` and `compileAll` passed, including the Scala 3.8 consumer.
- `testAll`: JVM core 712/712, Scala.js core 702/702, laws 54/54 on each platform.
- `parityTest`: 85/85; `interopBreezeTest`: 29/29.
- `docs/mdoc`: 28 files compiled, zero errors, five existing link warnings.

The normal formatter deliberately has a narrow allowlist. All seven new
Scala files were explicitly formatted and checked with the repository settings
and a temporary include list. The two existing files retain their original
layout outside the API changes. After narrowing this diff, the affected suites
passed 39/39 on both JVM and Scala.js. `docs/mdoc` also passed after the initial
formatting pass. This includes 18
new tests plus the existing 21 pseudo-inverse/numerical-contract tests.
Post-broad-gate changes contain identical non-comment, non-whitespace source
tokens; both pre-format and final source snapshots are retained separately.

New evidence covers analytically known wide minimum-norm coefficients, rotated
truncation, strict cutoff equality, default rank-policy disagreement, scaled
inputs, residual orthogonality, reusable strided matrix responses, zero/full
spaces, a pivot/SVD-rank counterexample, canonical null equations, complete wide
null spaces and seeded integer low-rank plants. This is numerical correctness
evidence for these contracts, not a measured performance claim.

Broad log: `/private/tmp/gale-cutoff-gates-01.log`.
Final focused/format log: `/private/tmp/gale-cutoff-final-narrow-check.log`.
The intermediate formatting/docs log is
`/private/tmp/gale-cutoff-final-format-check.log`. These have `.meta.json`
sidecars containing full commands and exit status zero.
Existing Scaladoc warnings and optional native-provider fallback warnings remain
in the broad log. `compileAll` publishes to the local Ivy cache for its consumer
check; no remote publication occurred.

## ScalaFIM helper-removal rehearsal

ScalaFIM source was archived from
`53097f3f43fa125cd166a99268f192f4a5b12311` into an isolated directory. The
[migration patch](cutoff-geometry-scalafim-migration.patch) replaces private
solve, row-space and residualizer kernels and deletes the private Gauss-Jordan
helper. Alias generation reuses retained factors and converts typed Gale
failures into ScalaFIM numerical failures. Application-specific alias display
cleanup and uncertainty metadata stay in ScalaFIM.

With `-Dscalafim.gale.build=/Users/bbuchsbaum/code/scala/gale`, the existing
`DesignDiagnosticsSuite`, `ContrastDiagnosticsSuite` and
`WeightedContrastDiagnosticsSuite` passed **30/30 on JVM and 30/30 on Scala.js**.
These include frozen independent diagnostics fixtures, rank-deficient and wide
designs, caller-selected rank policies, alias counts/residual bounds, F-row
orthonormalization and residualized-information oracles. The archive's lack of
Git metadata produces version-discovery diagnostics; compilation and tests exit
successfully. This is an isolated consumer rehearsal, not a live pin upgrade.

Consumer log: `/private/tmp/gale-scalafim-cutoff-check-01.log` and its
`.meta.json` sidecar. Patch SHA256:
`2009fd0088cea8d4935532d02c425b90764fd365534918efa3a468068b22787d`.

Immutable source snapshots and complete file-hash manifests are retained under
`/private/tmp/gale-scalafim-cutoff-fqgy3x1x/`:

- `gale-tested-source.tar.gz`, `gale-tested-source-manifest.json` bind the broad
  checks and consumer rehearsal.
- `gale-final-source.tar.gz`, `gale-final-source-manifest.json` bind the final
  formatting and focused checks.
- `source.tar` is the original ScalaFIM archive;
  `scalafim-migrated-source.tar.gz` is the tested migration candidate.

The live ScalaFIM checkout and its Gale pin remain unchanged. Adoption requires
a published Gale revision and the narrow consumer patch plus pin upgrade. The
Gale implementation is ready for review; these tickets are not claimed as a
completed live ScalaFIM migration.
