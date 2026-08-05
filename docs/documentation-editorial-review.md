# Gale documentation editorial review

> **Historical record:** This review describes the repository before the 0.1
> release-line decision. Current version, artifact, and compatibility claims
> live in `README.md`, `docs/release-policy.md`, and
> `docs/release-manifest.md`.

Evidence date: 2026-08-03. Audited revision: `e7c8579` on `main` before the
documentation changes described here. This is an internal editorial and API
review; `docs/user/` is the separate public site input.

## Research brief

### Readers and job

The primary reader is a Scala 3 developer who needs real-`Double` linear
algebra for scientific, statistical, or data software on the JVM and possibly
Scala.js. They know matrices and least squares but should not need to know
Gale's source layout, kernel architecture, or backend hierarchy.

Secondary readers are library authors migrating a bounded Breeze workload,
developers building matrix-free algorithms, and experts optimizing a measured
allocation or native-compute hot path.

Gale's main job is to provide one portable mathematical API for dense and
sparse values, factorizations, iterative solvers, and selected spectral
methods, while making numerical failure, convergence, ownership, and optional
acceleration visible.

### Three ideas to remember

1. Ordinary values are immutable-facing, and structural numerical failures are
   typed `LinAlgError` values on total entry points.
2. Matrices and `DoubleLinearOperator`s share solver and spectral contracts, so
   an algorithm need not materialize a matrix it only needs to apply.
3. The portable JVM/Scala.js route is the default. Mutable destinations,
   workspaces, Vector API, and native BLAS/LAPACK are explicit choices with
   narrower ownership and evidence boundaries.

### Core abstractions and consequences

- `DVec` and `DMat` are immutable-facing dense values. Transpose, row, column,
  and slice can share immutable storage; gathers and ordinary results own their
  data.
- `CSR`, `CSC`, and related sparse structures distinguish numeric updates from
  structural changes such as canonicalization or pruning.
- `DoubleLinearOperator` represents application without entry access. Some
  algorithms additionally require a transpose action, symmetry evidence, a
  positive-definite metric, a preconditioner, or a genuine solve.
- `LU`, `Cholesky`, and `QR` retain factors and diagnostics. Square solves and
  least squares are distinct capabilities.
- `DenseWorkspace`, builders, and sparse destinations are sequential,
  single-owner mutable resources. Returned ordinary results remain owned.
- Optional backends change eligible execution routes, not mathematical result
  types, error semantics, or ownership contracts.

### Canonical workflows and failures

The first complete workflow is a tall dense regression: construct a design and
observations, call `leastSquares`, inspect the owned coefficients, and keep the
`Either` until the application boundary. The ordinary follow-on workflows are
dense systems, sparse/operator iterative solves, and dense or partial spectral
analysis.

Expected failures include dimension mismatch, singularity, non-positive
definiteness, rank deficiency, unsupported selection, and non-finite operator
output. Iteration exhaustion is often a diagnostics-bearing result rather than
a structural `Left`; documentation must not flatten those cases together.

### Maturity, compatibility, and non-goals

The live build is `1.0.0-SNAPSHOT`, targets Scala 3.7.4, and has no configured
public binary publication. The intended `io.github.canardlapin` artifacts are
therefore names, not currently installable Maven Central coordinates. JDK 21 is
the core JVM line; finalized FFM modules require JDK 22. Node is the executed
Scala.js runtime in CI. Browser support is an intended public route, but the
current release policy correctly says browser CI is not yet a release gate.

Gale is not a source-compatible Breeze replacement, a generic scalar/tensor
library, a general complex matrix package, or a sparse-direct factorization
implementation.

## Reader question ladder

1. What is Gale, and does it fit my workload?
2. How can I use it before public artifacts exist?
3. Can I build a matrix, solve a real problem, and understand failure?
4. What is copied, what is a view, and when do dimensions fail?
5. Should I use a dense matrix, sparse matrix, or operator?
6. Which solve, factorization, or spectral route matches my problem?
7. How do I inspect convergence and numerical quality?
8. How do I reduce allocation without weakening result ownership?
9. When does an optional backend help, and what evidence scopes that claim?
10. How do I recover from a common compile, shape, rank, convergence, or
    runtime problem?

## Evidence ledger

| Claim | Status | Evidence | Scope or caveat | Documentation consequence |
| --- | --- | --- | --- | --- |
| Core dense, sparse, solver, and spectral APIs are shared by JVM and Scala.js | confirmed | `core` cross-project; `testAllFull` at `e7c8579` | Node executes JS tests; browser matrix is not a release gate | Say JVM and Scala.js; do not say browser-certified |
| Gale is not publicly released | confirmed | `1.0.0-SNAPSHOT`; no publish destination; no matching Central result in this audit | External state can change later | Teach source dependency and `publishLocal`; label future coordinates |
| Total solve/factorization entry points preserve typed failures | confirmed | public signatures and factorization/solver tests | primitive arithmetic validates by throwing | Explain the boundary, not “everything returns Either” |
| Ordinary returned dense values remain stable across workspace reuse | confirmed | ownership Scaladoc; workspace, QR, and mutable-builder suites | explicitly `unsafe` borrowed views are exceptions | Put preservation next to allocation-control examples |
| Matrix-free partial generalized symmetric solving is implemented | confirmed | LOBPCG/Lanczos source, focused suites, executable guide | only documented algebraic-end selections; convergence and extremality differ | Teach engine choice and diagnostics separately |
| Optional backends are universally faster | unverified and false as a general claim | dashboard contains enabled and rejected routes | results are machine, JDK, library, shape, and layout specific | Route readers to measured thresholds; avoid “fast” as an adjective |
| Current source is binary-compatible across future 1.x releases | unverified | no stable baseline or MiMa gate before 1.0.0 | policy is prospective | Call the current build pre-release; do not imply a live compatibility promise |

## Editorial findings

### Strengths

- The site already has the right technical stack: a dedicated
  `sbt-typelevel-site` project, mdoc examples, Laika navigation, JVM/JS
  Scaladoc, and a CI `docsCheck` gate.
- Public and internal documentation are physically separated: only
  `docs/user/` is rendered.
- The numerical, ownership, sparse-pattern, and generalized-eigensolver pages
  state unusually strong semantic boundaries.
- Tests expose real workflow and ownership contracts, including matrix RHS,
  operator diagnostics, strided inputs, workspace reuse, and typed failures.

### Problems found

- The README begins with kernel architecture and a module inventory before a
  meaningful result. A newcomer has no quick success on the first screen.
- `guides/examples.md` is a 600-line package tour with competing jobs. Several
  fences are illustrative rather than mdoc-checked.
- That page contains obsolete `io.gale` coordinates and incorrectly says
  `CgWorkspace.solution` is an alias; the live API returns a copy and names the
  borrowed route `unsafeSolutionView`.
- `advanced/numerical-contract.md` describes generalized block Lanczos as a
  future gated follow-up even though the public implementation and guide exist.
- `docs/release-policy.md` names the obsolete `io.gale` organization while the
  build and intended artifacts use `io.github.canardlapin`.
- Recent pivoted QR, matrix-RHS workspace solves, consuming builders, and
  row-scaled factorization APIs have strong tests but no coherent user guide.
- Troubleshooting is absent. Users must infer the difference between
  structural `Left`, thrown primitive precondition failures, non-converged
  iterative results, and unsupported native/runtime routes.
- The README links internal acceptance and release-evidence records as if they
  were normal user documentation.

## Documentation architecture

```text
README
Public guide
├── Overview
├── Getting started
├── Core concepts
├── Guides
│   ├── Choose a guide
│   ├── Worked-example map
│   ├── Dense systems and least squares
│   ├── Sparse matrices and operator solves
│   ├── Spectral analysis
│   ├── First-order composite optimization
│   └── Breeze migration
├── Advanced
│   ├── Numerical and backend contract
│   ├── Generalized operator eigenproblems
│   ├── Dense destinations
│   ├── Workspaces
│   ├── Sparse patterns and plans
│   └── Ownership
├── Reference
│   ├── Modules and platforms
│   └── Factorization capabilities
└── Troubleshooting
```

| Page | One question | Reader state | Proof or example | Next route |
| --- | --- | --- | --- | --- |
| README | Is Gale relevant, and what does using it feel like? | recognition | smallest least-squares workflow | getting started |
| Overview | Where should I begin? | recognition | site status and route map | getting started |
| Getting started | How do I reach one checked result? | first success | executable regression | core concepts |
| Core concepts | Which rules explain later behavior? | mental model | value/error/operator/backend consequences | task guides |
| Worked-example map | Which complete task matches mine? | ordinary work | links, not duplicated snippets | chosen guide |
| Dense systems and least squares | How do I solve dense and repeated fitted systems? | ordinary work | executable solve, pivoted QR, matrix RHS, row scaling | workspaces/numerical contract |
| Sparse matrices and operator solves | How do I solve without a dense matrix? | ordinary work | executable CSR and matrix-free CG | sparse advanced pages |
| Spectral analysis | Which dense or partial spectral route should I use? | ordinary work | executable symmetric eigen and SVD | generalized guide |
| Troubleshooting | Why did this operation fail or not converge? | recovery | symptom/cause/response table | relevant contract/reference |
| Modules and platforms | Which artifact/runtime would I choose? | reference | live build topology | API docs/release policy |

The large example inventory is split because dense solving, sparse/operator
work, and spectral interpretation have different prerequisites and failure
modes. Internal audits, architecture records, and benchmark receipts remain
outside the public renderer; public pages link only to the user consequence or
the qualified benchmark dashboard when needed.

## API friction ledger

| Priority | Workflow | Evidence | Reader expectation and current cost | Recommendation | Compatibility | Documentation now |
| --- | --- | --- | --- | --- | --- | --- |
| P1 | Install Gale | no public artifacts; source `ProjectRef` differs for JVM/JS | a normal dependency line; instead the reader must understand sbt source dependencies | configure and publish the intended artifacts before announcing 1.0 | release/process change, no source API change | show a complete source dependency and label future coordinates |
| P1 | Weighted/reweighted least squares | `qrScaledRows` and `solveLeastSquaresScaledRowsWith` require the same algebraic scales | a method named for weights; current API correctly accepts row multipliers, which may be signed | consider a separately named `weightedLeastSquares` façade using validated non-negative weights and documented square roots | additive API; semantic design required | call current values row scales, never statistical weights |
| P2 | Allocation-controlled QR solve | factor workspace and solve workspace share `DenseWorkspace` but have separate requirements | one obvious “workspace for this pipeline” constructor | consider a composed `qrPipelineRequirement(rows, cols, rhsCols, options)` convenience | additive API | show `alternative` composition and explain result ownership |
| P2 | Factor a filled transient matrix | `DMatBuilder.consumeQR` closes the builder, unlike `result()` followed by `qr` | a discoverable consuming operation near construction examples | retain the API; consider consistent `consume*` naming if more consuming factorizations appear | additive only if extended | state closure and storage transfer beside the call |
| P2 | Choose portable versus backend-routed QR | `qr`, `qr(options)`, and `qrWith` have deliberately different dispatch promises | overloads that differ only by options/workspace are expected to share routing | consider explicit Scaladoc cross-links or names if more routes appear; do not hide the distinction | renaming would be source-breaking | explain: default may route; explicit policy/workspace pins portable semantics |
| P3 | Discover API reference | site builds local JVM/JS Scaladoc but has no published API URL | a stable API link in site navigation | publish versioned Scaladoc with the first artifact release and configure `tlSiteApiUrl` | publication/site change | give local `docsCheck` route and source links; do not invent a live API URL |

## Final pre-edit assessment

Before this revision, Gale's documentation was technically rich but did not
feel announcement-ready. The release blockers were not a lack of pages; they
were an unpublished installation path, stale contradictory public claims, an
inventory-shaped main guide, no recovery path, and no coherent explanation of
the recently added QR allocation-control surface. Those are maturity gaps, with
public artifact availability and a stable API/reference URL remaining
release-blocking for a broad ecosystem announcement.

## Post-revision assessment

The revised documentation is ready for a pre-release, source-consumer audience.
The README now reaches a meaningful result before architecture detail; the
public site follows the reader-question ladder; ordinary dense, sparse/operator,
and spectral workflows are executable; and allocation-control features state
their ownership and dispatch consequences where readers encounter them.

It is not yet honest to present Gale as broadly 1.0-announcement-ready. The
remaining blockers are product and publication evidence rather than missing
prose:

- publish signed, versioned artifacts to a configured public repository;
- decide whether `gale-interop-ravel` belongs in the v1 artifact promise;
- publish stable versioned Scaladoc and wire the guide to it;
- promote browser execution from an intended route to a release gate if the
  announcement claims browser certification; and
- visually inspect the final site in an isolated browser before deployment.

Verification for this documentation revision:

- `sbt docs/tlSite`: 26 mdoc inputs compiled with no errors and 22 HTML pages
  rendered;
- rendered-file and HTML inspection: the new routes and anchors are present,
  and internal editorial/release records are not included as public pages;
- `sbt testAllFull`: passed for the current JVM and Scala.js project graph; and
- `git diff --check`: clean.

An isolated browser was unavailable in the authoring session, so screenshot-
level visual review is explicitly unverified. The user's existing Chrome
profile was not used as a substitute.
