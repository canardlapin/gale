# Release readiness repairs, 26 October 2026

Baseline: `f747d2aabd5ff21851982355f2093ac8b8f91ec5`. This record supplements
`docs/release-readiness-2026-10-26.md`; the original adverse evidence is retained.
No release tag or Central publication is authorized by this repair.

## Numerical changes

GMRES now caches the original residual at initialization and recomputes it after
every completed candidate/restart cycle. The preconditioned estimate only triggers
a check. Returned diagnostics describe the original system on success and failure.
Zero/nonfinite preconditioner output cannot certify an unsolved system. Happy
breakdown is scaled against the current Arnoldi vector, not the RHS magnitude.
Operator work is one initial application, one per Arnoldi iteration, and one
verification per completed cycle. There is no per-inner-step verification matvec.

Eight new regression tests cover both tolerance modes, scalar scaling, restart and
budget exits, unusable preconditioners, exact initial points and operator counts.
The initial six regressions and full solver selection passed 41 tests on JVM
and 41 on JS; large-RHS breakdown and late nonfinite-provider cases were then added.

## Digest diagnosis: fixtures, not changed kernel outputs

The old fixtures passed on ARM JDK 21.0.13 as well as the previous ARM JDK 25 run.
Running the unchanged fixtures on x86-64 JDK 21.0.13 under Rosetta reproduced the
hosted nonsymmetric eigen/Hessenberg and optimizer mismatch counts exactly:
47, 49 and 93 cases. QR passed in that short x86 run, but had failed in the warmed
hosted whole-suite run (175 cases).

A separate Java loop warmed `Math.pow(10, exponent)` and compared its raw bits
with `StrictMath.pow`. ARM had no differences over integer exponents -120..120;
x86 produced one-ulp differences, for example:

```text
10^-111 Math=28e33d4032c2c7f5 StrictMath=28e33d4032c2c7f4 deltaUlp=1.0
10^-106 Math=29ed5b561574765b StrictMath=29ed5b561574765c deltaUlp=1.0
```

These powers feed the graded QR and nonsymmetric matrices. The optimizer fixture
also calls platform-dependent exponentials/logarithms. Java explicitly permits
this variability in [Math](https://docs.oracle.com/en/java/javase/21/docs/api/java.base/java/lang/Math.html);
[StrictMath](https://docs.oracle.com/en/java/javase/21/docs/api/java.base/java/lang/StrictMath.html)
specifies the relevant portable algorithms.

Only fixture transcendental calls were changed to StrictMath. **No digest was
updated, removed, or relaxed:** all 13,112 QR digests, 462 optimizer digests and
3,614 nonsymmetric digests remained identical when exported with the deterministic
fixtures on ARM/JDK 21. Production kernels and their floating-point contracts were
unchanged. Independent residual and reconstruction suites remain mandatory.

## Release path changes

- `main` now requires the policy's eleven named CI checks, with strict/up-to-date
  checks and administrator enforcement; force pushes/deletions remain disabled.
- Release qualification calls the full CI workflow on the same source revision
  before staging. The separate four-coordinate consumer uses the staged Maven
  repository and executes JVM, fast JS and full-optimized JS tests.
- Bundle verification parses POMs, checks immutable dependencies, archive payloads,
  checksums and the exact file allowlist, and authenticates signatures against the
  expected primary fingerprint. Tests include an ephemeral signing key, wrong-key
  and fake-signature rejection. No owner private key is used by tests.
- Upload remains USER_MANAGED. The workflow waits for Portal VALIDATED and records
  its response; it never requests publication.
- The stability inventory now explicitly includes numeric and complex APIs and
  the implemented dense SVD and opt-in sparse Cholesky surfaces. The existing M1
  compatibility promise is preserved.

## Remaining owner prerequisite

The real signing key, Portal credentials and public fingerprint must be supplied
through secure configuration. A successful unsigned rehearsal or ephemeral-key
test does not establish ownership of `io.github.canardlapin` or validate a real
signed Portal deployment. No claim of public availability follows from this work.

## Receipts

Local raw logs, runtime archives and counterexamples are under
`.agent-work/release-fixes-20261026/`; they contain no release credentials.
The correction's final local/hosted and staged-consumer results are recorded below.

Local qualification initially passed 2,650 core/law tests and optimized links
on JDK 21.0.13/Node 22.23.4, then exhausted sbt's default 1 GB heap compiling
parity in the same process. Remaining gates are resumed with a 3 GB heap.
The x86 emulation run passed 23 focused tests but timed out in QR; this is not
a QR correctness verdict or a passing x86 qualification.

Automatic approval review rejected creating/committing/pushing a repair branch
as external publication without explicit user permission. No commit or branch
was created by that rejected command. Hosted CI and a hosted unsigned rehearsal
therefore await permission to push the reviewed local changes.

## Completed local qualification

On ARM macOS, JDK 21.0.13 and Node 22.23.4:

- Formatting, full core/law suite and both optimized JS test links passed before
  the combined-process heap exhaustion. All eight final GMRES regression tests
  then passed separately on each platform; the two added tests change no product code.
- Remaining court with `-J-Xmx3G` passed: Breeze parity **182**, Breeze goldens
  freshness, Breeze interop **29**, JVM/JS benchmark compilation, JVM/JS Scaladoc,
  executable guides/site, release preflight and four-artifact unsigned staging.
- Vector JDK 21 conformance passed **84** tests.
- Verifier tests passed **8**, including signer authentication and malformed bundle
  rejection. Credential-helper tests passed **3** (synthetic values only).
- The actual four-artifact unsigned bundle passed the hardened verifier. Signing
  the same payloads with a disposable test key passed authenticated verification.
  The temporary private key was destroyed; this is not owner-key qualification.
- The separate consumer proved both Gale dependencies on each platform resolve
  under `target/sona-staging`, then executed its solve/law/GMRES test on JVM,
  fast Scala.js and full-optimized Scala.js. All three passed.
- Workflow YAML parsed; `git diff --check` passed. The original malformed bundle
  now fails instead of being reported as verified.

The consumer's first resolver attempt was rejected by Coursier. The final build
uses a Maven file URL and project-scoped `externalResolvers`, with an explicit
resolved-path assertion excluding Ivy local and unintended cached artifacts.
The failed attempt is retained in `consumer.log`; the passing receipt is
`consumer-final.log`. A second bounded Rosetta QR diagnostic also timed out;
real hosted x86 qualification remains required, rather than increasing test limits.

The owner requested help configuring credentials. `docs/release-setup.md` and
`tools/release/configure_credentials.py` now provide hidden Terminal input,
dedicated encrypted-key setup, public-key distribution and account-checked GitHub
configuration. Its read-only check confirmed no key/secrets were configured yet.
The real setup has not been run by the agent because owner-only passphrase/token
input is required. No credentials were entered in chat or written to the repository.

The owner subsequently authorized committing, pushing and verifying the repair
branch. This supersedes the push permission blocker above; it does not authorize
merging main, creating a release tag, or publishing Gale.
