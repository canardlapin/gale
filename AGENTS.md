# AGENTS.md

## Cursor Cloud specific instructions

Gale is a Scala 3 linear algebra library (no GUI/server). Development is driven
entirely through `sbt`. The build cross-compiles every module to the JVM and to
Scala.js/Node. Command aliases live in `build.sbt`; the standard dev commands are
documented in `README.md` under "Development".

### Toolchain (already provisioned in the snapshot)
- JDK 21 is the default and is required for the core JVM/JS/laws/docs lanes.
- `sbt` (launcher 1.11.7, matching `project/build.properties`) is on `PATH`.
- Node 22 (via nvm) is available for the Scala.js/Node test lane.
- The startup update script only refreshes library dependencies
  (`coreJVM/coreJS/lawsJVM/lawsJS/parity/interopBreeze/docs` `update`). It does
  not compile, test, or build.

### Primary dev commands (JDK 21)
- Compile everything: `sbt compileAll`
- Lint / formatting gate: `sbt scalafmtCheckAll`
- Tests (JVM + JS, core + laws): `sbt testAll`
- Breeze differential parity + conversions: `sbt parityTest interopBreezeTest`
- Executable guides / API docs: `sbt docs/mdoc` (fast) or `sbt docsCheck` (full
  Scaladoc + Laika site; heavier).

### Non-obvious gotchas
- `sbt update` at the **root** fails: the aggregated `interop-ravel` module
  depends on an unpublished `ravel-core 1.0.0-SNAPSHOT`. Always scope
  dependency resolution / tasks to specific modules (e.g. `coreJVM/update`),
  which is why the update script does not use a bare `update`/`compile`.
- `sbt console` does **not** accept piped stdin here (no TTY); it just prints the
  banner and exits. To exercise the API non-interactively, run `sbt docs/mdoc`
  (executes the README/getting-started examples) or add a temporary MUnit test.
- The Vector API incubator lanes rely on `.jvmopts`
  (`--add-modules=jdk.incubator.vector`); leave that file in place.

### Out-of-scope lanes (not provisioned by default)
- FFM BLAS/LAPACK (`blasFfmBackendTest`, `nativeBackendTest`, `benchFfmCompile`)
  require JDK 22 with `--enable-native-access` and the OpenBLAS runtime
  (`apt-get install -y libopenblas-dev`). Neither JDK 22 nor OpenBLAS is
  installed.
- The experimental WebAssembly lane (`GALE_WASM=1`) needs Node 25+; the default
  Node 22 cannot instantiate its Wasm 3.0 bundles.
