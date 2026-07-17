import org.scalajs.linker.interface.{ESVersion, ModuleKind}
import org.scalajs.jsenv.nodejs.NodeJSEnv
import org.scalajs.sbtplugin.ScalaJSPlugin.autoImport.*
import pl.project13.scala.sbt.JmhPlugin
import sbtcrossproject.CrossPlugin.autoImport.*
import sbtcrossproject.CrossProject
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*

ThisBuild / organization := "io.gale"
ThisBuild / scalaVersion := "3.3.8"
ThisBuild / version      := "1.0.0-SNAPSHOT"

lazy val commonScalacOptions = Seq(
  "-deprecation",
  "-feature",
  "-unchecked",
  "-Xmax-inlines:64"
)

lazy val munitVersion = "1.3.0"

// Experimental WebAssembly output for Scala.js, toggled OFF by default so a plain
// `sbt testAll` produces exactly today's JavaScript build. Set GALE_WASM=1 (or
// true) to enable. The opt-in settings select ES2022 modules and pass Node's
// --experimental-wasm-exnref runtime flag. The env var is read once at load, so
// the default build never touches the linker or jsEnv config.
lazy val wasmEnabled: Boolean =
  sys.env.get("GALE_WASM").exists(v => v == "1" || v.equalsIgnoreCase("true"))

lazy val jsWasmSettings: Seq[Def.Setting[_]] =
  if (wasmEnabled)
    Seq(
      scalaJSLinkerConfig ~= {
        _.withESFeatures(
          _.withESVersion(ESVersion.ES2022).withUseWebAssembly(true)
        ).withModuleKind(ModuleKind.ESModule)
      },
      jsEnv := new NodeJSEnv(
        NodeJSEnv.Config().withArgs(List("--experimental-wasm-exnref"))
      )
    )
  else Seq.empty

lazy val commonSettings = Seq(
  scalacOptions ++= commonScalacOptions,
  Test / fork := false,
  libraryDependencies ++= Seq(
    "org.scalameta" %%% "munit"            % munitVersion % Test,
    "org.scalameta" %%% "munit-scalacheck" % munitVersion % Test
  )
)

lazy val core: CrossProject =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("core"))
    .settings(commonSettings)
    .settings(
      name := "gale-core"
    )
    .jsSettings(jsWasmSettings: _*)

// gale-laws: reusable, munit/scalacheck-backed law bundles built on the public
// core API. munit and scalacheck are MAIN dependencies here (the bundles are
// library code, not tests), so downstream suites can call them directly.
lazy val laws: CrossProject =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("laws"))
    .dependsOn(core)
    .settings(
      name := "gale-laws",
      scalacOptions ++= commonScalacOptions,
      Test / fork := false,
      libraryDependencies ++= Seq(
        "org.scalameta" %%% "munit"            % munitVersion,
        "org.scalameta" %%% "munit-scalacheck" % munitVersion
      )
    )
    .jsSettings(jsWasmSettings: _*)

lazy val lawsJS  = laws.js
lazy val lawsJVM = laws.jvm

lazy val coreJS  = core.js
lazy val coreJVM = core.jvm

lazy val breezeVersion = "2.1.0"

// gale-parity: a JVM-only correctness-parity harness that compares gale's public
// API against Scala Breeze on bit-identical random data. Breeze lives ONLY here,
// in test scope — gale-core stays 100% Breeze-free (PRD hard constraint), and no
// other module depends on this one.
//
// Breeze 2.1.0 is cross-published for Scala 3 (org.scalanlp:breeze_3:2.1.0, a
// native build carrying .tasty), so a plain `%%` resolves the Scala 3 artifact
// directly; no `CrossVersion.for3Use2_13` shim is needed (and its 2.13 variant is
// the more fragile path here). Its netlib backend (dev.ludovic.netlib) is a
// pure-Java reference implementation with a JVM fallback — it may log a one-time
// "native BLAS not found, using Java" notice, which is harmless.
lazy val parity =
  project
    .in(file("parity"))
    .dependsOn(coreJVM)
    .settings(
      name           := "gale-parity",
      publish / skip := true,
      Test / fork    := false,
      scalacOptions ++= commonScalacOptions,
      libraryDependencies ++= Seq(
        "org.scalameta" %% "munit"  % munitVersion  % Test,
        "org.scalanlp"  %% "breeze" % breezeVersion % Test
      )
    )

// gale-interop-breeze (PRD: "Breeze conversion helpers and migration aids"): the
// JVM-only bridge module where gale meets Breeze. Breeze is a COMPILE dependency
// here because conversion IS this module's purpose — but nothing in core/laws
// depends on it, so gale-core stays 100% Breeze-free. Same native Scala 3 artifact
// (breeze_3) used by the parity and benchmark modules.
lazy val interopBreeze =
  project
    .in(file("interop-breeze"))
    .dependsOn(coreJVM)
    .settings(
      name := "gale-interop-breeze",
      scalacOptions ++= commonScalacOptions,
      Test / fork := false,
      libraryDependencies ++= Seq(
        "org.scalanlp"  %% "breeze" % breezeVersion,
        "org.scalameta" %% "munit"  % munitVersion % Test
      )
    )

// gale-backend-jvm-vector: an OPTIONAL, JVM-only acceleration module supplying a
// `given Backend` whose dense `gemm` uses the JDK Vector API (jdk.incubator.vector)
// for SIMD. The incubator module must be resolvable at BOTH compile time (the
// in-process scalac in the sbt JVM — enabled by the repo-root `.jvmopts` carrying
// `--add-modules=jdk.incubator.vector`, which affects the whole sbt launch) and run
// time (the forked test JVM, via `Test / javaOptions` below). Nothing in core/laws
// depends on this module, so the pure build is untouched.
lazy val vectorBackend =
  project
    .in(file("backend-jvm-vector"))
    .dependsOn(coreJVM, lawsJVM % "test->compile")
    .settings(
      name := "gale-backend-jvm-vector",
      scalacOptions ++= commonScalacOptions,
      Test / fork := true,
      Test / javaOptions += "--add-modules=jdk.incubator.vector",
      libraryDependencies += "org.scalameta" %% "munit" % munitVersion % Test
    )

// JVM 22+ native storage. Kept separate so core and every Scala.js artifact stay
// free of java.lang.foreign references.
lazy val nativeBackend =
  project
    .in(file("backend-jvm-native"))
    .dependsOn(coreJVM)
    .settings(
      name := "gale-backend-jvm-native",
      scalacOptions ++= commonScalacOptions,
      Test / fork := true,
      Test / javaOptions += "--enable-native-access=ALL-UNNAMED",
      libraryDependencies += "org.scalameta" %% "munit" % munitVersion % Test
    )

// JVM 22+ FFM CBLAS provider. Loading is explicit and confined to this optional
// module; a core-only program never probes or loads a native library.
lazy val blasFfmBackend =
  project
    .in(file("backend-jvm-blas-ffm"))
    .dependsOn(nativeBackend, lawsJVM % "test->compile")
    .settings(
      name := "gale-backend-jvm-blas-ffm",
      scalacOptions ++= commonScalacOptions,
      Test / fork := true,
      Test / javaOptions += "--enable-native-access=ALL-UNNAMED",
      libraryDependencies += "org.scalameta" %% "munit" % munitVersion % Test
    )

lazy val benchmarkSettings = Seq(
  scalacOptions ++= Seq(
    "-deprecation",
    "-feature",
    "-unchecked",
    "-Xmax-inlines:64"
  ),
  publish / skip := true
)

lazy val benchmarksJVM =
  project
    .in(file("benchmarks/jvm"))
    .dependsOn(coreJVM, vectorBackend)
    .enablePlugins(JmhPlugin)
    .settings(benchmarkSettings)
    .settings(
      name := "gale-benchmarks-jvm",
      Jmh / javaOptions += "--add-modules=jdk.incubator.vector",
      // Breeze in COMPILE scope here (not test) so the paired gale-vs-Breeze JMH
      // benchmarks can call it. This module is publish-skipped and is never a
      // dependency of core/laws, so gale-core stays 100% Breeze-free. Same native
      // Scala 3 artifact (breeze_3) as the parity module — its netlib backend runs
      // the pure-Java F2J fallback here, which is deliberately the baseline the
      // benchmarks target (native-BLAS Breeze is a separate, deferred comparison).
      libraryDependencies += "org.scalanlp" %% "breeze" % breezeVersion
    )

// JDK 22+ copy-inclusive native crossover harness. Separate from benchmarksJVM
// so the Vector backend's JDK 21 compatibility job never compiles FFM sources.
lazy val benchmarksFfm =
  project
    .in(file("benchmarks/jvm-ffm"))
    .dependsOn(coreJVM, blasFfmBackend)
    .enablePlugins(JmhPlugin)
    .settings(benchmarkSettings)
    .settings(
      name := "gale-benchmarks-jvm-ffm",
      Jmh / javaOptions += "--enable-native-access=ALL-UNNAMED"
    )

lazy val benchmarksJS =
  project
    .in(file("benchmarks/js"))
    .enablePlugins(ScalaJSPlugin)
    .dependsOn(coreJS)
    .settings(benchmarkSettings)
    .settings(
      name := "gale-benchmarks-js",
      scalaJSUseMainModuleInitializer := true
    )
    .settings(jsWasmSettings: _*)

lazy val root =
  project
    .in(file("."))
    .aggregate(
      coreJS, coreJVM, lawsJS, lawsJVM, benchmarksJVM, benchmarksJS,
      parity, interopBreeze, vectorBackend
    )
    .settings(
      name := "gale",
      publish / skip := true
    )

addCommandAlias("compileAll", ";coreJVM/compile;coreJS/compile;lawsJVM/compile;lawsJS/compile")
addCommandAlias("testAll", ";coreJVM/test;coreJS/test;lawsJVM/test;lawsJS/test")
// Like testAll, then a full-optimizing Scala.js link of the JS test bundles as a
// stricter (Closure-level) check that fastLink-only builds can miss.
addCommandAlias("testAllFull", ";testAll;coreJS/Test/fullLinkJS;lawsJS/Test/fullLinkJS")
// Breeze parity harness (JVM-only correctness parity vs Scala Breeze 2.1.0).
addCommandAlias("parityTest", ";parity/test")
// Breeze interop module (conversions + migration aids).
addCommandAlias("interopBreezeTest", ";interopBreeze/test")
// JVM-only Vector-API (SIMD) acceleration backend.
addCommandAlias("vectorBackendTest", ";vectorBackend/test")
addCommandAlias("nativeBackendTest", ";nativeBackend/test")
addCommandAlias("blasFfmBackendTest", ";blasFfmBackend/test")
addCommandAlias("benchFfmCompile", ";benchmarksFfm/Jmh/compile")
addCommandAlias("benchCompile", ";benchmarksJVM/Jmh/compile;benchmarksJS/compile")
addCommandAlias("benchSmokeJS", ";benchmarksJS/run")
addCommandAlias("benchSmokeJSFull", ";set benchmarksJS/scalaJSStage := FullOptStage;benchmarksJS/run")
