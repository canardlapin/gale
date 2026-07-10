import org.scalajs.linker.interface.ModuleKind
import org.scalajs.sbtplugin.ScalaJSPlugin.autoImport.*
import pl.project13.scala.sbt.JmhPlugin
import sbtcrossproject.CrossPlugin.autoImport.*
import sbtcrossproject.CrossProject
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*

ThisBuild / organization := "io.gale"
ThisBuild / scalaVersion := "3.3.8"
ThisBuild / version      := "0.1.0-SNAPSHOT"

lazy val commonScalacOptions = Seq(
  "-deprecation",
  "-feature",
  "-unchecked",
  "-Xmax-inlines:64"
)

lazy val munitVersion = "1.3.0"

// Experimental WebAssembly output for Scala.js, toggled OFF by default so a plain
// `sbt testAll` produces exactly today's JavaScript build. Set GALE_WASM=1 (or
// true) to enable; running the tests then needs Node with an ES-module loader and
// --experimental-wasm-exnref (see README). The env var is read once at load, so
// the default build never touches the linker config.
lazy val wasmEnabled: Boolean =
  sys.env.get("GALE_WASM").exists(v => v == "1" || v.equalsIgnoreCase("true"))

lazy val jsWasmSettings: Seq[Def.Setting[_]] =
  if (wasmEnabled)
    Seq(
      scalaJSLinkerConfig ~= {
        _.withExperimentalUseWebAssembly(true).withModuleKind(ModuleKind.ESModule)
      }
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
    .dependsOn(coreJVM)
    .enablePlugins(JmhPlugin)
    .settings(benchmarkSettings)
    .settings(
      name := "gale-benchmarks-jvm",
      // Breeze in COMPILE scope here (not test) so the paired gale-vs-Breeze JMH
      // benchmarks can call it. This module is publish-skipped and is never a
      // dependency of core/laws, so gale-core stays 100% Breeze-free. Same native
      // Scala 3 artifact (breeze_3) as the parity module — its netlib backend runs
      // the pure-Java F2J fallback here, which is deliberately the baseline the
      // benchmarks target (native-BLAS Breeze is a separate, deferred comparison).
      libraryDependencies += "org.scalanlp" %% "breeze" % breezeVersion
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

lazy val root =
  project
    .in(file("."))
    .aggregate(coreJS, coreJVM, lawsJS, lawsJVM, benchmarksJVM, benchmarksJS, parity)
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
addCommandAlias("benchCompile", ";benchmarksJVM/Jmh/compile;benchmarksJS/compile")
addCommandAlias("benchSmokeJS", ";benchmarksJS/run")
addCommandAlias("benchSmokeJSFull", ";set benchmarksJS/scalaJSStage := FullOptStage;benchmarksJS/run")
