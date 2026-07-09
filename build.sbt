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
      name := "gale-benchmarks-jvm"
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
    .aggregate(coreJS, coreJVM, lawsJS, lawsJVM, benchmarksJVM, benchmarksJS)
    .settings(
      name := "gale",
      publish / skip := true
    )

addCommandAlias("compileAll", ";coreJVM/compile;coreJS/compile;lawsJVM/compile;lawsJS/compile")
addCommandAlias("testAll", ";coreJVM/test;coreJS/test;lawsJVM/test;lawsJS/test")
// Like testAll, then a full-optimizing Scala.js link of the JS test bundles as a
// stricter (Closure-level) check that fastLink-only builds can miss.
addCommandAlias("testAllFull", ";testAll;coreJS/Test/fullLinkJS;lawsJS/Test/fullLinkJS")
addCommandAlias("benchCompile", ";benchmarksJVM/Jmh/compile;benchmarksJS/compile")
addCommandAlias("benchSmokeJS", ";benchmarksJS/run")
addCommandAlias("benchSmokeJSFull", ";set benchmarksJS/scalaJSStage := FullOptStage;benchmarksJS/run")
