import org.scalajs.linker.interface.{ESVersion, ModuleKind}
import org.scalajs.jsenv.nodejs.NodeJSEnv
import org.scalajs.sbtplugin.ScalaJSPlugin.autoImport.*
import laika.ast.Path.Root
import laika.helium.config.{HeliumIcon, IconLink}
import pl.project13.scala.sbt.JmhPlugin
import sbtcrossproject.CrossPlugin.autoImport.*
import sbtcrossproject.CrossProject
import scalajscrossproject.ScalaJSCrossPlugin.autoImport.*

lazy val scalaBaselineVersion = "3.7.4"
lazy val scalaNextVersion = "3.8.4"

ThisBuild / organization := "io.github.canardlapin"
ThisBuild / scalaVersion := scalaBaselineVersion
ThisBuild / versionScheme := Some("early-semver")
ThisBuild / homepage     := Some(url("https://github.com/canardlapin/gale"))
ThisBuild / licenses     := Seq("Apache-2.0" -> url("https://www.apache.org/licenses/LICENSE-2.0.txt"))
ThisBuild / scmInfo := Some(
  ScmInfo(
    url("https://github.com/canardlapin/gale"),
    "scm:git:https://github.com/canardlapin/gale.git",
    Some("scm:git:git@github.com:canardlapin/gale.git")
  )
)
ThisBuild / developers := List(
  Developer(
    id = "canardlapin",
    name = "canardlapin",
    email = "307091466+canardlapin@users.noreply.github.com",
    url = url("https://github.com/canardlapin")
  )
)

// sbt-ci-release/sbt-dynver owns the release version. Gale is on the 0.1
// line: only an exact `v0.1.x`, `v0.1.x-Mn`, or `v0.1.x-RCn` tag yields a
// publishable version. Every other state remains a unique
// `0.1.0+...-SNAPSHOT` development build. 1.0 is a later compatibility freeze,
// not the current snapshot label.
def galeReleaseVersion(out: sbtdynver.GitDescribeOutput): String = {
  val taggedVersion = out.ref.value.stripPrefix("v")
  val allowedTag = taggedVersion.matches("0\\.1\\.[0-9]+(?:-(?:M|RC)[1-9][0-9]*)?")
  val exactTag = allowedTag && out.ref.value.startsWith("v") && out.commitSuffix.distance == 0
  val base = if (exactTag) taggedVersion else "0.1.0"
  val commit =
    if (out.commitSuffix.distance == 0) ""
    else s"+${out.commitSuffix.distance}-${out.commitSuffix.sha}"
  val dirty =
    if (out.dirtySuffix.value.isEmpty) ""
    else s"+${out.dirtySuffix.value.stripPrefix("+")}"
  val publishable = exactTag && dirty.isEmpty
  if (publishable) base else s"$base$commit$dirty-SNAPSHOT"
}

def galeReleaseFallbackVersion(date: java.util.Date): String = {
  s"0.1.0-SNAPSHOT-${sbtdynver.DynVer.timestamp(date)}"
}

inThisBuild(List(
  version := dynverGitDescribeOutput.value.mkVersion(
    galeReleaseVersion,
    galeReleaseFallbackVersion(dynverCurrentDate.value)
  ),
  // Own snapshot classification so a synthetic M1 dry-run version is not
  // treated as a snapshot and routed to Central's snapshot service.
  isSnapshot := version.value.endsWith("-SNAPSHOT"),
  dynver := {
    val date = new java.util.Date
    sbtdynver.DynVer
      .getGitDescribeOutput(date)
      .map(galeReleaseVersion)
      .getOrElse(galeReleaseFallbackVersion(date))
  }
))

lazy val commonScalacOptions = Seq(
  "-deprecation",
  "-feature",
  "-unchecked",
  "-Xmax-inlines:64",
  "-Werror"
)

// Release candidates must not silently pull an unpublished snapshot at
// compile or runtime. Keep this check attached to each admitted artifact so
// publication can fail before assembling a partial bundle.
lazy val releaseSnapshotCheck = taskKey[Unit](
  "Fail when an admitted artifact resolves a SNAPSHOT dependency"
)

lazy val releaseInternalModules = Set(
  "gale-core_3",
  "gale-core_sjs1_3",
  "gale-laws_3",
  "gale-laws_sjs1_3",
  "gale-interop-breeze_3",
  "gale-backend-jvm-vector_3",
  "gale-backend-jvm-native_3",
  "gale-backend-jvm-blas-ffm_3"
)

lazy val releaseSnapshotSettings = Seq(
  releaseSnapshotCheck := {
    val reports = Seq((Compile / update).value, (Runtime / update).value)
    val snapshots = reports
      .flatMap(_.configurations)
      .flatMap(_.modules)
      .map(_.module)
      .filter(_.revision.toUpperCase.contains("SNAPSHOT"))
      .filterNot(m =>
        m.organization == organization.value && releaseInternalModules.contains(m.name)
      )
      .distinct

    if (snapshots.nonEmpty) {
      val rendered = snapshots
        .map(m => s"${m.organization}:${m.name}:${m.revision}")
        .mkString(", ")
      sys.error(s"release artifact has prohibited SNAPSHOT dependencies: $rendered")
    }
  }
)

lazy val releaseVersionCheck = taskKey[Unit](
  "Require a clean, tag-derived 0.1 semantic version for publication"
)

lazy val releaseStagingCheck = taskKey[Unit](
  "Require a non-snapshot 0.1 candidate and sbt's local Central staging repository"
)

lazy val releaseVersionSettings = Seq(
  releaseVersionCheck := {
    val candidate = version.value
    val releasePattern = "0\\.1\\.[0-9]+(?:-(?:M|RC)[1-9][0-9]*)?"
    if (candidate.endsWith("-SNAPSHOT") || !candidate.matches(releasePattern))
      sys.error(
        s"0.1 publication requires a clean v0.1.x, v0.1.x-Mn, or v0.1.x-RCn tag; derived version was $candidate"
      )
  },
  releaseStagingCheck := {
    val candidate = version.value
    if (isSnapshot.value)
      sys.error(s"release candidate $candidate is still classified as a snapshot")
    publishTo.value match {
      case Some(destination) if destination.name == "local-staging" => ()
      case other =>
        sys.error(
          s"release candidate $candidate must publish to local-staging, found ${other.fold("no destination")(_.toString)}"
        )
    }
  }
)

lazy val munitVersion = "1.3.0"

// Experimental WebAssembly output for Scala.js, toggled OFF by default so a plain
// `sbt testAll` produces exactly today's JavaScript build. Set GALE_WASM=1 (or
// true) to enable. The opt-in settings select ES2022 modules. Scala.js 1.22's
// Wasm backend needs a Wasm 3.0 engine — Node.js 25+ (CI uses Node 25). The
// --experimental-wasm-exnref flag remains for engines that still gate exnref.
// The env var is read once at load, so the default build never touches the
// linker or jsEnv config.
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
    .settings(releaseSnapshotSettings)
    .settings(releaseVersionSettings)
    .settings(
      name := "gale-core",
      description := "Cross-platform linear algebra for Scala 3: dense and sparse matrices, factorizations, and solvers on shared strided kernels."
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
    .settings(releaseSnapshotSettings)
    .settings(releaseVersionSettings)
    .settings(
      name := "gale-laws",
      description := "Reusable munit/ScalaCheck law bundles for gale's public API.",
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

// Executable public guide site. User-facing Markdown lives under docs/user;
// internal design, audit, and release-evidence documents remain versioned under
// docs/ but are not renderer inputs. The sbt project lives in site/ so it does
// not collide with mdoc's input root. This is the standalone site plugin only:
// Gale does not opt into Typelevel branding or the full sbt-typelevel stack.
lazy val docs =
  project
    .in(file("site"))
    .dependsOn(coreJVM)
    .enablePlugins(TypelevelSitePlugin)
    .settings(
      name           := "gale-docs",
      publish / skip := true,
      scalacOptions ++= commonScalacOptions,
      mdocIn := file("docs/user"),
      tlSiteHelium := tlSiteHelium.value.site.topNavigationBar(
        homeLink = IconLink.internal(Root / "index.md", HeliumIcon.home)
      )
    )

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
        "org.scalameta" %% "munit"            % munitVersion  % Test,
        "org.scalameta" %% "munit-scalacheck" % munitVersion  % Test,
        "org.scalanlp"  %% "breeze"           % breezeVersion % Test
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
    .settings(releaseSnapshotSettings)
    .settings(releaseVersionSettings)
    .settings(
      name := "gale-interop-breeze",
      description := "Breeze interoperability for gale (JVM).",
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
    .settings(releaseSnapshotSettings)
    .settings(releaseVersionSettings)
    .settings(
      name := "gale-backend-jvm-vector",
      description := "Optional JDK Vector API GEMM backend for gale, with measured adaptive dispatch.",
      scalacOptions ++= commonScalacOptions,
      Test / fork := true,
      Test / javaOptions += "--add-modules=jdk.incubator.vector",
      libraryDependencies += "org.scalameta" %% "munit" % munitVersion % Test
    )

lazy val ravelVersion = "1.0.0-SNAPSHOT"

// gale-interop-ravel is the copy-only boundary between neutral dense Ravel
// storage and Gale's mathematical vector/matrix types. Neither core project
// depends on the other.
lazy val interopRavel: CrossProject =
  crossProject(JSPlatform, JVMPlatform)
    .crossType(CrossType.Full)
    .in(file("interop-ravel"))
    .dependsOn(core)
    .settings(commonSettings)
    .settings(
      name := "gale-interop-ravel",
      description := "Explicit copy conversions between Ravel arrays and Gale vectors and matrices.",
      // Ravel is still a development snapshot and is deliberately outside
      // the 0.1 milestone artifact set.
      publish / skip := true,
      libraryDependencies +=
        "io.github.canardlapin" %%% "ravel-core" % ravelVersion
    )
    .jsSettings(jsWasmSettings: _*)

lazy val interopRavelJVM = interopRavel.jvm
lazy val interopRavelJS  = interopRavel.js

// JVM 22+ native storage. Kept separate so core and every Scala.js artifact stay
// free of java.lang.foreign references.
lazy val nativeBackend =
  project
    .in(file("backend-jvm-native"))
    .dependsOn(coreJVM)
    .settings(releaseSnapshotSettings)
    .settings(releaseVersionSettings)
    .settings(
      name := "gale-backend-jvm-native",
      description := "Optional JDK 22+ off-heap matrix storage for gale over FFM MemorySegment.",
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
    .settings(releaseSnapshotSettings)
    .settings(releaseVersionSettings)
    .settings(
      name := "gale-backend-jvm-blas-ffm",
      description := "Optional JDK 22+ runtime-discovered BLAS/LAPACK backend for gale via FFM.",
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

// gale-demo: a browser-only PCA demo page (publish-skipped, never aggregated
// into the default build). Linked with ModuleKind.NoModule so the emitted
// script runs from a plain <script> tag over file:// — open demo/index.html
// directly after linking, no local server or bundler required. scalajs-dom is
// the only JS dependency.
lazy val demo =
  project
    .in(file("demo"))
    .enablePlugins(ScalaJSPlugin)
    .dependsOn(coreJS)
    .settings(
      name           := "gale-demo",
      publish / skip := true,
      scalacOptions ++= commonScalacOptions,
      scalaJSUseMainModuleInitializer := true,
      scalaJSLinkerConfig ~= { _.withModuleKind(ModuleKind.NoModule) },
      libraryDependencies += "org.scala-js" %%% "scalajs-dom" % "2.8.0"
    )

// Compile-only downstream source-consumption probe. The consumer uses the
// current Scala Next release against the locally published gale-core artifact
// compiled by this build's Scala 3.7.4 project. It therefore exercises the
// published TASTy/API boundary rather than compiling the probe against Gale's
// sibling source project.
lazy val scalaNextConsumer =
  project
    .in(file("compat/scala-next-consumer"))
    .settings(
      name           := "gale-scala-next-consumer-probe",
      scalaVersion   := scalaNextVersion,
      publish / skip := true,
      scalacOptions ++= commonScalacOptions,
      libraryDependencies += "io.github.canardlapin" %% "gale-core" % version.value
    )

// Compile-only downstream probe for the admitted Breeze interop artifact. It
// deliberately has no `.dependsOn` edge: the CI alias publishes Gale's local
// artifacts first, then resolves the coordinate through the generated POM.
// This catches missing packages, metadata, and transitive dependency mistakes
// that source-level interop tests cannot see.
lazy val publishedInteropConsumer =
  project
    .in(file("compat/published-interop-consumer"))
    .settings(
      name           := "gale-published-interop-consumer-probe",
      publish / skip := true,
      scalacOptions ++= commonScalacOptions,
      libraryDependencies += "io.github.canardlapin" %% "gale-interop-breeze" % version.value
    )

lazy val root =
  project
    .in(file("."))
    .aggregate(
      coreJS, coreJVM, lawsJS, lawsJVM, benchmarksJVM, benchmarksJS,
      parity, interopBreeze, interopRavelJVM, interopRavelJS, vectorBackend,
      nativeBackend, blasFfmBackend
    )
    .settings(
      name := "gale",
      publish / skip := true
    )

addCommandAlias("compileAll", ";coreJVM/compile;coreJS/compile;lawsJVM/compile;lawsJS/compile;coreJVM/publishLocal;scalaNextConsumer/compile")
addCommandAlias("testAll", ";coreJVM/test;coreJS/test;lawsJVM/test;lawsJS/test")
// Like testAll, then a full-optimizing Scala.js link of the JS test bundles as a
// stricter (Closure-level) check that fastLink-only builds can miss.
addCommandAlias("testAllFull", ";testAll;coreJS/Test/fullLinkJS;lawsJS/Test/fullLinkJS")
// Breeze parity harness (JVM-only correctness parity vs Scala Breeze 2.1.0).
addCommandAlias("parityTest", ";parity/test")
// Breeze interop module (conversions + migration aids).
addCommandAlias("interopBreezeTest", ";interopBreeze/test")
// Ravel interop module (copy-only dense vector/matrix conversions).
addCommandAlias("interopRavelTest", ";interopRavelJVM/test;interopRavelJS/test")
// Candidate modules are checked independently so an optional development
// integration cannot mask a release dependency failure. Provisional backends
// stay in the court; they are not part of the M1 publication aliases.
addCommandAlias(
  "releaseDependencyCheck",
  ";coreJVM/releaseSnapshotCheck;coreJS/releaseSnapshotCheck;lawsJVM/releaseSnapshotCheck;lawsJS/releaseSnapshotCheck;interopBreeze/releaseSnapshotCheck;vectorBackend/releaseSnapshotCheck;nativeBackend/releaseSnapshotCheck;blasFfmBackend/releaseSnapshotCheck"
)
addCommandAlias(
  "releaseM1Preflight",
  ";releaseDependencyCheck;coreJVM/releaseVersionCheck;coreJS/releaseVersionCheck;lawsJVM/releaseVersionCheck;lawsJS/releaseVersionCheck;coreJVM/releaseStagingCheck;coreJS/releaseStagingCheck;lawsJVM/releaseStagingCheck;lawsJS/releaseStagingCheck"
)
addCommandAlias(
  "releaseM1Unsigned",
  ";releaseM1Preflight;coreJVM/publish;coreJS/publish;lawsJVM/publish;lawsJS/publish"
)
addCommandAlias(
  "releaseM1Signed",
  ";releaseM1Preflight;coreJVM/publishSigned;coreJS/publishSigned;lawsJVM/publishSigned;lawsJS/publishSigned"
)
// JVM-only Vector-API (SIMD) acceleration backend.
addCommandAlias("vectorBackendTest", ";vectorBackend/test")
addCommandAlias("nativeBackendTest", ";nativeBackend/test")
addCommandAlias("blasFfmBackendTest", ";blasFfmBackend/test")
addCommandAlias("benchFfmCompile", ";benchmarksFfm/Jmh/compile")
addCommandAlias("benchCompile", ";benchmarksJVM/Jmh/compile;benchmarksJS/compile")
addCommandAlias("benchSmokeJS", ";benchmarksJS/run")
// Browser PCA demo: link, then open demo/index.html in a browser.
addCommandAlias("demoBuild", ";demo/fastLinkJS")
addCommandAlias("scalaNextConsumerProbe", ";coreJVM/publishLocal;scalaNextConsumer/compile")
addCommandAlias("publishedInteropProbe", ";coreJVM/publishLocal;interopBreeze/publishLocal;publishedInteropConsumer/compile")
addCommandAlias("benchSmokeJSFull", ";set benchmarksJS/scalaJSStage := FullOptStage;benchmarksJS/run")
// Compile API docs for both public platforms and execute/render the guide site.
addCommandAlias("docsCheck", ";coreJVM/doc;coreJS/doc;docs/tlSite")
