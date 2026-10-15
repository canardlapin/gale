// This branch prepares only the JVM dependency artifacts listed in docs/dependency-release.md.
// Publishing is an explicit local staging operation; Central upload is a separate step.
ThisBuild / isSnapshot := false
ThisBuild / version := "1.0.0-RC1"
ThisBuild / homepage := Some(url("https://github.com/canardlapin/gale"))
ThisBuild / scmInfo := Some(
  ScmInfo(
    url("https://github.com/canardlapin/gale"),
    "scm:git:https://github.com/canardlapin/gale.git"
  )
)
ThisBuild / publishMavenStyle := true
ThisBuild / publishTo := localStaging.value

addCommandAlias("stageDependencyRelease", ";coreJVM/publish")
