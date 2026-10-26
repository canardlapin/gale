ThisBuild / scalaVersion := "3.7.4"
val candidate = sys.props("gale.release.version")
val repository = file(sys.props("gale.release.repository")).getCanonicalFile
val fixture = file("src/test/scala").getCanonicalFile
val verifyStaged = taskKey[Unit]("Require both Gale JARs to come from the supplied staging repository")
lazy val consumerSettings = Seq(
  publish / skip := true,
  resolvers := Seq("gale-staged" at repository.toURI.toASCIIString, Resolver.mavenCentral),
  // Project scope overrides sbt's default inclusion of the Ivy local repository.
  externalResolvers := resolvers.value,
  Test / unmanagedSourceDirectories := Seq(fixture),
  verifyStaged := {
    val jars = update.value.allFiles.filter(_.getName.startsWith("gale-")).distinct
    require(jars.size == 2, s"expected core and laws JARs, found ${jars.toList}")
    require(jars.forall(_.getCanonicalPath.startsWith(repository.getPath + java.io.File.separator)),
      s"Gale dependency escaped staging: ${jars.toList}")
    streams.value.log.info("Verified staged Gale dependencies: " + jars.mkString(", "))
  }
)
lazy val jvm = project.in(file("jvm")).settings(consumerSettings).settings(
  libraryDependencies ++= Seq("io.github.canardlapin" %% "gale-core" % candidate, "io.github.canardlapin" %% "gale-laws" % candidate)
)
lazy val js = project.in(file("js")).enablePlugins(ScalaJSPlugin).settings(consumerSettings).settings(
  libraryDependencies ++= Seq("io.github.canardlapin" %%% "gale-core" % candidate, "io.github.canardlapin" %%% "gale-laws" % candidate)
)
lazy val root = project.in(file(".")).settings(publish / skip := true)
addCommandAlias("testRelease", ";jvm/verifyStaged;js/verifyStaged;jvm/test;js/test;set js / scalaJSStage := FullOptStage;js/test")
