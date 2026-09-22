ThisBuild / scalaVersion := "3.7.4"
ThisBuild / externalResolvers := Seq(
  "audited-staged-maven" at sys.props.getOrElse("gale.audit.maven", sys.error("Set -Dgale.audit.maven=file:///absolute/path/to/target/sona-staging/")),
  Resolver.mavenCentral
)
lazy val jvm = project.in(file("jvm")).settings(
  libraryDependencies += "io.github.canardlapin" %% "gale-laws" % "0.1.0-M1" % Test
)
lazy val js = project.in(file("js")).enablePlugins(ScalaJSPlugin).settings(
  libraryDependencies += "io.github.canardlapin" %%% "gale-laws" % "0.1.0-M1" % Test
)
