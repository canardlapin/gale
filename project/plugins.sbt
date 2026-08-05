addSbtPlugin("org.scala-js" % "sbt-scalajs" % "1.22.0")
addSbtPlugin("org.portable-scala" % "sbt-scalajs-crossproject" % "1.3.2")
addSbtPlugin("pl.project13.scala" % "sbt-jmh" % "0.4.8")
addSbtPlugin("org.typelevel" % "sbt-typelevel-site" % "0.8.7")
addSbtPlugin("org.scalameta" % "sbt-scalafmt" % "2.5.6")
addSbtPlugin("com.github.sbt" % "sbt-ci-release" % "1.11.2")
addSbtPlugin("org.scoverage" % "sbt-scoverage" % "2.4.4")
addSbtPlugin("ch.epfl.scala" % "sbt-version-policy" % "3.2.1")
addSbtPlugin("ch.epfl.scala" % "sbt-tasty-mima" % "1.4.0")

// Coursier 2.1.24 requests is-terminal 0.1.1, whose Java-22 multi-release
// class was accidentally emitted as Java 23 bytecode. 0.1.2 corrects the
// class version so a fresh Scalafmt runner cache works on the supported JDK 22.
dependencyOverrides += "io.github.alexarchambault" % "is-terminal" % "0.1.2"
