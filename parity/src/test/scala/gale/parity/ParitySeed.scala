package gale.parity

/** Replayable ScalaCheck seeds for the property-based parity suites.
  *
  * A suite passes its pinned default to [[initial]]. Setting the system property
  * `gale.parity.seed` to a seed printed by a failing run replays that run without
  * editing the suite:
  *
  * {{{
  * sbt -Dgale.parity.seed=<base64 seed> "parity/testOnly gale.parity.ReductionsNumericsParitySuite"
  * }}}
  *
  * The build forwards the property to the forked test JVM. munit-scalacheck
  * prints the failing seed with every property failure.
  */
object ParitySeed:
  val Property = "gale.parity.seed"

  def initial(default: String): String =
    sys.props.get(Property).map(_.trim).filter(_.nonEmpty).getOrElse(default)
