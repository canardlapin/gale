package gale.compat

import breeze.linalg.DenseMatrix
import gale.interop.breeze.BreezeMigration

/** A deliberately tiny external-style consumer for the published Breeze bridge.
  *
  * This has no source dependency on Gale. The release probe publishes the
  * admitted modules to the local repository, then compiles this project from
  * their coordinates and generated POMs.
  */
object PublishedInteropConsumer:
  def determinant(): Double =
    BreezeMigration.det(DenseMatrix((2.0, 1.0), (1.0, 3.0)))
