package gale.bench

import scala.io.Source
import java.nio.file.{Files, Path}

object OptimizationInputs:
  def read(relative: String = "docs/verification/optim-python/fixtures-v1.decimal.tsv"): String =
    val roots = Iterator.iterate(Path.of("").toAbsolutePath)(_.getParent).takeWhile(_ != null)
    val path = roots
      .map(_.resolve(relative))
      .find(Files.isRegularFile(_))
      .getOrElse(
        throw new IllegalArgumentException(s"cannot locate $relative above current directory")
      )
    val source = Source.fromFile(path.toFile)
    try source.mkString
    finally source.close()
