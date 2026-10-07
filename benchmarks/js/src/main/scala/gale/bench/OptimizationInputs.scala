package gale.bench

import scala.scalajs.js

object OptimizationInputs:
  def read(relative: String = "docs/verification/optim-python/fixtures-v1.decimal.tsv"): String =
    js.Dynamic.global
      .require("fs")
      .readFileSync(relative, "utf8")
      .asInstanceOf[String]
