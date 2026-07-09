package gale.bench

import gale.linalg.*
import scala.scalajs.js

object DenseKernelJsBench:
  def main(args: Array[String]): Unit =
    val x = Vec.tabulate(256)(i => ((i % 17) - 8).toDouble * 0.125)
    val y = Vec.tabulate(256)(i => ((i % 11) - 5).toDouble * 0.25)
    val a = Matrix.tabulate(32, 32)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    val b = Matrix.tabulate(32, 32)((i, j) => ((i * 7 - j * 2) % 23).toDouble * 0.015625)
    val gemvX = Vec.tabulate(32)(i => ((i % 7) - 3).toDouble * 0.5)

    time("dot", 20000) {
      var i = 0
      var acc = 0.0
      while i < 20000 do
        acc += x.dot(y)
        i += 1
      acc
    }

    time("axpy", 20000) {
      val out = y.mutableCopy
      var i = 0
      while i < 20000 do
        out += x
        i += 1
      out(0)
    }

    time("gemv", 2000) {
      val out = MutableVec.zeros(32)
      var i = 0
      var acc = 0.0
      while i < 2000 do
        a.mulInto(gemvX, out)
        acc += out(0)
        i += 1
      acc
    }

    time("gemm", 200) {
      var i = 0
      var acc = 0.0
      while i < 200 do
        acc += (a * b)(0, 0)
        i += 1
      acc
    }

  private def time(name: String, iterations: Int)(body: => Double): Unit =
    val start = js.Date.now()
    val result = body
    val elapsedMs = js.Date.now() - start
    val perSecond =
      if elapsedMs == 0.0 then Double.PositiveInfinity
      else iterations.toDouble / (elapsedMs / 1000.0)
    println(f"$name%-5s ${perSecond}%.2f ops/s result=$result%.6f")
