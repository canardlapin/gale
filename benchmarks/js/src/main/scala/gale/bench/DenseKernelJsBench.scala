package gale.bench

import gale.linalg.*
import scala.scalajs.js

object DenseKernelJsBench:
  def main(args: Array[String]): Unit =
    val x = Vec.tabulate(65536)(i => ((i % 17) - 8).toDouble * 0.125)
    val y = Vec.tabulate(65536)(i => ((i % 11) - 5).toDouble * 0.25)
    val a = Matrix.tabulate(512, 512)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    val gemvX = Vec.tabulate(512)(i => ((i % 7) - 3).toDouble * 0.5)
    val gemmA = Matrix.tabulate(128, 128)((i, j) => ((i * 3 + j * 5) % 19).toDouble * 0.03125)
    val gemmB = Matrix.tabulate(128, 128)((i, j) => ((i * 7 - j * 2) % 23).toDouble * 0.015625)

    profile("dot-65536", iterations = 500)(x.dot(y))

    val gemvOut = MutableVec.zeros(512)
    profile("gemv-512", iterations = 200) {
      a.mulInto(gemvX, gemvOut)
      gemvOut(0)
    }

    profile("gemm-128", iterations = 40)((gemmA * gemmB)(0, 0))

  /** A small wall-clock profile rather than a statistical benchmark framework:
    * five complete warmup batches, then the median of nine batches. Both JS and
    * Wasm run this exact linked program with the checksum consumed and printed.
    */
  private def profile(name: String, iterations: Int)(body: => Double): Unit =
    var warmup = 0
    var checksum = 0.0
    while warmup < 5 do
      var i = 0
      while i < iterations do
        checksum += body
        i += 1
      warmup += 1

    val samples = new Array[Double](9)
    var sample = 0
    while sample < samples.length do
      val start = js.Date.now()
      var i = 0
      while i < iterations do
        checksum += body
        i += 1
      samples(sample) = (js.Date.now() - start) * 1000000.0 / iterations.toDouble
      sample += 1
    scala.util.Sorting.quickSort(samples)
    val medianNs = samples(samples.length / 2)
    println(f"$name%-10s median=$medianNs%.1f ns/op samples=9 checksum=$checksum%.6f")
