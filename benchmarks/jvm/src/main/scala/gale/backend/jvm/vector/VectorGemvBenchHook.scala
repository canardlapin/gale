package gale.backend.jvm.vector

/** Bench-only access to the frozen pre-HC-4 row-major GEMV and the no-tree variant
  * of the live one ([[VectorGemvBaseline]]).
  */
object VectorGemvBenchHook:
  def gemvRowMajorBaseline(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    VectorGemvBaseline.gemvRowMajor(rows, cols, alpha, a, aOffset, aRowStride, x, xOffset, beta, y, yOffset)

  def gemvRowMajorNoTree(
      rows: Int,
      cols: Int,
      alpha: Double,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      x: Array[Double],
      xOffset: Int,
      beta: Double,
      y: Array[Double],
      yOffset: Int
  ): Unit =
    VectorGemvBaseline.gemvRowMajorNoTree(rows, cols, alpha, a, aOffset, aRowStride, x, xOffset, beta, y, yOffset)
