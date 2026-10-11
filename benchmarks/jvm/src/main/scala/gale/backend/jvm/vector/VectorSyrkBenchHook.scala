package gale.backend.jvm.vector

/** Benchmark-only access to the SIMD `syrk` path, bypassing the kernel's shape
  * guard, so `VectorSyrkJmh.simdForced` can time it at shapes the guard keeps
  * pure. Lives in the benchmarks module (never published) inside the backend's
  * package, so the kernel's `private[vector]` surface stays unchanged.
  */
object VectorSyrkBenchHook:
  def syrkSimdForced(
      m: Int,
      k: Int,
      a: Array[Double],
      aOffset: Int,
      aRowStride: Int,
      c: Array[Double],
      cOffset: Int,
      cRowStride: Int
  ): Unit =
    VectorDenseDoubleKernel.syrkPackedSimd(m, k, a, aOffset, aRowStride, c, cOffset, cRowStride)
