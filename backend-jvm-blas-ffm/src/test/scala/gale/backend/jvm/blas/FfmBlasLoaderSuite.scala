package gale.backend.jvm.blas

import gale.backend.{Backend, Capability}
import gale.backend.jvm.`native`.{Layout, NativeDMat}
import gale.linalg.Matrix

import java.util.List

class FfmBlasLoaderSuite extends munit.FunSuite:
  test("default loader reports the selected library and an honest capability set"):
    val backend = FfmBlasBackend.load().fold(throw _, identity)
    try
      assert(backend.libraryInfo.name.nonEmpty)
      assert(backend.capabilities.contains(Capability.NativeBlas))
      assert(!backend.capabilities.contains(Capability.NativeLapack))
      assert(backend.denseFactorizations.isEmpty)
      assert(backend.spectral.isEmpty)
    finally backend.close()

  test("given import resolves to an FFM BLAS backend"):
    import gale.backend.jvm.blas.given
    assert(summon[Backend].capabilities.contains(Capability.NativeBlas))

  test("zero threshold routes the public facade through native CBLAS"):
    val backend = FfmBlasBackend
      .load(thresholds = Some(FfmBlasThresholds(nativeGemmMinFlops = 0L)))
      .fold(throw _, identity)
    try
      val n = 64 // arithmetic intensity exceeds the copy-amortization guard
      val a = Matrix.tabulate(n, n)((i, j) => ((i * 7 + j * 3) % 13 - 6).toDouble / 7.0)
      val b = Matrix.tabulate(n, n)((i, j) => ((i * 5 - j * 11) % 17 - 8).toDouble / 9.0)
      val expected = a.*(b)(using gale.backend.PureBackend)
      val actual = a.*(b)(using backend)
      var i = 0
      while i < actual.rows do
        var j = 0
        while j < actual.cols do
          assertEqualsDouble(actual(i, j), expected(i, j), 1e-12)
          j += 1
        i += 1

      val expectedTranspose = a.*(b.t)(using gale.backend.PureBackend)
      val actualTranspose = a.*(b.t)(using backend)
      val expectedSyrk = a.t.*(a)(using gale.backend.PureBackend)
      val actualSyrk = a.t.*(a)(using backend)
      i = 0
      while i < n do
        var j = 0
        while j < n do
          assertEqualsDouble(actualTranspose(i, j), expectedTranspose(i, j), 1e-10)
          assertEqualsDouble(actualSyrk(i, j), expectedSyrk(i, j), 1e-10)
          j += 1
        i += 1
    finally backend.close()

  test("known optimized libraries get a conservative route; unknown BLAS stays disabled"):
    assertEquals(FfmBlasThresholds.forLibrary("Accelerate").nativeGemmMinFlops, 256L * 256L * 256L)
    assertEquals(FfmBlasThresholds.forLibrary("libopenblas.so").nativeGemmMinFlops, 256L * 256L * 256L)
    assertEquals(FfmBlasThresholds.forLibrary("libblas.so.3").nativeGemmMinFlops, Long.MaxValue)

  test("copy-free NativeDMat GEMM supports row, column, and mixed layouts"):
    val backend = FfmBlasBackend.load().fold(throw _, identity)
    try
      for (aLayout, bLayout, cLayout) <- Seq(
        (Layout.RowMajor, Layout.RowMajor, Layout.RowMajor),
        (Layout.ColMajor, Layout.ColMajor, Layout.ColMajor),
        (Layout.ColMajor, Layout.RowMajor, Layout.ColMajor),
        (Layout.RowMajor, Layout.ColMajor, Layout.RowMajor)
      ) do
        val a = NativeDMat.allocate(3, 4, aLayout)
        val b = NativeDMat.allocate(4, 2, bLayout)
        val c = NativeDMat.allocate(3, 2, cLayout)
        try
          var i = 0
          while i < 3 do
            var j = 0
            while j < 4 do
              a(i, j) = (i * 4 + j - 3).toDouble / 2.0
              j += 1
            i += 1
          i = 0
          while i < 4 do
            var j = 0
            while j < 2 do
              b(i, j) = (i * 3 - j * 2 + 1).toDouble / 3.0
              j += 1
            i += 1
          backend.gemm(a, b, c)
          val expected = a.toHeap.*(b.toHeap)(using gale.backend.PureBackend)
          val actual = c.toHeap
          i = 0
          while i < 3 do
            var j = 0
            while j < 2 do
              assertEqualsDouble(actual(i, j), expected(i, j), 1e-12)
              j += 1
            i += 1
        finally
          c.close(); b.close(); a.close()
    finally backend.close()

  test("skinny products decline native copying even when the facade threshold is forced"):
    val backend = FfmBlasBackend
      .load(thresholds = Some(FfmBlasThresholds(nativeGemmMinFlops = 0L)))
      .fold(throw _, identity)
    try
      val shared = 20000
      val a = Matrix.tabulate(1, shared)((_, j) => (j % 7 - 3).toDouble)
      val b = Matrix.tabulate(shared, 1)((i, _) => (i % 5 - 2).toDouble)
      val expected = a.*(b)(using gale.backend.PureBackend)
      val actual = a.*(b)(using backend)
      assertEqualsDouble(actual(0, 0), expected(0, 0), 1e-12)
    finally backend.close()

  test("explicit loader failure lists the attempted candidate"):
    val missing = "/definitely/not/a/cblas/library"
    val error = intercept[IllegalStateException](CblasBindings.load(List.of(missing)))
    assert(error.getMessage.contains(missing))
