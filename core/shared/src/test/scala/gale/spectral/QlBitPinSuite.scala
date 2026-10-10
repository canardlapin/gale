package gale.spectral

import gale.linalg.*

/** Raw-bit pins for the QL symmetric eigen path (`tred2` + implicit QL).
  *
  * The expected hashes fold `doubleToRawLongBits` of the eigenvalues and then
  * the row-major eigenvectors, and were generated from the kernel before the
  * blocked-`tred2` and paired-rotation rewrites (both bit-identical
  * refactors). Inputs are exact binary fractions, so the JVM and Scala.js see
  * the same matrix. A failure means the operation order changed: a deliberate
  * change must say so and regenerate these values, never relax them.
  *
  * The cases cover vectors at small orders (always QL), QL vectors at n = 64
  * through `divideAndConquer = false` (the block-Krylov projected path, which
  * runs QL at any order), values only above the divide-and-conquer order
  * (still QL), and a subnormal-scale tridiagonal. Divide-and-conquer pins at
  * n = 48, 77 and 101 (leaves of 24, 19/20, and 25 plus a 26-row block split
  * into 13/13) guard the merge, the secular solver and `multiplyRowMajor`;
  * they were generated after the power-of-two scaling change.
  *
  * '''Run this suite on the JS lane too.''' The claim that the symmetric
  * eigen kernels give identical bits on the JVM and Scala.js rests on these
  * hashes matching on both platforms, so `coreJS/test` must include it.
  * That last input took the `r == 0` underflow restart until QL gained its
  * tiny-scale lift (`2^600`); its pin was regenerated then (a deliberate bit
  * change for inputs below `2^-600` only), and on the lifted scale it no
  * longer underflows.
  */
class QlBitPinSuite extends munit.FunSuite:
  private def entry(i: Int, j: Int): Double =
    val (r, c) = if i >= j then (i, j) else (j, i)
    ((r + 1) * (c + 2) % 7 - 3).toDouble / 4.0 + (if r == c then r.toDouble / 2.0 else 0.0)

  private def hash(values: DVec, vectors: Option[DMat]): Long =
    var h = 17L
    var i = 0
    while i < values.length do
      h = h * 1000003L + java.lang.Double.doubleToRawLongBits(values(i))
      i += 1
    vectors.foreach: v =>
      var r = 0
      while r < v.rows do
        var c = 0
        while c < v.cols do
          h = h * 1000003L + java.lang.Double.doubleToRawLongBits(v(r, c))
          c += 1
        r += 1
    h

  private val vectorPins = Seq(
    5 -> -3930082448921340683L,
    6 -> 1897052786453825805L,
    7 -> 8183117621020379869L,
    9 -> 6056110939227575803L
  )

  for (n, expected) <- vectorPins do
    test(s"QL values and vectors keep their bits at n=$n (both routes)") {
      val a = Matrix.tabulate(n, n)(entry)
      val ordinary = DenseSpectralKernels.symmetricEigen(a, wantVectors = true).toOption.get
      assertEquals(hash(ordinary.values, ordinary.vectors), expected)
      val workspace =
        DenseSpectralKernels.symmetricEigenWith(a, wantVectors = true, DenseWorkspace.empty).toOption.get
      assertEquals(hash(workspace.values, workspace.vectors), expected)
    }

  test("QL values and vectors keep their bits at n=64 with divide and conquer disabled") {
    val a = Matrix.tabulate(64, 64)(entry)
    val ql = DenseSpectralKernels.symmetricEigen(a, wantVectors = true, divideAndConquer = false).toOption.get
    assertEquals(hash(ql.values, ql.vectors), -4215932567600665687L)
  }

  private val divideAndConquerPins = Seq(
    48 -> -4642693434710438725L,
    77 -> 6401054900483506815L,
    101 -> -3001978090945615853L
  )

  for (n, expected) <- divideAndConquerPins do
    test(s"divide-and-conquer values and vectors keep their bits at n=$n (both routes)") {
      assert(DenseSpectralKernels.usesDivideAndConquer(n))
      val a = Matrix.tabulate(n, n)(entry)
      val ordinary = DenseSpectralKernels.symmetricEigen(a, wantVectors = true).toOption.get
      assertEquals(hash(ordinary.values, ordinary.vectors), expected)
      val workspace =
        DenseSpectralKernels.symmetricEigenWith(a, wantVectors = true, DenseWorkspace.empty).toOption.get
      assertEquals(hash(workspace.values, workspace.vectors), expected)
    }

  private val valuePins = Seq(
    64 -> 2123037278644442257L,
    97 -> 5656159864977540248L
  )

  for (n, expected) <- valuePins do
    test(s"QL values-only keeps its bits at n=$n (both routes)") {
      val a = Matrix.tabulate(n, n)(entry)
      val ordinary = DenseSpectralKernels.symmetricEigen(a, wantVectors = false).toOption.get
      assertEquals(hash(ordinary.values, None), expected)
      val workspace =
        DenseSpectralKernels.symmetricEigenWith(a, wantVectors = false, DenseWorkspace.empty).toOption.get
      assertEquals(hash(workspace.values, None), expected)
    }

  test("QL keeps its bits on a subnormal-scale tridiagonal (solved on the lifted scale)") {
    val tiny = java.lang.Double.MIN_VALUE
    val diagonal = Vec(7 * tiny, -1e-310, 2 * tiny, 1e-310, -tiny, -2 * tiny)
    val offDiagonal = Vec(3 * tiny, 0.0, 5 * tiny, 2 * tiny, -2 * tiny)
    val result =
      DenseSpectralKernels.symmetricTridiagonalEigen(diagonal, offDiagonal, wantVectors = true).toOption.get
    assertEquals(hash(result.values, result.vectors), 8283310461317691503L)
  }
