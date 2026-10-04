package gale.linalg

import gale.spectral.FiniteInput
import gale.spectral.SpectralBackend
import gale.spectral.SvdCutoff
import gale.spectral.TruncatedSvd

final case class ReducedRowEchelon(matrix: DMat, pivotColumns: Vector[Int]):
  def rank: Int = pivotColumns.length

enum NullSpaceForm:
  /** Row-reduced null equations, transposed into basis columns. Sparse means interpretable elimination structure, not
    * globally minimal support.
    */
  case Sparse
  case Orthonormal

final case class NullSpaceBasis(basis: DMat, rank: Int, cutoff: Double):
  def nullity: Int = basis.cols

object RowReduction:
  /** Ordinary Gauss-Jordan reduction with an absolute pivot tolerance in the input's units. This tolerance is not an
    * SVD rank policy. Rows are partially pivoted; columns stay in input order; ties choose the first row. Normalized
    * coefficients are not pruned using the input's absolute tolerance.
    */
  def reducedRowEchelon(a: DMat, tolerance: Double): Either[LinAlgError, ReducedRowEchelon] =
    if !tolerance.isFinite || tolerance < 0.0 then
      Left(LinAlgError.InvalidArgument("pivot tolerance must be finite and non-negative"))
    else
      FiniteInput
        .matrix(a)
        .flatMap: _ =>
          val result = reduce(a, tolerance)
          FiniteInput.matrix(result.matrix).map(_ => result)

  /** RREF of the numerically retained row space, padded to the input shape. Spectral rank is selected before
    * elimination; it is never inferred from Gauss-Jordan pivot magnitudes. For deliberate truncation this is RREF(A_k).
    */
  def reducedRowEchelon(a: DMat, cutoff: SvdCutoff = SvdCutoff.Default)(using
      SpectralBackend
  ): Either[LinAlgError, ReducedRowEchelon] =
    TruncatedSvd
      .factor(a, cutoff)
      .flatMap: factor =>
        retainedReduction(factor).map: result =>
          ReducedRowEchelon(
            DMat.tabulate(a.rows, a.cols)((i, j) => if i < result.rank then result.matrix(i, j) else 0.0),
            result.pivotColumns
          )

  /** Complete numerical null space of A_k, including n-m missing directions for wide inputs. Sparse uses the canonical
    * RREF of null equations (the same exact-arithmetic convention as RREF(I-V_k V_k^T)); Orthonormal uses
    * twice-reorthogonalized Gram-Schmidt on the free-variable basis. No n x n projector is constructed. The basis
    * itself is n x (n-rank).
    */
  def nullSpaceBasis(
      a: DMat,
      cutoff: SvdCutoff = SvdCutoff.Default,
      form: NullSpaceForm = NullSpaceForm.Sparse
  )(using SpectralBackend): Either[LinAlgError, NullSpaceBasis] =
    TruncatedSvd.factor(a, cutoff).flatMap(nullSpaceBasis(_, form))

  /** Reuse the same rank selection as solves, projections and covariance. */
  def nullSpaceBasis(factor: TruncatedSvd, form: NullSpaceForm): Either[LinAlgError, NullSpaceBasis] =
    retainedReduction(factor).flatMap: reduced =>
      val n = factor.coefficientCount
      val free = (0 until n).filterNot(reduced.pivotColumns.contains).toVector
      val freeRows = DMat.tabulate(free.length, n): (i, j) =>
        val pivot = reduced.pivotColumns.indexOf(j)
        if pivot >= 0 then -reduced.matrix(pivot, free(i))
        else if j == free(i) then 1.0
        else 0.0
      FiniteInput
        .matrix(freeRows)
        .flatMap: _ =>
          val basisResult = form match
            case NullSpaceForm.Sparse =>
              val nullRows = reduce(freeRows, SvdCutoff.MachineEpsilon * math.max(1, n))
              if nullRows.rank != free.length then
                Left(LinAlgError.InvalidArgument("null-space basis lost directions during elimination"))
              else Right(nullRows.matrix.t)
            case NullSpaceForm.Orthonormal => Right(orthonormalize(freeRows))
          basisResult.flatMap: basis =>
            FiniteInput.matrix(basis).map(_ => NullSpaceBasis(basis, factor.rank, factor.cutoff))

  private def retainedReduction(factor: TruncatedSvd): Either[LinAlgError, ReducedRowEchelon] =
    val result = reduce(factor.vt, SvdCutoff.MachineEpsilon * math.max(1, factor.coefficientCount))
    if result.rank != factor.rank then
      Left(LinAlgError.InvalidArgument("retained row space lost rank during elimination"))
    else FiniteInput.matrix(result.matrix).map(_ => result)

  private def reduce(input: DMat, tolerance: Double): ReducedRowEchelon =
    val rows = input.rows
    val cols = input.cols
    val values = Array.tabulate(rows * cols)(i => input(i / cols, i % cols))
    val pivots = Vector.newBuilder[Int]
    var pivotRow = 0
    var column = 0
    while pivotRow < rows && column < cols do
      var best = pivotRow
      var r = pivotRow + 1
      while r < rows do
        if math.abs(values(r * cols + column)) > math.abs(values(best * cols + column)) then best = r
        r += 1
      if math.abs(values(best * cols + column)) <= tolerance then
        r = pivotRow
        while r < rows do
          values(r * cols + column) = 0.0
          r += 1
      else
        var c = 0
        while c < cols do
          val held = values(pivotRow * cols + c)
          values(pivotRow * cols + c) = values(best * cols + c)
          values(best * cols + c) = held
          c += 1
        val pivot = values(pivotRow * cols + column)
        c = 0
        while c < cols do
          values(pivotRow * cols + c) /= pivot
          c += 1
        values(pivotRow * cols + column) = 1.0
        r = 0
        while r < rows do
          if r != pivotRow then
            val multiplier = values(r * cols + column)
            c = 0
            while c < cols do
              values(r * cols + c) -= multiplier * values(pivotRow * cols + c)
              c += 1
            values(r * cols + column) = 0.0
          r += 1
        pivots += column
        pivotRow += 1
      column += 1
    ReducedRowEchelon(
      DMat.tabulate(rows, cols)((i, j) => if i < pivotRow then values(i * cols + j) else 0.0),
      pivots.result()
    )

  private def orthonormalize(rows: DMat): DMat =
    val n = rows.cols
    val q = rows.rows
    val columns = Array.ofDim[Double](q, n)
    var k = 0
    while k < q do
      // Scale first to avoid overflow for large free-variable coefficients.
      val scale = (0 until n).map(i => math.abs(rows(k, i))).max
      var i = 0
      while i < n do
        columns(k)(i) = rows(k, i) / scale
        i += 1
      var pass = 0
      while pass < 2 do
        var previous = 0
        while previous < k do
          var dot = 0.0
          i = 0
          while i < n do
            dot += columns(previous)(i) * columns(k)(i)
            i += 1
          i = 0
          while i < n do
            columns(k)(i) -= dot * columns(previous)(i)
            i += 1
          previous += 1
        pass += 1
      val length = math.sqrt(columns(k).map(x => x * x).sum)
      i = 0
      while i < n do
        columns(k)(i) /= length
        i += 1
      k += 1
    DMat.tabulate(n, q)((i, j) => columns(j)(i))
