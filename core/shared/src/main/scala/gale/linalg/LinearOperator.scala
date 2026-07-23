package gale.linalg

import gale.backend.Backend

trait LinearOperator[A]:
  def rows: Int
  def cols: Int

  /** Apply the operator: `into := this * x`. Total form; dimension and
    * representation failures are reported as values.
    */
  def applyTo(x: Vec[A], into: MutableVec[A]): Either[LinAlgError, Unit]

trait DoubleLinearOperator extends LinearOperator[Double]:
  /** Primitive apply: `into := this * x`. Throws `LinAlgError` on dimension
    * mismatch; the generic [[applyTo]] wraps this as a total form.
    */
  def applyTo(x: DVec, into: MutableDVec): Unit

  final override def applyTo(x: Vec[Double], into: MutableVec[Double]): Either[LinAlgError, Unit] =
    (x, into) match
      case (xd: DVec, yd: MutableDVec) =>
        try
          applyTo(xd, yd)
          Right(())
        catch case error: LinAlgError => Left(error)
      case _ =>
        Left(LinAlgError.UnsupportedRepresentation("DoubleLinearOperator requires DVec / MutableDVec arguments"))

  def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    throw LinAlgError.UnsupportedOperation("transposeApplyTo")

  /** The explicit real adjoint of this operator. Applying it delegates to
    * [[transposeApplyTo]]; its transpose delegates back to [[applyTo]].
    */
  def adjoint: DoubleLinearOperator =
    LinearOperator.adjoint(this)

  def apply(x: DVec): DVec =
    if cols != x.length then
      throw LinAlgError.VectorLengthMismatch(cols, x.length)
    val out = MutableDVec.zeros(rows)
    applyTo(x, out)
    // Implementations receive the destination and may retain it. Snapshot the
    // result so a user-defined operator cannot mutate an immutable return value.
    out.toVec

  def *(x: DVec)(using Backend): DVec =
    apply(x)

  /** Apply this operator to every column of `input` in one owned output
    * allocation. The result has shape `rows x input.cols`.
    */
  def applyTo(input: DMat): Either[LinAlgError, DMat] =
    if input.rows != cols then
      Left(
        LinAlgError.DimensionMismatch(
          Shape(Rows(cols), Cols(input.cols)),
          input.shape
        )
      )
    else
      try
        val output = DMat.newBuilder(rows, input.cols)
        var col = 0
        while col < input.cols do
          applyTo(input.col(col), output.mutableColumn(col))
          col += 1
        Right(output.result())
      catch case error: LinAlgError => Left(error)

  /** Apply the real adjoint to every column of `input` in one owned output
    * allocation. The result has shape `cols x input.cols`.
    */
  def transposeApplyTo(input: DMat): Either[LinAlgError, DMat] =
    if input.rows != rows then
      Left(
        LinAlgError.DimensionMismatch(
          Shape(Rows(rows), Cols(input.cols)),
          input.shape
        )
      )
    else
      try
        val output = DMat.newBuilder(cols, input.cols)
        var col = 0
        while col < input.cols do
          transposeApplyTo(input.col(col), output.mutableColumn(col))
          col += 1
        Right(output.result())
      catch case error: LinAlgError => Left(error)

  /** `this` after `before`, when the intermediate dimensions agree. */
  def compose(before: DoubleLinearOperator): Either[LinAlgError, DoubleLinearOperator] =
    LinearOperator.compose(this, before)

  /** `after` after `this`, when the intermediate dimensions agree. */
  def andThen(after: DoubleLinearOperator): Either[LinAlgError, DoubleLinearOperator] =
    LinearOperator.compose(after, this)

  def scaled(alpha: Double): DoubleLinearOperator =
    LinearOperator.scaled(this, alpha)

  def restrictRows(indices: IndexedSeq[Int]): Either[LinAlgError, DoubleLinearOperator] =
    LinearOperator.restrictRows(this, indices)

  def restrictColumns(indices: IndexedSeq[Int]): Either[LinAlgError, DoubleLinearOperator] =
    LinearOperator.restrictColumns(this, indices)

/** Matrix-free Kronecker product `left ⊗ right`.
  *
  * Vectors use the same block ordering as [[DMat.kron]]: the right factor's
  * coordinate varies fastest. Construction is total because the product shape
  * is checked before it is narrowed to `Int`.
  */
final class KroneckerLinearOperator private[linalg] (
    val left: DoubleLinearOperator,
    val right: DoubleLinearOperator,
    override val rows: Int,
    override val cols: Int
) extends DoubleLinearOperator:
  override def applyTo(x: DVec, into: MutableDVec): Unit =
    LinearOperator.applyKronecker(left, right, x, into)

  override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
    LinearOperator.applyKronecker(left.adjoint, right.adjoint, x, into)

  override def adjoint: KroneckerLinearOperator =
    new KroneckerLinearOperator(left.adjoint, right.adjoint, cols, rows)

object LinearOperator:
  def fromFunction(rowsValue: Int, colsValue: Int)(f: (DVec, MutableDVec) => Unit): DoubleLinearOperator =
    require(rowsValue >= 0 && colsValue >= 0, "operator shape must be non-negative")
    new DoubleLinearOperator:
      override def rows: Int =
        rowsValue

      override def cols: Int =
        colsValue

      override def applyTo(x: DVec, into: MutableDVec): Unit =
        if x.length != colsValue then
          throw LinAlgError.VectorLengthMismatch(colsValue, x.length)
        if into.length != rowsValue then
          throw LinAlgError.VectorLengthMismatch(rowsValue, into.length)
        f(x, into)

  def tabulate(rows: Int, cols: Int)(f: (Int, Int) => Double): DoubleLinearOperator =
    Matrix.tabulate(rows, cols)(f)

  /** Matrix-free operator with explicit forward and adjoint kernels. */
  def fromFunctions(rowsValue: Int, colsValue: Int)(
      forward: (DVec, MutableDVec) => Unit,
      transpose: (DVec, MutableDVec) => Unit
  ): DoubleLinearOperator =
    require(rowsValue >= 0 && colsValue >= 0, "operator shape must be non-negative")
    new DoubleLinearOperator:
      override def rows: Int =
        rowsValue

      override def cols: Int =
        colsValue

      override def applyTo(x: DVec, into: MutableDVec): Unit =
        requireApplicationShape(x, into, colsValue, rowsValue)
        forward(x, into)

      override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
        requireApplicationShape(x, into, rowsValue, colsValue)
        transpose(x, into)

  def adjoint(operator: DoubleLinearOperator): DoubleLinearOperator =
    new DoubleLinearOperator:
      override def rows: Int =
        operator.cols

      override def cols: Int =
        operator.rows

      override def applyTo(x: DVec, into: MutableDVec): Unit =
        operator.transposeApplyTo(x, into)

      override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
        operator.applyTo(x, into)

  def compose(
      after: DoubleLinearOperator,
      before: DoubleLinearOperator
  ): Either[LinAlgError, DoubleLinearOperator] =
    if before.rows != after.cols then
      Left(
        LinAlgError.DimensionMismatch(
          Shape(Rows(after.cols), Cols(1)),
          Shape(Rows(before.rows), Cols(1))
        )
      )
    else
      Right(
        new DoubleLinearOperator:
          override def rows: Int =
            after.rows

          override def cols: Int =
            before.cols

          override def applyTo(x: DVec, into: MutableDVec): Unit =
            requireApplicationShape(x, into, cols, rows)
            val intermediate = MutableDVec.zeros(before.rows)
            before.applyTo(x, intermediate)
            after.applyTo(intermediate.asVec, into)

          override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
            requireApplicationShape(x, into, rows, cols)
            val intermediate = MutableDVec.zeros(after.cols)
            after.transposeApplyTo(x, intermediate)
            before.transposeApplyTo(intermediate.asVec, into)
      )

  def scaled(operator: DoubleLinearOperator, alpha: Double): DoubleLinearOperator =
    new DoubleLinearOperator:
      override def rows: Int =
        operator.rows

      override def cols: Int =
        operator.cols

      override def applyTo(x: DVec, into: MutableDVec): Unit =
        operator.applyTo(x, into)
        into *= alpha

      override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
        operator.transposeApplyTo(x, into)
        into *= alpha

  /** Matrix-free Kronecker product `left ⊗ right`, ordered identically to
    * [[DMat.kron]]. Neither factor nor the product is materialized.
    */
  def kronecker(
      left: DoubleLinearOperator,
      right: DoubleLinearOperator
  ): Either[LinAlgError, KroneckerLinearOperator] =
    for
      rows <- checkedProduct(left.rows, right.rows, "Kronecker operator row count")
      cols <- checkedProduct(left.cols, right.cols, "Kronecker operator column count")
    yield new KroneckerLinearOperator(left, right, rows, cols)

  def restrictRows(
      operator: DoubleLinearOperator,
      indices: IndexedSeq[Int]
  ): Either[LinAlgError, DoubleLinearOperator] =
    validateSelection(indices, operator.rows).map: selected =>
      new DoubleLinearOperator:
        override def rows: Int =
          selected.length

        override def cols: Int =
          operator.cols

        override def applyTo(x: DVec, into: MutableDVec): Unit =
          requireApplicationShape(x, into, cols, rows)
          val full = MutableDVec.zeros(operator.rows)
          operator.applyTo(x, full)
          var i = 0
          while i < selected.length do
            into(i) = full(selected(i))
            i += 1

        override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
          requireApplicationShape(x, into, rows, cols)
          val full = MutableDVec.zeros(operator.rows)
          var i = 0
          while i < selected.length do
            full(selected(i)) = x(i)
            i += 1
          operator.transposeApplyTo(full.asVec, into)

  def restrictColumns(
      operator: DoubleLinearOperator,
      indices: IndexedSeq[Int]
  ): Either[LinAlgError, DoubleLinearOperator] =
    validateSelection(indices, operator.cols).map: selected =>
      new DoubleLinearOperator:
        override def rows: Int =
          operator.rows

        override def cols: Int =
          selected.length

        override def applyTo(x: DVec, into: MutableDVec): Unit =
          requireApplicationShape(x, into, cols, rows)
          val full = MutableDVec.zeros(operator.cols)
          var i = 0
          while i < selected.length do
            full(selected(i)) = x(i)
            i += 1
          operator.applyTo(full.asVec, into)

        override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
          requireApplicationShape(x, into, rows, cols)
          val full = MutableDVec.zeros(operator.cols)
          operator.transposeApplyTo(x, full)
          var i = 0
          while i < selected.length do
            into(i) = full(selected(i))
            i += 1

  def blockDiagonal(
      operators: IndexedSeq[DoubleLinearOperator]
  ): Either[LinAlgError, DoubleLinearOperator] =
    if operators.isEmpty then
      Left(LinAlgError.InvalidArgument("block diagonal operator list must be non-empty"))
    else
      val rowOffsets = prefixOffsets(operators.map(_.rows))
      val colOffsets = prefixOffsets(operators.map(_.cols))
      Right(
        new DoubleLinearOperator:
          override def rows: Int =
            rowOffsets.last

          override def cols: Int =
            colOffsets.last

          override def applyTo(x: DVec, into: MutableDVec): Unit =
            requireApplicationShape(x, into, cols, rows)
            var block = 0
            while block < operators.length do
              val op = operators(block)
              val input = x.slice(colOffsets(block), colOffsets(block + 1))
              val output = MutableDVec.zeros(op.rows)
              op.applyTo(input, output)
              copyInto(output.asVec, into, rowOffsets(block))
              block += 1

          override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
            requireApplicationShape(x, into, rows, cols)
            var block = 0
            while block < operators.length do
              val op = operators(block)
              val input = x.slice(rowOffsets(block), rowOffsets(block + 1))
              val output = MutableDVec.zeros(op.cols)
              op.transposeApplyTo(input, output)
              copyInto(output.asVec, into, colOffsets(block))
              block += 1
      )

  /** Assemble a rectangular block operator. Every operator in a block row must
    * have the same row count; every operator in a block column must have the
    * same column count.
    */
  def block(
      blocks: IndexedSeq[IndexedSeq[DoubleLinearOperator]]
  ): Either[LinAlgError, DoubleLinearOperator] =
    validateBlocks(blocks).map: layout =>
      new DoubleLinearOperator:
        override def rows: Int =
          layout.rowOffsets.last

        override def cols: Int =
          layout.colOffsets.last

        override def applyTo(x: DVec, into: MutableDVec): Unit =
          requireApplicationShape(x, into, cols, rows)
          var blockRow = 0
          while blockRow < blocks.length do
            val accumulator = MutableDVec.zeros(layout.rowSizes(blockRow))
            var blockCol = 0
            while blockCol < blocks(blockRow).length do
              val op = blocks(blockRow)(blockCol)
              val input = x.slice(layout.colOffsets(blockCol), layout.colOffsets(blockCol + 1))
              val contribution = MutableDVec.zeros(op.rows)
              op.applyTo(input, contribution)
              accumulator += contribution.asVec
              blockCol += 1
            copyInto(accumulator.asVec, into, layout.rowOffsets(blockRow))
            blockRow += 1

        override def transposeApplyTo(x: DVec, into: MutableDVec): Unit =
          requireApplicationShape(x, into, rows, cols)
          into.clear()
          var blockRow = 0
          while blockRow < blocks.length do
            val input = x.slice(layout.rowOffsets(blockRow), layout.rowOffsets(blockRow + 1))
            var blockCol = 0
            while blockCol < blocks(blockRow).length do
              val op = blocks(blockRow)(blockCol)
              val contribution = MutableDVec.zeros(op.cols)
              op.transposeApplyTo(input, contribution)
              addInto(contribution.asVec, into, layout.colOffsets(blockCol))
              blockCol += 1
            blockRow += 1

  private final case class BlockLayout(
      rowSizes: IndexedSeq[Int],
      colSizes: IndexedSeq[Int],
      rowOffsets: IndexedSeq[Int],
      colOffsets: IndexedSeq[Int]
  )

  private def validateBlocks(
      blocks: IndexedSeq[IndexedSeq[DoubleLinearOperator]]
  ): Either[LinAlgError, BlockLayout] =
    if blocks.isEmpty || blocks.headOption.forall(_.isEmpty) then
      Left(LinAlgError.InvalidArgument("block operator layout must be non-empty"))
    else
      val blockCols = blocks.head.length
      if blocks.exists(_.length != blockCols) then
        Left(LinAlgError.InvalidArgument("block operator rows must have equal length"))
      else
        val rowSizes = blocks.map(_.head.rows)
        val colSizes = blocks.head.map(_.cols)
        var row = 0
        var error = Option.empty[LinAlgError]
        while row < blocks.length && error.isEmpty do
          var col = 0
          while col < blockCols && error.isEmpty do
            val op = blocks(row)(col)
            if op.rows != rowSizes(row) then
              error = Some(LinAlgError.InvalidArgument(s"block row $row has inconsistent heights"))
            else if op.cols != colSizes(col) then
              error = Some(LinAlgError.InvalidArgument(s"block column $col has inconsistent widths"))
            col += 1
          row += 1
        error match
          case Some(value) => Left(value)
          case None =>
            Right(
              BlockLayout(
                rowSizes,
                colSizes,
                prefixOffsets(rowSizes),
                prefixOffsets(colSizes)
              )
            )

  private def validateSelection(
      indices: IndexedSeq[Int],
      bound: Int
  ): Either[LinAlgError, IndexedSeq[Int]] =
    if indices.isEmpty then
      Left(LinAlgError.InvalidArgument("operator restriction must be non-empty"))
    else
      val seen = scala.collection.mutable.HashSet.empty[Int]
      var i = 0
      var error = Option.empty[LinAlgError]
      while i < indices.length && error.isEmpty do
        val index = indices(i)
        if index < 0 || index >= bound then
          error = Some(LinAlgError.IndexOutOfBounds(index, bound))
        else if seen.contains(index) then
          error = Some(LinAlgError.InvalidArgument(s"duplicate operator restriction index $index"))
        else
          seen += index
        i += 1
      error match
        case Some(value) => Left(value)
        case None        => Right(indices.toVector)

  private def requireApplicationShape(
      input: DVec,
      output: MutableDVec,
      expectedInput: Int,
      expectedOutput: Int
  ): Unit =
    if input.length != expectedInput then
      throw LinAlgError.VectorLengthMismatch(expectedInput, input.length)
    if output.length != expectedOutput then
      throw LinAlgError.VectorLengthMismatch(expectedOutput, output.length)

  private[linalg] def applyKronecker(
      left: DoubleLinearOperator,
      right: DoubleLinearOperator,
      input: DVec,
      output: MutableDVec
  ): Unit =
    val inputOuter = left.cols
    val inputInner = right.cols
    val outputOuter = left.rows
    val outputInner = right.rows
    requireApplicationShape(input, output, inputOuter * inputInner, outputOuter * outputInner)
    output.clear()

    val innerResults = Array.fill(inputOuter)(MutableDVec.zeros(outputInner))
    var outer = 0
    while outer < inputOuter do
      val segment = input.slice(outer * inputInner, (outer + 1) * inputInner)
      right.applyTo(segment, innerResults(outer))
      outer += 1

    val acrossOuter = MutableDVec.zeros(inputOuter)
    val transformed = MutableDVec.zeros(outputOuter)
    var inner = 0
    while inner < outputInner do
      outer = 0
      while outer < inputOuter do
        acrossOuter(outer) = innerResults(outer)(inner)
        outer += 1
      left.applyTo(acrossOuter.asVec, transformed)
      outer = 0
      while outer < outputOuter do
        output(outer * outputInner + inner) = transformed(outer)
        outer += 1
      inner += 1

  private def checkedProduct(left: Int, right: Int, label: String): Either[LinAlgError, Int] =
    if left < 0 || right < 0 then
      Left(LinAlgError.InvalidArgument(s"$label requires non-negative factors, got $left and $right"))
    else
      val product = left.toLong * right.toLong
      if product > Int.MaxValue.toLong then
        Left(LinAlgError.InvalidArgument(s"$label $product exceeds ${Int.MaxValue}"))
      else Right(product.toInt)

  private def prefixOffsets(sizes: IndexedSeq[Int]): IndexedSeq[Int] =
    val out = Array.ofDim[Int](sizes.length + 1)
    var i = 0
    while i < sizes.length do
      out(i + 1) = out(i) + sizes(i)
      i += 1
    out.toIndexedSeq

  private def copyInto(source: DVec, destination: MutableDVec, offset: Int): Unit =
    var i = 0
    while i < source.length do
      destination(offset + i) = source(i)
      i += 1

  private def addInto(source: DVec, destination: MutableDVec, offset: Int): Unit =
    var i = 0
    while i < source.length do
      destination(offset + i) = destination(offset + i) + source(i)
      i += 1
