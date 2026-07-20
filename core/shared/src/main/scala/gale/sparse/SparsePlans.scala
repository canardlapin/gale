package gale.sparse

import gale.linalg.LinAlgError
import gale.platform.DoubleArray
import gale.platform.DoubleArray.*
import gale.platform.IndexArray
import gale.platform.IndexArray.*
import scala.collection.mutable.ArrayBuffer

/** Reusable mutable numeric storage tied to one immutable [[CSRPattern]].
  *
  * No backing array is exposed or adopted. Symbolic plans overwrite the stored
  * values in place; [[snapshot]] makes an independent immutable CSR when a value
  * must outlive later evaluations.
  */
final class CSRValuesDestination private[sparse] (
    val pattern: CSRPattern,
    private[sparse] val values: DoubleArray
):
  def length: Int = pattern.nnz

  def apply(storedIndex: Int): Double =
    checkStoredIndex(storedIndex)
    values(storedIndex)

  def update(storedIndex: Int, value: Double): Unit =
    checkStoredIndex(storedIndex)
    values(storedIndex) = value

  /** Stored value at a matrix coordinate, or zero when the canonical result
    * pattern has no entry there.
    */
  def apply(row: Int, col: Int): Double =
    if row < 0 || row >= pattern.rows then throw LinAlgError.IndexOutOfBounds(row, pattern.rows)
    if col < 0 || col >= pattern.cols then throw LinAlgError.IndexOutOfBounds(col, pattern.cols)
    val storage = pattern.storage
    var p = storage.laneStart(row)
    val end = storage.laneEnd(row)
    while p < end do
      if storage.minorAtUnchecked(p) == col then return values(p)
      p += 1
    0.0

  def fill(value: Double): Unit =
    var p = 0
    while p < length do
      values(p) = value
      p += 1

  /** Allocation-free traversal in deterministic CSR storage order. */
  def foreachStoredEntry(consumer: SparseEntryConsumer): Unit =
    val storage = pattern.storage
    var row = 0
    while row < pattern.rows do
      var p = storage.laneStart(row)
      val end = storage.laneEnd(row)
      while p < end do
        consumer(row, storage.minorAtUnchecked(p), values(p))
        p += 1
      row += 1

  /** Independent immutable CSR snapshot. Explicit zeros remain stored. */
  def snapshot(): CSR =
    pattern.bindOwned(DoubleArray.copy(values))

  private[sparse] def matches(expected: CSRPattern): Boolean =
    pattern.sharesStorageWith(expected)

  private def checkStoredIndex(index: Int): Unit =
    if index < 0 || index >= length then throw LinAlgError.IndexOutOfBounds(index, length)

/** Symbolic union of two exact canonical CSR patterns.
  *
  * Analysis merges structure once and records, for every result position, the
  * corresponding left/right stored indices (`-1` for a structural absence).
  * Numeric replay therefore performs one fixed pass with no index search.
  */
final class CSRUnionPlan private (
    val leftPattern: CSRPattern,
    val rightPattern: CSRPattern,
    val resultPattern: CSRPattern,
    private val leftByResult: IndexArray,
    private val rightByResult: IndexArray
):
  def newDestination(initialValue: Double = 0.0): CSRValuesDestination =
    resultPattern.valuesDestination(initialValue)

  /** Allocate one owned numeric result while sharing the analyzed result
    * pattern. Numeric cancellation remains an explicit stored zero.
    */
  def evaluate(
      left: CSR,
      right: CSR,
      leftScale: Double = 1.0,
      rightScale: Double = 1.0
  ): Either[LinAlgError, CSR] =
    validateInputs(left, right) match
      case Left(error) => Left(error)
      case Right(()) =>
        val out = DoubleArray.alloc(resultPattern.nnz)
        evaluateValues(left.values, right.values, out, leftScale, rightScale)
        Right(resultPattern.bindOwned(out))

  /** Replay into caller-owned numeric storage. On success no result structure or
    * numeric array is allocated; callers snapshot only when immutable ownership
    * is required.
    */
  def evaluateInto(
      left: CSR,
      right: CSR,
      destination: CSRValuesDestination,
      leftScale: Double = 1.0,
      rightScale: Double = 1.0
  ): Either[LinAlgError, Unit] =
    validateInputs(left, right) match
      case Left(error) => Left(error)
      case Right(()) =>
        if !destination.matches(resultPattern) then
          Left(LinAlgError.InvalidArgument("union destination pattern differs from the analyzed result pattern"))
        else
          evaluateValues(left.values, right.values, destination.values, leftScale, rightScale)
          Right(())

  private def validateInputs(left: CSR, right: CSR): Either[LinAlgError, Unit] =
    if !left.sharesPatternStorage(leftPattern) then
      Left(LinAlgError.InvalidArgument("union left input pattern differs from the pattern used during analysis"))
    else if !right.sharesPatternStorage(rightPattern) then
      Left(LinAlgError.InvalidArgument("union right input pattern differs from the pattern used during analysis"))
    else Right(())

  private def evaluateValues(
      leftValues: DoubleArray,
      rightValues: DoubleArray,
      out: DoubleArray,
      leftScale: Double,
      rightScale: Double
  ): Unit =
    var p = 0
    while p < out.length do
      var value = 0.0
      val leftIndex = leftByResult(p)
      if leftIndex >= 0 && leftScale != 0.0 then value += leftScale * leftValues(leftIndex)
      val rightIndex = rightByResult(p)
      if rightIndex >= 0 && rightScale != 0.0 then value += rightScale * rightValues(rightIndex)
      out(p) = value
      p += 1

object CSRUnionPlan:
  /** Checked symbolic analysis. Inputs must be same-shaped canonical patterns;
    * the result pattern is deterministic and canonical.
    */
  def analyze(left: CSRPattern, right: CSRPattern): Either[LinAlgError, CSRUnionPlan] =
    if left.rows != right.rows || left.cols != right.cols then
      Left(
        LinAlgError.InvalidArgument(
          s"CSR union requires equal shapes, got ${left.rows}x${left.cols} and ${right.rows}x${right.cols}"
        )
      )
    else if !left.hasCanonicalFormat then
      Left(LinAlgError.InvalidArgument("CSR union left pattern must be canonical"))
    else if !right.hasCanonicalFormat then
      Left(LinAlgError.InvalidArgument("CSR union right pattern must be canonical"))
    else
      val offsets = new Array[Int](left.rows + 1)
      val columns = ArrayBuffer.empty[Int]
      val leftMap = ArrayBuffer.empty[Int]
      val rightMap = ArrayBuffer.empty[Int]
      val leftStorage = left.storage
      val rightStorage = right.storage
      var row = 0
      while row < left.rows do
        offsets(row) = columns.length
        var pa = leftStorage.laneStart(row)
        val aEnd = leftStorage.laneEnd(row)
        var pb = rightStorage.laneStart(row)
        val bEnd = rightStorage.laneEnd(row)
        while pa < aEnd || pb < bEnd do
          if columns.length == Int.MaxValue then
            return Left(
              LinAlgError.InvalidArgument(
                s"CSR union result pattern exceeds ${Int.MaxValue} addressable entries"
              )
            )
          if pb >= bEnd || (pa < aEnd && leftStorage.minorAtUnchecked(pa) < rightStorage.minorAtUnchecked(pb)) then
            columns += leftStorage.minorAtUnchecked(pa)
            leftMap += pa
            rightMap += -1
            pa += 1
          else if pa >= aEnd || rightStorage.minorAtUnchecked(pb) < leftStorage.minorAtUnchecked(pa) then
            columns += rightStorage.minorAtUnchecked(pb)
            leftMap += -1
            rightMap += pb
            pb += 1
          else
            columns += leftStorage.minorAtUnchecked(pa)
            leftMap += pa
            rightMap += pb
            pa += 1
            pb += 1
        row += 1
      offsets(left.rows) = columns.length
      CSRPattern.checked(left.rows, left.cols, offsets, columns.toArray).map: result =>
        new CSRUnionPlan(
          left,
          right,
          result,
          IndexArray.fromArray(leftMap.toArray),
          IndexArray.fromArray(rightMap.toArray)
        )

/** Symbolic sparse matrix-product plan for exact canonical CSR patterns.
  *
  * Analysis computes the Boolean product pattern and stores every contributing
  * pair of input positions grouped by output position. Numeric replay is a
  * deterministic sequence of dot products with no sparse lookup or scratch.
  */
final class CSRProductPlan private (
    val leftPattern: CSRPattern,
    val rightPattern: CSRPattern,
    val resultPattern: CSRPattern,
    private val termOffsets: IndexArray,
    private val leftTerms: IndexArray,
    private val rightTerms: IndexArray
):
  def contributionCount: Int = leftTerms.length

  def newDestination(initialValue: Double = 0.0): CSRValuesDestination =
    resultPattern.valuesDestination(initialValue)

  /** Allocate one owned numeric product while sharing the analyzed result
    * pattern. Numeric cancellation remains an explicit stored zero.
    */
  def evaluate(
      left: CSR,
      right: CSR,
      scale: Double = 1.0
  ): Either[LinAlgError, CSR] =
    validateInputs(left, right) match
      case Left(error) => Left(error)
      case Right(()) =>
        val out = DoubleArray.alloc(resultPattern.nnz)
        evaluateValues(left.values, right.values, out, scale)
        Right(resultPattern.bindOwned(out))

  /** Replay into caller-owned numeric storage without allocating result
    * structure, numeric arrays, or algorithm scratch on the success path.
    */
  def evaluateInto(
      left: CSR,
      right: CSR,
      destination: CSRValuesDestination,
      scale: Double = 1.0
  ): Either[LinAlgError, Unit] =
    validateInputs(left, right) match
      case Left(error) => Left(error)
      case Right(()) =>
        if !destination.matches(resultPattern) then
          Left(LinAlgError.InvalidArgument("product destination pattern differs from the analyzed result pattern"))
        else
          evaluateValues(left.values, right.values, destination.values, scale)
          Right(())

  private def validateInputs(left: CSR, right: CSR): Either[LinAlgError, Unit] =
    if !left.sharesPatternStorage(leftPattern) then
      Left(LinAlgError.InvalidArgument("product left input pattern differs from the pattern used during analysis"))
    else if !right.sharesPatternStorage(rightPattern) then
      Left(LinAlgError.InvalidArgument("product right input pattern differs from the pattern used during analysis"))
    else Right(())

  private def evaluateValues(
      leftValues: DoubleArray,
      rightValues: DoubleArray,
      out: DoubleArray,
      scale: Double
  ): Unit =
    var output = 0
    while output < out.length do
      var sum = 0.0
      if scale != 0.0 then
        var term = termOffsets(output)
        val end = termOffsets(output + 1)
        while term < end do
          sum += leftValues(leftTerms(term)) * rightValues(rightTerms(term))
          term += 1
        sum *= scale
      out(output) = sum
      output += 1

object CSRProductPlan:
  /** Checked symbolic analysis of the Boolean CSR product. Inputs must be
    * canonical and inner dimensions compatible. Contribution mappings are
    * rejected if their addressable count would overflow an `Int`.
    */
  def analyze(left: CSRPattern, right: CSRPattern): Either[LinAlgError, CSRProductPlan] =
    if left.cols != right.rows then
      Left(
        LinAlgError.InvalidArgument(
          s"CSR product requires left.cols == right.rows, got ${left.rows}x${left.cols} and ${right.rows}x${right.cols}"
        )
      )
    else if !left.hasCanonicalFormat then
      Left(LinAlgError.InvalidArgument("CSR product left pattern must be canonical"))
    else if !right.hasCanonicalFormat then
      Left(LinAlgError.InvalidArgument("CSR product right pattern must be canonical"))
    else
      analyzeCanonical(left, right)

  private def analyzeCanonical(left: CSRPattern, right: CSRPattern): Either[LinAlgError, CSRProductPlan] =
    val leftStorage = left.storage
    val rightStorage = right.storage
    val resultOffsets = new Array[Int](left.rows + 1)
    val resultColumns = ArrayBuffer.empty[Int]
    val seenInRow = Array.fill(right.cols)(-1)
    val touched = ArrayBuffer.empty[Int]
    var totalTerms = 0L
    var row = 0
    while row < left.rows do
      resultOffsets(row) = resultColumns.length
      touched.clear()
      var pa = leftStorage.laneStart(row)
      val aEnd = leftStorage.laneEnd(row)
      while pa < aEnd do
        val inner = leftStorage.minorAtUnchecked(pa)
        var pb = rightStorage.laneStart(inner)
        val bEnd = rightStorage.laneEnd(inner)
        while pb < bEnd do
          totalTerms += 1L
          if totalTerms > Int.MaxValue.toLong then
            return Left(
              LinAlgError.InvalidArgument(
                s"CSR product symbolic contribution count exceeds ${Int.MaxValue} addressable terms"
              )
            )
          val col = rightStorage.minorAtUnchecked(pb)
          if seenInRow(col) != row then
            seenInRow(col) = row
            touched += col
          pb += 1
        pa += 1
      val sorted = touched.toArray
      scala.util.Sorting.quickSort(sorted)
      var i = 0
      while i < sorted.length do
        if resultColumns.length == Int.MaxValue then
          return Left(
            LinAlgError.InvalidArgument(
              s"CSR product result pattern exceeds ${Int.MaxValue} addressable entries"
            )
          )
        resultColumns += sorted(i)
        i += 1
      row += 1
    resultOffsets(left.rows) = resultColumns.length

    val columns = resultColumns.toArray
    CSRPattern.checked(left.rows, right.cols, resultOffsets, columns) match
      case Left(error) => Left(error)
      case Right(resultPattern) =>
        val outputCount = columns.length
        val termCounts = new Array[Int](outputCount)
        val positionByColumn = Array.fill(right.cols)(-1)

        row = 0
        while row < left.rows do
          var output = resultOffsets(row)
          val outputEnd = resultOffsets(row + 1)
          while output < outputEnd do
            positionByColumn(columns(output)) = output
            output += 1
          var pa = leftStorage.laneStart(row)
          val aEnd = leftStorage.laneEnd(row)
          while pa < aEnd do
            val inner = leftStorage.minorAtUnchecked(pa)
            var pb = rightStorage.laneStart(inner)
            val bEnd = rightStorage.laneEnd(inner)
            while pb < bEnd do
              val resultIndex = positionByColumn(rightStorage.minorAtUnchecked(pb))
              termCounts(resultIndex) += 1
              pb += 1
            pa += 1
          row += 1

        val offsets = new Array[Int](outputCount + 1)
        var sum = 0
        var output = 0
        while output < outputCount do
          offsets(output) = sum
          sum += termCounts(output)
          output += 1
        offsets(outputCount) = sum
        val next = offsets.clone()
        val leftTerms = IndexArray.alloc(sum)
        val rightTerms = IndexArray.alloc(sum)

        row = 0
        while row < left.rows do
          output = resultOffsets(row)
          val outputEnd = resultOffsets(row + 1)
          while output < outputEnd do
            positionByColumn(columns(output)) = output
            output += 1
          var pa = leftStorage.laneStart(row)
          val aEnd = leftStorage.laneEnd(row)
          while pa < aEnd do
            val inner = leftStorage.minorAtUnchecked(pa)
            var pb = rightStorage.laneStart(inner)
            val bEnd = rightStorage.laneEnd(inner)
            while pb < bEnd do
              val resultIndex = positionByColumn(rightStorage.minorAtUnchecked(pb))
              val termIndex = next(resultIndex)
              leftTerms(termIndex) = pa
              rightTerms(termIndex) = pb
              next(resultIndex) += 1
              pb += 1
            pa += 1
          row += 1

        Right(
          new CSRProductPlan(
            left,
            right,
            resultPattern,
            IndexArray.fromArray(offsets),
            leftTerms,
            rightTerms
          )
        )
