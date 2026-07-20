package gale.sparse

import gale.linalg.LinAlgError

object MatrixMarket:
  private val Header = "%%MatrixMarket matrix coordinate real general"

  def writeCoordinate(A: CSR): String =
    val builder = new StringBuilder
    builder.append(Header).append('\n')
    builder.append(A.rows).append(' ').append(A.cols).append(' ').append(A.nnz).append('\n')
    A.foreachStoredEntry { (row, col, value) =>
      builder
        .append(row + 1)
        .append(' ')
        .append(col + 1)
        .append(' ')
        .append(value)
        .append('\n')
    }
    builder.toString

  def readCoordinate(text: String): CSR =
    val lines = text
      .split('\n')
      .map(_.trim)
      .filter(_.nonEmpty)
      .toVector
    require(lines.nonEmpty, "missing Matrix Market header")
    validateHeader(lines.head)
    val data = lines.tail.filterNot(_.startsWith("%"))
    require(data.nonEmpty, "missing Matrix Market shape line")
    val shape = fields(data.head)
    require(shape.length == 3, "Matrix Market coordinate shape must have rows cols nnz")
    val rows = shape(0).toInt
    val cols = shape(1).toInt
    val expectedNnz = shape(2).toInt
    val builder = Sparse.coo(rows, cols)
    data.tail.foreach { line =>
      val entry = fields(line)
      require(entry.length == 3, s"invalid Matrix Market entry: $line")
      builder.add(entry(0).toInt - 1, entry(1).toInt - 1, entry(2).toDouble)
    }
    require(builder.nnz == expectedNnz, s"expected $expectedNnz entries, got ${builder.nnz}")
    builder.toCSR()

  /** Accept only `matrix coordinate real|integer general`. Any other object,
    * storage format, field, or symmetry is rejected rather than mis-parsed.
    */
  private def validateHeader(line: String): Unit =
    val tokens = fields(line)
    if tokens.length != 5 || tokens(0) != "%%MatrixMarket" then
      throw LinAlgError.UnsupportedRepresentation(s"invalid Matrix Market banner: $line")
    val obj = tokens(1).toLowerCase
    val format = tokens(2).toLowerCase
    val field = tokens(3).toLowerCase
    val symmetry = tokens(4).toLowerCase
    if obj != "matrix" then
      throw LinAlgError.UnsupportedRepresentation(s"unsupported Matrix Market object '$obj' (only 'matrix')")
    if format != "coordinate" then
      throw LinAlgError.UnsupportedRepresentation(s"unsupported Matrix Market format '$format' (only 'coordinate')")
    if field != "real" && field != "integer" then
      throw LinAlgError.UnsupportedRepresentation(s"unsupported Matrix Market field '$field' (only 'real' or 'integer')")
    if symmetry != "general" then
      throw LinAlgError.UnsupportedRepresentation(s"unsupported Matrix Market symmetry '$symmetry' (only 'general')")

  private def fields(line: String): Array[String] =
    line.split("\\s+").filter(_.nonEmpty)
