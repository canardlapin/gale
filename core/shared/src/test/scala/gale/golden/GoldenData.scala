package gale.golden

/** One named, row-major array of a Breeze golden case: an input (`isInput`) or a Breeze reference result.
  *
  * Values are stored as the 16-hex-digit IEEE-754 bit patterns of each double, so a golden replays exactly the
  * inputs and references the generator saw on every platform. Decoding is lazy and cached.
  */
final class GoldenArray(val key: String, val rows: Int, val cols: Int, val isInput: Boolean, hex: String):
  require(hex.length == 16 * rows * cols, s"golden array $key: ${hex.length} hex digits for ${rows}x$cols values")

  lazy val values: Array[Double] =
    val out = new Array[Double](rows * cols)
    var i = 0
    while i < out.length do
      val hi = GoldenArray.parseHex8(hex, 16 * i)
      val lo = GoldenArray.parseHex8(hex, 16 * i + 8)
      out(i) = java.lang.Double.longBitsToDouble((hi.toLong << 32) | (lo.toLong & 0xffffffffL))
      i += 1
    out

  def size: Int = rows * cols

  def apply(i: Int, j: Int): Double = values(i * cols + j)

object GoldenArray:
  private def parseHex8(s: String, from: Int): Int =
    var acc = 0
    var k = from
    while k < from + 8 do
      val c = s.charAt(k)
      val d =
        if c >= '0' && c <= '9' then c - '0'
        else if c >= 'a' && c <= 'f' then c - 'a' + 10
        else throw new IllegalArgumentException(s"bad hex digit '$c'")
      acc = (acc << 4) | d
      k += 1
    acc

/** A Breeze golden case: a family (the operation group the replay dispatches on), a unique name, and its named input
  * and reference arrays.
  */
final case class GoldenCase(family: String, name: String, arrays: IndexedSeq[GoldenArray]):
  private lazy val byKey: Map[String, GoldenArray] = arrays.map(a => a.key -> a).toMap

  def has(key: String): Boolean = byKey.contains(key)

  def array(key: String): GoldenArray =
    byKey.getOrElse(key, throw new NoSuchElementException(s"golden case $name has no array '$key'"))

  def values(key: String): Array[Double] = array(key).values

  def scalar(key: String): Double =
    val a = array(key)
    require(a.size == 1, s"golden case $name: '$key' is not a scalar")
    a.values(0)

/** Constructors used by the generated `BreezeGoldens` source. */
object GoldenSyntax:
  /** An input array, built by the generator from seeds with host-independent arithmetic. */
  def in(key: String, rows: Int, cols: Int, chunks: String*): GoldenArray =
    new GoldenArray(key, rows, cols, isInput = true, chunks.mkString)

  /** A Breeze reference result. */
  def ref(key: String, rows: Int, cols: Int, chunks: String*): GoldenArray =
    new GoldenArray(key, rows, cols, isInput = false, chunks.mkString)
