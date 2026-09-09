package gale.linalg

import gale.platform.DoubleArray.*

/** Scaling of the columns of a DCT-II basis matrix. */
enum DctNormalization:
  case Orthonormal, Unscaled

/** Explicit DCT-II basis columns. This constructs a dense matrix; it does not perform a fast cosine transform. */
object DctBasis:

  /** Construct columns `first` through `first + count - 1` of the `samples`-point DCT-II basis.
    *
    * Entry `(row, col)` is `cos(pi * (row + 0.5) * (first + col) / samples)` times the selected scale. Orthonormal
    * scaling is `1 / sqrt(samples)` for component zero and `sqrt(2 / samples)` otherwise. Unscaled columns have scale
    * one. Rows index samples and columns index components, both starting at zero.
    *
    * `samples` must be positive and the selected range must lie in `[0, samples]`. An empty range is valid, including
    * `first == samples`, and returns a `samples` by zero matrix without traversing its rows. Invalid ranges and shapes
    * exceeding `Int.MaxValue` elements return [[LinAlgError.InvalidArgument]]. Only the requested columns are
    * allocated: time and storage are O(samples * count).
    */
  def columns(
      samples: Int,
      first: Int,
      count: Int,
      normalization: DctNormalization = DctNormalization.Orthonormal
  ): Either[LinAlgError, DMat] =
    if samples <= 0 then Left(LinAlgError.InvalidArgument(s"DCT samples must be positive, got $samples"))
    else if first < 0 || count < 0 || first.toLong + count.toLong > samples.toLong then
      Left(LinAlgError.InvalidArgument(s"DCT component range first=$first count=$count must lie in [0, $samples]"))
    else if samples.toLong * count.toLong > Int.MaxValue.toLong then
      Left(LinAlgError.InvalidArgument(s"DCT matrix ${samples}x${count} exceeds ${Int.MaxValue} storable elements"))
    else
      val out = DMat.zeros(samples, count)
      if count > 0 then
        val dcScale = normalization match
          case DctNormalization.Orthonormal => 1.0 / math.sqrt(samples.toDouble)
          case DctNormalization.Unscaled    => 1.0
        val otherScale = normalization match
          case DctNormalization.Orthonormal => math.sqrt(2.0 / samples.toDouble)
          case DctNormalization.Unscaled    => 1.0
        var row = 0
        while row < samples do
          var col = 0
          while col < count do
            val component = first + col
            val angle = math.Pi * (row.toDouble + 0.5) * component.toDouble / samples.toDouble
            val scale = if component == 0 then dcScale else otherScale
            out.data(row * count + col) = math.cos(angle) * scale
            col += 1
          row += 1
      Right(out)
