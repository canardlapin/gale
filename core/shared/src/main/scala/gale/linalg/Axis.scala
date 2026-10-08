package gale.linalg

/** The direction of a per-axis matrix reduction, named after '''what the result
  * is indexed by''':
  *
  *   - [[Axis.Rows]] reduces '''within each row''' and returns one value per
  *     row (a vector of length `rows`). It corresponds to Breeze's
  *     `sum(A(*, ::))` and NumPy's `A.sum(axis=1)`.
  *   - [[Axis.Cols]] reduces '''within each column''' and returns one value per
  *     column (a vector of length `cols`). It corresponds to Breeze's
  *     `sum(A(::, *))` (which Breeze returns as a transposed row vector) and
  *     NumPy's `A.sum(axis=0)`.
  *
  * The same reading holds for shape-preserving operations: a softmax with
  * `Axis.Rows` normalizes each row independently.
  */
enum Axis:
  /** One result per row: reduce across the columns of each row. */
  case Rows

  /** One result per column: reduce down the rows of each column. */
  case Cols
