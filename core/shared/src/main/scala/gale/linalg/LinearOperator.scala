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
  /** Primitive apply: `into := this * x`. Throws [[LinAlgError]] on dimension
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

  def apply(x: DVec): DVec =
    if cols != x.length then
      throw LinAlgError.VectorLengthMismatch(cols, x.length)
    val out = MutableDVec.zeros(rows)
    applyTo(x, out)
    out.asVec

  def *(x: DVec)(using Backend): DVec =
    apply(x)

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
