package gale.interop.ravel

import gale.linalg.{DMat, DVec}
import ravel.{Array1, Array2, BorrowedNDArray, NDArray, Rank}
import scala.annotation.targetName

/** Copy a rank-one Ravel array into an independently owned Gale vector. */
@targetName("fromOwnedRavelVectorCopy")
def fromRavelCopy(array: Array1[Double]): DVec =
  val builder = DVec.newBuilder(array.size)
  var index = 0
  while index < array.size do
    builder(index) = array(index)
    index += 1
  builder.result()

/** Copy a borrowed rank-one Ravel array without retaining its external alias. */
@targetName("fromBorrowedRavelVectorCopy")
def fromRavelCopy(array: BorrowedNDArray[Double, Rank[1]]): DVec =
  val builder = DVec.newBuilder(array.size)
  var index = 0
  while index < array.size do
    builder(index) = array(index)
    index += 1
  builder.result()

/** Copy a rank-two Ravel array into an independently owned row-major Gale matrix. */
@targetName("fromOwnedRavelMatrixCopy")
def fromRavelCopy(array: Array2[Double]): DMat =
  val builder = DMat.newBuilder(array.shape(0), array.shape(1))
  var row = 0
  while row < array.shape(0) do
    var column = 0
    while column < array.shape(1) do
      builder(row, column) = array(row, column)
      column += 1
    row += 1
  builder.result()

/** Copy a borrowed rank-two Ravel array without retaining its external alias. */
@targetName("fromBorrowedRavelMatrixCopy")
def fromRavelCopy(array: BorrowedNDArray[Double, Rank[2]]): DMat =
  val builder = DMat.newBuilder(array.shape(0), array.shape(1))
  var row = 0
  while row < array.shape(0) do
    var column = 0
    while column < array.shape(1) do
      builder(row, column) = array(row, column)
      column += 1
    row += 1
  builder.result()

/** Copy a Gale vector into an independently owned rank-one Ravel array. */
def toRavelCopy(vector: DVec): Array1[Double] =
  NDArray.tabulate(vector.length)(vector.apply)

/** Copy a Gale matrix view into an independently owned rank-two Ravel array. */
def toRavelCopy(matrix: DMat): Array2[Double] =
  NDArray.tabulate(matrix.rows, matrix.cols)(matrix.apply)
