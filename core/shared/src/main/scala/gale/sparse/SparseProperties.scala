package gale.sparse

import gale.linalg.LinAlgError

/** Zero-cost evidence that a sparse value is in its canonical representation. */
opaque type CanonicalSparse[A] <: A = A

object CanonicalSparse:
  private[gale] inline def unsafe[A](value: A): CanonicalSparse[A] = value

  extension [A](wrapped: CanonicalSparse[A])
    inline def value: A = wrapped

extension (matrix: COO)
  /** Assert canonical COO structure without inspecting the matrix. */
  inline def assumeCanonicalSparse: CanonicalSparse[COO] =
    CanonicalSparse.unsafe(matrix)

  /** Verify the COO canonical-format contract. */
  def verifyCanonicalSparse: Either[LinAlgError, CanonicalSparse[COO]] =
    verifyCanonical(matrix, matrix.hasCanonicalFormat, "COO")

extension (matrix: CSR)
  /** Assert canonical CSR structure without inspecting the matrix. */
  inline def assumeCanonicalSparse: CanonicalSparse[CSR] =
    CanonicalSparse.unsafe(matrix)

  /** Verify sorted unique columns and the absence of explicit zeros. */
  def verifyCanonicalSparse: Either[LinAlgError, CanonicalSparse[CSR]] =
    verifyCanonical(matrix, matrix.hasCanonicalFormat, "CSR")

extension (matrix: CSC)
  /** Assert canonical CSC structure without inspecting the matrix. */
  inline def assumeCanonicalSparse: CanonicalSparse[CSC] =
    CanonicalSparse.unsafe(matrix)

  /** Verify sorted unique rows and the absence of explicit zeros. */
  def verifyCanonicalSparse: Either[LinAlgError, CanonicalSparse[CSC]] =
    verifyCanonical(matrix, matrix.hasCanonicalFormat, "CSC")

private def verifyCanonical[A](matrix: A, canonical: Boolean, representation: String)
    : Either[LinAlgError, CanonicalSparse[A]] =
  if canonical then Right(CanonicalSparse.unsafe(matrix))
  else Left(LinAlgError.InvalidArgument(s"$representation matrix is not in canonical sparse format"))
