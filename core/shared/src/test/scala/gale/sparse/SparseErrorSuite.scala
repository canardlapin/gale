package gale.sparse

import gale.linalg.*

class SparseErrorSuite extends munit.FunSuite:
  // Item 13: argument-validation failures in the sparse API must surface as
  // LinAlgError, consistent with the rest of the library, rather than a bare
  // java IllegalArgumentException.
  test("DuplicatePolicy.Error reports a LinAlgError on duplicate entries") {
    intercept[LinAlgError.InvalidArgument] {
      Sparse.coo(2, 2).add(0, 0, 1.0).add(0, 0, 2.0).toCOO(DuplicatePolicy.Error)
    }
  }

  test("Sparse.permutation reports a LinAlgError on duplicate targets") {
    intercept[LinAlgError.InvalidArgument] {
      Sparse.permutation(0, 0)
    }
  }
