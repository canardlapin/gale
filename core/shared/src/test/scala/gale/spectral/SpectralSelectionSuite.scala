package gale.spectral

import gale.linalg.Matrix
import gale.linalg.assumePositiveDefinite
import gale.linalg.orThrow
import gale.solvers.SolverConfig

class SpectralSelectionSuite extends munit.FunSuite:
  test("EigenOrder ships the six directional cases plus BothEnds, no imaginary") {
    assertEquals(EigenOrder.values.length, 7)
    assert(EigenOrder.values.contains(EigenOrder.BothEnds))
    assert(!EigenOrder.values.exists(_.toString.contains("Imaginary")))
  }

  test("SingularOrder is Largest / Smallest") {
    assertEquals(SingularOrder.values.toSet, Set(SingularOrder.Largest, SingularOrder.Smallest))
  }

  test("EigenSelection carries selection data without validating it") {
    EigenSelection.Count(3, EigenOrder.LargestAlgebraic) match
      case EigenSelection.Count(k, order) =>
        assertEquals(k, 3)
        assertEquals(order, EigenOrder.LargestAlgebraic)
      case other => fail(s"expected Count, got $other")

    // No constructor invariants: k <= 0 and reversed ranges are the solver's job.
    assertEquals(EigenSelection.Count(0, EigenOrder.LargestMagnitude), EigenSelection.Count(0, EigenOrder.LargestMagnitude))
    assertEquals(EigenSelection.IndexRange(4, 1), EigenSelection.IndexRange(4, 1))

    EigenSelection.ValueInterval(-1.0, 2.0) match
      case EigenSelection.ValueInterval(lo, hi) =>
        assertEqualsDouble(lo, -1.0, 0.0)
        assertEqualsDouble(hi, 2.0, 0.0)
      case other => fail(s"expected ValueInterval, got $other")

    assertEquals(EigenSelection.All, EigenSelection.All)
  }

  test("SingularSelection covers All and Count") {
    assertEquals(
      SingularSelection.Count(5, SingularOrder.Smallest),
      SingularSelection.Count(5, SingularOrder.Smallest)
    )
    assertEquals(SingularSelection.All, SingularSelection.All)
  }

  test("SpectralTarget.ShiftInvert names an explicit real sigma and executable solve plan") {
    val system = Matrix.eye(3).assumePositiveDefinite
    val plan =
      LinearSolvePlan
        .iterative(system, SolverConfig(tolerance = 1e-8))
        .orThrow

    SpectralTarget.ShiftInvert(1.5, plan) match
      case SpectralTarget.ShiftInvert(sigma, LinearSolvePlan.Use(solver)) =>
        assertEqualsDouble(sigma, 1.5, 0.0)
        assertEquals(solver.kind, LinearSolveKind.ConjugateGradient)
      case other => fail(s"expected ShiftInvert / Use, got $other")

    assertEquals(LinearSolvePlan.Backend, LinearSolvePlan.Backend)

    SpectralTarget.Around(0.25) match
      case SpectralTarget.Around(v) => assertEqualsDouble(v, 0.25, 0.0)
      case other                    => fail(s"expected Around, got $other")
  }

  test("SpectralOptions defaults match the documented contract") {
    val o = SpectralOptions()
    assertEqualsDouble(o.tolerance, 1e-10, 0.0)
    assertEquals(o.maxIterations, 1000)
    assertEquals(o.subspaceDimension, None)
    assertEquals(o.startVector, None)
    assertEquals(o.returnVectors, EigenVectors.Right)
  }

  test("GeneralizedSpectralOptions carries a block initial subspace") {
    val defaults = GeneralizedSpectralOptions()
    assertEqualsDouble(defaults.tolerance, 1e-10, 0.0)
    assertEquals(defaults.maxIterations, 200)
    assertEquals(defaults.initialSubspace, None)
    assertEquals(defaults.returnVectors, EigenVectors.Right)

    val initial = Matrix.eye(3)
    val configured = GeneralizedSpectralOptions(
      tolerance = 1e-8,
      maxIterations = 40,
      initialSubspace = Some(initial),
      returnVectors = EigenVectors.ValuesOnly
    )
    assert(configured.initialSubspace.contains(initial))
    assertEquals(configured.returnVectors, EigenVectors.ValuesOnly)
  }

  test("EigenVectors enumerates all four cases") {
    assertEquals(
      EigenVectors.values.toSet,
      Set(EigenVectors.ValuesOnly, EigenVectors.Right, EigenVectors.Left, EigenVectors.LeftAndRight)
    )
  }
