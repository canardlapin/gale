package gale.solvers

import gale.linalg.*

/** Tests for the P5 solver rewrites: Givens-based GMRES with left
  * preconditioning and preconditioned BiCGSTAB (K6), plus the relative
  * tolerance mode (K7).
  */
class PreconditionedSolverSuite extends munit.FunSuite:
  /** 2D 5-point Laplacian on a `g x g` grid (SPD, size `g*g`). */
  private def laplacian2D(g: Int): DMat =
    val n = g * g
    Matrix.tabulate(n, n) { (a, b) =>
      val ra = a / g
      val ca = a % g
      val rb = b / g
      val cb = b % g
      if a == b then 4.0
      else if (ra == rb && math.abs(ca - cb) == 1) || (ca == cb && math.abs(ra - rb) == 1) then -1.0
      else 0.0
    }

  /** Nonsymmetric, strongly diagonally dominant matrix whose diagonal spans four
    * orders of magnitude, so Jacobi (diagonal) preconditioning clusters the
    * spectrum of `M^{-1}A` near 1.
    */
  private def illConditioned(n: Int): DMat =
    Matrix.tabulate(n, n) { (i, j) =>
      if i == j then math.pow(10.0, 4.0 * i / (n - 1)) + 5.0
      else if math.abs(i - j) == 1 then 1.0
      else 0.0
    }

  // ---- K6: GMRES preconditioning ------------------------------------------

  test("preconditioned GMRES converges in fewer iterations than unpreconditioned") {
    val n = 40
    val A = illConditioned(n)
    val truth = Vec.tabulate(n)(i => (i + 1).toDouble)
    val b = A * truth
    val config = SolverConfig(tolerance = 1e-8, maxIterations = n, restart = n)

    val plain = gmres(A, b, config)
    val prec = gmres(A, b, config, Preconditioner.Jacobi(A))

    assert(plain.converged, s"unpreconditioned residual=${plain.residual}")
    assert(prec.converged, s"preconditioned residual=${prec.residual}")
    assert(
      prec.iterations < plain.iterations,
      s"preconditioned=${prec.iterations} unpreconditioned=${plain.iterations}"
    )
    assert(norm(prec.x - truth) < 1e-6)
  }

  test("Givens GMRES still solves a random nonsymmetric system") {
    val rng = new scala.util.Random(5L)
    val n = 12
    val A = Matrix.tabulate(n, n) { (i, j) =>
      if i == j then n.toDouble + rng.nextDouble()
      else rng.nextDouble() * 2.0 - 1.0
    }
    val truth = Vec.tabulate(n)(i => rng.nextDouble() * 4.0 - 2.0)
    val b = A * truth
    val result = gmres(A, b, SolverConfig(tolerance = 1e-11, maxIterations = 200, restart = n))
    assert(result.converged, s"residual=${result.residual}")
    assert(norm(result.x - truth) < 1e-7)
  }

  // ---- K6: BiCGSTAB preconditioning ---------------------------------------

  test("preconditioned BiCGSTAB solves an SPD system") {
    val A = laplacian2D(6) // 36 x 36
    val truth = Vec.tabulate(36)(i => (i % 5 + 1).toDouble)
    val b = A * truth
    val config = SolverConfig(tolerance = 1e-10, maxIterations = 500)

    val prec = bicgstab(A, b, config, Preconditioner.Jacobi(A))
    assert(prec.converged, s"residual=${prec.residual}")
    assert(norm(prec.x - truth) < 1e-7)

    // Identity preconditioning reproduces the unpreconditioned recurrence.
    val plain = bicgstab(A, b, config)
    val plainIdentity = bicgstab(A, b, config, Preconditioner.Identity)
    assertEquals(plainIdentity.iterations, plain.iterations)
    assert(math.abs(plainIdentity.residual - plain.residual) < 1e-14)
  }

  // ---- K7: tolerance modes -------------------------------------------------

  test("RelativeToRhs stops at a looser absolute residual than Absolute") {
    val A = laplacian2D(8) // 64 x 64
    // A high-frequency random RHS makes CG converge gradually, so the looser
    // relative threshold is crossed several iterations before the absolute one.
    val rng = new scala.util.Random(3L)
    val b = DVec.fromSeq(Seq.fill(64)((rng.nextDouble() * 2.0 - 1.0) * 1000.0))
    val bNorm = b.norm2
    val tol = 1e-6
    val config = SolverConfig(tolerance = tol, maxIterations = 500)

    val absolute = cg(A, b, config, toleranceMode = ToleranceMode.Absolute)
    val relative = cg(A, b, config, toleranceMode = ToleranceMode.RelativeToRhs)

    assert(absolute.converged)
    assert(relative.converged)

    // Each mode respects its own threshold.
    assert(absolute.residual <= tol, s"absolute residual=${absolute.residual}")
    assert(relative.residual <= tol * bNorm, s"relative residual=${relative.residual}")

    // Relative uses the looser bar (tol * ||b|| >> tol), so it stops earlier with
    // a strictly larger residual and does strictly less work.
    assert(relative.residual > tol, s"relative residual=${relative.residual}")
    assert(
      relative.iterations < absolute.iterations,
      s"relative=${relative.iterations} absolute=${absolute.iterations}"
    )
  }

  test("Absolute mode is the default and matches an explicit Absolute request") {
    val A = laplacian2D(6)
    val b = DVec.fill(36)(3.0)
    val config = SolverConfig(tolerance = 1e-9, maxIterations = 500)

    val default = cg(A, b, config)
    val explicit = cg(A, b, config, toleranceMode = ToleranceMode.Absolute)
    assertEquals(default.iterations, explicit.iterations)
    assert(math.abs(default.residual - explicit.residual) < 1e-15)
  }
