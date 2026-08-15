#!/usr/bin/env python3
"""Generate checked-in NumPy/SciPy reference fixtures for Gale parity tests.

Breeze has no honest public reference for generalized eigen, GSVD, QZ,
sparse-direct factorization, near-cutoff rank/pinv/cond, or Krylov algorithm
diagnostics. This script records NumPy/SciPy answers (the R counterparts are
geigen, svd/gsvd via geigen or pracma, kappa, Matrix::lu, and the same
Krylov family) so `sbt parityTest` does not need Python or R at runtime.

Regenerate from the repository root:

    python3 parity/scripts/generate_numpy_references.py
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh, pinv
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import bicgstab, cg, gmres, lsqr, splu

OUT = Path(__file__).resolve().parents[1] / "src/test/scala/gale/parity/NumpyScipyFixtures.scala"
EPS = float(np.finfo(np.float64).eps)


def fmt(x: float) -> str:
    x = float(x)
    if math.isnan(x):
        return "Double.NaN"
    if math.isinf(x):
        return "Double.PositiveInfinity" if x > 0.0 else "Double.NegativeInfinity"
    return f"{x:.16e}"


def emit_vec(xs: np.ndarray, indent: str = "        ") -> str:
    body = ", ".join(fmt(float(v)) for v in np.asarray(xs, dtype=np.float64).ravel())
    return f"Array({body})"


def emit_mat(a: np.ndarray, indent: str = "        ") -> str:
    rows = [emit_vec(row) for row in np.asarray(a, dtype=np.float64)]
    inner = f",\n{indent}  ".join(rows)
    return f"Array(\n{indent}  {inner}\n{indent})"


def orthonormal(n: int, rng: np.random.Generator) -> np.ndarray:
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    # Fix the sign of each diagonal of R-equivalent so regeneration is stable.
    for j in range(n):
        if q[0, j] < 0.0:
            q[:, j] *= -1.0
    return q


def symmetric(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    s = rng.uniform(-1.0, 1.0, size=(n, n))
    return 0.5 * (s + s.T)


def spd(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    b = rng.uniform(-1.0, 1.0, size=(n, n))
    return b @ b.T + n * np.eye(n)


def diag_dom(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = rng.uniform(-1.0, 1.0, size=(n, n))
    for i in range(n):
        a[i, i] = np.sum(np.abs(a[i])) - np.abs(a[i, i]) + 1.0
    return a


def random_matrix(rows: int, cols: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(-1.0, 1.0, size=(rows, cols))


def with_singular_values(sigmas: list[float], m: int, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    u = orthonormal(m, rng)
    v = orthonormal(n, rng)
    s = np.zeros((m, n))
    for i, sigma in enumerate(sigmas):
        s[i, i] = sigma
    return u @ s @ v.T


def gsvd_ratios(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Generalized singular values via the Gram pencil (B full column rank)."""
    w = eigh(a.T @ a, b.T @ b, eigvals_only=True)
    ratios = np.sqrt(np.maximum(w, 0.0))
    return np.sort(ratios)[::-1]


def count_iters(callback_holder: list[int]):
    def cb(*_args, **_kwargs):
        callback_holder[0] += 1

    return cb


def generalized_eigen_refs() -> list[str]:
    cases = [
        ("n4_well", 4, 11, 12),
        ("n6_well", 6, 21, 22),
        ("n8_well", 8, 31, 32),
        ("n5_moderate_kappa_B", 5, 41, 42),
    ]
    blocks = []
    for name, n, a_seed, b_seed in cases:
        a = symmetric(n, a_seed)
        if name.endswith("kappa_B"):
            # Moderately conditioned metric: eigenvalues of B span ~1e4.
            rng = np.random.default_rng(b_seed)
            q = orthonormal(n, rng)
            spectrum = np.geomspace(1.0, 1.0e4, n)
            b = q @ np.diag(spectrum) @ q.T
            b = 0.5 * (b + b.T)
        else:
            b = spd(n, b_seed)
        vals = eigh(a, b, eigvals_only=True)
        blocks.append(
            "    GeneralizedEigenRef(\n"
            f'      "{name}",\n'
            f"      {emit_mat(a)},\n"
            f"      {emit_mat(b)},\n"
            f"      {emit_vec(vals)}\n"
            "    )"
        )
    return blocks


def gsvd_refs() -> list[str]:
    blocks = []
    cases = [
        ("rand_5x3_4x3", random_matrix(5, 3, 101), random_matrix(4, 3, 102) + 0.75 * np.eye(4, 3)),
        ("rand_6x4_5x4", random_matrix(6, 4, 103), random_matrix(5, 4, 104) + 0.75 * np.eye(5, 4)),
        ("rand_7x5_6x5", random_matrix(7, 5, 105), random_matrix(6, 5, 106) + 0.5 * np.eye(6, 5)),
    ]
    a_bi = random_matrix(6, 4, 107)
    cases.append(("b_identity_6x4", a_bi, np.eye(4)))
    for name, a, b in cases:
        ratios = gsvd_ratios(a, b)
        svd_a = np.linalg.svd(a, compute_uv=False) if name.startswith("b_identity") else None
        extra = ""
        if svd_a is not None:
            extra = f",\n      ordinarySvd = Some({emit_vec(svd_a)})"
        blocks.append(
            "    GsvdRef(\n"
            f'      "{name}",\n'
            f"      {emit_mat(a)},\n"
            f"      {emit_mat(b)},\n"
            f"      {emit_vec(ratios)}{extra}\n"
            "    )"
        )
    # Analytic Infinite / Finite(4/3) / Zero pencil (full column rank).
    a = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 0.8, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    b = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.6, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
        ]
    )
    blocks.append(
        "    GsvdRef(\n"
        '      "analytic_infinite_finite_zero",\n'
        f"      {emit_mat(a)},\n"
        f"      {emit_mat(b)},\n"
        f"      {emit_vec(np.array([np.inf, 4.0 / 3.0, 0.0]))}\n"
        "    )"
    )
    return blocks


def near_cutoff_refs() -> list[str]:
    blocks = []

    def add(name: str, a: np.ndarray, kind: str, definite_rank: int, near_count: int) -> None:
        m, n = a.shape
        rcond = max(m, n) * EPS
        plus = pinv(a, atol=0.0, rtol=rcond)
        svd_rank = int(np.linalg.matrix_rank(a, rtol=rcond))
        if m == n:
            cond1 = float(np.linalg.cond(a, 1))
        else:
            cond1 = float("nan")
        blocks.append(
            "    NearCutoffRef(\n"
            f'      "{name}",\n'
            f'      "{kind}",\n'
            f"      {emit_mat(a)},\n"
            f"      {emit_mat(plus)},\n"
            f"      numpySvdRank = {svd_rank},\n"
            f"      definiteRank = {definite_rank},\n"
            f"      nearCutoffCount = {near_count},\n"
            f"      cond1 = {fmt(cond1)}\n"
            "    )"
        )

    add("diag_clear", np.diag([0.5, 1.0, 4.0, 25.0]), "clear", 4, 0)
    add("diag_exact_zero", np.diag([1.0, 0.0, 2.0, -4.0]), "clear_deficient", 3, 0)
    add(
        "prescribed_clear_full",
        with_singular_values([2.0, 1.1, 0.7, 0.3, 0.15, 0.08], 8, 6, 201),
        "clear",
        6,
        0,
    )
    add(
        "prescribed_exact_zero_col",
        with_singular_values([1.5, 0.9, 0.4, 0.2, 0.1, 0.0], 8, 6, 202),
        "clear_deficient",
        5,
        0,
    )
    m, n = 8, 6
    cutoff = max(m, n) * EPS * 1.0
    add(
        "prescribed_near_pinv_cutoff",
        with_singular_values([1.0, 0.4, 0.2, 1.0e-4, 2.0 * cutoff, 0.5 * cutoff], m, n, 203),
        "near_cutoff",
        4,
        2,
    )
    add(
        "square_near_singular_diag",
        np.diag([1.0, 1.0, 1.0, 1.0e-12]),
        "near_singular",
        4,
        0,
    )
    add("square_singular", np.array([[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [3.0, 6.0, 9.0]]), "singular", 1, 0)
    return blocks


def iterative_refs() -> list[str]:
    blocks = []
    rtol = 1e-10

    def add_square(name: str, algorithm: str, a: np.ndarray, seed: int) -> None:
        rng = np.random.default_rng(seed)
        x_true = rng.uniform(-1.0, 1.0, size=a.shape[0])
        b = a @ x_true
        ncalls = [0]
        cb = count_iters(ncalls)
        if algorithm == "cg":
            x, info = cg(a, b, rtol=rtol, atol=0.0, maxiter=500, callback=cb)
        elif algorithm == "bicgstab":
            x, info = bicgstab(a, b, rtol=rtol, atol=0.0, maxiter=500, callback=cb)
        elif algorithm == "gmres":
            x, info = gmres(
                a,
                b,
                rtol=rtol,
                atol=0.0,
                restart=40,
                maxiter=500,
                callback=cb,
                callback_type="legacy",
            )
        else:
            raise ValueError(algorithm)
        residual = float(np.linalg.norm(a @ x - b))
        blocks.append(
            "    IterativeRef(\n"
            f'      "{name}",\n'
            f'      "{algorithm}",\n'
            f"      {emit_mat(a)},\n"
            f"      {emit_vec(b)},\n"
            f"      {emit_vec(x)},\n"
            f"      iterations = {int(ncalls[0])},\n"
            f"      residual = {fmt(residual)},\n"
            f"      rtol = {fmt(rtol)},\n"
            f"      converged = {str(info == 0).lower()}\n"
            "    )"
        )

    add_square("cg_spd_n8", "cg", spd(8, 301), 311)
    add_square("cg_spd_n16", "cg", spd(16, 302), 312)
    add_square("bicgstab_dd_n10_s1", "bicgstab", diag_dom(10, 321), 331)
    add_square("bicgstab_dd_n10_s2", "bicgstab", diag_dom(10, 322), 332)
    add_square("gmres_dd_n8_s1", "gmres", diag_dom(8, 341), 351)
    add_square("gmres_dd_n12_s2", "gmres", diag_dom(12, 342), 352)

    rng = np.random.default_rng(361)
    a = rng.uniform(-1.0, 1.0, size=(20, 6))
    for i in range(6):
        a[i, i] += 3.0
    x_true = rng.uniform(-1.0, 1.0, size=6)
    b = a @ x_true
    res = lsqr(a, b, atol=rtol, btol=rtol, iter_lim=500)
    x, istop, itn, arnorm = res[0], res[1], res[2], res[7]
    blocks.append(
        "    IterativeRef(\n"
        '      "lsqr_20x6",\n'
        '      "lsqr",\n'
        f"      {emit_mat(a)},\n"
        f"      {emit_vec(b)},\n"
        f"      {emit_vec(x)},\n"
        f"      iterations = {int(itn)},\n"
        f"      residual = {fmt(float(arnorm))},\n"
        f"      rtol = {fmt(rtol)},\n"
        f"      converged = {str(istop in (1, 2)).lower()}\n"
        "    )"
    )
    return blocks


def sparse_lu_refs() -> list[str]:
    blocks = []

    def add(name: str, a: np.ndarray, seed: int) -> None:
        rng = np.random.default_rng(seed)
        x_true = rng.uniform(-1.0, 1.0, size=a.shape[0])
        b = a @ x_true
        x = splu(csc_matrix(a)).solve(b)
        blocks.append(
            "    SparseLuRef(\n"
            f'      "{name}",\n'
            f"      {emit_mat(a)},\n"
            f"      {emit_vec(b)},\n"
            f"      {emit_vec(x)}\n"
            "    )"
        )

    # Sparse-ish diagonally dominant: tridiagonal plus a few fill-ins.
    n = 6
    a = np.zeros((n, n))
    for i in range(n):
        a[i, i] = 4.0
        if i > 0:
            a[i, i - 1] = -1.0
        if i + 1 < n:
            a[i, i + 1] = -1.0
    a[0, n - 1] = 0.3
    a[n - 1, 0] = 0.3
    add("tridiag_plus_corners", a, 401)

    add("spd_dense_stored", spd(5, 402), 412)
    return blocks


def main() -> None:
    parts = []
    parts.append(
        """package gale.parity

/** NumPy / SciPy reference fixtures for operations with no honest Breeze counterpart.
  *
  * Generated by `parity/scripts/generate_numpy_references.py`. Do not edit by
  * hand. Inputs are stored explicitly so the Scala tests reconstruct the same
  * matrices the references were computed from.
  */
object NumpyScipyFixtures:

  final case class GeneralizedEigenRef(
      name: String,
      a: Array[Array[Double]],
      b: Array[Array[Double]],
      eigenvalues: Array[Double]
  )

  final case class GsvdRef(
      name: String,
      a: Array[Array[Double]],
      b: Array[Array[Double]],
      ratios: Array[Double],
      ordinarySvd: Option[Array[Double]] = None
  )

  final case class NearCutoffRef(
      name: String,
      kind: String,
      a: Array[Array[Double]],
      pinv: Array[Array[Double]],
      numpySvdRank: Int,
      definiteRank: Int,
      nearCutoffCount: Int,
      cond1: Double
  )

  final case class IterativeRef(
      name: String,
      algorithm: String,
      a: Array[Array[Double]],
      b: Array[Double],
      x: Array[Double],
      iterations: Int,
      residual: Double,
      rtol: Double,
      converged: Boolean
  )

  final case class SparseLuRef(
      name: String,
      a: Array[Array[Double]],
      b: Array[Double],
      x: Array[Double]
  )
"""
    )
    parts.append("  val generalizedEigen: IndexedSeq[GeneralizedEigenRef] =\n    IndexedSeq(\n")
    parts.append(",\n".join(generalized_eigen_refs()))
    parts.append("\n    )\n\n")
    parts.append("  val gsvd: IndexedSeq[GsvdRef] =\n    IndexedSeq(\n")
    parts.append(",\n".join(gsvd_refs()))
    parts.append("\n    )\n\n")
    parts.append("  val nearCutoff: IndexedSeq[NearCutoffRef] =\n    IndexedSeq(\n")
    parts.append(",\n".join(near_cutoff_refs()))
    parts.append("\n    )\n\n")
    parts.append("  val iterative: IndexedSeq[IterativeRef] =\n    IndexedSeq(\n")
    parts.append(",\n".join(iterative_refs()))
    parts.append("\n    )\n\n")
    parts.append("  val sparseLu: IndexedSeq[SparseLuRef] =\n    IndexedSeq(\n")
    parts.append(",\n".join(sparse_lu_refs()))
    parts.append("\n    )\n")

    OUT.write_text("".join(parts))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
