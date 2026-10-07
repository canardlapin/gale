#!/usr/bin/env python3
"""Generate portable fixtures and independently run Python optimization baselines.

The generated fixture is deliberately plain JSON: Scala/JVM and Scala.js can
parse its materialized IEEE-754 decimal inputs without implementing NumPy's RNG.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

# Set these before NumPy/SciPy import when invoked directly.  The benchmark
# recipe also sets them in its shell environment, so timings are single-threaded.
for _thread_env in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_thread_env, "1")

import numpy as np
import scipy
from scipy.optimize import minimize

try:
    from sklearn.linear_model import Lasso
    from threadpoolctl import threadpool_info
except ImportError:  # The protocol remains runnable without the optional comparator.
    Lasso = None
    threadpool_info = None


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURES = ROOT / "docs/verification/optim-python/fixtures-v1.json"
DEFAULT_FIXTURES_FLAT = ROOT / "docs/verification/optim-python/fixtures-v1.flat.tsv"
DEFAULT_FIXTURES_DECIMAL = ROOT / "docs/verification/optim-python/fixtures-v1.decimal.tsv"
DEFAULT_REFERENCES = ROOT / "docs/verification/optim-python/reference-results-v1.json"
DEFAULT_REFERENCES_FLAT = ROOT / "docs/verification/optim-python/reference-results-v1.flat.tsv"
DEFAULT_REFERENCES_DECIMAL = ROOT / "docs/verification/optim-python/reference-results-v1.decimal.tsv"
DEFAULT_BENCHMARKS = ROOT / "docs/verification/optim-python/python-baseline-v1.json"
DEFAULT_TSV = ROOT / "docs/verification/optim-python/python-baseline-v1.tsv"
SCHEMA_VERSION = 1


def vector_center(n: int) -> np.ndarray:
    i = np.arange(1, n + 1, dtype=np.float64)
    return np.sin(0.17 * i) + 0.25 * np.cos(0.11 * i)


def diagonal_case(case_id: str, n: int, bounded: bool = False) -> dict[str, Any]:
    d = np.power(10.0, 3.0 * np.arange(n, dtype=np.float64) / (n - 1))
    result: dict[str, Any] = {
        "id": case_id,
        "family": "diagonal_quadratic",
        "dimensions": {"n": n},
        "dtype": "float64",
        "x0": [0.0] * n,
        "data": {"diagonal": d.tolist(), "center": vector_center(n).tolist()},
        "objective": {"kind": "half_weighted_squared_distance", "formula": "0.5 * sum(d_i * (x_i-c_i)^2)"},
        "accuracy": {"stationarity_inf": 1e-6, "reference_stationarity_inf": 1e-9, "objective_abs_delta": 1e-8, "objective_rel_delta": 1e-8},
        "budget": {"max_iter": 20000, "max_eval": 100000},
    }
    if bounded:
        result["bounds"] = {"lower": [-0.5] * n, "upper": [0.5] * n}
        result["objective"]["constraint"] = "-0.5 <= x_i <= 0.5"
    return result


def rotated_case() -> dict[str, Any]:
    n = 32
    rng = np.random.default_rng(2026100701)
    q, r = np.linalg.qr(rng.standard_normal((n, n)))
    q *= np.sign(np.diag(r))  # Canonicalize QR signs before serializing material data.
    d = np.power(10.0, 3.0 * np.arange(n, dtype=np.float64) / (n - 1))
    h = (q * d) @ q.T
    return {
        "id": "rotated-quadratic-n32-k1e3",
        "family": "rotated_quadratic",
        "dimensions": {"n": n}, "dtype": "float64", "x0": [0.0] * n,
        "data": {"hessian": h.tolist(), "center": vector_center(n).tolist()},
        "objective": {"kind": "half_quadratic", "formula": "0.5 * (x-c)^T H (x-c)"},
        "accuracy": {"stationarity_inf": 1e-6, "reference_stationarity_inf": 1e-9, "objective_abs_delta": 1e-8, "objective_rel_delta": 1e-8},
        "budget": {"max_iter": 30000, "max_eval": 150000},
    }


def rosenbrock_case(n: int) -> dict[str, Any]:
    x0 = [-1.2 if i % 2 == 0 else 1.0 for i in range(n)]
    return {
        "id": f"rosenbrock-n{n}", "family": "rosenbrock", "dimensions": {"n": n},
        "dtype": "float64", "x0": x0, "data": {},
        "objective": {"kind": "rosenbrock", "formula": "sum_i 100*(x_(i+1)-x_i^2)^2 + (1-x_i)^2"},
        "accuracy": {"stationarity_inf": 1e-6, "reference_stationarity_inf": 1e-9, "objective_abs_delta": 1e-8, "objective_rel_delta": 1e-8},
        "budget": {"max_iter": 50000, "max_eval": 250000},
    }


def logistic_case(m: int, n: int) -> dict[str, Any]:
    rng = np.random.default_rng(2026101000 + 100 * m + n)
    x = rng.standard_normal((m, n)) / np.sqrt(n)
    beta = vector_center(n)
    noise = 0.20 * rng.standard_normal(m)
    y = np.where(x @ beta + noise >= 0.0, 1.0, -1.0)
    return {
        "id": f"logistic-l2-m{m}-n{n}", "family": "l2_logistic", "dimensions": {"m": m, "n": n},
        "dtype": "float64", "x0": [0.0] * n,
        "data": {"x": x.tolist(), "y": y.tolist(), "lambda": 0.1},
        "objective": {"kind": "mean_logistic_plus_l2", "formula": "mean(logaddexp(0,-y*(Xx))) + lambda/2*||x||_2^2"},
        "accuracy": {"stationarity_inf": 1e-6, "reference_stationarity_inf": 1e-9, "objective_abs_delta": 1e-8, "objective_rel_delta": 1e-8},
        "budget": {"max_iter": 30000, "max_eval": 150000},
    }


def lasso_case(m: int, n: int) -> dict[str, Any]:
    rng = np.random.default_rng(2026102000 + 100 * m + n)
    x = rng.standard_normal((m, n))
    beta = vector_center(n) * (np.arange(n) % 3 == 0)
    y = x @ beta + 0.02 * rng.standard_normal(m)
    return {
        "id": f"lasso-m{m}-n{n}", "family": "lasso", "dimensions": {"m": m, "n": n},
        "dtype": "float64", "x0": [0.0] * n,
        "data": {"x": x.tolist(), "y": y.tolist(), "lambda": 0.05},
        "objective": {"kind": "least_squares_plus_l1", "formula": "0.5/m*||Xx-y||_2^2 + lambda*||x||_1"},
        "accuracy": {"kkt_inf": 1e-6, "reference_kkt_inf": 1e-9, "objective_abs_delta": 1e-8, "objective_rel_delta": 1e-8},
        "budget": {"max_iter": 100000, "max_eval": 100000},
    }


def fixtures() -> dict[str, Any]:
    return {"schema_version": SCHEMA_VERSION, "generator": {"numpy_rng": "PCG64", "materialized": True}, "cases": [
        diagonal_case("diag-quadratic-n16-k1e3", 16), diagonal_case("diag-quadratic-n512-k1e3", 512),
        diagonal_case("box-diag-quadratic-n16-k1e3", 16, bounded=True), rotated_case(),
        rosenbrock_case(2), rosenbrock_case(32), logistic_case(1024, 16), lasso_case(256, 32), lasso_case(1024, 128),
    ]}


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def numeric_line(name: str, value: Any, decimal: bool = False) -> str:
    """Flatten numbers as C99/Java hex floats, preserving every binary64 bit."""
    values = np.asarray(value, dtype=np.float64).reshape(-1)
    render = (lambda value: format(float(value), ".17g")) if decimal else (lambda value: float(value).hex())
    return name + "\t" + "\t".join(render(v) for v in values)


def write_flat_fixtures(path: Path, doc: dict[str, Any], decimal: bool = False) -> None:
    lines = [f"GALE_OPTIM_FIXTURES_V{SCHEMA_VERSION}"]
    for case in doc["cases"]:
        metadata = {k: v for k, v in case.items() if k not in ("x0", "data", "bounds")}
        scalar = lambda value: format(float(value), ".17g") if decimal else float(value).hex()
        lines.extend((f"CASE\t{case['id']}", "META\t" + json.dumps(metadata, separators=(",", ":"), allow_nan=False), numeric_line("X0", case["x0"], decimal)))
        for name, value in case["data"].items():
            lines.append(numeric_line(name.upper(), value, decimal) if isinstance(value, list) else name.upper() + "\t" + scalar(value))
        if "bounds" in case:
            lines.extend((numeric_line("LOWER", case["bounds"]["lower"], decimal), numeric_line("UPPER", case["bounds"]["upper"], decimal)))
        lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_flat_references(path: Path, fixtures_doc: dict[str, Any], references_doc: dict[str, Any], decimal: bool = False) -> None:
    cases = {case["id"]: case for case in fixtures_doc["cases"]}
    lines = [f"GALE_OPTIM_REFERENCES_V{SCHEMA_VERSION}"]
    for result in references_doc["results"]:
        expected_result = result["expected"]
        lines.extend((f"CASE\t{result['case_id']}", "META\t" + json.dumps({k: v for k, v in result.items() if k != "expected"}, separators=(",", ":"), allow_nan=False)))
        if expected_result.get("x") is not None:
            solution = np.asarray(expected_result["x"], dtype=np.float64)
            case = cases[result["case_id"]]
            lines.append(numeric_line("EXPECTED_X", solution, decimal))
            if case["family"] == "lasso":
                _, kkt = lasso_metrics(case, solution)
                lines.append("EXPECTED_KKT_INF\t" + (format(float(kkt), ".17g") if decimal else float(kkt).hex()))
            else:
                _, gradient = objective(case)[0](solution)
                stationarity = smooth_stationarity(case, solution, gradient)
                lines.extend((numeric_line("EXPECTED_GRADIENT", gradient, decimal), "EXPECTED_STATIONARITY_INF\t" + (format(stationarity, ".17g") if decimal else stationarity.hex())))
        lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def objective(case: dict[str, Any]) -> tuple[Callable[[np.ndarray], tuple[float, np.ndarray]], Callable[[np.ndarray], float]]:
    family, data = case["family"], case["data"]
    if family == "diagonal_quadratic":
        d, c = np.asarray(data["diagonal"]), np.asarray(data["center"])
        return lambda z: (float(0.5 * np.dot(d * (z-c), z-c)), d * (z-c)), lambda z: float(0.5 * np.dot(d * (z-c), z-c))
    if family == "rotated_quadratic":
        h, c = np.asarray(data["hessian"]), np.asarray(data["center"])
        return lambda z: (float(0.5 * np.dot(z-c, h @ (z-c))), h @ (z-c)), lambda z: float(0.5 * np.dot(z-c, h @ (z-c)))
    if family == "rosenbrock":
        def fg(z: np.ndarray) -> tuple[float, np.ndarray]:
            r = z[1:] - z[:-1] ** 2
            g = np.zeros_like(z)
            g[:-1] = -400.0 * z[:-1] * r - 2.0 * (1.0-z[:-1])
            g[1:] += 200.0 * r
            return float(np.sum(100.0*r*r + (1.0-z[:-1])**2)), g
        return fg, lambda z: fg(z)[0]
    if family == "l2_logistic":
        x, y, lam = np.asarray(data["x"]), np.asarray(data["y"]), float(data["lambda"])
        def fg(z: np.ndarray) -> tuple[float, np.ndarray]:
            yz = y * (x @ z)
            p = 1.0 / (1.0 + np.exp(np.clip(yz, -745.0, 745.0)))
            return float(np.mean(np.logaddexp(0.0, -yz)) + 0.5*lam*np.dot(z,z)), -(x.T @ (y*p))/x.shape[0] + lam*z
        return fg, lambda z: fg(z)[0]
    raise ValueError(f"smooth objective unavailable for {family}")


def lasso_metrics(case: dict[str, Any], z: np.ndarray) -> tuple[float, float]:
    x, y, lam = np.asarray(case["data"]["x"]), np.asarray(case["data"]["y"]), float(case["data"]["lambda"])
    grad = x.T @ (x @ z - y) / x.shape[0]
    nonzero = np.abs(z) > 1e-10
    kkt = np.empty_like(z)
    kkt[nonzero] = np.abs(grad[nonzero] + lam * np.sign(z[nonzero]))
    kkt[~nonzero] = np.maximum(np.abs(grad[~nonzero]) - lam, 0.0)
    return float(0.5*np.dot(x@z-y, x@z-y)/x.shape[0] + lam*np.sum(np.abs(z))), float(np.max(kkt))


def smooth_stationarity(case: dict[str, Any], x: np.ndarray, gradient: np.ndarray) -> float:
    if "bounds" not in case:
        return float(np.max(np.abs(gradient)))
    lo, hi = np.asarray(case["bounds"]["lower"]), np.asarray(case["bounds"]["upper"])
    projected = np.where((x <= lo + 1e-10) & (gradient > 0), 0.0, np.where((x >= hi - 1e-10) & (gradient < 0), 0.0, gradient))
    return float(np.max(np.abs(projected)))


def prepare(case: dict[str, Any]) -> dict[str, Any]:
    prepared = dict(case)
    prepared["x0"] = np.asarray(case["x0"], dtype=np.float64)
    prepared["data"] = {key: np.asarray(value, dtype=np.float64) if isinstance(value, list) else value for key, value in case["data"].items()}
    return prepared


def solve(case: dict[str, Any], reference: bool = True, cached_objective=None) -> dict[str, Any]:
    start = np.asarray(case["x0"], dtype=np.float64)
    if case["family"] == "lasso":
        if Lasso is None:
            return {"solver": "scikit-learn Lasso", "available": False, "reason": "scikit-learn is not installed"}
        d = case["data"]
        model = Lasso(alpha=float(d["lambda"]), fit_intercept=False, max_iter=case["budget"]["max_iter"], tol=1e-12 if reference else 1e-8, selection="cyclic")
        model.fit(np.asarray(d["x"]), np.asarray(d["y"]))
        z = model.coef_.astype(np.float64, copy=False)
        value, kkt = lasso_metrics(case, z)
        return {"solver": "scikit-learn Lasso coordinate descent", "available": True, "x": z.tolist(), "objective": value, "kkt_inf": kkt, "iterations": int(model.n_iter_), "dual_gap": float(model.dual_gap_), "success": kkt <= case["accuracy"]["reference_kkt_inf"]}
    fg, fun = cached_objective if cached_objective is not None else objective(case)
    bounds = None
    if "bounds" in case:
        bounds = list(zip(case["bounds"]["lower"], case["bounds"]["upper"], strict=True))
    res = minimize(fg, start, jac=True, method="L-BFGS-B", bounds=bounds, options={"gtol": 1e-12 if reference else 1e-7, "ftol": 1e-15 if reference else 0.0, "maxiter": case["budget"]["max_iter"], "maxfun": case["budget"]["max_eval"], "maxls": 50})
    _, grad = fg(res.x)
    stationarity = smooth_stationarity(case, res.x, grad)
    return {"solver": "SciPy minimize L-BFGS-B" + (" bounded" if bounds is not None else " unbounded"), "available": True, "x": res.x.tolist(), "objective": fun(res.x), "stationarity_inf": stationarity, "iterations": int(res.nit), "nfev": int(res.nfev), "njev": int(res.njev), "status": int(res.status), "message": str(res.message), "success": stationarity <= case["accuracy"]["reference_stationarity_inf"]}


def environment() -> dict[str, Any]:
    return {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "scikit_learn": None if Lasso is None else __import__("sklearn").__version__, "threads": {key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS")}, "threadpools": [] if threadpool_info is None else threadpool_info()}


def expected(case: dict[str, Any]) -> dict[str, Any]:
    """Use an analytic minimizer where one exists; otherwise retain a strict external solve."""
    family = case["family"]
    if family == "diagonal_quadratic":
        x = np.asarray(case["data"]["center"], dtype=np.float64)
        if "bounds" in case:
            x = np.clip(x, np.asarray(case["bounds"]["lower"]), np.asarray(case["bounds"]["upper"]))
    elif family == "rotated_quadratic":
        x = np.asarray(case["data"]["center"], dtype=np.float64)
    elif family == "rosenbrock":
        x = np.ones(case["dimensions"]["n"], dtype=np.float64)
    else:
        result = solve(case)
        return {"source": result["solver"], "x": result.get("x"), "objective": result.get("objective"), "metric": result.get("kkt_inf", result.get("stationarity_inf")), "reference_success": result.get("success", False)}
    if family == "lasso":
        value, metric = lasso_metrics(case, x)
    else:
        fg, _ = objective(case)
        value, grad = fg(x)
        metric = smooth_stationarity(case, x, grad)
    return {"source": "analytic minimizer", "x": x.tolist(), "objective": value, "metric": metric, "reference_success": metric <= case["accuracy"].get("reference_kkt_inf", case["accuracy"].get("reference_stationarity_inf"))}


def references(fixtures_doc: dict[str, Any]) -> dict[str, Any]:
    results = []
    for case in fixtures_doc["cases"]:
        results.append({"case_id": case["id"], "standard_comparator": solve(case), "expected": expected(case)})
    return {"schema_version": SCHEMA_VERSION, "purpose": "Standard comparator results plus stronger analytic/strict expected endpoints; comparator failures are retained.", "environment": environment(), "results": results}


def benchmark(fixtures_doc: dict[str, Any], warmup: int, repeats: int) -> dict[str, Any]:
    rows = []
    for material in fixtures_doc["cases"]:
        case = prepare(material)
        cached = None if case["family"] == "lasso" else objective(case)
        for _ in range(warmup):
            solve(case, reference=False, cached_objective=cached)
        times, final = [], None
        for _ in range(repeats):
            began = time.perf_counter_ns()
            final = solve(case, reference=False, cached_objective=cached)
            times.append((time.perf_counter_ns() - began) / 1e6)
        assert final is not None
        metric = final.get("kkt_inf", final.get("stationarity_inf"))
        target = case["accuracy"].get("kkt_inf", case["accuracy"].get("stationarity_inf"))
        rows.append({"case_id": case["id"], "solver": final.get("solver"), "available": final.get("available"), "median_ms": float(np.median(times)), "min_ms": float(np.min(times)), "max_ms": float(np.max(times)), "repeats": repeats, "warmup": warmup, "times_ms": times, "metric": metric, "target": target, "timed_endpoint_valid": isinstance(metric, float) and metric <= target, "result": final})
    return {"schema_version": SCHEMA_VERSION, "purpose": "Solve-only wall time; fixture parse/construction and process startup are excluded and reported separately.", "environment": environment(), "warmup": warmup, "repeats": repeats, "rows": rows}


def write_tsv(path: Path, report: dict[str, Any]) -> None:
    lines = ["case_id\tsolver\tavailable\tmedian_ms\tmetric\ttarget\ttimed_endpoint_valid"]
    for row in report["rows"]:
        lines.append("\t".join(str(row.get(k, "")) for k in ("case_id", "solver", "available", "median_ms", "metric", "target", "timed_endpoint_valid")))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("generate", "reference", "benchmark", "all"))
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--fixtures-flat", type=Path, default=DEFAULT_FIXTURES_FLAT)
    parser.add_argument("--fixtures-decimal", type=Path, default=DEFAULT_FIXTURES_DECIMAL)
    parser.add_argument("--references", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--references-flat", type=Path, default=DEFAULT_REFERENCES_FLAT)
    parser.add_argument("--references-decimal", type=Path, default=DEFAULT_REFERENCES_DECIMAL)
    parser.add_argument("--benchmarks", type=Path, default=DEFAULT_BENCHMARKS)
    parser.add_argument("--tsv", type=Path, default=DEFAULT_TSV)
    parser.add_argument("--warmup", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.command in ("generate", "all"):
        generated = fixtures()
        write_json(args.fixtures, generated)
        write_flat_fixtures(args.fixtures_flat, generated)
        write_flat_fixtures(args.fixtures_decimal, generated, decimal=True)
    data = json.loads(args.fixtures.read_text(encoding="utf-8"))
    if args.command in ("reference", "all"):
        refs = references(data)
        write_json(args.references, refs)
        write_flat_references(args.references_flat, data, refs)
        write_flat_references(args.references_decimal, data, refs, decimal=True)
    if args.command in ("benchmark", "all"):
        report = benchmark(data, args.warmup, args.repeats)
        write_json(args.benchmarks, report)
        write_tsv(args.tsv, report)


if __name__ == "__main__":
    main()
