#!/usr/bin/env python3
"""Materialize shared fixtures, time prepared SciPy solves, and independently check endpoints.

Use the pinned environment from tools/optim-python/requirements.txt. Thread limits must
be set before Python starts. NumPy callbacks are timed; scalar Python verification is not.
"""
import argparse
import csv
import json
import math
import platform
import statistics
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/verification/optim-extensions"


def read_flat(path):
    cases = []
    for line in path.read_text().splitlines():
        cells = line.split("\t")
        if cells[0] == "CASE":
            case = {"id": cells[1], "data": {}}
        elif cells[0] == "END":
            cases.append(case)
        elif cells[0] != "META" and len(cells) > 1:
            case["data"][cells[0]] = [float(x) for x in cells[1:]]
    return cases


def fixtures():
    import numpy as np
    old = {c["id"]: c["data"] for c in read_flat(ROOT / "docs/verification/optim-python/fixtures-v1.decimal.tsv")}
    cases = []
    for name, m, offset, noise in [("lm-exp-64", 64, False, 0.0), ("lm-exp-offset-256", 256, True, 0.01)]:
        t = np.linspace(0.0, 8.0, m)
        y = 2.0 * np.exp(-0.35 * t) + (0.2 if offset else 0.0) + noise * np.sin(np.arange(m) * 1.7)
        cases.append({"id": name, "data": {"X0": [1.0, 0.1, 0.0] if offset else [1.0, 0.1], "T": t.tolist(), "Y": y.tolist()}})
    cases.append({"id": "lm-rosenbrock", "data": {"X0": [-1.2, 1.0]}})
    diagonal = np.geomspace(1e-4, 1e4, 8)
    a = np.vstack([np.diag(diagonal), 1.3 * np.diag(diagonal)])
    cases.append({"id": "lm-scaled-linear-8", "data": {"X0": [0.0] * 8, "X": a.ravel().tolist(), "Y": (a @ (1.0 / diagonal)).tolist()}})
    t = np.linspace(0.1, 5.0, 64)
    a = np.column_stack([np.sin(t), 2 * np.sin(t), np.cos(t)])
    y = a @ np.array([1.0, 1.0, 0.4]) + 0.01 * np.cos(5 * t)
    cases.append({"id": "lm-rank-deficient-3", "data": {"X0": [0.0] * 3, "X": a.ravel().tolist(), "Y": y.tolist()}})
    for n in (16, 512):
        data = dict(old[f"diag-quadratic-n{n}-k1e3"])
        data.update(LOWER=[-0.5] * n, UPPER=[0.5] * n)
        cases.append({"id": f"box-diagonal-{n}", "data": data})
    data = dict(old["rotated-quadratic-n32-k1e3"])
    data.update(LOWER=[-0.5] * 32, UPPER=[0.5] * 32)
    cases.append({"id": "box-rotated-32", "data": data})
    cases.append({"id": "box-rosenbrock-2", "data": {"X0": [-1.2, 1.0], "LOWER": [-2.0, -1.0], "UPPER": [0.8, 2.0]}})
    data = dict(old["logistic-l2-m1024-n16"])
    data.update(LOWER=[-0.15] * 16, UPPER=[0.15] * 16)
    cases.append({"id": "box-logistic-1024-16", "data": data})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "fixtures.json").write_text(json.dumps(cases, indent=2) + "\n")
    lines = []
    for case in cases:
        lines.append("CASE\t" + case["id"])
        for key, values in case["data"].items():
            lines.append(key + "\t" + "\t".join(format(x, ".17g") for x in values))
        lines.append("END")
    (OUT / "fixtures.tsv").write_text("\n".join(lines) + "\n")
    print(f"Materialized {len(cases)} fixtures")


def prepared(case):
    import numpy as np
    from scipy.special import expit
    name, d = case["id"], {k: np.asarray(v, dtype=np.float64) for k, v in case["data"].items()}
    n = len(d["X0"])
    if name.startswith("lm-exp"):
        t, y = d["T"], d["Y"]
        def residual(x):
            return x[0] * np.exp(-x[1] * t) + (x[2] if n == 3 else 0.0) - y
        def jacobian(x):
            exp = np.exp(-x[1] * t)
            columns = [exp, -x[0] * t * exp]
            if n == 3:
                columns.append(np.ones_like(t))
            return np.column_stack(columns)
        return d, residual, jacobian
    if name == "lm-rosenbrock":
        return d, lambda x: np.array([10 * (x[1] - x[0] ** 2), 1 - x[0]]), lambda x: np.array([[-20 * x[0], 10.0], [-1.0, 0.0]])
    if name.startswith("lm-"):
        a, y = d["X"].reshape(-1, n), d["Y"]
        return d, lambda x: a @ x - y, lambda x: a
    if "DIAGONAL" in d:
        def objective(x):
            delta = x - d["CENTER"]
            g = d["DIAGONAL"] * delta
            return 0.5 * np.dot(delta, g), g
    elif "HESSIAN" in d:
        h = d["HESSIAN"].reshape(n, n)
        def objective(x):
            delta = x - d["CENTER"]
            g = h @ delta
            return 0.5 * np.dot(delta, g), g
    elif name.startswith("box-rosenbrock"):
        def objective(x):
            r = x[1] - x[0] * x[0]
            return 100 * r * r + (1 - x[0]) ** 2, np.array([-400 * x[0] * r - 2 * (1 - x[0]), 200 * r])
    else:
        a, y, penalty = d["X"].reshape(-1, n), d["Y"], d["LAMBDA"][0]
        def objective(x):
            margin = y * (a @ x)
            value = np.logaddexp(0, -margin).mean() + 0.5 * penalty * np.dot(x, x)
            gradient = a.T @ (-y * expit(-margin)) / len(y) + penalty * x
            return value, gradient
    return d, objective, None


def independent(case, x):
    """Scalar-loop formulas: no NumPy callbacks or solver-reported residuals."""
    d, name, n = case["data"], case["id"], len(x)
    if len(x) != len(d["X0"]) or not all(math.isfinite(z) for z in x):
        raise ValueError("invalid endpoint")
    if name.startswith("lm-"):
        if name.startswith("lm-exp"):
            residual, jacobian = [], []
            for t, y in zip(d["T"], d["Y"], strict=True):
                e = math.exp(-x[1] * t)
                residual.append(x[0] * e + (x[2] if n == 3 else 0.0) - y)
                jacobian.append([e, -x[0] * t * e] + ([1.0] if n == 3 else []))
        elif name == "lm-rosenbrock":
            residual = [10 * (x[1] - x[0] ** 2), 1 - x[0]]
            jacobian = [[-20 * x[0], 10.0], [-1.0, 0.0]]
        else:
            jacobian = [d["X"][i * n:(i + 1) * n] for i in range(len(d["Y"]))]
            residual = [math.fsum(a * b for a, b in zip(row, x, strict=True)) - y for row, y in zip(jacobian, d["Y"], strict=True)]
        value = 0.5 * math.fsum(r * r for r in residual)
        gradient = [math.fsum(row[j] * r for row, r in zip(jacobian, residual, strict=True)) for j in range(n)]
        return value, max(map(abs, gradient)), 0.0
    if "DIAGONAL" in d:
        delta = [x[i] - d["CENTER"][i] for i in range(n)]
        gradient = [d["DIAGONAL"][i] * delta[i] for i in range(n)]
        value = 0.5 * math.fsum(delta[i] * gradient[i] for i in range(n))
    elif "HESSIAN" in d:
        delta = [x[i] - d["CENTER"][i] for i in range(n)]
        gradient = [math.fsum(d["HESSIAN"][i * n + j] * delta[j] for j in range(n)) for i in range(n)]
        value = 0.5 * math.fsum(delta[i] * gradient[i] for i in range(n))
    elif name.startswith("box-rosenbrock"):
        r = x[1] - x[0] ** 2
        value = 100 * r * r + (1 - x[0]) ** 2
        gradient = [-400 * x[0] * r - 2 * (1 - x[0]), 200 * r]
    else:
        m, a, y, penalty = len(d["Y"]), d["X"], d["Y"], d["LAMBDA"][0]
        losses, residual = [], []
        for i in range(m):
            z = y[i] * math.fsum(a[i * n + j] * x[j] for j in range(n))
            tail = math.exp(-abs(z))
            losses.append(max(0.0, -z) + math.log1p(tail))
            residual.append(-y[i] * (tail if z >= 0 else 1.0) / (1 + tail))
        value = math.fsum(losses) / m + 0.5 * penalty * math.fsum(z * z for z in x)
        gradient = [math.fsum(a[i * n + j] * residual[i] for i in range(m)) / m + penalty * x[j] for j in range(n)]
    feasibility = max(max(l - z, z - u, 0.0) for z, l, u in zip(x, d["LOWER"], d["UPPER"], strict=True))
    pg = [min(g, z - l) if g > 0 else max(g, z - u) for z, g, l, u in zip(x, gradient, d["LOWER"], d["UPPER"], strict=True)]
    return value, max(map(abs, pg)), feasibility


def benchmark(output):
    import numpy as np
    import scipy
    from scipy.optimize import least_squares, minimize
    from threadpoolctl import threadpool_info
    cases = json.loads((OUT / "fixtures.json").read_text())
    rows = []
    for case in cases:
        d, fun, jac = prepared(case)
        bounds = list(zip(d["LOWER"], d["UPPER"])) if "LOWER" in d else None
        def solve(strict=False):
            if case["id"].startswith("lm-"):
                return least_squares(fun, d["X0"], jac=jac, method="lm", max_nfev=20000,
                                     ftol=1e-14 if strict else 1e-12, xtol=1e-14 if strict else 1e-12,
                                     gtol=1e-13 if strict else 1e-8)
            return minimize(fun, d["X0"], jac=True, method="L-BFGS-B", bounds=bounds,
                            options={"maxiter": 50000, "maxfun": 250000, "maxls": 40, "maxcor": 10,
                                     "gtol": 1e-11 if strict else 1e-8, "ftol": 0.0})
        reference = solve(True)
        refvalue, refmetric, reffeas = independent(case, reference.x.tolist())
        reference_source = "strict-scipy"
        reference_x = reference.x.copy()
        if "DIAGONAL" in d:
            reference_x = np.clip(d["CENTER"], d["LOWER"], d["UPPER"])
            reference_source = "analytic-clipped-center"
        elif "HESSIAN" in d:
            # Independently solve the linear KKT equations for the candidate active set.
            h = d["HESSIAN"].reshape(len(reference_x), -1)
            active = (reference_x == d["LOWER"]) | (reference_x == d["UPPER"])
            free = ~active
            delta = reference_x - d["CENTER"]
            delta[free] = np.linalg.solve(h[np.ix_(free, free)], -h[np.ix_(free, active)] @ delta[active])
            reference_x[free] = d["CENTER"][free] + delta[free]
            reference_source = "active-set-linear-kkt"
        refvalue, refmetric, reffeas = independent(case, reference_x.tolist())
        if refmetric > 1e-6 or reffeas != 0.0:
            raise ValueError(f"reference fails independent checks: {case['id']} {refmetric}")
        for _ in range(10):
            solve()
        times = []
        for _ in range(11):
            before = time.perf_counter_ns()
            result = solve()
            times.append((time.perf_counter_ns() - before) / 1e6)
        value, metric, feasibility = independent(case, result.x.tolist())
        valid = metric <= 1e-6 and feasibility == 0 and abs(value - refvalue) <= 1e-7 * (1 + abs(refvalue))
        rows.append({"id": case["id"], "x": result.x.tolist(), "objective": value, "metric": metric,
                     "valid": valid, "success": bool(result.success), "message": str(result.message),
                     "nfev": int(result.nfev), "njev": int(result.njev), "median_ms": statistics.median(times),
                     "times_ms": times, "reference_objective": refvalue, "reference_metric": refmetric,
                     "reference_success": bool(reference.success), "reference_source": reference_source,
                     "strict_scipy_metric": independent(case, reference.x.tolist())[1], "reference_x": reference_x.tolist()})
        print(case["id"], "valid=", valid, "ms=", round(statistics.median(times), 4), "metric=", metric)
    output.write_text(json.dumps({"cases": rows, "environment": {"python": platform.python_version(),
        "platform": platform.platform(), "numpy": np.__version__, "scipy": scipy.__version__, "threadpools": threadpool_info()}}, indent=2) + "\n")
    print("Timed Python qualification failures:", [row["id"] for row in rows if not row["valid"]])


def verify(harness, baseline, output):
    cases = {c["id"]: c for c in json.loads((OUT / "fixtures.json").read_text())}
    python = {c["id"]: c for c in json.loads(baseline.read_text())["cases"]}
    rows = list(csv.DictReader(harness.open(), delimiter="\t"))
    if len(rows) != len(cases) or {r["case_id"] for r in rows} != set(cases):
        raise ValueError("missing, duplicate or unknown harness cases")
    report = []
    for row in rows:
        key = row["case_id"]
        x = list(map(float, row["x"].split(",")))
        value, metric, feasible = independent(cases[key], x)
        reported = float(row["objective"])
        reported_residual = float(row["residual"])
        times = list(map(float, row["times_ms"].split(",")))
        median = float(row["median_ms"])
        reference = python[key]["reference_objective"]
        certificate_valid = (row["status"] in {"Converged", "IterationLimit", "NumericalStagnation", "EvaluationLimit", "Cancelled", "LineSearchFailed"} and
                             math.isfinite(reported_residual) and reported_residual >= 0 and
                             abs(reported_residual - metric) <= 1e-10 * max(1.0, metric) and
                             (row["status"] != "Converged" or (reported_residual <= 1e-8 and metric <= 1.01e-8)))
        valid = (all(math.isfinite(v) for v in [value, metric, feasible, reported, median]) and certificate_valid and
                 metric <= 1e-6 and feasible == 0.0 and abs(value - reference) <= 1e-7 * (1 + abs(reference)) and
                 abs(reported - value) <= 1e-10 * (1 + abs(value)) and len(times) == 11 and
                 all(math.isfinite(t) and t > 0 for t in times) and
                 abs(statistics.median(times) - median) <= 1e-9)
        report.append({"id": key, "valid": valid, "certificate_valid": certificate_valid, "python_valid": python[key]["valid"],
                       "objective_error": abs(value - reference), "metric": metric,
                       "status": row["status"], "gale_ms": median, "python_ms": python[key]["median_ms"],
                       "ratio": median / python[key]["median_ms"], "iterations": int(row["iterations"]),
                       "callbacks": int(row["callbacks"]), "python_nfev": python[key]["nfev"]})
    output.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    lines = ["# SciPy comparison", "", "Independent scalar-loop endpoint checks; time ratio is Gale/Python.", "",
             "| Case | Gale valid | Python valid | Status | Metric | Gale ms | Python ms | Ratio |", "|---|:---:|:---:|---|---:|---:|---:|---:|"]
    for r in report:
        ratio = f"{r['ratio']:.3f}" if r['valid'] and r['python_valid'] else "unqualified"
        lines.append(f"| {r['id']} | {r['valid']} | {r['python_valid']} | {r['status']} | {r['metric']:.3e} | {r['gale_ms']:.3f} | {r['python_ms']:.3f} | {ratio} |")
    output.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"{sum(r['valid'] for r in report)}/{len(report)} pass; {sum(r['valid'] and r['python_valid'] and r['ratio'] < 1 for r in report)} qualified lower Gale medians")
    if not all(r["valid"] for r in report):
        raise ValueError("Gale endpoint qualification failed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["fixtures", "python", "verify"])
    parser.add_argument("--harness", type=Path)
    parser.add_argument("--baseline", type=Path, default=OUT / "python.json")
    parser.add_argument("--output", type=Path, default=OUT / "comparison")
    args = parser.parse_args()
    if args.command == "fixtures":
        fixtures()
    elif args.command == "python":
        benchmark(args.baseline)
    else:
        verify(args.harness, args.baseline, args.output)
