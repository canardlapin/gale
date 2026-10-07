#!/usr/bin/env python3
"""Independently verify Gale smooth-optimization harness output.

The verifier intentionally evaluates fixture formulae with scalar Python loops
instead of calling ``optim_protocol.objective``: a shared objective bug must not
turn a Gale comparison into self-confirmation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURES = ROOT / "docs/verification/optim-python/fixtures-v1.json"
DEFAULT_REFERENCES = ROOT / "docs/verification/optim-python/reference-results-v1.json"
DEFAULT_PYTHON_BASELINE = ROOT / "docs/verification/optim-python/python-baseline-v1.json"
DEFAULT_JSON = ROOT / "docs/verification/optim-python/gale-comparison-v1.json"
DEFAULT_MARKDOWN = ROOT / "docs/verification/optim-python/gale-comparison-v1.md"
ACCURACY = 1e-6
OBJECTIVE_FACTOR = 1e-7
ZERO_COORDINATE = 1e-10


def finite(values: list[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def json_number(value: float) -> float | None:
    return value if math.isfinite(value) else None


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right, strict=True))


def logistic_loss(value: float) -> float:
    # log(1 + exp(-value)), without overflowing at a bad trial endpoint.
    return math.log1p(math.exp(-value)) if value >= 0.0 else -value + math.log1p(math.exp(value))


def evaluate(case: dict[str, Any], point: list[float]) -> tuple[float, float, bool]:
    """Return objective, the prescribed stationarity/KKT infinity norm, feasibility."""
    family, data = case["family"], case["data"]
    n = case["dimensions"]["n"]
    if len(point) != n:
        raise ValueError(f"expected {n} coordinates, got {len(point)}")
    if family == "diagonal_quadratic":
        diagonal, center = data["diagonal"], data["center"]
        gradient = [diagonal[i] * (point[i] - center[i]) for i in range(n)]
        value = 0.5 * sum(diagonal[i] * (point[i] - center[i]) ** 2 for i in range(n))
        if "bounds" not in case:
            return value, max(abs(item) for item in gradient), True
        lower, upper = case["bounds"]["lower"], case["bounds"]["upper"]
        feasibility = all(lower[i] - 1e-12 <= point[i] <= upper[i] + 1e-12 for i in range(n))
        # Unit projected-gradient residual: ||x - P_[l,u](x - grad)||_inf.
        projected = [min(upper[i], max(lower[i], point[i] - gradient[i])) for i in range(n)]
        return value, max(abs(point[i] - projected[i]) for i in range(n)), feasibility
    if family == "rotated_quadratic":
        center, h = data["center"], data["hessian"]
        difference = [point[i] - center[i] for i in range(n)]
        gradient = [sum(h[row][column] * difference[column] for column in range(n)) for row in range(n)]
        return 0.5 * dot(difference, gradient), max(abs(item) for item in gradient), True
    if family == "rosenbrock":
        gradient = [0.0] * n
        value = 0.0
        for i in range(n - 1):
            residual = point[i + 1] - point[i] * point[i]
            value += 100.0 * residual * residual + (1.0 - point[i]) ** 2
            gradient[i] += -400.0 * point[i] * residual - 2.0 * (1.0 - point[i])
            gradient[i + 1] += 200.0 * residual
        return value, max(abs(item) for item in gradient), True
    if family == "l2_logistic":
        matrix, labels, lam = data["x"], data["y"], data["lambda"]
        m = case["dimensions"]["m"]
        gradient = [lam * coordinate for coordinate in point]
        value = 0.5 * lam * dot(point, point)
        for row, label in zip(matrix, labels, strict=True):
            margin = label * dot(row, point)
            value += logistic_loss(margin) / m
            # sigmoid(-margin), in an overflow-safe form.
            probability = math.exp(-margin) / (1.0 + math.exp(-margin)) if margin >= 0.0 else 1.0 / (1.0 + math.exp(margin))
            for j in range(n):
                gradient[j] -= row[j] * label * probability / m
        return value, max(abs(item) for item in gradient), True
    if family == "lasso":
        matrix, response, lam = data["x"], data["y"], data["lambda"]
        m = case["dimensions"]["m"]
        gradient = [0.0] * n
        squared = 0.0
        for row, observed in zip(matrix, response, strict=True):
            residual = dot(row, point) - observed
            squared += residual * residual
            for j in range(n):
                gradient[j] += row[j] * residual / m
        kkt = []
        for coordinate, derivative in zip(point, gradient, strict=True):
            if abs(coordinate) <= ZERO_COORDINATE:
                kkt.append(max(abs(derivative) - lam, 0.0))
            else:
                kkt.append(abs(derivative + lam * math.copysign(1.0, coordinate)))
        return 0.5 * squared / m + lam * sum(abs(value) for value in point), max(kkt), True
    raise ValueError(f"unsupported fixture family {family}")


def parse_point(text: str) -> list[float]:
    if not text:
        return []
    return [float(value) for value in text.split(",")]


def parse_harness(path: Path) -> list[dict[str, str]]:
    expected = {"case_id", "method", "status", "objective", "residual", "iterations", "callbacks", "median_ms", "times_ms", "x"}
    with path.open(newline="", encoding="utf-8") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if not rows or set(rows[0]) != expected:
        actual = set(rows[0]) if rows else set()
        raise ValueError(f"harness columns must be exactly {sorted(expected)}, got {sorted(actual)}")
    return rows


def markdown(report: dict[str, Any]) -> str:
    lines = ["# Gale versus Python optimization comparison", "", report["summary"], "", "Timing ratio is Gale/Python; below 1 favors Gale. Near ties need more evidence.", "", "| case | method | status | objective error | stationarity/KKT | valid | Gale ms | Python ms | ratio |", "|---|---|---|---:|---:|:---:|---:|---:|---:|"]
    for row in report["rows"]:
        ratio = "" if row["timing_ratio"] is None else f"{row['timing_ratio']:.3f}"
        python_ms = "" if row["python_median_ms"] is None else f"{row['python_median_ms']:.3f}"
        objective_error = "" if row["objective_error"] is None else f"{row['objective_error']:.3e}"
        metric = "" if row["metric"] is None else f"{row['metric']:.3e}"
        median = "" if row["median_ms"] is None else f"{row['median_ms']:.3f}"
        lines.append(f"| {row['case_id']} | {row['method']} | {row['status']} | {objective_error} | {metric} | {row['accuracy_valid']} | {median} | {python_ms} | {ratio} |")
    lines.extend(("", "## Environments", "", "```json", json.dumps(report["environments"], indent=2, sort_keys=True), "```", ""))
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("harness", type=Path, help="tab-separated Gale harness result")
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--references", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--python-baseline", type=Path, default=DEFAULT_PYTHON_BASELINE)
    parser.add_argument("--json-out", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--markdown-out", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args()

    fixtures = json.loads(args.fixtures.read_text(encoding="utf-8"))
    references = json.loads(args.references.read_text(encoding="utf-8"))
    baseline = json.loads(args.python_baseline.read_text(encoding="utf-8"))
    cases = {case["id"]: case for case in fixtures["cases"]}
    expected = {row["case_id"]: row["expected"] for row in references["results"]}
    baseline_rows = {row["case_id"]: row for row in baseline["rows"]}
    rows = []
    harness = parse_harness(args.harness)
    if len(harness) != len(cases) or {raw["case_id"] for raw in harness} != set(cases):
        raise ValueError("harness must contain each of the nine fixtures exactly once")
    for raw in harness:
        case_id = raw["case_id"]
        try:
            case = cases[case_id]
            point = parse_point(raw["x"])
            reported_objective = float(raw["objective"])
            reported_residual = float(raw["residual"])
            median_ms = float(raw["median_ms"])
            times_ms = parse_point(raw["times_ms"])
            value, metric, feasible = evaluate(case, point)
            reference = float(expected[case_id]["objective"])
            objective_error = abs(value - reference)
            objective_limit = OBJECTIVE_FACTOR * (1.0 + abs(reference))
            all_finite = bool(times_ms) and finite(point + times_ms + [reported_objective, reported_residual, median_ms, value, metric])
            accuracy_valid = (all_finite and feasible and metric <= ACCURACY and objective_error <= objective_limit
                              and abs(reported_objective-value) <= 1e-12*(1+abs(value))
                              and median_ms > 0 and min(times_ms) > 0
                              and math.isclose(median_ms, statistics.median(times_ms), rel_tol=1e-10))
            python = baseline_rows.get(case_id)
            python_ms = python["median_ms"] if python and python.get("timed_endpoint_valid") else None
            ratio = median_ms / python_ms if accuracy_valid and python_ms and python_ms > 0.0 else None
            rows.append({"case_id": case_id, "method": raw["method"], "status": raw["status"], "objective": json_number(value), "reported_objective": json_number(reported_objective), "reported_residual": json_number(reported_residual), "objective_error": json_number(objective_error), "objective_limit": objective_limit, "metric": json_number(metric), "feasible": feasible, "all_finite": all_finite, "accuracy_valid": accuracy_valid, "iterations": int(raw["iterations"]), "callbacks": int(raw["callbacks"]), "median_ms": json_number(median_ms), "times_ms": [json_number(value) for value in times_ms], "python_median_ms": python_ms, "timing_ratio": ratio, "x": [json_number(value) for value in point]})
        except (KeyError, ValueError, OverflowError) as error:
            rows.append({"case_id": case_id, "method": raw.get("method", ""), "status": raw.get("status", ""), "objective_error": None, "metric": None, "median_ms": None, "python_median_ms": None, "timing_ratio": None, "accuracy_valid": False, "error": str(error)})
    failures = sum(not row["accuracy_valid"] for row in rows)
    report = {"schema_version": 1, "accuracy_target": ACCURACY, "objective_factor": OBJECTIVE_FACTOR, "summary": f"{len(rows) - failures}/{len(rows)} Gale harness rows meet independently recomputed accuracy and feasibility checks.", "environments": {"python_baseline": baseline.get("environment"), "fixtures": str(args.fixtures), "references": str(args.references), "harness": str(args.harness)}, "rows": rows}
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    args.markdown_out.write_text(markdown(report), encoding="utf-8")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
