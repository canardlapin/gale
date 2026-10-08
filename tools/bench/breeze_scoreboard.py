#!/usr/bin/env python3
"""Build a gale-vs-Breeze scoreboard from JMH JSON results and the netlib sidecar.

Usage::

    python3 tools/bench/breeze_scoreboard.py --lane A \\
        --netlib benchmarks/jvm/target/breeze-netlib.jsonl \\
        target/laneA.json [more.json ...] [-o scoreboard.md]

Inputs
------
* One or more JMH result files written with ``-rf json`` (the ``breezeLaneA`` and
  ``breezeLaneB`` sbt aliases set ``-rf json``; pass ``-rff <path>``).
* One or more netlib sidecars (JSON Lines). Every paired Breeze bench appends one
  record per trial from ``BreezeBenchData.recordNetlib``: lane, benchmark, params,
  JDK, and the ``dev.ludovic.netlib`` BLAS/LAPACK classes Breeze resolved. The
  default path is ``target/breeze-netlib.jsonl`` relative to the JMH fork's working
  directory, i.e. ``benchmarks/jvm/target/breeze-netlib.jsonl``.

Pairing convention
------------------
A paired benchmark is named ``<package>.<Class>.<lib><Op>`` where ``<lib>`` is
``gale`` or ``breeze`` and ``<Op>`` starts with an upper-case letter (for example
``gale.bench.BlasL3BreezeJmh.galeGemm`` and ``...breezeGemm``). A gale row pairs
with the Breeze row that has the same class, the same ``<Op>``, and the same JMH
params once the gale-only ``backend`` param is removed. Each gale backend variant
(``pure``/``vector``) therefore gets its own row against the one Breeze twin.
Benchmarks whose method name does not follow the convention are ignored and
counted; a gale or Breeze row without a twin is listed as unpaired.

Ratio and verdict
-----------------
``ratio`` is expressed as gale speed over Breeze speed, so ``> 1`` means gale is
faster: ``gale/breeze`` for throughput modes and ``breeze/gale`` for time-per-op
modes (``avgt``, ``sample``, ``ss``). The verdict compares JMH's 99.9% confidence
intervals: disjoint intervals give ``ahead`` or ``behind``; overlapping intervals
give ``tie``; a missing interval (a single measurement iteration) gives ``n/a``.

Receipt checks (errors, exit status 2)
--------------------------------------
* Lane A (out-of-box scalar): any forked JVM run with ``jdk.incubator.vector``, a
  Breeze BLAS class that is ``VectorBLAS`` or native (``JNIBLAS``/``NativeBLAS``),
  a native LAPACK class, or a gale row whose backend is not ``pure``.
* Lane B (SIMD): a Breeze BLAS class other than ``VectorBLAS``.
* Either lane: an empty sidecar, a sidecar record from another lane, more than one
  BLAS or LAPACK class, a paired Breeze row with no sidecar record, or one
  benchmark/params combination reported twice.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

NAME_RE = re.compile(
    r"^(?P<prefix>.*)\.(?P<cls>[A-Za-z_]\w*)\.(?P<lib>gale|breeze)(?P<op>[A-Z]\w*)$"
)
TIME_MODES = {"avgt", "sample", "ss"}
VECTOR_BLAS = "dev.ludovic.netlib.blas.VectorBLAS"
NATIVE_MARKERS = ("JNIBLAS", "NativeBLAS", "JNILAPACK", "NativeLAPACK")


class ReceiptError(Exception):
    """A receipt that must not be published as a scoreboard."""


@dataclass(frozen=True)
class Row:
    benchmark: str
    cls: str
    lib: str
    op: str
    params: dict
    mode: str
    unit: str
    score: float
    error: float
    lo: float
    hi: float
    jdk: str
    vm: str
    jvm_args: tuple

    @property
    def pair_params(self) -> tuple:
        return tuple(sorted((k, v) for k, v in self.params.items() if k != "backend"))

    @property
    def backend(self) -> str:
        return self.params.get("backend", "-")


def _num(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def load_results(paths: list[Path]) -> tuple[list[Row], int]:
    rows: list[Row] = []
    ignored = 0
    seen: set[tuple] = set()
    for path in paths:
        data = json.loads(path.read_text())
        if not isinstance(data, list):
            raise ReceiptError(f"{path}: expected a JMH -rf json array")
        for entry in data:
            name = entry["benchmark"]
            match = NAME_RE.match(name)
            if not match:
                ignored += 1
                continue
            params = {str(k): str(v) for k, v in (entry.get("params") or {}).items()}
            key = (name, tuple(sorted(params.items())))
            if key in seen:
                raise ReceiptError(
                    f"{name} {params} appears more than once across the result files"
                )
            seen.add(key)
            metric = entry["primaryMetric"]
            lo, hi = (_num(x) for x in metric.get("scoreConfidence", ["NaN", "NaN"]))
            rows.append(
                Row(
                    benchmark=name,
                    cls=match["cls"],
                    lib=match["lib"],
                    op=match["op"],
                    params=params,
                    mode=entry["mode"],
                    unit=metric["scoreUnit"],
                    score=_num(metric["score"]),
                    error=_num(metric.get("scoreError")),
                    lo=lo,
                    hi=hi,
                    jdk=str(entry.get("jdkVersion", "?")),
                    vm=str(entry.get("vmName", "?")),
                    jvm_args=tuple(entry.get("jvmArgs") or ()),
                )
            )
    return rows, ignored


def load_sidecars(paths: list[Path]) -> list[dict]:
    records: list[dict] = []
    for path in paths:
        for number, line in enumerate(path.read_text().splitlines(), start=1):
            if line.strip():
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError as exc:
                    raise ReceiptError(
                        f"{path}:{number}: malformed sidecar line ({exc})"
                    ) from exc
    return records


def validate(lane: str, rows: list[Row], sidecar: list[dict]) -> tuple[str, str]:
    if not sidecar:
        raise ReceiptError(
            "netlib sidecar has no records; the receipt cannot name Breeze's BLAS"
        )
    lanes = {r.get("lane") for r in sidecar}
    if lanes != {lane}:
        raise ReceiptError(
            f"sidecar lanes {sorted(map(str, lanes))} do not match --lane {lane}"
        )
    blas = {r.get("blas") for r in sidecar}
    lapack = {r.get("lapack") for r in sidecar}
    if len(blas) != 1 or len(lapack) != 1:
        raise ReceiptError(
            f"mixed netlib implementations in one receipt: blas={sorted(blas)} lapack={sorted(lapack)}"
        )
    blas_cls, lapack_cls = blas.pop(), lapack.pop()

    recorded = {r.get("benchmark") for r in sidecar}
    missing = sorted({row.benchmark for row in rows if row.lib == "breeze"} - recorded)
    if missing:
        raise ReceiptError(
            f"no netlib sidecar record for Breeze benchmarks: {', '.join(missing)}"
        )

    if lane == "A":
        if blas_cls == VECTOR_BLAS or any(m in blas_cls for m in NATIVE_MARKERS):
            raise ReceiptError(
                f"lane A receipt rejected: Breeze BLAS is {blas_cls} (must be scalar Java BLAS)"
            )
        if any(m in lapack_cls for m in NATIVE_MARKERS):
            raise ReceiptError(
                f"lane A receipt rejected: Breeze LAPACK is native ({lapack_cls})"
            )
        vector_jvms = sorted(
            {
                r.benchmark
                for r in rows
                if any("jdk.incubator.vector" in a for a in r.jvm_args)
            }
        )
        if vector_jvms:
            raise ReceiptError(
                f"lane A receipt rejected: forks ran with jdk.incubator.vector: {', '.join(vector_jvms)}"
            )
        non_pure = sorted(
            {
                f"{r.benchmark}[{r.backend}]"
                for r in rows
                if r.lib == "gale" and r.backend not in ("pure", "-")
            }
        )
        if non_pure:
            raise ReceiptError(
                f"lane A receipt rejected: gale must run pure: {', '.join(non_pure)}"
            )
    elif blas_cls != VECTOR_BLAS:
        raise ReceiptError(
            f"lane B receipt rejected: Breeze BLAS is {blas_cls}, expected {VECTOR_BLAS}"
        )
    return blas_cls, lapack_cls


def verdict(gale: Row, breeze: Row) -> tuple[float, str]:
    time_mode = gale.mode in TIME_MODES
    ratio = breeze.score / gale.score if time_mode else gale.score / breeze.score
    if any(math.isnan(x) for x in (gale.lo, gale.hi, breeze.lo, breeze.hi)):
        return ratio, "n/a"
    if gale.lo > breeze.hi:
        return ratio, "behind" if time_mode else "ahead"
    if gale.hi < breeze.lo:
        return ratio, "ahead" if time_mode else "behind"
    return ratio, "tie"


def pair(rows: list[Row]) -> tuple[list[tuple[Row, Row]], list[Row]]:
    breeze = {(r.cls, r.op, r.pair_params): r for r in rows if r.lib == "breeze"}
    pairs: list[tuple[Row, Row]] = []
    used: set[tuple] = set()
    unpaired: list[Row] = []
    for g in (r for r in rows if r.lib == "gale"):
        key = (g.cls, g.op, g.pair_params)
        b = breeze.get(key)
        if b is None:
            unpaired.append(g)
            continue
        if (g.mode, g.unit) != (b.mode, b.unit):
            raise ReceiptError(
                f"{g.benchmark} vs {b.benchmark}: mode/unit differ ({g.mode} {g.unit} vs {b.mode} {b.unit})"
            )
        pairs.append((g, b))
        used.add(key)
    unpaired += [b for k, b in breeze.items() if k not in used]
    return pairs, unpaired


def _size_key(params: tuple) -> tuple:
    return tuple((k, int(v)) if v.lstrip("-").isdigit() else (k, v) for k, v in params)


def _fmt(score: float, error: float) -> str:
    err = "?" if math.isnan(error) else f"{error:.4g}"
    return f"{score:.6g} ± {err}"


def git_commit(repo: Path) -> str:
    try:
        sha = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def render(
    lane: str,
    rows: list[Row],
    blas: str,
    lapack: str,
    commit: str,
    machine: str | None,
    ignored: int,
) -> str:
    pairs, unpaired = pair(rows)
    pairs.sort(
        key=lambda p: (p[0].cls, p[0].op, _size_key(p[0].pair_params), p[0].backend)
    )
    jdks = sorted({f"{r.jdk} ({r.vm})" for r in rows})
    jvm_args = sorted({" ".join(r.jvm_args) or "(none)" for r in rows})
    counts = {"ahead": 0, "tie": 0, "behind": 0, "n/a": 0}
    body = []
    for g, b in pairs:
        ratio, v = verdict(g, b)
        counts[v] += 1
        params = ", ".join(f"{k}={val}" for k, val in g.pair_params) or "-"
        body.append(
            f"| {g.cls} | {g.op} | {params} | {g.backend} | {g.mode} | {_fmt(g.score, g.error)} | "
            f"{_fmt(b.score, b.error)} | {g.unit} | {ratio:.2f}x | {v} |"
        )
    title = "out-of-box scalar" if lane == "A" else "SIMD"
    out = [
        f"# gale vs Breeze scoreboard — lane {lane} ({title})",
        "",
        f"- Lane: {lane}",
        f"- JDK: {'; '.join(jdks) or 'unknown'}",
        f"- Fork JVM args: {'; '.join(jvm_args)}",
        f"- Breeze netlib BLAS: `{blas}`",
        f"- Breeze netlib LAPACK: `{lapack}`",
        f"- Commit: `{commit}`",
    ]
    if machine:
        out.append(f"- Machine: {machine}")
    out += [
        "",
        f"**{counts['ahead']} ahead, {counts['tie']} tie, {counts['behind']} behind, {counts['n/a']} n/a** "
        f"of {len(pairs)} pairs. Ratio > 1 means gale is faster; verdicts use non-overlapping 99.9% CIs.",
        "",
        "| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict |",
        "|---|---|---|---|---|---:|---:|---|---:|---|",
        *body,
    ]
    if unpaired:
        out += [
            "",
            "Unpaired benchmarks: "
            + ", ".join(sorted(f"`{r.benchmark}` {r.params}" for r in unpaired)),
        ]
    if ignored:
        out += [
            "",
            f"{ignored} result(s) outside the gale/breeze naming convention were ignored.",
        ]
    return "\n".join(out) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "results", nargs="+", type=Path, help="JMH -rf json result files"
    )
    parser.add_argument("--lane", required=True, choices=["A", "B"])
    parser.add_argument(
        "--netlib",
        action="append",
        required=True,
        type=Path,
        help="netlib sidecar (repeatable)",
    )
    parser.add_argument(
        "--commit", help="source commit (default: git rev-parse HEAD of this checkout)"
    )
    parser.add_argument(
        "--machine", help="free-text machine description for the header"
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="write markdown here instead of stdout"
    )
    args = parser.parse_args(argv)
    try:
        rows, ignored = load_results(args.results)
        blas, lapack = validate(args.lane, rows, load_sidecars(args.netlib))
        commit = args.commit or git_commit(Path(__file__).resolve().parent)
        text = render(args.lane, rows, blas, lapack, commit, args.machine, ignored)
    except ReceiptError as exc:
        print(f"breeze_scoreboard: error: {exc}", file=sys.stderr)
        return 2
    if args.output:
        args.output.write_text(text)
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
