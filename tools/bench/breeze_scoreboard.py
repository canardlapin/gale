#!/usr/bin/env python3
"""Build a gale-vs-Breeze scoreboard from JMH JSON results and the netlib sidecar.

Usage::

    python3 tools/bench/breeze_scoreboard.py --lane A \\
        --netlib benchmarks/jvm/target/breeze-netlib.jsonl \\
        target/laneA.json [more.json ...] [--strict] [-o scoreboard.md]

Inputs
------
* One or more JMH result files written with ``-rf json`` (the ``breezeLaneA`` and
  ``breezeLaneB`` sbt aliases set ``-rf json``; pass ``-rff <path>``).
* One or more netlib sidecars (JSON Lines). Every paired Breeze bench appends one
  record per trial from ``BreezeBenchData.recordNetlib``: lane, benchmark, params,
  JDK, whether ``jdk.incubator.vector`` was resolved, and the ``dev.ludovic.netlib``
  BLAS/LAPACK classes Breeze resolved. The default path is
  ``benchmarks/jvm/target/breeze-netlib.jsonl``; both lane aliases delete it first,
  so each run writes a fresh sidecar. Keep it with the JSON receipt.

Pairing convention
------------------
A paired benchmark is named ``<package>.<Class>.<lib><Op>`` where ``<lib>`` is
``gale`` or ``breeze`` and ``<Op>`` starts with an upper-case letter (for example
``gale.bench.BlasL3BreezeJmh.galeGemm`` and ``...breezeGemm``). A gale row pairs
with the Breeze row that has the same class, the same ``<Op>``, and the same JMH
params once the gale-only ``backend`` param is removed. Each gale backend variant
(``pure``/``vector``) therefore gets its own row against the one Breeze twin.
Benchmarks whose method name does not follow the convention are ignored and
counted; a gale or Breeze row without a twin is unpaired (an error with
``--strict``).

The gale backend column reads:

* ``pure`` / ``vector`` — the ``GaleBackendState`` the gale method ran with;
* ``vector (gemm-routed only)`` — factorizations and least squares, which reach
  the Vector backend only where they route a product through its gemm;
* ``backend-insensitive`` — gale methods without a ``backend`` param, because the
  operation takes no ``Backend`` (L1 ``dot``/``norm2``/``axpyInPlace``, symmetric
  eigen). These rows measure pure gale in both lanes.

Ratio and verdict
-----------------
``ratio`` is expressed as gale speed over Breeze speed, so ``> 1`` means gale is
faster: ``gale/breeze`` for throughput modes and ``breeze/gale`` for time-per-op
modes (``avgt``, ``sample``, ``ss``). The verdict compares JMH's 99.9% confidence
intervals: disjoint intervals give ``ahead`` or ``behind``; overlapping intervals
give ``tie``; a missing interval (a single measurement iteration) gives ``n/a``.

Caveats
-------
``CAVEATS`` maps ``(class, op)`` or ``(class, op, params)`` to a note shown in the
``note`` column. A caveat marked ``withhold`` replaces the verdict with
``withheld``: the two benchmarks do not measure the same work (for example a gale
CSR product against Breeze's CSC product), so no ahead/behind claim is made. Other
caveats keep their verdict and name a known asymmetry (allocation count, flop
count, storage order) that a reader must weigh.

Receipt checks (errors, exit status 2)
--------------------------------------
* Lane A (out-of-box scalar): any forked JVM run with ``jdk.incubator.vector`` in
  its JMH args, any sidecar record with ``vectorModule`` true (this also catches the
  module arriving through ``JDK_JAVA_OPTIONS``), a Breeze BLAS class that is
  ``VectorBLAS`` or native (``JNIBLAS``/``NativeBLAS``), a native LAPACK class, or a
  gale row whose backend is not ``pure``.
* Lane B (SIMD): a sidecar record without the Vector module, or a Breeze BLAS class
  other than ``VectorBLAS``.
* Either lane: an empty sidecar; a sidecar record from another lane; more than one
  BLAS or LAPACK class; sidecar JDK versions that differ from the results'; a paired
  result with no sidecar record, or a sidecar record with no result (a stale or
  foreign sidecar); one benchmark/params combination reported twice; a score that
  is zero, negative or not a number; mismatched modes or units within a pair; and,
  with ``--strict``, any unpaired row.
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
GEMM_ROUTED = {
    ("FactorizationBreezeJmh", "Lu"),
    ("FactorizationBreezeJmh", "Chol"),
    ("FactorizationBreezeJmh", "Solve"),
    ("FactorizationBreezeJmh", "Qr"),
    ("LeastSquaresBreezeJmh", "Lstsq"),
    ("FactorizationLargeBreezeJmh", "Lu"),
    ("FactorizationLargeBreezeJmh", "Chol"),
    ("FactorizationLargeBreezeJmh", "Solve"),
    ("FactorizationLargeBreezeJmh", "Qr"),
    ("FactorizationLargeBreezeJmh", "Lstsq"),
    ("MultiRhsBreezeJmh", "LuSolve"),
    ("MultiRhsBreezeJmh", "CholSolve"),
    ("DenseDecompositionBreezeJmh", "Inv"),
    ("DenseDecompositionBreezeJmh", "Det"),
    ("SmallDenseBreezeJmh", "Solve"),
    ("SmallDenseBreezeJmh", "Inv"),
    ("SmallDenseBreezeJmh", "Det"),
}
INSENSITIVE = "backend-insensitive"
WITHHELD = "withheld"


@dataclass(frozen=True)
class Caveat:
    """A known asymmetry in a pairing. ``withhold`` suppresses the verdict
    because the two benchmarks do not measure the same work."""

    note: str
    withhold: bool = False


_NO_CSR = Caveat(
    "not like-for-like: Breeze has no CSR, its twin runs the CSC product", withhold=True
)
_OWN_TOLERANCE = Caveat(
    "not like-for-like: each library stops on its own convergence test", withhold=True
)
_SOFTMAX = "Breeze idiom `exp(x - softmax(x))` makes 2 allocations and an extra pass"
_INV = Caveat(
    "gale `solve(I)` is ~8/3 n^3 flops vs Breeze `dgetri` ~2 n^3 (disfavours gale)"
)
_SUM = Caveat("gale multi-accumulator sum vs Breeze's single-accumulator loop")
_ROWS = "gale row-major: rows contiguous for gale, strided for Breeze (column-major)"
_COLS = "gale row-major: columns strided for gale, contiguous for Breeze (column-major)"

# Keyed by (class, op) or (class, op, ((param, value), ...)); a param-qualified key
# applies when its params are a subset of the pair's params and wins over (class, op).
CAVEATS: dict[tuple, Caveat] = {
    ("SparseMatrixBreezeJmh", "CsrMatvec"): _NO_CSR,
    ("SparseMatrixBreezeJmh", "CsrMatmul"): _NO_CSR,
    ("LbfgsBreezeJmh", "Rosenbrock", (("budget", "tolerance"),)): _OWN_TOLERANCE,
    ("LbfgsBreezeJmh", "Logistic", (("budget", "tolerance"),)): _OWN_TOLERANCE,
    ("ReductionBreezeJmh", "Softmax"): Caveat(_SOFTMAX),
    ("MatrixReductionBreezeJmh", "SoftmaxRows"): Caveat(f"{_SOFTMAX}; {_ROWS}"),
    ("DenseDecompositionBreezeJmh", "Inv"): _INV,
    ("SmallDenseBreezeJmh", "Inv"): _INV,
    ("ReductionBreezeJmh", "Mean"): Caveat(
        "Breeze `stats.mean` is a running mean (a division per element); gale is sum/n"
    ),
    ("ReductionBreezeJmh", "Sum"): _SUM,
    ("MatrixReductionBreezeJmh", "Sum"): _SUM,
    ("MatrixReductionBreezeJmh", "SumRows"): Caveat(_ROWS),
    ("MatrixReductionBreezeJmh", "MaxRows"): Caveat(_ROWS),
    ("MatrixReductionBreezeJmh", "LogSumExpRows"): Caveat(_ROWS),
    ("MatrixReductionBreezeJmh", "SumCols"): Caveat(_COLS),
    ("MatrixReductionBreezeJmh", "MaxCols"): Caveat(_COLS),
}


def caveat(cls: str, op: str, params: tuple) -> Caveat | None:
    for key, value in CAVEATS.items():
        if len(key) == 3 and key[:2] == (cls, op) and set(key[2]) <= set(params):
            return value
    return CAVEATS.get((cls, op))


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

    @property
    def backend_label(self) -> str:
        if "backend" not in self.params:
            return INSENSITIVE
        if self.backend == "vector" and (self.cls, self.op) in GEMM_ROUTED:
            return "vector (gemm-routed only)"
        return self.backend


def _num(value) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _params_key(benchmark: str, params: dict) -> tuple:
    return (benchmark, tuple(sorted(params.items())))


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
            key = _params_key(name, params)
            if key in seen:
                raise ReceiptError(
                    f"{name} {params} appears more than once across the result files"
                )
            seen.add(key)
            metric = entry["primaryMetric"]
            score = _num(metric["score"])
            if not score > 0 or math.isinf(score):
                raise ReceiptError(
                    f"{name} {params}: score {metric['score']!r} is not a positive finite number"
                )
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
                    score=score,
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


def _sidecar_params(text: str) -> dict:
    return dict(
        part.split("=", 1) for part in str(text or "").split(",") if "=" in part
    )


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

    sidecar_jdks = {str(r.get("jdk")) for r in sidecar}
    result_jdks = {r.jdk for r in rows}
    if sidecar_jdks != result_jdks:
        raise ReceiptError(
            f"sidecar JDK {sorted(sidecar_jdks)} does not match results JDK {sorted(result_jdks)}"
        )

    recorded = {
        _params_key(str(r.get("benchmark")), _sidecar_params(r.get("params")))
        for r in sidecar
    }
    measured = {_params_key(r.benchmark, r.params) for r in rows}
    missing = sorted(f"{b} {dict(p)}" for b, p in measured - recorded)
    if missing:
        raise ReceiptError(f"no netlib sidecar record for: {', '.join(missing)}")
    stale = sorted(f"{b} {dict(p)}" for b, p in recorded - measured)
    if stale:
        raise ReceiptError(
            f"sidecar records with no matching result (stale or foreign sidecar): {', '.join(stale)}"
        )

    vector_modules = {r.get("vectorModule") for r in sidecar}
    if lane == "A":
        if True in vector_modules:
            raise ReceiptError(
                "lane A receipt rejected: a fork resolved jdk.incubator.vector (sidecar vectorModule=true)"
            )
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
    else:
        if vector_modules != {True}:
            raise ReceiptError(
                "lane B receipt rejected: a fork ran without jdk.incubator.vector (sidecar vectorModule)"
            )
        if blas_cls != VECTOR_BLAS:
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
    strict: bool = False,
) -> str:
    pairs, unpaired = pair(rows)
    if strict and unpaired:
        raise ReceiptError(
            "--strict: unpaired benchmarks: "
            + ", ".join(sorted(f"{r.benchmark} {r.params}" for r in unpaired))
        )
    pairs.sort(
        key=lambda p: (p[0].cls, p[0].op, _size_key(p[0].pair_params), p[0].backend)
    )
    jdks = sorted({f"{r.jdk} ({r.vm})" for r in rows})
    jvm_args = sorted({" ".join(r.jvm_args) or "(none)" for r in rows})
    counts = {"ahead": 0, "tie": 0, "behind": 0, "n/a": 0, WITHHELD: 0}
    body = []
    for g, b in pairs:
        ratio, v = verdict(g, b)
        note = caveat(g.cls, g.op, g.pair_params)
        if note is not None and note.withhold:
            v = WITHHELD
        counts[v] += 1
        params = ", ".join(f"{k}={val}" for k, val in g.pair_params) or "-"
        body.append(
            f"| {g.cls} | {g.op} | {params} | {g.backend_label} | {g.mode} | {_fmt(g.score, g.error)} | "
            f"{_fmt(b.score, b.error)} | {g.unit} | {ratio:.2f}x | {v} | {note.note if note else ''} |"
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
        f"**{counts['ahead']} ahead, {counts['tie']} tie, {counts['behind']} behind, {counts['n/a']} n/a, "
        f"{counts[WITHHELD]} withheld** of {len(pairs)} pairs; {len(unpaired)} unpaired. Ratio > 1 means gale "
        "is faster; verdicts use non-overlapping 99.9% CIs. `backend-insensitive` rows run pure gale in every "
        "lane. A `withheld` pair is not like-for-like (see its note) and carries no verdict; other notes name "
        "a known asymmetry behind a verdict.",
        "",
        "| class | op | params | gale backend | mode | gale | breeze | unit | gale speedup | verdict | note |",
        "|---|---|---|---|---|---:|---:|---|---:|---|---|",
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
        "--strict",
        action="store_true",
        help="fail when any gale or Breeze row is unpaired",
    )
    parser.add_argument(
        "-o", "--output", type=Path, help="write markdown here instead of stdout"
    )
    args = parser.parse_args(argv)
    try:
        rows, ignored = load_results(args.results)
        blas, lapack = validate(args.lane, rows, load_sidecars(args.netlib))
        commit = args.commit or git_commit(Path(__file__).resolve().parent)
        text = render(
            args.lane, rows, blas, lapack, commit, args.machine, ignored, args.strict
        )
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
