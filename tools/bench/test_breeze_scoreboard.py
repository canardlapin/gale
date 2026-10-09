#!/usr/bin/env python3
"""Self-test for breeze_scoreboard.py on the synthetic fixtures in testdata/.

Run isolated from the repository: ``python3 -I tools/bench/test_breeze_scoreboard.py``.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "testdata"

sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location(
    "breeze_scoreboard", HERE / "breeze_scoreboard.py"
)
scoreboard = importlib.util.module_from_spec(_spec)
sys.modules["breeze_scoreboard"] = scoreboard
_spec.loader.exec_module(scoreboard)


def run(*args: str) -> tuple[int, str, str]:
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = scoreboard.main([*args, "--commit", "0123abc"])
    return code, out.getvalue(), err.getvalue()


def lane(name: str, sidecar: str, results: str, *extra: str) -> tuple[int, str, str]:
    return run(
        "--lane", name, "--netlib", str(DATA / sidecar), str(DATA / results), *extra
    )


def table_row(
    text: str, cls: str, op: str, params: str, backend: str, note: bool = False
) -> list[str]:
    """The row's cells; the trailing note column is dropped unless ``note``."""
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if cells[:4] == [cls, op, params, backend]:
            return cells if note else cells[:-1]
    raise AssertionError(f"no row for {cls}.{op} {params} {backend}:\n{text}")


class LaneA(unittest.TestCase):
    def test_scoreboard(self) -> None:
        code, out, err = lane("A", "laneA-netlib.jsonl", "laneA-results.json")
        self.assertEqual(code, 0, err)
        self.assertIn("Breeze netlib BLAS: `dev.ludovic.netlib.blas.Java11BLAS`", out)
        self.assertIn("Commit: `0123abc`", out)
        self.assertIn("JDK: 25.0.1 (OpenJDK 64-Bit Server VM)", out)
        self.assertIn(
            "**2 ahead, 2 tie, 2 behind, 0 n/a, 0 withheld** of 6 pairs; 1 unpaired.", out
        )
        self.assertEqual(
            table_row(out, "BlasL3BreezeJmh", "Gemm", "n=16", "pure")[-2:],
            ["2.00x", "ahead"],
        )
        self.assertEqual(
            table_row(out, "BlasL3BreezeJmh", "Gemm", "n=256", "pure")[-2:],
            ["0.80x", "behind"],
        )
        self.assertEqual(
            table_row(out, "BlasL1BreezeJmh", "Dot", "n=65536", "backend-insensitive")[
                -1
            ],
            "tie",
        )
        self.assertIn("Unpaired benchmarks: `gale.bench.BlasL3BreezeJmh.galeAtA`", out)
        self.assertIn(
            "1 result(s) outside the gale/breeze naming convention were ignored", out
        )

    def test_average_time_is_inverted(self) -> None:
        code, out, err = lane("A", "laneA-netlib.jsonl", "laneA-results.json")
        self.assertEqual(code, 0, err)
        # avgt: lower time is better, so the speedup is breeze/gale.
        faster = table_row(
            out, "TimedBreezeJmh", "Solve", "n=64", "pure"
        )  # gale 2 us, breeze 4 us
        self.assertEqual(faster[4], "avgt")
        self.assertEqual(faster[-2:], ["2.00x", "ahead"])
        slower = table_row(
            out, "TimedBreezeJmh", "Solve", "n=128", "pure"
        )  # gale 8 us, breeze 4 us
        self.assertEqual(slower[-2:], ["0.50x", "behind"])
        overlap = table_row(
            out, "TimedBreezeJmh", "Solve", "n=256", "pure"
        )  # CIs overlap
        self.assertEqual(overlap[-1], "tie")

    def test_strict_rejects_unpaired(self) -> None:
        code, _, err = lane("A", "laneA-netlib.jsonl", "laneA-results.json", "--strict")
        self.assertEqual(code, 2)
        self.assertIn(
            "--strict: unpaired benchmarks: gale.bench.BlasL3BreezeJmh.galeAtA", err
        )

    def test_rejects_vector_blas(self) -> None:
        code, _, err = lane("A", "laneA-netlib-vectorblas.jsonl", "laneA-results.json")
        self.assertEqual(code, 2)
        self.assertIn(
            "lane A receipt rejected: Breeze BLAS is dev.ludovic.netlib.blas.VectorBLAS",
            err,
        )

    def test_rejects_native_blas(self) -> None:
        code, _, err = lane("A", "laneA-netlib-jniblas.jsonl", "laneA-results.json")
        self.assertEqual(code, 2)
        self.assertIn("JNIBLAS", err)

    def test_rejects_vector_module_from_environment(self) -> None:
        # e.g. --add-modules arriving through JDK_JAVA_OPTIONS: absent from the JMH jvmArgs.
        code, _, err = lane(
            "A", "laneA-netlib-vectormodule.jsonl", "laneA-results.json"
        )
        self.assertEqual(code, 2)
        self.assertIn("vectorModule=true", err)

    def test_rejects_vector_fork(self) -> None:
        code, _, err = lane(
            "A", "laneA-netlib-vector-fork.jsonl", "laneA-results-vector-fork.json"
        )
        self.assertEqual(code, 2)
        self.assertIn("forks ran with jdk.incubator.vector", err)

    def test_rejects_jdk_mismatch(self) -> None:
        code, _, err = lane("A", "laneA-netlib-jdk24.jsonl", "laneA-results.json")
        self.assertEqual(code, 2)
        self.assertIn("does not match results JDK", err)

    def test_rejects_stale_sidecar(self) -> None:
        code, _, err = lane("A", "laneA-netlib-stale.jsonl", "laneA-results.json")
        self.assertEqual(code, 2)
        self.assertIn("stale or foreign sidecar", err)

    def test_rejects_missing_sidecar_record(self) -> None:
        code, _, err = lane("A", "laneA-netlib-vector-fork.jsonl", "laneA-results.json")
        self.assertEqual(code, 2)
        self.assertIn("no netlib sidecar record for", err)

    def test_rejects_zero_score(self) -> None:
        code, _, err = lane(
            "A", "laneA-netlib-vector-fork.jsonl", "laneA-results-zero.json"
        )
        self.assertEqual(code, 2)
        self.assertIn("is not a positive finite number", err)

    def test_rejects_lane_mismatch(self) -> None:
        code, _, err = lane("B", "laneA-netlib.jsonl", "laneA-results.json")
        self.assertEqual(code, 2)
        self.assertIn("do not match --lane B", err)

    def test_rejects_duplicate_results(self) -> None:
        results = str(DATA / "laneA-results.json")
        code, _, err = run(
            "--lane",
            "A",
            "--netlib",
            str(DATA / "laneA-netlib.jsonl"),
            results,
            results,
        )
        self.assertEqual(code, 2)
        self.assertIn("appears more than once", err)


class Caveats(unittest.TestCase):
    def setUp(self) -> None:
        code, self.out, err = lane(
            "A", "laneA-caveats-netlib.jsonl", "laneA-caveats-results.json", "--strict"
        )
        self.assertEqual(code, 0, err)

    def test_not_like_for_like_withholds_verdict(self) -> None:
        csr = table_row(
            self.out, "SparseMatrixBreezeJmh", "CsrMatvec", "density=0.01, n=1000",
            "backend-insensitive", note=True,
        )
        self.assertEqual(csr[-3:-1], ["2.00x", "withheld"])
        self.assertIn("not like-for-like", csr[-1])
        tol = table_row(
            self.out, "LbfgsBreezeJmh", "Rosenbrock", "budget=tolerance",
            "backend-insensitive", note=True,
        )
        self.assertEqual(tol[-2], "withheld")
        self.assertIn("own convergence test", tol[-1])
        self.assertIn("**2 ahead, 0 tie, 1 behind, 0 n/a, 2 withheld** of 5 pairs", self.out)

    def test_caveat_keeps_verdict_and_adds_note(self) -> None:
        inv = table_row(
            self.out, "DenseDecompositionBreezeJmh", "Inv", "n=64", "pure", note=True,
        )
        self.assertEqual(inv[-3:-1], ["0.50x", "behind"])
        self.assertIn("dgetri", inv[-1])

    def test_unqualified_pairs_have_no_note(self) -> None:
        csc = table_row(
            self.out, "SparseMatrixBreezeJmh", "CscMatvec", "density=0.01, n=1000",
            "backend-insensitive", note=True,
        )
        self.assertEqual(csc[-3:], ["2.00x", "ahead", ""])
        fixed = table_row(
            self.out, "LbfgsBreezeJmh", "Rosenbrock", "budget=fixed",
            "backend-insensitive", note=True,
        )
        self.assertEqual(fixed[-3:], ["2.00x", "ahead", ""])


class LaneB(unittest.TestCase):
    def test_backends_pair_with_one_breeze_row(self) -> None:
        code, out, err = lane(
            "B", "laneB-netlib.jsonl", "laneB-results.json", "--strict"
        )
        self.assertEqual(code, 0, err)
        self.assertIn("Breeze netlib BLAS: `dev.ludovic.netlib.blas.VectorBLAS`", out)
        self.assertIn("of 4 pairs; 0 unpaired.", out)
        self.assertEqual(
            table_row(out, "BlasL3BreezeJmh", "Gemm", "n=64", "pure")[-2:],
            ["0.60x", "behind"],
        )
        self.assertEqual(
            table_row(out, "BlasL3BreezeJmh", "Gemm", "n=64", "vector")[-2:],
            ["1.20x", "ahead"],
        )
        chol = table_row(
            out, "FactorizationBreezeJmh", "Chol", "n=64", "vector (gemm-routed only)"
        )
        self.assertEqual(chol[-2:], ["1.80x", "ahead"])
        dot = table_row(out, "BlasL1BreezeJmh", "Dot", "n=65536", "backend-insensitive")
        self.assertEqual(dot[-2:], ["0.75x", "behind"])

    def test_rejects_scalar_blas(self) -> None:
        code, _, err = lane("B", "laneB-netlib-scalar.jsonl", "laneB-results.json")
        self.assertEqual(code, 2)
        self.assertIn("lane B receipt rejected: Breeze BLAS is dev.ludovic.netlib.blas.Java11BLAS", err)

    def test_rejects_missing_vector_module(self) -> None:
        code, _, err = lane("B", "laneB-netlib-novector.jsonl", "laneB-results.json")
        self.assertEqual(code, 2)
        self.assertIn("ran without jdk.incubator.vector", err)


if __name__ == "__main__":
    unittest.main()
