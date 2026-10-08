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


def table_row(text: str, cls: str, op: str, params: str, backend: str) -> list[str]:
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip("|").split("|")]
        if cells[:4] == [cls, op, params, backend]:
            return cells
    raise AssertionError(f"no row for {cls}.{op} {params} {backend}:\n{text}")


class LaneA(unittest.TestCase):
    def test_scoreboard(self) -> None:
        code, out, err = run(
            "--lane",
            "A",
            "--netlib",
            str(DATA / "laneA-netlib.jsonl"),
            str(DATA / "laneA-results.json"),
        )
        self.assertEqual(code, 0, err)
        self.assertIn("Breeze netlib BLAS: `dev.ludovic.netlib.blas.Java11BLAS`", out)
        self.assertIn("Commit: `0123abc`", out)
        self.assertIn("JDK: 25.0.1 (OpenJDK 64-Bit Server VM)", out)
        self.assertIn("**2 ahead, 1 tie, 1 behind, 0 n/a** of 4 pairs", out)
        ahead = table_row(out, "BlasL3BreezeJmh", "Gemm", "n=16", "pure")
        self.assertEqual(ahead[-2:], ["2.00x", "ahead"])
        behind = table_row(out, "BlasL3BreezeJmh", "Gemm", "n=256", "pure")
        self.assertEqual(behind[-2:], ["0.80x", "behind"])
        tie = table_row(out, "BlasL1BreezeJmh", "Dot", "n=65536", "pure")
        self.assertEqual(tie[-1], "tie")
        timed = table_row(out, "TimedBreezeJmh", "Solve", "n=64", "pure")
        self.assertEqual(timed[-2:], ["2.00x", "ahead"])  # avgt: lower time is better
        self.assertIn("Unpaired benchmarks: `gale.bench.BlasL3BreezeJmh.galeAtA`", out)
        self.assertIn(
            "1 result(s) outside the gale/breeze naming convention were ignored", out
        )

    def test_rejects_vector_blas(self) -> None:
        code, _, err = run(
            "--lane",
            "A",
            "--netlib",
            str(DATA / "laneA-netlib-vectorblas.jsonl"),
            str(DATA / "laneA-results.json"),
        )
        self.assertEqual(code, 2)
        self.assertIn(
            "lane A receipt rejected: Breeze BLAS is dev.ludovic.netlib.blas.VectorBLAS",
            err,
        )

    def test_rejects_native_blas(self) -> None:
        code, _, err = run(
            "--lane",
            "A",
            "--netlib",
            str(DATA / "laneA-netlib-jniblas.jsonl"),
            str(DATA / "laneA-results.json"),
        )
        self.assertEqual(code, 2)
        self.assertIn("JNIBLAS", err)

    def test_rejects_vector_fork(self) -> None:
        code, _, err = run(
            "--lane",
            "A",
            "--netlib",
            str(DATA / "laneA-netlib.jsonl"),
            str(DATA / "laneA-results-vector-fork.json"),
        )
        self.assertEqual(code, 2)
        self.assertIn("forks ran with jdk.incubator.vector", err)

    def test_rejects_lane_mismatch(self) -> None:
        code, _, err = run(
            "--lane",
            "B",
            "--netlib",
            str(DATA / "laneA-netlib.jsonl"),
            str(DATA / "laneA-results.json"),
        )
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


class LaneB(unittest.TestCase):
    def test_backends_pair_with_one_breeze_row(self) -> None:
        code, out, err = run(
            "--lane",
            "B",
            "--netlib",
            str(DATA / "laneB-netlib.jsonl"),
            str(DATA / "laneB-results.json"),
        )
        self.assertEqual(code, 0, err)
        self.assertIn("Breeze netlib BLAS: `dev.ludovic.netlib.blas.VectorBLAS`", out)
        self.assertEqual(
            table_row(out, "BlasL3BreezeJmh", "Gemm", "n=64", "pure")[-2:],
            ["0.60x", "behind"],
        )
        self.assertEqual(
            table_row(out, "BlasL3BreezeJmh", "Gemm", "n=64", "vector")[-2:],
            ["1.20x", "ahead"],
        )

    def test_rejects_scalar_blas(self) -> None:
        code, _, err = run(
            "--lane",
            "B",
            "--netlib",
            str(DATA / "laneA-netlib.jsonl"),
            str(DATA / "laneB-results.json"),
        )
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
