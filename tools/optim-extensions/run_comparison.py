#!/usr/bin/env python3
"""Reproduce JVM/optimized Node/SciPy comparisons using existing JDK21/Node22/sbt and pinned Python.

Inputs remain the materialized repository fixtures. --output preserves checked-in results.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def run(argv, log, env):
    with log.open("w") as stream:
        result = subprocess.run(argv, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"exit {result.returncode}: see {log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/verification/optim-extensions")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    run(["sbt", "benchmarksJVM/runMain gale.bench.OptimizationExtensionsBench 100 11",
         'set benchmarksJS / Compile / mainClass := Some("gale.bench.OptimizationExtensionsBench")',
         "set benchmarksJS / scalaJSStage := FullOptStage", "benchmarksJS/run"], out / "gale-run.log", env)
    lines = (out / "gale-run.log").read_text().splitlines()
    blocks = [lines[i:i+11] for i, line in enumerate(lines) if line.startswith("case_id\t")]
    if len(blocks) != 2 or any(len(block) != 11 for block in blocks):
        raise RuntimeError("expected two ten-fixture tables")
    for name, block in zip(("jvm", "node"), blocks, strict=True):
        (out / f"{name}.tsv").write_text("\n".join(block) + "\n")
    run([args.python, "tools/optim-extensions/compare.py", "python", "--baseline", str(out / "python.json")],
        out / "python-run.log", env)
    for name in ("jvm", "node"):
        run([args.python, "tools/optim-extensions/compare.py", "verify", "--harness", str(out / f"{name}.tsv"),
             "--baseline", str(out / "python.json"), "--output", str(out / f"{name}-comparison")],
            out / f"verify-{name}.log", env)
    print(f"Both ten-case Gale comparisons pass independent endpoint/certificate checks: {out}")


if __name__ == "__main__":
    main()
