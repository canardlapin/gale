#!/usr/bin/env python3
"""Re-run portable Gale/Python solve comparisons from the repository root.

Use JDK 21, Node 22, sbt, and a Python environment containing requirements.txt.
No package installation, global configuration, or publication is performed.
"""
from pathlib import Path
import argparse
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def run(argv, log, env):
    with log.open("w", encoding="utf-8") as output:
        result = subprocess.run(argv, cwd=ROOT, env=env, stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}); see {log}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable, help="Python with the pinned comparator packages")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/verification/optim-python")
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", OMP_NUM_THREADS="1")
    log = out / "gale-run.log"
    run(["sbt", "benchmarksJVM/runMain gale.bench.OptimizationBench 100 11",
         'set benchmarksJS / Compile / mainClass := Some("gale.bench.OptimizationBench")',
         "set benchmarksJS / scalaJSStage := FullOptStage", "benchmarksJS/run"], log, env)
    lines = log.read_text().splitlines()
    blocks = [lines[i:i+10] for i, line in enumerate(lines) if line.startswith("case_id\t")]
    if len(blocks) != 2 or any(len(block) != 10 for block in blocks):
        raise RuntimeError("expected exactly two nine-fixture result blocks")
    for platform, block in zip(("jvm", "node"), blocks, strict=True):
        (out / f"gale-{platform}-v1.tsv").write_text("\n".join(block) + "\n")
    run([args.python, "tools/optim-python/optim_protocol.py", "benchmark", "--warmup", "10", "--repeats", "11",
         "--benchmarks", str(out / "python-baseline-v1.json"), "--tsv", str(out / "python-baseline-v1.tsv")], out / "python-run.log", env)
    for platform in ("jvm", "node"):
        run([sys.executable, "tools/optim-python/compare_gale.py", str(out / f"gale-{platform}-v1.tsv"),
             "--python-baseline", str(out / "python-baseline-v1.json"),
             "--json-out", str(out / f"gale-{platform}-comparison-v1.json"),
             "--markdown-out", str(out / f"gale-{platform}-comparison-v1.md")], out / f"verify-{platform}.log", env)
    print(f"Both nine-case comparisons pass; reports in {out}")


if __name__ == "__main__":
    main()
