#!/usr/bin/env python3
"""Run portable Gale/pycma solves and independent endpoint/update checks with existing pinned tools."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def run(argv, log, env):
    with log.open('w') as stream:
        result = subprocess.run(argv, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'exit {result.returncode}; see {log}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', default=sys.executable)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seeds', type=int, default=20)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--warmups', type=int, default=5)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    counts = f'{args.warmups} {args.seeds} {args.repeats}'
    initializer = ('set benchmarksJS / Compile / scalaJSMainModuleInitializer := Some('
                   'org.scalajs.linker.interface.ModuleInitializer.mainMethodWithArgs('
                   '"gale.bench.CMAESBench", "main", List('
                   + ', '.join('"' + value + '"' for value in counts.split()) + ')))')
    run(['sbt', f'benchmarksJVM/runMain gale.bench.CMAESBench {counts}',
         'set benchmarksJS / Compile / mainClass := Some("gale.bench.CMAESBench")',
         'set benchmarksJS / scalaJSStage := FullOptStage', initializer, 'benchmarksJS/run'], out/'gale-run.log', env)
    common = [args.python, 'tools/optim-cmaes/compare.py']
    parameters = ['--output', str(out), '--seeds', str(args.seeds), '--repeats', str(args.repeats), '--warmups', str(args.warmups)]
    run(common + ['parse'] + parameters + ['--log', str(out/'gale-run.log')], out/'parse.log', env)
    run(common + ['python'] + parameters, out/'python-run.log', env)
    run(common + ['verify'] + parameters, out/'verify.log', env)
    print(f'All endpoint validity, work, repeatability and update-reference checks pass: {out}')
    print((out/'comparison.md').read_text())


if __name__ == '__main__':
    main()
