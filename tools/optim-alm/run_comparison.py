#!/usr/bin/env python3
"""Run sequential JVM, full-optimized Node, and calibrated SciPy constrained benchmarks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def run(argv, log, env):
    with log.open('w') as stream:
        result = subprocess.run(argv, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'exit {result.returncode}: see {log}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--python', default=sys.executable)
    p.add_argument('--warmups', type=int, default=50)
    p.add_argument('--repeats', type=int, default=15)
    a = p.parse_args()
    if a.warmups < 0 or a.repeats < 1: p.error('invalid run counts')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    compare = [a.python, 'tools/optim-alm/compare.py']
    common = ['--output', str(out), '--warmups', str(a.warmups), '--repeats', str(a.repeats)]
    run(compare + ['fixtures'] + common, out / 'fixtures.log', env)
    path = str(out / 'fixtures.tsv')
    if any(c in path for c in ('"', '\\', '\n', '\r')):
        raise ValueError('output path contains unsupported sbt quoting characters')
    parameters = f'{a.warmups} {a.repeats} "{path}" JVM'
    initializer = ('set benchmarksJS / Compile / scalaJSMainModuleInitializer := Some('
        'org.scalajs.linker.interface.ModuleInitializer.mainMethodWithArgs('
        '"gale.bench.OptimizationALMBench", "main", List('
        + ', '.join(json.dumps(s) for s in (str(a.warmups), str(a.repeats), path, 'Node')) + ')))')
    run(['sbt', f'benchmarksJVM/runMain gale.bench.OptimizationALMBench {parameters}',
         'set benchmarksJS / Compile / mainClass := Some("gale.bench.OptimizationALMBench")',
         'set benchmarksJS / scalaJSStage := FullOptStage', initializer, 'benchmarksJS/run'], out / 'gale-run.log', env)
    run(compare + ['parse'] + common + ['--log', str(out / 'gale-run.log')], out / 'parse.log', env)
    run(compare + ['python'] + common, out / 'python-run.log', env)
    run(compare + ['verify'] + common, out / 'verify.log', env)
    paths = sorted(set(list((ROOT / 'core/shared/src/main/scala/gale/optim').glob('*.scala')) +
        list((ROOT / 'core/shared/src/test/scala/gale/optim').glob('*.scala')) +
        list((ROOT / 'tools/optim-alm').glob('*.py')) + [ROOT / 'tools/optim-alm/requirements.txt',
        ROOT / 'benchmarks/shared/src/main/scala/gale/bench/OptimizationALMBench.scala', ROOT / 'build.sbt']))
    manifest = dict(git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        platform=platform.platform(), warmups=a.warmups, repeats=a.repeats,
        java=subprocess.check_output([str(Path(env['JAVA_HOME']) / 'bin/java'), '-version'], stderr=subprocess.STDOUT, text=True),
        node=subprocess.check_output(['node', '--version'], text=True).strip())
    (out / 'source-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print((out / 'comparison.md').read_text())


if __name__ == '__main__': main()
