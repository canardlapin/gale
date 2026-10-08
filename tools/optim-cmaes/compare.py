#!/usr/bin/env python3
"""Pinned pycma solves, scalar endpoint checks, and an independent cmaes update oracle.

The update oracle zeroes negative weights and aligns the initial generation-counter origin; it does not replace update code.
Its even population avoids the packages' different odd-population weight conventions.
"""
import argparse
import csv
import hashlib
import json
import inspect
from importlib.metadata import version
import math
import os
import platform
import sys
from pathlib import Path
import statistics
import time
import warnings

import numpy as np
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import cma
from cmaes import CMA
from threadpoolctl import threadpool_info

ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "docs/verification/optim-cmaes/fixtures.tsv"


def fixtures():
    with FIXTURES.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    for f in rows:
        f.update(dimension=int(f['dimension']), sigma=float(f['sigma']), budget=int(f['budget']),
                 target=float(f['target']), initial=np.array(list(map(float, f['initial'].split(',')))))
    return rows


def numpy_objective(f):
    n, name = f['dimension'], f['case_id']
    v = np.arange(1., n + 1)
    scales = np.power(1e4, np.arange(n) / (n - 1))
    center = np.where(np.arange(n) % 2 == 0, .7, -.6)
    if name == 'rosenbrock-4':
        return lambda x: float(np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2))
    if name == 'rotated-ellipsoid-8':
        def rotated(x):
            y = x - 2 * v * np.dot(v, x) / np.dot(v, v)
            return float(np.dot(scales, y * y))
        return rotated
    if name == 'nonsmooth-8':
        return lambda x: float(np.dot(v, np.abs(x)))
    if name == 'rastrigin-4':
        return lambda x: float(np.sum(x * x + 10 * (1 - np.cos(2 * np.pi * x))))
    if name == 'box-sphere-8':
        return lambda x: float(np.dot(x - center, x - center))
    if name == 'box-boundary-4':
        return lambda x: float(np.dot(x - 2, x - 2))
    return lambda x: float(np.dot(x, x))


def scalar_objective(f, x):
    """Independent scalar endpoint oracle; never used for timed Python evaluations."""
    name, n = f['case_id'], f['dimension']
    if name == 'rosenbrock-4':
        return math.fsum(100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2 for i in range(n - 1))
    if name == 'rotated-ellipsoid-8':
        p = math.fsum((i + 1) * x[i] for i in range(n))
        d = sum((i + 1) ** 2 for i in range(n))
        return math.fsum(1e4 ** (i / (n - 1)) * (x[i] - 2 * (i + 1) * p / d) ** 2 for i in range(n))
    if name == 'nonsmooth-8':
        return math.fsum((i + 1) * abs(x[i]) for i in range(n))
    if name == 'rastrigin-4':
        return math.fsum(v * v + 10 * (1 - math.cos(2 * math.pi * v)) for v in x)
    if name == 'box-sphere-8':
        return math.fsum((v - (.7 if i % 2 == 0 else -.6)) ** 2 for i, v in enumerate(x))
    if name == 'box-boundary-4':
        return math.fsum((v - 2) ** 2 for v in x)
    return math.fsum(v * v for v in x)


def py_solve(f, seed):
    n, budget = f['dimension'], f['budget']
    population = 4 + int(3 * math.log(n))
    mu = population // 2
    positive = np.log(mu + .5) - np.log(np.arange(1, mu + 1))
    weights = list(positive / positive.sum()) + [0.] * (population - mu)
    clean = numpy_objective(f)
    noise = np.random.RandomState(seed ^ 0x54321)
    objective = (lambda x: clean(x) * math.exp(.1 * noise.standard_normal() - .005)) if f['case_id'] == 'noisy-sphere-8' else clean
    bounded = f['case_id'].startswith('box-')
    opts = dict(seed=seed, popsize=population, CMA_active=False, CMA_mirrors=0,
                CMA_diagonal_decoding=0, CMA_recombination_weights=weights,
                conditioncov_alleviate=[math.inf, math.inf], updatecovwait=1, verbose=-9, verb_log=0)
    begin = time.perf_counter_ns()
    es = cma.CMAEvolutionStrategy(f['initial'].copy(), f['sigma'], opts)
    evaluations = generations = sampled = rejected = 0
    best_value, best_x = math.inf, None
    status = 'EvaluationLimit'
    while evaluations < budget and generations < 10000:
        count = min(population, budget - evaluations)
        points = []
        attempts = 0
        while len(points) < count and attempts < 100000:
            proposals = es.ask(number=min(count - len(points), 100000 - attempts))
            attempts += len(proposals)
            sampled += len(proposals)
            for x in proposals:
                if not bounded or np.all((x >= -1) & (x <= 1)):
                    points.append(x)
                else:
                    rejected += 1
        if len(points) < count:
            status = 'SamplingLimit'
            break
        if all(np.array_equal(x, es.mean) for x in points):
            status = 'NoRepresentableStep'
            break
        values = [objective(x) for x in points]
        evaluations += len(values)
        index = int(np.argmin(values))
        if values[index] < best_value:
            best_value, best_x = values[index], points[index].copy()
        if count == population:
            generations += 1
        if best_value <= f['target']:
            status = 'TargetReached'
            break
        if count == population:
            es.tell(points, values)
            es.sm.update_now()  # check the current covariance, as Gale does after each refresh
            if es.sm.condition_number > 1e14:
                status = 'ConditionLimit'
                break
        if generations == 10000:
            status = 'GenerationLimit'
    milliseconds = (time.perf_counter_ns() - begin) / 1e6
    assert best_x is not None
    return dict(case_id=f['case_id'], seed=seed, observed=best_value, clean=clean(best_x),
                status=status, evaluations=evaluations, generations=generations,
                sampled=sampled, rejected=rejected, ms=milliseconds, x=best_x.tolist())


def python_baseline(out, warmups, seeds, repeats):
    rows = []
    for f in fixtures():
        for seed in range(10000, 10000 + warmups):
            py_solve(f, seed)
        for seed in range(1, seeds + 1):
            for repeat in range(repeats):
                row = py_solve(f, seed)
                row['repeat'] = repeat
                rows.append(row)
        print(f"Python {f['case_id']}: {seeds} seeds x {repeats} runs", flush=True)
    out.write_text(json.dumps(dict(cma=cma.__version__, numpy=np.__version__, python=sys.version,
        platform=platform.platform(), threadpools=threadpool_info(),
        threads={k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
        settings=dict(active=False, mirrors=0, diagonal_decoding=0, covariance_alleviation=False,
                      eigen_refresh=1, weights='positive log(mu+0.5)-log(rank)', restarts=0), rows=rows), indent=2) + '\n')


def parse_gale(log, out):
    blocks, refs = [], []
    for line in log.read_text().splitlines():
        if line.startswith('case_id\tseed\t'):
            blocks.append([line])
            refs.append([])
        elif line.startswith('REF\t'):
            if not refs:
                refs.append([])
            refs[-1].append(line)
        elif blocks and line.split('\t', 1)[0] in {f['case_id'] for f in fixtures()}:
            blocks[-1].append(line)
    assert len(blocks) == 2, f"expected JVM and Node tables, got {len(blocks)}"
    for platform, block, ref in zip(('jvm', 'node'), blocks, refs, strict=True):
        (out / f'{platform}.tsv').write_text('\n'.join(block) + '\n')
        (out / f'{platform}-reference.tsv').write_text('\n'.join(ref) + '\n')


def reference_check(path):
    assert version("cmaes") == "0.12.0", "the recurrence adapter requires pinned cmaes 0.12.0"
    instances, records = {}, []
    for line in path.read_text().splitlines():
        if not line.startswith('REF\t'):
            continue
        _, seed, generation, positions, fitness, mean, sigma, covariance = line.split('\t')
        seed, generation = int(seed), int(generation)
        points = np.array([list(map(float, point.split(','))) for point in positions.split(';')])
        values = np.fromstring(fitness, sep=',')
        if seed not in instances:
            instances[seed] = CMA(np.array([1., -2., .5]), .3, population_size=6)
            instances[seed]._weights[instances[seed]._mu:] = 0.0  # positive-weight algorithm
            instances[seed]._g = -1  # tell increments before h-sigma uses (g + 1); canonical first exponent is 2
        oracle = instances[seed]
        n = 3
        w = np.log(3.5) - np.log(np.arange(1, 4))
        w /= w.sum()
        eff = 1 / np.dot(w, w)
        c1 = 2 / ((n + 1.3)**2 + eff)
        cmu = min(1-c1, 2*(eff-2+1/eff)/((n+2)**2+eff))
        cs = (eff+2)/(n+eff+5)
        expected_parameters = [eff, (4+eff/n)/(n+4+2*eff/n), c1, cmu, cs,
                               1+2*max(0, math.sqrt((eff-1)/(n+1))-1)+cs]
        np.testing.assert_allclose([oracle._mu_eff, oracle._cc, oracle._c1, oracle._cmu,
                                    oracle._c_sigma, oracle._d_sigma], expected_parameters, rtol=2e-15, atol=2e-16)
        assert cmu < 1-c1-1e-8  # external cap variant is inactive
        assert np.linalg.eigvalsh(oracle._C).min() > 1e-10  # eigenvalue repair is inactive
        assert oracle._sigma < 1e10  # external sigma cap is inactive
        assert generation == oracle.generation + 2
        oracle.tell([(p.copy(), float(v)) for p, v in zip(points, values, strict=True)])
        # Verify the adapter aligns the canonical generation exponent, including h-sigma=false.
        norm = np.linalg.norm(oracle._p_sigma)
        bound = (1.4 + 2 / 4) * oracle._chi_n
        canonical = norm / math.sqrt(1 - (1 - oracle._c_sigma) ** (2 * generation)) < bound
        external = norm / math.sqrt(1 - (1 - oracle._c_sigma) ** (2 * (oracle.generation + 1))) < bound
        assert canonical == external, 'reference h-sigma variants disagree; this fixture is not an identical recurrence'
        expected = dict(mean=oracle.mean.tolist(), sigma=float(oracle._sigma), covariance=oracle._C.ravel().tolist())
        observed = dict(mean=np.fromstring(mean, sep=','), sigma=float(sigma), covariance=np.fromstring(covariance, sep=','))
        errors = {}
        for field in expected:
            np.testing.assert_allclose(observed[field], expected[field], rtol=2e-11, atol=2e-12,
                                       err_msg=f'seed={seed} generation={generation} {field}')
            errors[field] = float(np.max(np.abs(np.asarray(observed[field]) - expected[field])))
        records.append(dict(seed=seed, generation=generation, h_sigma=bool(canonical), errors=errors, expected=expected))
    assert len(records) == 12 and set(instances) == {7, 42}, 'missing or duplicate reference rows'
    return records


def validate_rows(rows, seeds, repeats):
    by_id = {f['case_id']: f for f in fixtures()}
    seen = set()
    for row in rows:
        name = row['case_id']
        assert name in by_id
        f = by_id[name]
        key = (name, int(row['seed']), int(row['repeat']))
        assert key not in seen
        seen.add(key)
        assert key[1] in range(1, seeds + 1) and key[2] in range(repeats)
        x = row['x']
        if isinstance(x, str):
            x = list(map(float, x.split(',')))
        assert len(x) == f['dimension'] and all(math.isfinite(v) for v in x)
        assert math.isfinite(float(row['observed'])) and math.isfinite(float(row['clean']))
        assert math.isfinite(float(row['ms'])) and float(row['ms']) > 0
        assert 0 < int(row['evaluations']) <= f['budget']
        assert 0 <= int(row['generations']) <= 10000
        assert 0 <= int(row['rejected']) <= int(row['sampled'])
        assert row['status'] in {'TargetReached', 'EvaluationLimit', 'GenerationLimit', 'NoRepresentableStep', 'ConditionLimit', 'StepTolerance', 'FitnessStagnation', 'SamplingLimit'}
        if name.startswith('box-'):
            assert all(-1 <= v <= 1 for v in x)
        truth = scalar_objective(f, x)
        assert abs(truth - float(row['clean'])) <= 2e-12 * (1 + abs(truth))
        if name != 'noisy-sphere-8':
            assert abs(truth - float(row['observed'])) <= 2e-12 * (1 + abs(truth))
        if row['status'] == 'TargetReached':
            assert float(row['observed']) <= f['target']
        row['success'] = truth <= (4 + 1e-8 if name == 'box-boundary-4' else 1e-8)
        row['independent_value'] = truth
    assert len(seen) == len(by_id) * seeds * repeats
    # Runtime-specific seeds must reproduce coordinates, receipts, and success across repetitions.
    for name in by_id:
        for seed in range(1, seeds + 1):
            selected = [r for r in rows if r['case_id'] == name and int(r['seed']) == seed]
            for field in ('x', 'evaluations', 'generations', 'observed', 'success', 'status'):
                assert all(r[field] == selected[0][field] for r in selected), f'non-reproducible {name}/{seed}/{field}'
    return rows


def verify(out, seeds, repeats, warmups):
    platforms = {}
    reference = {}
    for platform in ('jvm', 'node'):
        with (out / f'{platform}.tsv').open() as stream:
            platforms[platform] = validate_rows(list(csv.DictReader(stream, delimiter='\t')), seeds, repeats)
        reference[platform] = reference_check(out / f'{platform}-reference.tsv')
    platforms['python'] = validate_rows(json.loads((out / 'python.json').read_text())['rows'], seeds, repeats)
    summary = []
    for f in fixtures():
        name = f['case_id']
        counts, evaluations, times = {}, {}, {}
        for platform, rows in platforms.items():
            selected = [r for r in rows if r['case_id'] == name]
            successful = [r for r in selected if r['success']]
            counts[platform] = sum(int(r['repeat']) == 0 for r in successful)
            evaluations[platform] = statistics.median(int(r['evaluations']) for r in successful) if successful else None
            times[platform] = statistics.median(float(r['ms']) for r in successful) if successful else None
        summary.append(dict(case_id=name, successes=counts, median_success_evaluations=evaluations, median_success_ms=times))
    result = dict(seeds=seeds, repeats=repeats, warmups=warmups, summary=summary, update_reference=reference,
                  reference_adapter=dict(package='cmaes==0.12.0', negative_weights='zeroed', initial_generation=-1,
                      population=6, learning_rate_adaptation=False, eigen_update_period=1,
                      package_source_sha256=hashlib.sha256(Path(inspect.getsourcefile(CMA)).read_bytes()).hexdigest()),
                  source_hashes={name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
                      'core/shared/src/main/scala/gale/optim/CMAES.scala', 'core/shared/src/main/scala/gale/optim/CMAESContracts.scala',
                      'benchmarks/shared/src/main/scala/gale/bench/CMAESBench.scala', 'tools/optim-cmaes/compare.py',
                      'docs/verification/optim-cmaes/fixtures.tsv')})
    (out / 'comparison.json').write_text(json.dumps(result, indent=2) + '\n')
    lines = ['# CMA-ES multi-seed comparison', '', f'{seeds} seeds, {repeats} timed repetitions per seed; {warmups} warmup solves per case.', '',
             'Success is independently checked clean objective error at most 1e-8, including the noisy fixture. Bounds must hold exactly.',
             'Times and evaluation counts below condition on success; differing success sets prevent treating their ratios as paired speedups.', '',
             '| Case | JVM successes | Node successes | Python successes | JVM evals | Node evals | Python evals | JVM ms | Node ms | Python ms |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for row in summary:
        values = [row['case_id']] + [f"{row['successes'][p]}/{seeds}" for p in ('jvm', 'node', 'python')]
        values += [('—' if row[k][p] is None else f"{row[k][p]:.3f}") for k in ('median_success_evaluations', 'median_success_ms') for p in ('jvm', 'node', 'python')]
        lines.append('| ' + ' | '.join(values) + ' |')
    lines += ['', 'All endpoint validity, budget, repeatability and external update checks passed. Failures to reach a target remain in the raw tables.']
    (out / 'comparison.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['python', 'parse', 'verify', 'reference'])
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--log', type=Path)
    parser.add_argument('--warmups', type=int, default=5)
    parser.add_argument('--seeds', type=int, default=20)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.action == 'python':
        python_baseline(args.output / 'python.json', args.warmups, args.seeds, args.repeats)
    elif args.action == 'parse':
        parse_gale(args.log, args.output)
    elif args.action == 'reference':
        result = reference_check(args.log)
        (args.output/'reference.json').write_text(json.dumps(result, indent=2)+'\n')
        print(f'{len(result)} update states match independent cmaes recurrences')
    else:
        verify(args.output, args.seeds, args.repeats, args.warmups)


if __name__ == '__main__':
    main()
