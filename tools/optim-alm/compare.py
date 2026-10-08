#!/usr/bin/env python3
"""Reproducible constrained fixtures, SciPy timings, and independent original-problem verification."""
import argparse
import hashlib
import json
import math
import os
import platform
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FEAS = 1e-8
STAT = 1e-6
COMP = 1e-8
POINT = 1e-5
BUDGET = 50000


def fixtures(out):
    import numpy as np
    families = []
    for n, m, condition, dense in [(8, 1, 10, False), (32, 4, 100, True), (128, 8, 10, True)]:
        diagonal = np.geomspace(1, condition, n)
        v = np.arange(1, n + 1, dtype=float)
        h = np.eye(n) - 2 * np.outer(v, v) / (v @ v)
        q = (h * diagonal) @ h if dense else np.diag(diagonal)
        a = np.array([[math.cos(math.pi * (j + .5) * (i + 1) / n) / math.sqrt(n / 2)
                       for j in range(n)] for i in range(m)])
        target = np.array([.3 * math.sin(j + .2) for j in range(n)])
        multipliers = np.linspace(.5, 1.5, m)
        center = target + np.linalg.solve(q, a.T @ multipliers)
        families.append((f'equality-{"dense" if dense else "diagonal"}-{n}', dict(
            KIND=[0], EQUALITIES=[m], INEQUALITIES=[0], Q=q.ravel().tolist(), A=a.ravel().tolist(),
            B=(a @ target).tolist(), CENTER=center.tolist(), TARGET=target.tolist(), MULTIPLIERS=multipliers.tolist())))
    families.append(('mixed-box-3', dict(KIND=[0], EQUALITIES=[1], INEQUALITIES=[2], Q=np.eye(3).ravel().tolist(),
        CENTER=[3., 0., -1.], A=[1., 1., 0., 0., -1., 0., 0., 0., 1.], B=[1., 0., .5],
        LOWER=[-10., -10., 0.], UPPER=[10., 10., 10.], TARGET=[1., 0., 0.], MULTIPLIERS=[2., 2., 0.])))
    for n in (8, 64):
        center = np.linspace(.5, 1.5, n)
        length = np.linalg.norm(center)
        families.append((f'ball-{n}', dict(KIND=[1], EQUALITIES=[0], INEQUALITIES=[1], Q=np.eye(n).ravel().tolist(),
            CENTER=center.tolist(), TARGET=(center / length).tolist(), MULTIPLIERS=[(length - 1) / 2])))
    lo, hi = 1., 2.
    for _ in range(80):
        mid = (lo + hi) / 2
        if 2 * mid ** 3 - mid - 2 > 0: hi = mid
        else: lo = mid
    root = (lo + hi) / 2
    for n in (2, 16):
        families.append((f'parabola-{n}', dict(KIND=[2], EQUALITIES=[n // 2], INEQUALITIES=[0],
            Q=np.eye(n).ravel().tolist(), CENTER=[2., 1.] * (n // 2),
            TARGET=[root, root * root] * (n // 2), MULTIPLIERS=[1 - root * root] * (n // 2))))
    cases = []
    for name, data in families:
        n = len(data['CENTER'])
        for start in range(2):
            initial = [0.] * n if start == 0 else [.3 * math.sin(i + .4) for i in range(n)]
            if 'LOWER' in data:
                initial = [max(l, min(u, x)) for x, l, u in zip(initial, data['LOWER'], data['UPPER'], strict=True)]
            cases.append(dict(id=f'{name}-s{start}', data=dict(data, X0=initial)))
    out.mkdir(parents=True, exist_ok=True)
    (out / 'fixtures.json').write_text(json.dumps(cases, indent=2, allow_nan=False) + '\n')
    lines = []
    for c in cases:
        lines.append('CASE\t' + c['id'])
        lines += [key + '\t' + '\t'.join(format(x, '.17g') for x in values) for key, values in c['data'].items()]
        lines.append('END')
    (out / 'fixtures.tsv').write_text('\n'.join(lines) + '\n')


def prepared(case):
    import numpy as np
    d = {k: np.asarray(v, dtype=float) for k, v in case['data'].items()}
    n = len(d['X0'])
    q = d['Q'].reshape(n, n)
    eq, ineq = int(d['EQUALITIES'][0]), int(d['INEQUALITIES'][0])
    a = d['A'].reshape(eq + ineq, n) if 'A' in d else None
    def objective(x):
        delta = x - d['CENTER']
        g = q @ delta
        return .5 * (delta @ g), g
    def constraints(x):
        if d['KIND'][0] == 0: return a @ x - d['B'], a
        if d['KIND'][0] == 1: return np.array([x @ x - 1]), (2 * x)[None, :]
        j = np.zeros((eq, n))
        j[np.arange(eq), 2 * np.arange(eq)] = -2 * x[::2]
        j[np.arange(eq), 2 * np.arange(eq) + 1] = 1
        return x[1::2] - x[::2] ** 2, j
    return d, eq, ineq, objective, constraints


def benchmark(out, warmups, repeats):
    import numpy as np
    import scipy
    from scipy.optimize import Bounds, NonlinearConstraint, minimize
    from threadpoolctl import threadpool_info
    cases = json.loads((out / 'fixtures.json').read_text())
    rows = []
    calibration = []
    for case in cases:
        d, eq, ineq, f, c = prepared(case)
        bounds = Bounds(d['LOWER'], d['UPPER']) if 'LOWER' in d else None
        for method in ('SLSQP', 'trust-constr'):
            def run(options):
                # Cache both values and Jacobian together, just like Gale's fused constraint callback.
                calls = dict(objective=0, constraints=0)
                last_x = None
                last_c = None
                def charge(kind):
                    if sum(calls.values()) >= BUDGET: raise RuntimeError('global callback budget exceeded')
                    calls[kind] += 1
                def fun(x):
                    charge('objective')
                    return f(x)
                def con(x):
                    nonlocal last_x, last_c
                    if last_x is None or not np.array_equal(last_x, x):
                        charge('constraints')
                        last_c = c(x)
                        last_x = x.copy()
                    return last_c
                if method == 'SLSQP':
                    restrictions = []
                    if eq: restrictions.append(dict(type='eq', fun=lambda x: con(x)[0][:eq], jac=lambda x: con(x)[1][:eq]))
                    if ineq: restrictions.append(dict(type='ineq', fun=lambda x: -con(x)[0][eq:], jac=lambda x: -con(x)[1][eq:]))
                else:
                    lb = np.r_[np.zeros(eq), np.full(ineq, -np.inf)]
                    restrictions = NonlinearConstraint(lambda x: con(x)[0], lb, np.zeros(eq + ineq), jac=lambda x: con(x)[1])
                # Solver wrapper/cache/constraint construction is included; prebuilt data and base callbacks are not.
                result = minimize(fun, d['X0'], jac=True, method=method, bounds=bounds, constraints=restrictions, options=options)
                return result, calls
            # Give SciPy the fastest accurate setting from a declared sweep, rather than requiring arbitrary extra digits.
            candidates = ([dict(ftol=tol, maxiter=1000) for tol in (1e-6, 1e-8, 1e-10, 1e-12, 1e-14)] if method == 'SLSQP' else
                [dict(gtol=tol, xtol=1e-14, barrier_tol=1e-12, maxiter=1000,
                      initial_barrier_parameter=barrier, initial_barrier_tolerance=barrier)
                 for barrier in (.1, 1e-9) for tol in (1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-16)])
            qualified = []
            for options in candidates:
                run(options)  # one untimed calibration warmup
                trials, accurate = [], True
                for _ in range(3):
                    before = time.perf_counter_ns()
                    r, _ = run(options)
                    trials.append((time.perf_counter_ns() - before) / 1e6)
                    metrics = independent(case, r.x.tolist(), case['data']['MULTIPLIERS'])
                    accurate &= accurate_endpoint(case, metrics)
                calibration.append(dict(case=case['id'], runtime=method, options=options, accurate=bool(accurate),
                    median_ms=statistics.median(trials), metrics=metrics))
                if accurate: qualified.append((statistics.median(trials), options))
            options = min(qualified, key=lambda x: x[0])[1] if qualified else candidates[-1]
            for _ in range(warmups): run(options)
            for repeat in range(repeats):
                before = time.perf_counter_ns()
                r, counts = run(options)
                elapsed = (time.perf_counter_ns() - before) / 1e6
                rows.append(dict(runtime=method, case=case['id'], repeat=repeat, ms=elapsed, status=str(r.message),
                    solver_success=bool(r.success), options=options, objective=float(r.fun), x=r.x.tolist(), iterations=int(r.nit),
                    objective_calls=counts['objective'], constraint_calls=counts['constraints'], callbacks=sum(counts.values())))
    (out / 'python.json').write_text(json.dumps(dict(rows=rows, calibration=calibration, environment=dict(
        python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform(),
        machine=platform.machine(), processor=platform.processor(), threadpools=threadpool_info(),
        thread_env={k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
        warmups=warmups, repeats=repeats)), indent=2, allow_nan=False) + '\n')


def parse(out, log):
    rows = []
    for line in log.read_text().splitlines():
        if not line.startswith('ALM\t'): continue
        c = line.split('\t')
        rows.append(dict(runtime=c[1], case=c[2], repeat=int(c[3]), ms=float(c[4]), status=c[5],
            solver_success=c[5] == 'Converged', objective=float(c[6]), iterations=int(c[7]), inner_iterations=int(c[8]),
            callbacks=int(c[9]), objective_calls=int(c[10]), constraint_calls=int(c[11]),
            reported=dict(feasibility=float(c[12]), scaled_feasibility=float(c[13]), stationarity=float(c[14]), complementarity=float(c[15])),
            penalty=float(c[16]), x=list(map(float, c[17].split(','))), multipliers=list(map(float, c[18].split(',')))))
    if not rows: raise ValueError('no Gale benchmark rows')
    (out / 'gale.json').write_text(json.dumps(rows, indent=2, allow_nan=False) + '\n')


def independent(case, x, multipliers):
    """Scalar math.fsum implementation, independent of timed NumPy/Gale kernels and solver diagnostics."""
    d = case['data']
    n, eq = len(d['X0']), int(d['EQUALITIES'][0])
    m = eq + int(d['INEQUALITIES'][0])
    if len(x) != n or len(multipliers) != m or not all(math.isfinite(v) for v in x + multipliers):
        raise ValueError('invalid endpoint or multipliers')
    delta = [x[i] - d['CENTER'][i] for i in range(n)]
    gradient = [math.fsum(d['Q'][i * n + j] * delta[j] for j in range(n)) for i in range(n)]
    objective = .5 * math.fsum(delta[i] * gradient[i] for i in range(n))
    if d['KIND'][0] == 0:
        jac = [d['A'][i * n:(i + 1) * n] for i in range(m)]
        values = [math.fsum(jac[i][j] * x[j] for j in range(n)) - d['B'][i] for i in range(m)]
    elif d['KIND'][0] == 1:
        values, jac = [math.fsum(v * v for v in x) - 1], [[2 * v for v in x]]
    else:
        values = [x[2 * i + 1] - x[2 * i] ** 2 for i in range(eq)]
        jac = [[-2 * x[j] if j == 2 * i else 1. if j == 2 * i + 1 else 0. for j in range(n)] for i in range(eq)]
    lagrangian = [math.fsum([gradient[j]] + [multipliers[i] * jac[i][j] for i in range(m)]) for j in range(n)]
    lower = d.get('LOWER', [-math.inf] * n)
    upper = d.get('UPPER', [math.inf] * n)
    projected = [min(g, z - l) if g > 0 else max(g, z - u) for z, g, l, u in zip(x, lagrangian, lower, upper, strict=True)]
    bound_violation = max([0.] + [max(l - z, z - u, 0.) for z, l, u in zip(x, lower, upper, strict=True)])
    feasibility = max([bound_violation] + [abs(v) if i < eq else max(v, 0.) for i, v in enumerate(values)])
    complementarity = max([0.] + [abs(multipliers[i] * values[i]) for i in range(eq, m)])
    return dict(objective=objective, feasibility=feasibility, scaled_feasibility=feasibility,
        stationarity=max(map(abs, projected)), complementarity=complementarity,
        dual_feasibility=max([0.] + [-v for v in multipliers[eq:]]),
        point_error=max(abs(a - b) for a, b in zip(x, d['TARGET'], strict=True)), bound_violation=bound_violation)


def accurate_endpoint(case, metrics):
    target_value = independent(case, case['data']['TARGET'], case['data']['MULTIPLIERS'])['objective']
    metrics['objective_error'] = abs(metrics['objective'] - target_value)
    return (metrics['feasibility'] <= FEAS and metrics['stationarity'] <= STAT and metrics['complementarity'] <= COMP
        and metrics['point_error'] <= POINT and metrics['objective_error'] <= 1e-7 * max(1., abs(target_value))
        and metrics['bound_violation'] == 0. and metrics['dual_feasibility'] == 0.)


def verify(out, repeats):
    cases = {c['id']: c for c in json.loads((out / 'fixtures.json').read_text())}
    rows = json.loads((out / 'gale.json').read_text()) + json.loads((out / 'python.json').read_text())['rows']
    keys = [(r['runtime'], r['case'], r['repeat']) for r in rows]
    expected = {(runtime, case, repeat) for runtime in ('JVM', 'Node', 'SLSQP', 'trust-constr') for case in cases for repeat in range(repeats)}
    if len(keys) != len(set(keys)) or set(keys) != expected: raise ValueError('missing, duplicate or unexpected benchmark rows')
    for row in rows:
        c = cases[row['case']]
        metrics = independent(c, row['x'], c['data']['MULTIPLIERS'])
        if not math.isfinite(row['ms']) or row['ms'] <= 0: raise ValueError('invalid timing')
        if not math.isfinite(row['objective']) or abs(row['objective'] - metrics['objective']) > 1e-10 * max(1., abs(metrics['objective'])):
            raise ValueError('incorrect reported objective')
        if row['callbacks'] != row['objective_calls'] + row['constraint_calls'] or not 0 < row['callbacks'] <= BUDGET:
            raise ValueError('incorrect work counts')
        if row['runtime'] in ('JVM', 'Node'):
            if row['status'] not in ('Converged', 'IterationLimit', 'InnerIterationLimit', 'LineSearchFailed', 'EvaluationLimit', 'Cancelled', 'PenaltyLimit', 'NumericalStagnation'):
                raise ValueError('unknown Gale termination')
            if not isinstance(row['solver_success'], bool) or row['solver_success'] != (row['status'] == 'Converged'):
                raise ValueError('inconsistent Gale convergence fields')
            own = independent(c, row['x'], row['multipliers'])
            for name, reported in row['reported'].items():
                if not math.isfinite(reported) or abs(reported - own[name]) > 1e-10 * max(1., abs(own[name])):
                    raise ValueError(f'incorrect {name} diagnostic: {row["case"]}')
            if row['status'] == 'Converged' and not (own['feasibility'] <= FEAS * (1 + 1e-5) and own['stationarity'] <= 1e-7 * (1 + 1e-5)
                and own['complementarity'] <= COMP * (1 + 1e-5) and own['dual_feasibility'] == 0):
                raise ValueError('false Gale convergence')
        row['metrics'] = metrics
        row['accurate'] = accurate_endpoint(c, metrics)
    summary = []
    for case in cases:
        record = dict(case=case)
        for runtime in ('JVM', 'Node', 'SLSQP', 'trust-constr'):
            subset = [r for r in rows if r['case'] == case and r['runtime'] == runtime]
            accurate = [r for r in subset if r['accurate']]
            record[runtime] = dict(accurate=len(accurate), runs=len(subset), solver_success=sum(r['solver_success'] for r in subset),
                median_ms=statistics.median(r['ms'] for r in subset), min_ms=min(r['ms'] for r in subset), max_ms=max(r['ms'] for r in subset),
                objective_calls=statistics.median(r['objective_calls'] for r in subset), constraint_calls=statistics.median(r['constraint_calls'] for r in subset),
                feasibility=max(r['metrics']['feasibility'] for r in subset), stationarity=max(r['metrics']['stationarity'] for r in subset),
                complementarity=max(r['metrics']['complementarity'] for r in subset), statuses=sorted(set(r['status'] for r in subset)))
        summary.append(record)
    (out / 'comparison.json').write_text(json.dumps(dict(summary=summary, rows=rows,
        thresholds=dict(feasibility=FEAS, stationarity=STAT, complementarity=COMP, point_error=POINT, objective_relative=1e-7)), indent=2, allow_nan=False) + '\n')
    lines = ['# Augmented Lagrangian comparison', '', 'Median full-solve milliseconds. **FAIL** means at least one endpoint failed the common accuracy gate; its time is not a qualified speed result.', '',
             '| Case | JVM | Node | SciPy SLSQP | SciPy trust-constr |', '|---|---:|---:|---:|---:|']
    for r in summary:
        cells = []
        for runtime in ('JVM', 'Node', 'SLSQP', 'trust-constr'):
            v = r[runtime]
            cells.append(f'{v["median_ms"]:.3f}' + ('' if v['accurate'] == v['runs'] else f' **FAIL {v["accurate"]}/{v["runs"]}**'))
        lines.append('| ' + r['case'] + ' | ' + ' | '.join(cells) + ' |')
    (out / 'comparison.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['fixtures', 'python', 'parse', 'verify'])
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--log', type=Path)
    p.add_argument('--warmups', type=int, default=20)
    p.add_argument('--repeats', type=int, default=11)
    a = p.parse_args()
    if a.mode == 'fixtures': fixtures(a.output)
    elif a.mode == 'python': benchmark(a.output, a.warmups, a.repeats)
    elif a.mode == 'parse': parse(a.output, a.log)
    else: verify(a.output, a.repeats)


if __name__ == '__main__': main()
