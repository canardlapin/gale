#!/usr/bin/env python3
"""Independently qualify instrumented ALM policy ablations; report work counts, not trace-distorted timings."""
import argparse
import json
import math
from pathlib import Path
import tempfile

import compare


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixtures', type=Path, required=True)
    p.add_argument('--log', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--require-all-converged', action='store_true')
    a = p.parse_args()
    cases = {c['id']: c for c in json.loads(a.fixtures.read_text())}
    with tempfile.TemporaryDirectory(prefix='gale-alm-ablation-') as temp:
        out = Path(temp)
        compare.parse(out, a.log)
        rows = json.loads((out / 'gale.json').read_text())
    variants = ('fixed', 'adaptive', 'scale', 'combined')
    keys = [(r['runtime'], r['case']) for r in rows]
    assert len(keys) == len(set(keys)) and set(keys) == {(v, c) for v in variants for c in cases}
    for r in rows:
        case = cases[r['case']]
        own = compare.independent(case, r['x'], r['multipliers'])
        ref = compare.independent(case, r['x'], case['data']['MULTIPLIERS'])
        assert r['solver_success'] == (r['status'] == 'Converged')
        valid_kkt = (own['feasibility'] <= compare.FEAS * (1 + 1e-5)
            and own['stationarity'] <= 1e-7 * (1 + 1e-5)
            and own['complementarity'] <= compare.COMP * (1 + 1e-5)
            and own['dual_feasibility'] == 0 and own['bound_violation'] == 0)
        r['accurate'] = compare.accurate_endpoint(case, ref)
        if r['status'] == 'Converged': assert valid_kkt and r['accurate']
        if a.require_all_converged or r['runtime'] == 'combined': assert r['status'] == 'Converged' and r['accurate']
        for key, value in r['reported'].items():
            assert math.isfinite(value) and abs(value - own[key]) <= 1e-10 * max(1., abs(own[key]))
        assert abs(r['objective'] - own['objective']) <= 1e-10 * max(1., abs(own['objective']))
        assert r['callbacks'] == r['objective_calls'] + r['constraint_calls'] <= compare.BUDGET
        del r['ms']
        r['independent'] = own
    result = dict(rows=rows, totals={v: {field: sum(r[field] for r in rows if r['runtime'] == v)
        for field in ('objective_calls', 'constraint_calls', 'inner_iterations', 'iterations')} for v in variants})
    result['converged'] = {v: sum(r['solver_success'] for r in rows if r['runtime'] == v) for v in variants}
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))


if __name__ == '__main__': main()
