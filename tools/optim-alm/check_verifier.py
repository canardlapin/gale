#!/usr/bin/env python3
"""Mutation probes ensure the benchmark verifier rejects corrupted evidence."""
import argparse
import contextlib
import copy
import io
import json
from pathlib import Path
import tempfile

import compare


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('--repeats', type=int, default=15)
    a = p.parse_args()
    original = json.loads((a.output / 'gale.json').read_text())
    def nan_point(rows): rows[0]['x'][0] = float('nan')
    def wrong_value(rows): rows[0]['objective'] += 1
    def wrong_diagnostic(rows): rows[0]['reported']['stationarity'] += 1
    def wrong_budget(rows): rows[0]['callbacks'] = compare.BUDGET + 1
    def wrong_status(rows): rows[0]['status'] = 'ProbablyConverged'
    def duplicate(rows): rows.append(copy.deepcopy(rows[0]))
    cases = {c['id']: c for c in json.loads((a.output / 'fixtures.json').read_text())}
    def inconsistent_convergence(rows):
        row = rows[0]
        row['status'] = 'Converged'
        row['solver_success'] = False
        row['multipliers'] = [0.] * len(row['multipliers'])
        metrics = compare.independent(cases[row['case']], row['x'], row['multipliers'])
        assert metrics['stationarity'] > 1e-7
        row['reported'] = {key: metrics[key] for key in row['reported']}
    def inconsistent_failure(rows):
        rows[0]['status'] = 'IterationLimit'
        rows[0]['solver_success'] = True
    probes = [nan_point, wrong_value, wrong_diagnostic, wrong_budget, wrong_status, duplicate,
              inconsistent_convergence, inconsistent_failure]
    with tempfile.TemporaryDirectory(prefix='gale-alm-verifier-') as temp:
        out = Path(temp)
        for name in ('fixtures.json', 'python.json'):
            (out / name).write_bytes((a.output / name).read_bytes())
        for probe in probes:
            rows = copy.deepcopy(original)
            probe(rows)
            (out / 'gale.json').write_text(json.dumps(rows))
            try:
                with contextlib.redirect_stdout(io.StringIO()): compare.verify(out, a.repeats)
            except ValueError:
                print(f'PASS rejected {probe.__name__}')
            else:
                raise AssertionError(f'verifier accepted {probe.__name__}')


if __name__ == '__main__': main()
