#!/usr/bin/env python3
"""Mutation probes ensure invalid receipts and numerical reference states fail qualification."""
import argparse
import copy
import json
from pathlib import Path
import tempfile

from compare import validate_rows, reference_check

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('evidence', type=Path)
parser.add_argument('--seeds', type=int, default=20)
parser.add_argument('--repeats', type=int, default=3)
args = parser.parse_args()
original = json.loads((args.evidence/'python.json').read_text())['rows']
validate_rows(copy.deepcopy(original), args.seeds, args.repeats)
mutations = [lambda r: r[0].update(clean=float('nan')),
             lambda r: r[0].update(evaluations=15001),
             lambda r: r[0].update(status='Converged'),
             lambda r: r[0].update(x=[999.] * len(r[0]['x'])),
             lambda r: r.append(copy.deepcopy(r[0]))]
for mutate in mutations:
    rows = copy.deepcopy(original)
    mutate(rows)
    try:
        validate_rows(rows, args.seeds, args.repeats)
    except AssertionError:
        continue
    raise RuntimeError('invalid endpoint mutation was accepted')
with tempfile.TemporaryDirectory() as directory:
    lines = (args.evidence/'jvm-reference.tsv').read_text().splitlines()
    parts = lines[0].split('\t')
    parts[6] = str(float(parts[6]) * 1.1)
    lines[0] = '\t'.join(parts)
    path = Path(directory)/'mutated.tsv'
    path.write_text('\n'.join(lines)+'\n')
    try:
        reference_check(path)
    except AssertionError:
        pass
    else:
        raise RuntimeError('mutated step-size update was accepted')
print('Six endpoint/update verifier mutation probes rejected as expected.')
