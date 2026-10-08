from pathlib import Path
import csv, hashlib, json, math, re
out=Path('docs/verification/optim-alm-review')
new=json.loads((out/'comparison.json').read_text())
old=json.loads(Path('docs/verification/optim-alm/comparison.json').read_text())
assert len(new['rows'])==960 and all(r['accurate'] for r in new['rows'])
assert (out/'fixtures.json').read_bytes()==Path('docs/verification/optim-alm/fixtures.json').read_bytes()
checks=[]
for runtime in ['jvm','node']:
 for lane,before,after in [
  ('first-order',Path(f'docs/verification/optim-python/followup/gale-{runtime}-v1.tsv'),out/f'regressions/first-order/gale-{runtime}-v1.tsv'),
  ('extensions',Path(f'docs/verification/optim-extensions/{runtime}.tsv'),out/f'regressions/extensions/{runtime}.tsv')]:
  previous={r['case_id']:r for r in csv.DictReader(before.open(),delimiter='\t')}
  current={r['case_id']:r for r in csv.DictReader(after.open(),delimiter='\t')}
  assert previous.keys()==current.keys()
  for name,a in previous.items():
   b=current[name]
   assert (a['status'],a['callbacks'],a['iterations'])==(b['status'],b['callbacks'],b['iterations']),(runtime,name,a['status'],b['status'],a['callbacks'],b['callbacks'])
   checks.append(dict(runtime=runtime,lane=lane,case=name,status=b['status'],callbacks=int(b['callbacks']),iterations=int(b['iterations']),baseline=str(before)))
(out/'checks/standalone.json').write_text(json.dumps(checks,indent=2)+'\n')
ratios=[]
for rt in ['JVM','Node']:
 for ref in ['SLSQP','trust-constr']:
  r=[x[ref]['median_ms']/x[rt]['median_ms'] for x in new['summary']]
  ratios.append(dict(runtime=rt,reference=ref,wins=sum(v>1 for v in r),geomean=math.exp(sum(map(math.log,r))/len(r))))
oldcalls=sum(r['JVM']['objective_calls'] for r in old['summary'])
newcalls=sum(r['JVM']['objective_calls'] for r in new['summary'])
probe=(out/'checks/comparisons.log').read_text()
allocation=float(re.search(r'bytes_per_zero_solve=([\d.]+)',probe).group(1))
summary=dict(accurate_endpoints=960, ratios=ratios, old_objective_calls=oldcalls,new_objective_calls=newcalls,
 callback_reduction=1-newcalls/oldcalls, allocation_before=169659.96,allocation_after=allocation,
 allocation_reduction=1-allocation/169659.96,standalone_status_and_work_matches=len(checks))
(out/'review-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
for line in probe.splitlines():
 if line.startswith(('ALLOCATION','rho=')):print(line)
