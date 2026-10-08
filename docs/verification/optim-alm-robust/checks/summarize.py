from pathlib import Path
import csv,hashlib,json,math,re
out=Path('docs/verification/optim-alm-robust')
new=json.loads((out/'comparison.json').read_text())
old=json.loads(Path('docs/verification/optim-alm-review/comparison.json').read_text())
assert len(new['rows'])==960 and all(r['accurate'] for r in new['rows'])
assert all(r['status']=='Converged' for r in new['rows'] if r['runtime'] in ['JVM','Node'])
assert (out/'fixtures.json').read_bytes()==Path('docs/verification/optim-alm-review/fixtures.json').read_bytes()
standalone=[]
for rt in ['jvm','node']:
 for lane,prev,curr in [
  ('first-order',Path(f'docs/verification/optim-alm-review/regressions/first-order/gale-{rt}-v1.tsv'),out/f'regressions/first-order/gale-{rt}-v1.tsv'),
  ('extensions',Path(f'docs/verification/optim-alm-review/regressions/extensions/{rt}.tsv'),out/f'regressions/extensions/{rt}.tsv')]:
  a={r['case_id']:r for r in csv.DictReader(prev.open(),delimiter='\t')}
  b={r['case_id']:r for r in csv.DictReader(curr.open(),delimiter='\t')}
  assert a.keys()==b.keys()
  for key in a:
   assert all(a[key][f]==b[key][f] for f in ['status','callbacks','iterations']),(rt,key)
   standalone.append(dict(runtime=rt,case=key,status=b[key]['status'],callbacks=int(b[key]['callbacks']),iterations=int(b[key]['iterations'])))
(out/'checks/standalone.json').write_text(json.dumps(standalone,indent=2)+'\n')
ratio=[]
for rt in ['JVM','Node']:
 for ref in ['SLSQP','trust-constr']:
  vals=[r[ref]['median_ms']/r[rt]['median_ms'] for r in new['summary']]
  ratio.append(dict(runtime=rt,reference=ref,wins=sum(v>1 for v in vals),geomean=math.exp(sum(map(math.log,vals))/len(vals))))
policies={}
for rt in ['jvm','node']:
 a=json.loads((out/f'{rt}-ablations.json').read_text())
 assert all(v==16 for v in a['converged'].values())
 counts={v:0 for v in a['converged']}
 cases={}
 for line in (out/f'{rt}-ablations.log').read_text().splitlines():
  if line.startswith('INNER\t'):
   fields=line.split('\t'); n=int(fields[-1].rsplit(',',1)[1].rstrip(')'))
   counts[fields[2]]+=n
   cases[(fields[2],fields[1])]=cases.get((fields[2],fields[1]),0)+n
 policies[rt]=dict(converged=a['converged'],totals=a['totals'],approximate_steps=counts,
  case_approximate_steps=[dict(policy=k[0],case=k[1],steps=v) for k,v in cases.items() if v])
summary=dict(accurate_endpoints=960,standalone_status_work_matches=len(standalone),ratios=ratio,policies=policies,
 old_default_objective_calls=sum(r['JVM']['objective_calls'] for r in old['summary']),
 default_objective_calls=sum(r['JVM']['objective_calls'] for r in new['summary']))
(out/'robustness-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='policies'},indent=2))
for rt,info in policies.items(): print(rt,info['converged'],info['approximate_steps'])
manifest=json.loads((out/'source-manifest.json').read_text())
for f,h in manifest['source_sha256'].items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h,f
print('All measured source hashes match')
