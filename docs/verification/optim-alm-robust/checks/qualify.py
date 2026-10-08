import json, os, subprocess, sys
from pathlib import Path
out=Path('docs/verification/optim-alm-robust')
checks=out/'checks'
env=os.environ.copy()
env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
def run(args, log):
 print('START',log,flush=True)
 with Path(log).open('w') as f:r=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=env)
 print('END',log,'exit',r.returncode,flush=True)
 if r.returncode:raise SystemExit(r.returncode)
run(['sbt','scalafmtAll','scalafmtCheckAll','compileAll','testAll','parityTest','interopBreezeTest','docsCheck','benchCompile',
 'set coreJS / scalaJSStage := FullOptStage',
 'coreJS/testOnly gale.optim.StrongWolfeSuite gale.optim.SmoothOptimizationSuite gale.optim.LBFGSBSuite gale.optim.AugmentedLagrangianSuite'],checks/'gates.log')
run([sys.executable,'tools/optim-alm/run_comparison.py','--python',sys.executable,'--output',str(out),'--warmups','50','--repeats','15'],checks/'comparison-driver.log')
run([sys.executable,'tools/optim-alm/check_verifier.py',str(out)],checks/'verifier-probes.log')
variants=['fixed','adaptive','scale','combined'];fixture=str(out/'fixtures.tsv')
run(['sbt']+[f'benchmarksJVM/runMain gale.bench.OptimizationALMBench 5 1 {fixture} {v} {v} trace' for v in variants],out/'jvm-ablations.log')
commands=['sbt','set benchmarksJS / Compile / mainClass := Some("gale.bench.OptimizationALMBench")','set benchmarksJS / scalaJSStage := FullOptStage']
for v in variants:
 args=['5','1',fixture,v,v,'trace']
 commands += ['set benchmarksJS / Compile / scalaJSMainModuleInitializer := Some(org.scalajs.linker.interface.ModuleInitializer.mainMethodWithArgs("gale.bench.OptimizationALMBench", "main", List('+','.join(json.dumps(x) for x in args)+')))', 'benchmarksJS/run']
run(commands,out/'node-ablations.log')
for rt in ['jvm','node']:
 run([sys.executable,'tools/optim-alm/check_ablations.py','--fixtures',str(out/'fixtures.json'),'--log',str(out/f'{rt}-ablations.log'),'--output',str(out/f'{rt}-ablations.json'),'--require-all-converged'],checks/f'{rt}-ablations-verify.log')
run([sys.executable,'tools/optim-python/run_comparison.py','--python',sys.executable,'--output',str(out/'regressions/first-order')],checks/'first-order-driver.log')
run([sys.executable,'tools/optim-extensions/run_comparison.py','--python',sys.executable,'--output',str(out/'regressions/extensions')],checks/'extensions-driver.log')
