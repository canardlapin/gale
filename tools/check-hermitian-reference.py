"""Frozen independent NumPy Hermitian oracle; never import the candidate."""
import json, pathlib, numpy as np, hashlib
assert np.__version__ == '2.3.5'
root=pathlib.Path(__file__).resolve().parents[1]
p=root/'core/shared/src/test/resources'
rng=np.random.Generator(np.random.PCG64(2026101402))
def split(a):return {'real':a.real.ravel().tolist(),'imaginary':a.imag.ravel().tolist()}
cases=[]
for n,r,real in [(2,2,False),(3,2,False),(8,8,False),(16,7,False),(4,4,True)]:
    a=rng.normal(size=(n,r))+1j*(0 if real else rng.normal(size=(n,r)))
    h=a@a.conj().T
    if n==2: h=np.array([[2,1j],[-1j,2]],dtype=np.complex128)
    b=rng.normal(size=(n,2))+1j*rng.normal(size=(n,2))
    w,v=np.linalg.eigh(h);keep=w>max(0,w[-1]*1e-10)
    pinv=(v[:,keep]/w[keep])@v[:,keep].conj().T
    cases.append(dict(n=n,rank=int(keep.sum()),real_only=real,matrix=split(h),rhs=split(b),eigenvalues=w.tolist(),projector=split(v[:,keep]@v[:,keep].conj().T),inverse=split(pinv),solution=split(pinv@b)))
d=dict(schema='hermitian-reference-v1',numpy=np.__version__,seed=2026101402,relative_cutoff=1e-10,tolerances=dict(eigenvalues_abs=1e-10,projector_abs=1e-9,inverse_abs=1e-9,solution_abs=1e-9,residual_certificate=1e-9),cases=cases)
encoded=(json.dumps(d,indent=2,allow_nan=False)+'\n').encode()
target=p/'hermitian-reference-v1.json'
if target.exists(): assert target.read_bytes()==encoded
else:target.write_bytes(encoded)
print(hashlib.sha256(encoded).hexdigest())

# The compiled shared fixtures must contain exactly the independently generated numbers.
import re
source=(root/'core/shared/src/test/scala/gale/complex/HermitianReference.scala').read_text()
for index,c in enumerate(cases):
    body=source.split('private def case'+str(index),1)[1].split('HermitianReferenceCase(',1)[1]
    body=body.split('private def case',1)[0].split('val cases',1)[0]
    numbers=[float(v) for v in re.findall(r'(?<![A-Za-z_])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?',body)]
    expected=[c['n'],c['rank']]+[v for key in ['matrix','rhs'] for part in ['real','imaginary'] for v in c[key][part]]+c['eigenvalues']+[v for key in ['projector','inverse','solution'] for part in ['real','imaginary'] for v in c[key][part]]
    assert numbers==expected,('compiled fixture differs',index)
print('all compiled fixtures match independent NumPy reference')
