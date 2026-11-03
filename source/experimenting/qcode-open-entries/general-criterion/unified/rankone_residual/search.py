"""Exact task enumeration for the residual q1/q2 column (interface for runners).

expected_slices(plan,row) yields EVERY (branch, N-orbit, distinguished H,
identification-orbit representative) of a suffix class, with its first-label
chunk count. Coverage rationale is in PROOF.md. Nothing here excludes alone.
"""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from fixed_kernel.symmetry import suffix_slice
from rankone_residual.common import CHUNK,slice_setup,task_name,digest,run_task,fingerprints,validate
from mixed_frames.stage3.common import atomic_json

def expected_slices(plan,row,*,branches=None,profiles=None):
 for branch in row['branches']:
  q,t=branch['q'],branch['t']
  if branches is not None and (q,t) not in branches:continue
  for oi,orbit in enumerate(branch['N_orbits']):
   for prefix in [p for p in plan['prefix'] if (p['q'],p['t'])==(q,t)]:
    if profiles is not None and prefix['profile_index'] not in profiles:continue
    model,data,H,order,contexts,io=slice_setup(tuple(row['C']),tuple(orbit['N']),tuple(prefix['H']),prefix['p'],q,t)
    for o in io['orbits']:
     ci=o['representative'];ss=suffix_slice(model,orbit['N'],data,contexts[ci])
     n=max(1,-(-len(ss['first_representatives'])//CHUNK))
     yield dict(branch=branch,q=q,t=t,oi=oi,prefix=prefix,ci=ci,members=o['members'],chunks=n,orbit_digest=digest(ss),
                names=[task_name(row['cls'],q,t,oi,prefix['profile_index'],ci,k) for k in range(n)])

def run_slice(plan,row,s,seconds,hashes,folder,*,rewrite=False,log=print,backend='compact'):
 """Run every chunk of one slice; returns list of statuses. Stops on FOUND."""
 out=[]
 for k,name in enumerate(s['names']):
  path=Path(folder)/name
  if path.exists() and not rewrite:
   old=json.loads(path.read_text())
   if old['hashes']!=hashes:raise RuntimeError(f'semantic sources changed since {path}')
   if old['status']=='FOUND':raise RuntimeError(f'Inspect target candidate {path}')
   if old['status']=='EXCLUDED_SLICE_CHUNK':out.append(old['status']);continue
  r,ss=run_task(row,s['branch'],s['oi'],s['prefix'],s['ci'],k,seconds,hashes,backend=backend)
  assert r['orbit_digest']==s['orbit_digest'] and r['chunks']==s['chunks']
  atomic_json(path,r);out.append(r['status'])
  log(json.dumps(dict(task=name,status=r['status'],nodes=r['z'].get('nodes'),leaves=r['z']['leaves'],seconds=round(r['seconds'],3))))
  if r['status']=='FOUND':raise RuntimeError(f'Inspect target candidate {path}')
 return out
