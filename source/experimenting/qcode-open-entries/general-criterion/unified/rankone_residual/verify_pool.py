"""Two process workers must reproduce sequential complete chunk evidence."""
import json,sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from rankone_residual.common import validate,fingerprints
from rankone_residual.search import expected_slices
from rankone_residual.pool import work
from rankone_residual import verify_aggregate
from mixed_frames.stage3.common import atomic_json,sha

def scrub(x):
 if isinstance(x,dict):return {k:scrub(v) for k,v in x.items() if k!='seconds'}
 if isinstance(x,list):return [scrub(v) for v in x]
 return x

def main():
 verify_aggregate.main();plan=validate(require_controls=True);hashes=fingerprints();row=next(r for r in plan['records'] if r['cls']==5)
 slices=[s for s in expected_slices(plan,row,profiles=[4]) if s['oi'] in (0,1)]
 assert len(slices)==2
 total=nodes=0
 for backend in ['compact','joint']:
  folders=[HERE/f'pool_control_{backend}_sequential',HERE/f'pool_control_{backend}_parallel']
  for folder in folders:folder.mkdir(exist_ok=True)
  sequential=[(row,s,k,hashes,None,str(folders[0]),backend) for s in slices for k in range(s['chunks'])]
  parallel=[(row,s,k,hashes,None,str(folders[1]),backend) for s in slices for k in range(s['chunks'])]
  a=[work(task) for task in sequential]
  with ProcessPoolExecutor(max_workers=2) as pool:b=list(pool.map(work,parallel))
  assert all(r['status']=='EXCLUDED_SLICE_CHUNK' for r in a+b)
  names=sorted(r['task'] for r in a);assert names==sorted(r['task'] for r in b)
  for name in names:assert scrub(json.loads((folders[0]/name).read_text()))==scrub(json.loads((folders[1]/name).read_text()))
  total+=len(names);nodes+=sum(r['nodes'] for r in a)
 report=dict(status='PASS',complete_chunks=total,backends=['compact','joint'],nodes=nodes,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'pool.py',HERE/'controls.json',HERE/'aggregate.py',HERE/'aggregation_verification.json',HERE/'common.py',HERE/'search.py',HERE/'PROOF.md']})
 atomic_json(HERE/'pool_verification.json',report);print(json.dumps(report),flush=True)

if __name__=='__main__':main()
