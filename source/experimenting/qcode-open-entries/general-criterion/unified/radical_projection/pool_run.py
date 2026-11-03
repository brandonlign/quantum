"""Dynamic exact-N workers using the unchanged certified search implementation."""
import argparse,json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from radical_projection.run import (target,search,fingerprints,quotient_input,independent_quotient,syndrome_orbits,sha,atomic_json)

def pool_hashes():return {str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}

def task(item):
 row,i,expected,provenance=item;N=row['orbits'][i]['representative'];model=target(tuple(row['C']));start=time.monotonic();data=quotient_input(model,N);check=independent_quotient(model,N)
 assert all(data[k]==check[k] for k in ('units','labels','distances','Gamma'))
 orbits=syndrome_orbits(model,data,N);result=search(model,N,data,orbits['first_representatives'],None)
 result.update(cls=row['cls'],orbit_index=i,N=N,orbit_data=orbits,hashes=expected,pool_hashes=provenance,seconds=time.monotonic()-start,scope='One exact-N orbit in hull-five exclusive rank-one column, conditional on suffix classification')
 path=HERE/'runs'/f'{row["cls"]}_{i}.json';atomic_json(path,result);return {k:result[k] for k in ('cls','orbit_index','status','seconds','z','tally')}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=1);parser.add_argument('--max-tasks',type=int);args=parser.parse_args();assert 1<=args.workers<=8 and __debug__
 expected=fingerprints();provenance=pool_hashes();verified=json.loads((HERE/'verification.json').read_text());assert verified['status']=='PASS'
 for p,h in verified['hashes'].items():assert sha(ROOT/p)==h
 plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text());screen=json.loads((ROOT/'rankone_recenter/column_screen.json').read_text());lookup={r['cls']:r for r in screen['records']};tasks=[]
 assert plan['source_sha256']==sha(ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json') and plan['preparation_sha256']==sha(ROOT/'rankone_recenter/prepare.py')
 for row in plan['records']:
  for i,o in enumerate(row['orbits']):
   if lookup[row['cls']]['orbits'][i]['reason'] is not None:continue
   path=HERE/'runs'/f'{row["cls"]}_{i}.json'
   if path.exists():
    old=json.loads(path.read_text());assert old['hashes']==expected
    if 'pool_hashes' in old:assert old['pool_hashes']==provenance
    if old['status']=='EXCLUDED_FIXED_N':continue
    if old['status']=='FOUND':raise RuntimeError('Inspect saved target candidate')
   tasks.append((row,i,expected,provenance))
 if args.max_tasks is not None:tasks=tasks[:args.max_tasks]
 print(json.dumps(dict(tasks=len(tasks),workers=args.workers)),flush=True)
 # Submit a bounded window; only this manager owns the selected tasks.
 with ProcessPoolExecutor(max_workers=args.workers) as pool:
  iterator=iter(tasks);pending={}
  for _ in range(args.workers):
   item=next(iterator,None)
   if item is not None:pending[pool.submit(task,item)]=item
  while pending:
   future=next(as_completed(pending));pending.pop(future);result=future.result();print(json.dumps(result),flush=True)
   if result['status']=='FOUND':
    for f in pending:f.cancel()
    raise RuntimeError('A physical target candidate was found; stop dispatch and inspect it')
   item=next(iterator,None)
   if item is not None:pending[pool.submit(task,item)]=item
 atomic_json(HERE/'pool_finished.json',dict(status='FINISHED_REQUIRES_AGGREGATION',hashes=expected,pool_hashes=provenance,selected_tasks=len(tasks)))
if __name__=='__main__':main()
