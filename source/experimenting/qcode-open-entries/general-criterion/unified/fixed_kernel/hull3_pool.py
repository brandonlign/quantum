"""Bounded dynamic process pool, using the exact validated five-H worker."""
import argparse,json,os,sys
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from functools import lru_cache
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.hull3_run import validate,fingerprints,run_one
from mixed_frames.stage3.common import atomic_json,sha

def pool_hashes():return {str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'hull3_pool_verification.json']}

@lru_cache(maxsize=8)
def target(C):
 model=ExtensionModel(17,3,6,8,[255<<8],list(C));assert model.hull==3 and model.C_distance>=6
 return model

def work(task):
 row,oi,prefixes,hashes,ph,folder,seconds=task;N=row['orbits'][oi]['N'];model=target(tuple(row['C']));data=quotient_input(model,N);out=[]
 for p in prefixes:
  path=Path(folder)/f'{row["cls"]}_{oi}_{p["profile_index"]}.json'
  if path.exists():
   old=json.loads(path.read_text());assert old['hashes']==hashes
   if 'pool_hashes' in old:assert old['pool_hashes']==ph
   if old['status']=='EXCLUDED_FIXED_H_N':continue
   if old['status']=='FOUND':raise RuntimeError(f'Inspect candidate {path}')
  r=run_one(model,N,p,data,seconds,hashes,row['cls'],oi)
  if ph is not None:r['pool_hashes']=ph
  atomic_json(path,r);out.append({k:r[k] for k in ('cls','orbit_index','profile_index','status','seconds','z','tally')})
  if r['status']=='FOUND':break
 return out

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--workers',type=int,default=1);parser.add_argument('--max-tasks',type=int);parser.add_argument('--seconds',type=float);args=parser.parse_args();assert __debug__ and 1<=args.workers<=7
 proof=json.loads((HERE/'hull3_pool_verification.json').read_text());assert proof['status']=='PASS'
 for path,digest in proof['hashes'].items():assert sha(ROOT/path)==digest,path
 plan=validate();hashes=fingerprints();ph=pool_hashes();folder=HERE/'hull3_runs';folder.mkdir(exist_ok=True)
 def tasks():
  count=0
  for row in plan['records']:
   for oi,o in enumerate(row['orbits']):
    if o['screen'] is not None:continue
    complete=True
    for p in plan['prefix']:
     path=folder/f'{row["cls"]}_{oi}_{p["profile_index"]}.json'
     if not path.exists():complete=False;break
     r=json.loads(path.read_text());assert r['hashes']==hashes
     if 'pool_hashes' in r:assert r['pool_hashes']==ph
     if r['status']=='FOUND':raise RuntimeError(f'Inspect candidate {path}')
     if r['status']!='EXCLUDED_FIXED_H_N':complete=False;break
    if complete:continue
    yield row,oi,plan['prefix'],hashes,ph,str(folder),args.seconds
    count+=1
    if args.max_tasks is not None and count>=args.max_tasks:return
 iterator=iter(tasks())
 with ProcessPoolExecutor(max_workers=args.workers) as pool:
  pending=set()
  def submit():
   try:task=next(iterator)
   except StopIteration:return False
   pending.add(pool.submit(work,task));return True
  for _ in range(2*args.workers):
   if not submit():break
  while pending:
   done,pending=wait(pending,return_when=FIRST_COMPLETED)
   for f in done:
    rows=f.result()
    for r in rows:print(json.dumps(r),flush=True)
    if any(r['status']=='FOUND' for r in rows):raise RuntimeError('Inspect target candidate; stop census')
    submit()
 atomic_json(HERE/'hull3_finished.json',dict(status='SELECTED_BATCH_FINISHED' if args.max_tasks is not None else 'FINISHED_REQUIRES_STRICT_AGGREGATION',hashes=hashes,pool_hashes=ph))
if __name__=='__main__':main()
