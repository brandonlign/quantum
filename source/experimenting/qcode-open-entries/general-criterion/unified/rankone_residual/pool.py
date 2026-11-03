"""Bounded dynamic first-label-chunk pool; saved LIMITs remain unresolved."""
import argparse,json,sys
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from rankone_residual.common import validate,fingerprints,run_task
from rankone_residual.search import expected_slices
from mixed_frames.stage3.common import atomic_json,sha

def pool_hashes():return {str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'pool_verification.json']}

def work(task):
 row,s,k,hashes,provenance,folder,*choice=task;backend=choice[0] if choice else 'compact';folder=Path(folder);path=folder/s['names'][k]
 r,_=run_task(row,s['branch'],s['oi'],s['prefix'],s['ci'],k,None,hashes,candidate_dir=folder/'candidates',backend=backend)
 if provenance is not None:r['pool_hashes']=provenance
 assert r['chunks']==s['chunks'] and r['orbit_digest']==s['orbit_digest']
 atomic_json(path,r)
 return dict(task=path.name,status=r['status'],backend=backend,seconds=r['seconds'],nodes=r['z']['nodes'],leaves=r['z']['leaves'])

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=1);ap.add_argument('--classes',type=int,nargs='+');ap.add_argument('--branch',choices=['1,0','2,2','2,0']);ap.add_argument('--profiles',type=int,nargs='+');ap.add_argument('--max-tasks',type=int);ap.add_argument('--backend',choices=['compact','joint-auto'],default='compact');args=ap.parse_args();assert __debug__ and 1<=args.workers<=8
 proof=json.loads((HERE/'pool_verification.json').read_text());assert proof['status']=='PASS'
 for path,digest in proof['hashes'].items():assert sha(ROOT/path)==digest,path
 plan=validate(require_controls=True);hashes=fingerprints();provenance=pool_hashes();folder=HERE/'runs';folder.mkdir(exist_ok=True)
 selected=None if args.classes is None else set(args.classes);branches=None if args.branch is None else {tuple(map(int,args.branch.split(',')))}
 if selected is not None:assert selected<={r['cls'] for r in plan['records']}
 def tasks():
  count=0
  for row in plan['records']:
   if selected is not None and row['cls'] not in selected:continue
   for s in expected_slices(plan,row,branches=branches,profiles=args.profiles):
    for k,name in enumerate(s['names']):
     path=folder/name
     if path.exists():
      old=json.loads(path.read_text());assert old['hashes']==hashes,path
      if 'pool_hashes' in old:assert old['pool_hashes']==provenance,path
      if old['status']=='FOUND':raise RuntimeError(f'Inspect target candidate {path}')
      if old['status']=='EXCLUDED_SLICE_CHUNK':continue
     backend='joint' if args.backend=='joint-auto' and s['prefix']['profile_index'] in (0,8) else 'compact'
     yield row,s,k,hashes,provenance,str(folder),backend;count+=1
     if args.max_tasks is not None and count>=args.max_tasks:return
 iterator=iter(tasks())
 with ProcessPoolExecutor(max_workers=args.workers) as executor:
  pending=set()
  def submit():
   try:task=next(iterator)
   except StopIteration:return False
   pending.add(executor.submit(work,task));return True
  for _ in range(2*args.workers):
   if not submit():break
  while pending:
   done,pending=wait(pending,return_when=FIRST_COMPLETED)
   for future in done:
    result=future.result();print(json.dumps(result),flush=True)
    if result['status']=='FOUND':
     for f in pending:f.cancel()
     raise RuntimeError('Target candidate persisted; stop dispatch and inspect it')
    assert result['status']=='EXCLUDED_SLICE_CHUNK','Unlimited chunk did not finish; preserve unresolved record'
    submit()
 atomic_json(HERE/'pool_finished.json',dict(status='BATCH_FINISHED_REQUIRES_STRICT_AGGREGATION',args=vars(args),hashes=hashes,pool_hashes=provenance))

if __name__=='__main__':main()
