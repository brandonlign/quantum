"""Eight-core continuation of every open rank-one branch, with strict gates."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scheduler.handoff_fast_frames import log,python_manager,stop_tree
from scheduler.rankone_after_888 import verify_ready
from mixed_frames.stage3.common import sha,atomic_json

JOBS={
 'hull5':dict(module='radical_projection.pool_run',aggregate='radical_projection.aggregate',report='radical_projection/aggregate.json',complete='CONDITIONALLY_COMPLETE',log='radical_projection/pool.log',pid='radical_projection/pool.pid',cap=8),
 'hull3':dict(module='fixed_kernel.hull3_pool',aggregate='fixed_kernel.hull3_aggregate',report='fixed_kernel/hull3_aggregate.json',complete='CONDITIONALLY_COMPLETE',log='fixed_kernel/hull3_pool.log',pid='fixed_kernel/hull3_pool.pid',cap=7),
 'residual':dict(module='rankone_residual.pool',aggregate='rankone_residual.aggregate',report='rankone_residual/aggregate.json',complete='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE',log='rankone_residual/pool.log',pid='rankone_residual/pool.pid',cap=8),
}

def validate():
 verify_ready();r=json.loads((ROOT/'rankone_residual/READY.json').read_text());assert r['status']=='READY_FOR_RESIDUAL_SWEEP'
 for path,digest in r['hashes'].items():assert sha(ROOT/path)==digest,path
 from rankone_residual.common import validate as check
 check(require_controls=True)

def allocation(keys):
 keys=set(keys)
 if 'residual' in keys:
  q3=keys-{'residual'}
  if len(q3)==2:return dict(hull5=4,hull3=3,residual=1)
  if q3:
   key=next(iter(q3));return {key:min(7,JOBS[key]['cap']),'residual':1}
  return dict(residual=8)
 if len(keys)==2:return dict(hull5=4,hull3=4)
 return {key:min(8,JOBS[key]['cap']) for key in keys}

def launch(key,workers):
 job=JOBS[key];assert not python_manager(job['module'])
 with (ROOT/job['log']).open('a') as f:
  args=['python3','-B','-m',job['module'],'--workers',str(workers)]
  if key=='residual':args+=['--backend',json.loads((ROOT/'rankone_residual/READY.json').read_text())['backend_policy']]
  p=subprocess.Popen(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 (ROOT/job['pid']).write_text(str(p.pid)+'\n');log(f'started {key} manager {p.pid}, {workers} workers');return p

def main():
 validate()
 # The caller retires earlier supervisors alone, preserving children until
 # this controller takes ownership. Completed atomic checkpoints are reused.
 for job in JOBS.values():
  for pid in python_manager(job['module']):stop_tree(pid)
 active=set(JOBS);running={};workers={}
 while active:
  want=allocation(active);assert sum(want.values())<=8
  for key in list(running):
   if workers[key]!=want[key]:stop_tree(running[key].pid);del running[key]
  for key,count in want.items():
   if key not in running:running[key]=launch(key,count)
  workers=want
  atomic_json(ROOT/'scheduler/remaining_state.json',dict(status='RUNNING',managers={k:dict(pid=p.pid,workers=workers[k]) for k,p in running.items()}))
  finished=[]
  while not finished:
   time.sleep(10)
   finished=[key for key,p in running.items() if p.poll() is not None]
  for key in finished:
   p=running.pop(key);assert p.returncode==0,f'{key} exited {p.returncode}; inspect candidates/logs, do not infer exclusion'
   job=JOBS[key];subprocess.run(['python3','-B','-m',job['aggregate']],cwd=ROOT,check=True)
   r=json.loads((ROOT/job['report']).read_text());assert r['status']==job['complete'],f'{key} remains unresolved'
   active.remove(key);log(f'{key} strict aggregation complete conditionally')
 atomic_json(ROOT/'scheduler/remaining_state.json',dict(status='ALL_SEARCH_BRANCHES_CONDITIONALLY_CLOSED_REQUIRES_FINAL_PROOF_REVIEW'))
 log('All rank-one search branches conditionally closed. Assemble and independently review the final proof; no automatic novelty claim.')

if __name__=='__main__':main()
