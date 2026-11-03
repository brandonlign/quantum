"""Use seven rank-one cores after certified 688 completion, preserving 888."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scheduler.handoff_fast_frames import log,python_manager,stop_tree

def launch(module,workers,logfile,pidfile):
 with (ROOT/logfile).open('a') as f:
  p=subprocess.Popen(['python3','-B','-m',module,'--workers',str(workers)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 (ROOT/pidfile).write_text(str(p.pid)+'\n');log(f'started {module} {p.pid}, {workers} workers');return p

def aggregate(module,path):
 subprocess.run(['python3','-B','-m',module],cwd=ROOT,check=True)
 r=json.loads((ROOT/path).read_text());assert r['status']=='CONDITIONALLY_COMPLETE',path

def ready():
 p=ROOT/'fixed_kernel/READY.json'
 if not p.exists():return False
 from mixed_frames.stage3.common import sha
 r=json.loads(p.read_text());assert r['status']=='READY_FOR_Q3_HULL3_SWEEP'
 for path,digest in r['hashes'].items():assert sha(ROOT/path)==digest,path
 return True

def main():
 log('waiting for strict fast688 completion and fixed-kernel validation; current budget remains 6+1+1')
 while True:
  r=json.loads((ROOT/'speedup_review/aggregate.json').read_text())
  if r['status']=='COMPLETE' and r['verified'] and ready():break
  time.sleep(15)
 while python_manager('speedup_review.run'):time.sleep(2)
 from speedup_review.common import validate_manifest
 from fixed_kernel.hull3_run import validate
 validate_manifest();validate()
 for module in ['radical_projection.pool_run','fixed_kernel.hull3_pool','fixed_kernel.hull3_run']:
  if python_manager(module):raise RuntimeError(f'Unexpected existing manager {module}')
 for pid in python_manager('radical_projection.run'):stop_tree(pid)
 a=launch('radical_projection.pool_run',3,'radical_projection/pool.log','radical_projection/pool.pid')
 b=launch('fixed_kernel.hull3_pool',4,'fixed_kernel/hull3_pool.log','fixed_kernel/hull3_pool.pid')
 log('3 hull5 + 4 q3 hull3 + 1 existing888 = 8 cores')
 while a.poll() is None and b.poll() is None:time.sleep(10)
 if a.poll() is not None:
  if a.returncode!=0:raise RuntimeError(f'hull5 pool exited {a.returncode}; inspect log and candidates')
  aggregate('radical_projection.aggregate','radical_projection/aggregate.json')
  if b.poll() is None:
   stop_tree(b.pid);b=launch('fixed_kernel.hull3_pool',7,'fixed_kernel/hull3_pool.log','fixed_kernel/hull3_pool.pid')
  if b.wait()!=0:raise RuntimeError('q3 hull3 pool failed; inspect log and candidates')
  aggregate('fixed_kernel.hull3_aggregate','fixed_kernel/hull3_aggregate.json')
 else:
  if b.returncode!=0:raise RuntimeError(f'q3 hull3 pool exited {b.returncode}; inspect log and candidates')
  aggregate('fixed_kernel.hull3_aggregate','fixed_kernel/hull3_aggregate.json')
  if a.poll() is None:
   stop_tree(a.pid);a=launch('radical_projection.pool_run',7,'radical_projection/pool.log','radical_projection/pool.pid')
  if a.wait()!=0:raise RuntimeError('hull5 pool failed; inspect log and candidates')
  aggregate('radical_projection.aggregate','radical_projection/aggregate.json')
 log('Both rank-one branches strictly closed conditionally; q1/q2 and 888 still need their own complete work')
if __name__=='__main__':main()
