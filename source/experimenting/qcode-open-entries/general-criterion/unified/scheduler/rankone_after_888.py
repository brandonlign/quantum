"""Use the retired 888 core for rank one; keep completed evidence intact."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scheduler.handoff_fast_frames import log,python_manager,stop_tree

def verify_ready():
 from mixed_frames.stage3.common import sha
 r=json.loads((ROOT/'fixed_kernel/READY.json').read_text())
 assert r['status']=='READY_FOR_Q3_HULL3_SWEEP'
 for path,digest in r['hashes'].items():assert sha(ROOT/path)==digest,path
 r=json.loads((ROOT/'rank2_kernel/combined_888.json').read_text())
 assert r['status']=='CONDITIONALLY_COMPLETE_888_COLUMN' and r['classes']==56
 for path,digest in r['hashes'].items():assert sha(ROOT/path)==digest,path

def launch(module,workers,logfile,pidfile):
 assert not python_manager(module),module
 with (ROOT/logfile).open('a') as f:
  p=subprocess.Popen(['python3','-B','-m',module,'--workers',str(workers)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 (ROOT/pidfile).write_text(str(p.pid)+'\n');log(f'started {module} {p.pid}, {workers} workers');return p

def aggregate(module,path):
 subprocess.run(['python3','-B','-m',module],cwd=ROOT,check=True)
 r=json.loads((ROOT/path).read_text());assert r['status']=='CONDITIONALLY_COMPLETE',path

def main():
 verify_ready()
 assert not python_manager('speedup_review.run'),'688 manager still active'
 r=json.loads((ROOT/'speedup_review/aggregate.json').read_text());assert r['status']=='COMPLETE' and r['verified']
 # The old supervisor is retired by the caller without killing its children.
 # This replacement owns all pool managers before any dispatch changes.
 for module in ['mixed_frames.stage4.run','radical_projection.run','radical_projection.pool_run','fixed_kernel.hull3_pool']:
  for pid in python_manager(module):stop_tree(pid)
 a=launch('radical_projection.pool_run',4,'radical_projection/pool.log','radical_projection/pool.pid')
 b=launch('fixed_kernel.hull3_pool',4,'fixed_kernel/hull3_pool.log','fixed_kernel/hull3_pool.pid')
 log('888 column closed conditionally: 4 hull5 + 4 q3 hull3 = 8 search cores')
 while a.poll() is None and b.poll() is None:time.sleep(10)
 if a.poll() is not None:
  if a.returncode!=0:raise RuntimeError(f'hull5 pool exited {a.returncode}; inspect evidence')
  aggregate('radical_projection.aggregate','radical_projection/aggregate.json')
  if b.poll() is None:
   stop_tree(b.pid);b=launch('fixed_kernel.hull3_pool',7,'fixed_kernel/hull3_pool.log','fixed_kernel/hull3_pool.pid')
  if b.wait()!=0:raise RuntimeError('q3 hull3 pool failed; inspect evidence')
  aggregate('fixed_kernel.hull3_aggregate','fixed_kernel/hull3_aggregate.json')
 else:
  if b.returncode!=0:raise RuntimeError(f'q3 hull3 pool exited {b.returncode}; inspect evidence')
  aggregate('fixed_kernel.hull3_aggregate','fixed_kernel/hull3_aggregate.json')
  if a.poll() is None:
   stop_tree(a.pid);a=launch('radical_projection.pool_run',7,'radical_projection/pool.log','radical_projection/pool.pid')
  if a.wait()!=0:raise RuntimeError('hull5 pool failed; inspect evidence')
  aggregate('radical_projection.aggregate','radical_projection/aggregate.json')
 log('Both q3 rank-one branches closed conditionally; residual q1/q2 require separate complete coverage')

if __name__=='__main__':main()
