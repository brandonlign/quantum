"""Give rank-one seven cores only after certified fast688 finishes."""
import json,sys,time,subprocess,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scheduler.handoff_fast_frames import log,python_manager,stop_tree

def main():
 log('waiting for full fast688 completion before rank-one fanout')
 while True:
  p=ROOT/'speedup_review/aggregate.json'
  if p.exists():
   r=json.loads(p.read_text())
   if r['status']=='COMPLETE' and r['verified']:break
  time.sleep(30)
 if python_manager('speedup_review.run'):raise RuntimeError('fast688 aggregate complete but manager still live; inspect before allocation')
 from speedup_review.common import validate_manifest
 validate_manifest()
 for pid in python_manager('radical_projection.run'):stop_tree(pid)
 if python_manager('radical_projection.pool_run'):raise RuntimeError('Rank-one pool already active')
 with (ROOT/'radical_projection/pool.log').open('a') as f:
  p=subprocess.Popen(['python3','-B','-m','radical_projection.pool_run','--workers','7'],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 (ROOT/'radical_projection/pool.pid').write_text(str(p.pid)+'\n');log(f'rank-one pool {p.pid}: 7 workers, existing888 keeps1')
 code=p.wait()
 if code!=0:raise RuntimeError(f'Rank-one pool exited {code}')
 subprocess.run(['python3','-B','-m','radical_projection.aggregate'],cwd=ROOT,check=True)
 log('rank-one pool finished; hull-one/three, classification and final proof still require completion')
if __name__=='__main__':main()
