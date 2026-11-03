"""Borrow the live 888 core for a bounded full-H sweep, then resume it."""
import json,os,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from mixed_frames.stage3.common import atomic_json
from rank2_kernel.run import validate,fingerprints

def main():
 validate();hashes=fingerprints();old=int((ROOT/'mixed_frames/stage4/manager.pid').read_text());child=None;paused=False
 def resume():
  nonlocal paused
  if paused:
   try:os.killpg(old,signal.SIGCONT)
   except ProcessLookupError:pass
   paused=False
 def stop(signum=None,frame=None):
  if child is not None and child.poll() is None:
   try:os.killpg(child.pid,signal.SIGTERM)
   except ProcessLookupError:pass
  resume()
  if signum is not None:raise SystemExit(128+signum)
 signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
 status='ERROR';started=time.time()
 try:
  command=subprocess.check_output(['ps','-p',str(old),'-o','command='],text=True)
  assert '-m mixed_frames.stage4.run' in command and os.getpgid(old)==old
  os.killpg(old,signal.SIGSTOP);paused=True
  with (HERE/'run.log').open('a') as log:
   child=subprocess.Popen(['python3','-B','-m','rank2_kernel.run'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
   atomic_json(HERE/'lease.json',dict(status='RUNNING',wrapper=os.getpid(),borrowed_group=old,child=child.pid,deadline_seconds=600,hashes=hashes))
   try:code=child.wait(timeout=600);status='FINISHED' if code==0 else 'ERROR'
   except subprocess.TimeoutExpired:stop();child.wait();status='PARTIAL_TIMEBOX'
 finally:
  stop();atomic_json(HERE/'lease.json',dict(status=status,borrowed_group=old,resumed=True,seconds=time.time()-started,hashes=hashes))
 print(status,flush=True)
if __name__=='__main__':main()
