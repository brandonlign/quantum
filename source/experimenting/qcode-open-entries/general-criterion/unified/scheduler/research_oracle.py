"""Eight Mac workers plus two Oracle workers, immutable disjoint classes."""
import json
import os
import subprocess
import sys
import time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from mixed_frames.stage3.common import atomic_json
from rankone_reuse.pool import validate
from oracle_residual.common import assignment, fingerprints
from oracle_residual.bridge import sync
from scheduler.handoff_fast_frames import stop_tree


def main():
    validate()
    fingerprints()
    a=assignment()
    for f in ['radical_projection/aggregate.json','fixed_kernel/hull3_aggregate.json']:
        assert json.loads((ROOT/f).read_text())['status']=='CONDITIONALLY_COMPLETE'
    with (ROOT/'rankone_reuse/pool.log').open('a') as log:
        p=subprocess.Popen(['python3','-B','-m','rankone_reuse.pool','--workers','8',
            '--backend','joint-auto','--classes',*map(str,a['local'])],cwd=ROOT,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
    (ROOT/'rankone_reuse/pool.pid').write_text(str(p.pid)+'\n')
    remote=None
    failures=0
    while True:
        try:
            remote=sync()
            failures=0
        except (subprocess.SubprocessError, OSError, ValueError) as e:
            failures+=1
            print(json.dumps(dict(event='SYNC_RETRY',count=failures,error=str(e))),flush=True)
        local_code=p.poll()
        state=dict(status='RUNNING',managers=dict(residual=dict(pid=p.pid,
            workers=8,returncode=local_code)),oracle=remote,sync_failures=failures)
        atomic_json(ROOT/'scheduler/remaining_state.json',state)
        if (ROOT/'scheduler/STOP_REQUESTED').exists():
            stop_tree(p.pid)
            subprocess.run(['ssh','-o','BatchMode=yes','oracle-vm',
                'touch ~/qcode1736-oracle-v1/project/scheduler/STOP_REQUESTED'],timeout=20,check=True)
            raise RuntimeError('User stop flag: local stopped, remote dispatch stopped')
        if remote is not None and remote['status']=='ERROR':
            stop_tree(p.pid)
            raise RuntimeError('Oracle error/candidate: checkpoints collected; inspect before resuming')
        if local_code is not None and local_code!=0:
            subprocess.run(['ssh','-o','BatchMode=yes','oracle-vm',
                'touch ~/qcode1736-oracle-v1/project/scheduler/STOP_REQUESTED'],timeout=20,check=True)
            raise RuntimeError(f'Local error/candidate {local_code}: remote dispatch stopped')
        if local_code==0 and remote is not None and remote['status']=='BATCH_FINISHED':
            break
        time.sleep(45)
    subprocess.run(['python3','-B','-m','oracle_residual.aggregate'],cwd=ROOT,check=True)
    r=json.loads((ROOT/'oracle_residual/aggregate.json').read_text())
    assert r['status']=='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE'
    atomic_json(ROOT/'scheduler/remaining_state.json',dict(
        status='ALL_SEARCH_BRANCHES_CONDITIONALLY_CLOSED_REQUIRES_FINAL_PROOF_REVIEW'))


if __name__=='__main__':
    main()
