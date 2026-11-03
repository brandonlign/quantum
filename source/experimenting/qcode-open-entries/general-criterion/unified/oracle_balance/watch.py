"""Avoid a static-shard bottleneck without changing either admitted producer."""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from oracle_residual.common import assignment, fingerprints
from oracle_residual.bridge import sync
from scheduler.handoff_fast_frames import stop_tree
from mixed_frames.stage3.common import atomic_json, sha


def launch(classes):
    with (ROOT/'rankone_reuse/pool.log').open('a') as log:
        p=subprocess.Popen(['python3','-B','-m','rankone_reuse.pool','--workers','8',
            '--backend','joint-auto','--classes',*map(str,classes)],cwd=ROOT,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
            env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
    (ROOT/'rankone_reuse/pool.pid').write_text(str(p.pid)+'\n')
    return p


def stop_vm():
    subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10','oracle-vm',
        'cd ~/qcode1736-oracle-v1/project && ../venv/bin/python -B -m oracle_balance.stop_remote'],
        check=True,timeout=25)


def frontier(folder,selected,op):
    values=[]
    for path in folder.glob('*.json'):
        try:cls=int(path.name.split('_')[0])
        except ValueError:continue
        if cls in selected:values.append(cls)
    return op(values) if values else None


def main():
    ready=json.loads((ROOT/'oracle_balance/READY.json').read_text())
    assert ready['status']=='READY_OPTIONAL_ORACLE_BALANCING'
    for path,digest in ready['hashes'].items():assert sha(ROOT/path)==digest,path
    fingerprints()
    a=assignment()
    (ROOT/'oracle_balance/watch.pid').write_text(str(os.getpid())+'\n')
    with open(os.devnull,'w') as out:
        subprocess.Popen(['caffeinate','-i','-w',str(os.getpid())],stdout=out,stderr=out)
    # Existing supervisor retains first-phase ownership and candidate handling.
    while True:
        state=json.loads((ROOT/'scheduler/remaining_state.json').read_text())
        if state['status']=='ALL_SEARCH_BRANCHES_CONDITIONALLY_CLOSED_REQUIRES_FINAL_PROOF_REVIEW':
            return
        if (ROOT/'scheduler/STOP_REQUESTED').exists():
            return
        local=state.get('managers',{}).get('residual',{})
        remote=state.get('oracle') or {}
        if local.get('returncode')==0 and remote.get('status')=='RUNNING':
            break
        time.sleep(45)
    owner=int((ROOT/'scheduler/research_oracle.pid').read_text())
    command=subprocess.check_output(['ps','-p',str(owner),'-o','command='],text=True)
    assert 'scheduler/research_oracle.py' in command
    os.kill(owner,signal.SIGTERM)
    print('First local partition finished; Mac assists Oracle from the opposite end.',flush=True)
    selected=set(a['remote'])
    p=launch(a['remote'])
    while True:
        remote=sync()
        code=p.poll()
        state=dict(status='BALANCING_ORACLE_SHARD',managers=dict(residual=dict(
            pid=p.pid,workers=8,returncode=code)),oracle=remote)
        atomic_json(ROOT/'scheduler/remaining_state.json',state)
        if (ROOT/'scheduler/STOP_REQUESTED').exists():
            stop_tree(p.pid);stop_vm();return
        if remote['status']=='ERROR' or (code is not None and code!=0):
            stop_tree(p.pid);stop_vm()
            raise RuntimeError('Candidate/error during balancing: preserve and inspect')
        hi=frontier(ROOT/'rankone_reuse/runs',selected,max)
        lo=frontier(ROOT/'oracle_residual/runs',selected,min)
        if code==0 or remote['status']=='BATCH_FINISHED' or (hi is not None and lo is not None and hi>=lo):
            if code is None:stop_tree(p.pid)
            if remote['status']=='RUNNING':stop_vm()
            sync()
            break
        time.sleep(45)
    # Frontiers only trigger verification. They never establish exclusion.
    # Re-run all still-open classes on the original Mac producer; all partial
    # evidence remains available and the full strict union gates final closure.
    from oracle_balance.aggregate import aggregate
    result=aggregate()
    assert not result['errors'] and not result['orphans'], 'Strict audit error: inspect'
    pending=[r['cls'] for r in result['records'] if r['status']!='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES']
    if pending:
        p=launch(pending)
        atomic_json(ROOT/'scheduler/remaining_state.json',dict(status='FINAL_RESIDUAL_GAPS',
            managers=dict(residual=dict(pid=p.pid,workers=8)),classes=pending))
        assert p.wait()==0,'Final-gap search error/candidate: inspect'
        result=aggregate()
    assert result['status']=='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE'
    atomic_json(ROOT/'scheduler/remaining_state.json',dict(
        status='ALL_SEARCH_BRANCHES_CONDITIONALLY_CLOSED_REQUIRES_FINAL_PROOF_REVIEW'))


if __name__=='__main__':
    main()
