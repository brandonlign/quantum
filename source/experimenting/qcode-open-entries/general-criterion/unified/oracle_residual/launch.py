"""Detach-safe remote launch and explicit exit state, with no callbacks."""
import json
import os
import subprocess
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from mixed_frames.stage3.common import atomic_json
HERE = ROOT/'oracle_residual'


def main():
    (HERE/'runs').mkdir(exist_ok=True)
    (HERE/'candidates').mkdir(exist_ok=True)
    with (HERE/'pool.log').open('a') as log:
        p = subprocess.Popen([sys.executable,'-B','-m','oracle_residual.pool'],
            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT,
            env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'))
        (HERE/'pool.pid').write_text(str(p.pid)+'\n')
        atomic_json(HERE/'pool_state.json', dict(status='RUNNING',pid=p.pid,workers=2))
        code = p.wait()
    atomic_json(HERE/'pool_state.json', dict(status='BATCH_FINISHED' if code==0 else 'ERROR',
        pid=p.pid, workers=2, returncode=code))


if __name__=='__main__':
    main()
