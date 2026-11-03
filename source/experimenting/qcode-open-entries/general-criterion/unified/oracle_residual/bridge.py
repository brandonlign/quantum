"""Pull atomic remote checkpoints; no model invocation or destructive sync."""
import json
import subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'oracle_residual'
REMOTE = 'qcode1736-oracle-v1/project'


def sync():
    (HERE/'runs').mkdir(exist_ok=True)
    subprocess.run(['rsync', '-az', '--timeout=35', '-e',
        'ssh -o BatchMode=yes -o ConnectTimeout=10',
        f'oracle-vm:{REMOTE}/oracle_residual/runs/', str(HERE/'runs')+'/'],
        check=True, timeout=50)
    p = subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=10',
        'oracle-vm', f'cat ~/{REMOTE}/oracle_residual/pool_state.json'],
        capture_output=True, text=True, check=True, timeout=20)
    state = json.loads(p.stdout)
    subprocess.run(['rsync','-az','--timeout=15','-e',
        'ssh -o BatchMode=yes -o ConnectTimeout=10',
        f'oracle-vm:{REMOTE}/oracle_residual/pool.log', str(HERE/'remote_pool.log')],
        check=True, timeout=25)
    if state['status'] != 'RUNNING':
        subprocess.run(['rsync','-az','--timeout=15','-e',
            'ssh -o BatchMode=yes -o ConnectTimeout=10',
            f'oracle-vm:{REMOTE}/oracle_residual/candidates/', str(HERE/'candidates')+'/'],
            check=True, timeout=25)
    return state


if __name__ == '__main__':
    print(json.dumps(sync()), flush=True)
