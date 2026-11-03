"""Stop only this search's remote process tree, preserving all checkpoints."""
import os
import signal
import subprocess
import time
from pathlib import Path


def main():
    root=Path(__file__).resolve().parents[1]
    pid=int((root/'oracle_residual/pool.pid').read_text())
    lines=subprocess.check_output(['ps','-eo','pid=,ppid=,args='],text=True).splitlines()
    rows=[x.strip().split(None,2) for x in lines]
    rows=[x for x in rows if len(x)==3]
    own=[x for x in rows if len(x)==3 and int(x[0])==pid]
    if not own:
        return
    assert '-m oracle_residual.pool' in own[0][2], 'PID no longer belongs to our search'
    tree={pid}
    while True:
        new={int(p) for p,pp,c in rows if int(pp) in tree}-tree
        if not new:
            break
        tree|=new
    for p in [pid,*sorted(tree-{pid})]:
        try:os.kill(p,signal.SIGTERM)
        except ProcessLookupError:pass
    # Native enumerators handle SIGTERM by persisting LIMIT; neither native
    # nor Python partial work is promoted into a completed record here.
    time.sleep(1)


if __name__=='__main__':
    main()
