"""Second check for m=9: phisearch3 (forward checking + proven symmetry breaking) on coset labels
from sat_check.own_dmin (independent of full.syndromes, used by run_c9.log)."""
import os as _os
CELLS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'census', 'cells.json')
import json, sys, subprocess, time
from concurrent.futures import ProcessPoolExecutor
from sat_check import own_dmin
from phi_input import make_input
cells = json.load(open(CELLS))
def job(i):
    dmin, r = own_dmin(cells['7,4'][i], 7); assert r == 10; t = time.time()
    out = subprocess.run(['./phisearch3'], input=make_input(dmin, 9), capture_output=True, text=True).stdout.strip()
    return i, out, time.time() - t
if __name__ == '__main__':
    ids = range(43) if sys.argv[1] == 'all' else [int(x) for x in sys.argv[1].split(',')]
    W = int(sys.argv[2]) if len(sys.argv) > 2 else 2; surv = []
    with ProcessPoolExecutor(W) as ex:
        for i, out, dt in ex.map(job, ids):
            print(9, i, out, f'{dt:.1f}s', flush=True)
            if not out.startswith('NONE'): surv.append(i)   # FOUND, or LIMIT (unresolved)
    print('m=9 survivors', surv, flush=True)
