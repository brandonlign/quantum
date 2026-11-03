"""Run phisearch on every class of the given cells, in parallel."""
import os as _os
CELLS = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'census', 'cells.json')
import json, sys, subprocess, time
from concurrent.futures import ProcessPoolExecutor
from full import syndromes
from phi_input import make_input
import os; PROG = [os.environ.get('PHIPROG', './phisearch2'), os.environ.get('PHILIMIT', '0')]
cells = json.load(open(CELLS))
KEY = {6: '10,10', 7: '9,8', 8: '8,6', 9: '7,4'}
def job(a):
    m, i = a
    dmin, syn, r = syndromes(cells[KEY[m]][i], 16 - m)
    assert r == 10
    t = time.time()
    out = subprocess.run(PROG, input=make_input(dmin, m), capture_output=True, text=True).stdout.strip()
    return m, i, out, time.time() - t
if __name__ == '__main__':
    ms = [int(x) for x in sys.argv[1:]]
    jobs = [(m, i) for m in ms for i in range(len(cells[KEY[m]]))]
    surv = {m: [] for m in ms}
    with ProcessPoolExecutor(7) as ex:
        for m, i, out, dt in ex.map(job, jobs):
            print(m, i, out, f'{dt:.1f}s', flush=True)
            if not out.startswith('NONE'): surv[m].append(i)   # FOUND, or LIMIT (unresolved)
    for m in ms: print(f'm={m}: {len(cells[KEY[m]])} classes, survivors {surv[m]}', flush=True)
