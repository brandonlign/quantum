"""[[17,6,5]] prefix-map search, r=1 cases, reusing the [[16,5,5]] cells: C is (17-m, 2^(24-2m)), W/C dim 10."""
import json, sys, subprocess, time, os
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search')
os.chdir('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/search')
from full import syndromes
from phi_input import make_input
from concurrent.futures import ProcessPoolExecutor
cells = json.load(open('../census/cells.json'))
KEY = {7: '10,10', 8: '9,8', 9: '8,6', 10: '7,4'}
PROG = os.environ.get('PHIPROG', './phisearch2')
def job(a):
    m, i = a
    dmin, syn, r = syndromes(cells[KEY[m]][i], 17 - m)
    assert r == 10
    t = time.time()
    out = subprocess.run([PROG, os.environ.get('PHILIMIT', '0')], input=make_input(dmin, m), capture_output=True, text=True).stdout.strip()
    return m, i, out, time.time() - t
if __name__ == '__main__':
    m = int(sys.argv[1]); ids = range(len(cells[KEY[m]])) if len(sys.argv) < 3 else map(int, sys.argv[2].split(','))
    surv = []
    with ProcessPoolExecutor(int(os.environ.get('W', '5'))) as ex:
        for m_, i, out, dt in ex.map(job, [(m, i) for i in ids]):
            print(m_, i, out, f'{dt:.1f}s', flush=True)
            if not out.startswith('NONE'): surv.append(i)
    print(f'm={m} survivors {surv}', flush=True)
