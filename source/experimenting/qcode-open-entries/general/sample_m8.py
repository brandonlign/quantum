"""Difficulty probe: phisearch3 with a node limit on sampled (9,2^5,6) classes for [[17,3,6]] m=8."""
import sys, pickle, subprocess, time, random
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, 'general')
from prefix_test import to_xz, coset_dmin, prefix_input
raw = pickle.load(open('census6/cells2_d6_9_5.pkl', 'rb'))['9,5']
def job(i):
    g = [to_xz(w, 9) for w in raw[i]]; dmin, codim = coset_dmin(g, 9); t = time.time()
    out = subprocess.run(['shorten/phisearch3', '300000000'], input=prefix_input(dmin, 8, 1, 6), capture_output=True, text=True).stdout.strip()
    return i, out, time.time() - t
if __name__ == '__main__':
  random.seed(1); ids = sorted(random.sample(range(len(raw)), 12))
  with ProcessPoolExecutor(2) as ex:
    for i, out, dt in ex.map(job, ids): print(i, out, f'{dt:.0f}s', flush=True)
