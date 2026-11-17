"""Survey of the 32,865 dimension-five shortened codes: does the sixteen-qubit
clique/circuit test transfer?

For each class C (rank-one local stabilizer R = <Z^8>, d = 6):
  * coset minimum weights delta on W/C = F_2^13 by breadth-first search;
  * the deep-hole graph: vertices u with delta(u) >= 5 (possible images
    xi_i = M(Z_i)), edges u ~ v with delta(u + v) >= 4;
  * a search, capped at TIME_CAP seconds, for one 8-circuit satisfying every
    Z-only subset condition delta(sum_A xi) >= 6 - min(|A|, 8-|A|) and
    sum xi = 0, without the Gram/isotropy condition.  A circuit found means
    the sixteen-qubit style circuit test does not exclude the class.
Results go to survey_rankone.json.  Usage:
    python3 survey_rankone.py [--workers 6] [--limit N]
"""
import argparse, json, os, sys, time
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rankone_exact as rx

ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
CLASSES = os.path.join(ROOT, 'experimenting/qcode-open-entries/general-criterion/unified/'
                       'opus_parallel/suffix9_dimension5/suffix9_dim5.json')
TIME_CAP = 5.0


def work(item):
    idx, gens = item
    t0 = time.time()
    S = rx.Suffix(gens, d=6)
    deep = np.nonzero(S.delta >= 5)[0]
    pair = S.delta[deep[:, None] ^ deep[None, :]] >= 4
    edges = int((pair.sum() - 0) // 2)
    hist = np.bincount(S.delta, minlength=7).tolist()
    circuit = None
    t1 = time.time()
    try:
        for xi in rx.enumerate_Z(S, use_isotropy=False, deadline=t1 + TIME_CAP):
            circuit = list(xi)
            break
    except TimeoutError:
        circuit = None
    if circuit is not None:          # direct re-check of all 254 proper subsets
        for A in range(1, 255):
            s = 0
            for i in range(8):
                if A >> i & 1:
                    s ^= circuit[i]
            a = A.bit_count()
            assert S.delta[s] >= 6 - min(a, 8 - a), (idx, A)
        assert len(set(circuit)) == 8
    return dict(index=idx, hull=S.hull, delta_hist=hist, deep=len(deep), edges=edges,
                density=round(2 * edges / (len(deep) * (len(deep) - 1)), 4),
                circuit=circuit, circuit_seconds=round(time.time() - t1, 3),
                seconds=round(time.time() - t0, 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--limit', type=int, default=None)
    a = ap.parse_args()
    cls = json.load(open(CLASSES))['classes']
    assert len(cls) == 32865
    items = [(c['class_index'], c['representative']['generators']) for c in cls]
    if a.limit:
        items = items[:a.limit]
    t0 = time.time()
    out = []
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(work, items, chunksize=20):
            out.append(r)
            if len(out) % 2000 == 0:
                print(len(out), 'classes, %.0fs' % (time.time() - t0), flush=True)
    out.sort(key=lambda r: r['index'])
    no_circuit = [r['index'] for r in out if r['circuit'] is None]  # not found within TIME_CAP
    summary = dict(classes=len(out), hull_counts={h: sum(1 for r in out if r['hull'] == h) for h in (1, 3, 5)},
                   deep_min=min(r['deep'] for r in out), deep_max=max(r['deep'] for r in out),
                   deep_mean=round(sum(r['deep'] for r in out) / len(out), 1),
                   density_min=min(r['density'] for r in out), density_max=max(r['density'] for r in out),
                   classes_without_circuit=no_circuit, max_circuit_seconds=max(r['circuit_seconds'] for r in out),
                   wall_seconds=round(time.time() - t0, 1))
    print(json.dumps(summary, indent=1))
    name = 'survey_rankone.json' if not a.limit else 'survey_rankone_sample.json'
    json.dump(dict(summary=summary, classes=out), open(os.path.join(HERE, name), 'w'))


if __name__ == '__main__':
    main()
