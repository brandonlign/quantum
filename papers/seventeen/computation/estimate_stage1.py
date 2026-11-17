"""Random-probe (Knuth) estimate of the Stage-1 search-tree size of
rankone_exact.enumerate_Z (d = 6, no kernel restriction) on a sample of
shortened codes.  Each probe walks one random root-to-leaf path through the
same tree that enumerate_Z traverses (same subset-sum filters, same isotropy
rank test, same tight-radical restriction) and multiplies the branching
factors; the average over probes is an unbiased estimate of the number of
nodes at each depth.

Usage: python3 estimate_stage1.py [--per-hull 1:60,3:30,5:10] [--probes 200]
Writes estimate_stage1.json.
"""
import argparse, json, os, random, sys, time
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rankone_exact as rx
from survey_rankone import CLASSES


def estimate(item):
    idx, gens, probes = item
    S = rx.Suffix(gens, d=6)
    TH = rx.thresholds_Z(6)
    OK = {t: S.delta >= t for t in range(8)}
    root = np.nonzero(OK[TH[1]])[0]
    rng = random.Random(idx)
    est = [0.0] * 8
    t0 = time.time()
    for _ in range(probes):
        sums, cand, Ub, w = {}, root, [1 << j for j in range(13)], 1.0
        for k in range(7):
            kids = []
            for v in cand:
                v = int(v)
                Ub2 = rx.restrict_annihilator(Ub, v)
                if Ub2 is None:
                    continue
                ok, rad = rx.isotropy_state(S, Ub2, k + 1)
                if ok:
                    kids.append((v, Ub2, rad))
            if not kids:
                break
            w *= len(kids)
            est[k + 1] += w
            v, Ub, rad = kids[rng.randrange(len(kids))]
            if k == 6:
                break
            c2 = cand[cand > v]
            for s, t in [(1, v)] + [(s + 1, x ^ v) for s, l in sums.items() for x in l]:
                c2 = c2[OK[TH[s + 1]][c2 ^ t]]
            if rad:
                for r in rad:
                    c2 = c2[(np.bitwise_count(c2 & r) & 1) == 0]
            s2 = {s: list(l) for s, l in sums.items()}
            s2.setdefault(1, []).append(v)
            for s, l in sums.items():
                s2.setdefault(s + 1, []).extend(x ^ v for x in l)
            sums, cand = s2, c2
    return dict(index=idx, hull=S.hull, est_nodes=[e / probes for e in est[1:]],
                est_total=sum(est[1:]) / probes, seconds=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-hull', default='1:60,3:30,5:10')
    ap.add_argument('--probes', type=int, default=200)
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    cls = json.load(open(CLASSES))['classes']
    rng = random.Random(17)
    items = []
    for part in a.per_hull.split(','):
        h, m = map(int, part.split(':'))
        pool = [c for c in cls if c['representative']['hull'] == h]
        for c in rng.sample(pool, m):
            items.append((c['class_index'], c['representative']['generators'], a.probes))
    with Pool(a.workers) as p:
        res = sorted(p.map(estimate, items), key=lambda r: r['index'])
    summ = {}
    for h in (1, 3, 5):
        r = [x for x in res if x['hull'] == h]
        if r:
            tot = sorted(x['est_total'] for x in r)
            summ[h] = dict(sample=len(r), mean_nodes=float(np.mean(tot)), median_nodes=tot[len(tot) // 2],
                           max_nodes=tot[-1])
    print(json.dumps(summ, indent=1))
    json.dump(dict(probes=a.probes, summary=summ, classes=res), open(os.path.join(HERE, 'estimate_stage1.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
