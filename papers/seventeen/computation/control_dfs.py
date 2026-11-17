"""Control for the Stage-1 depth-first search of rankone_exact.py.

For a few shortened codes at d = 5 (where Stage 1 has many solutions), take
a valid Z-tuple xi, force its first five entries, and compare the set of
completions returned by enumerate_Z (with its subset filters, partial
isotropy rank test and tight-radical restriction) with a brute-force list:
all pairs xi_5 < xi_6 of syndromes such that xi_7 = xi_0 + ... + xi_6 exceeds
xi_6, every proper subset A of the eight images satisfies
delta(sum_A) >= d - min(|A|, 8-|A|), and the 6-dimensional annihilator of the
images is Gamma-isotropic.  The two sets must be equal.
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rankone_exact as rx
from survey_rankone import CLASSES


def brute(S, prefix):
    d = S.d
    TH = rx.thresholds_Z(d)
    delta = S.delta
    out = []
    # subset sums of the first five with sizes
    base = [(0, 0)]
    for v in prefix:
        base += [(s ^ v, c + 1) for s, c in base]
    allv = np.arange(8192)
    for v5 in range(prefix[-1] + 1, 8192):
        # conditions involving v5 but not v6
        if any(delta[s ^ v5] < TH[c + 1] for s, c in base):
            continue
        cand = allv[allv > v5]
        for s, c in base:
            # subsets containing v6 only, and both v5 and v6
            cand = cand[delta[cand ^ s] >= TH[c + 1]]
            cand = cand[delta[cand ^ s ^ v5] >= TH[c + 2]]
        tot = 0
        for v in prefix:
            tot ^= v
        for v6 in cand:
            v6 = int(v6)
            x7 = tot ^ v5 ^ v6
            if x7 <= v6:
                continue
            xi = list(prefix) + [v5, v6]
            Ub = [1 << j for j in range(13)]
            for v in xi:
                Ub = rx.restrict_annihilator(Ub, v)
            if all(S.gamma(a, b) == 0 for a in Ub for b in Ub):
                out.append(tuple(xi) + (x7,))
    return sorted(out)


cls = json.load(open(CLASSES))['classes']
report = []
for ci in (0, 1, 5, 150):
    S = rx.Suffix(cls[ci]['representative']['generators'], d=5)
    xi = next(rx.enumerate_Z(S))
    t0 = time.time()
    dfs = sorted(rx.enumerate_Z(S, forced=xi[:5]))
    t1 = time.time()
    bf = brute(S, xi[:5])
    t2 = time.time()
    assert dfs == bf, (ci, len(dfs), len(bf))
    report.append(dict(cls=ci, hull=S.hull, prefix=list(xi[:5]), completions=len(dfs),
                       dfs_seconds=round(t1 - t0, 1), brute_seconds=round(t2 - t1, 1)))
    print('class', ci, 'hull', S.hull, 'completions', len(dfs), 'DFS = brute force', flush=True)
json.dump(report, open(os.path.join(HERE, 'control_dfs.json'), 'w'), indent=1)
print('PASS')
