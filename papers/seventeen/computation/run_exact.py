"""Run the exact rank-one test (rankone_exact.py) on chosen shortened codes.

    python3 run_exact.py CLASS [CLASS ...] [--ndim k] [--nindex i,j,...]

Without --ndim the search is over all of W/C: an empty result means that the
class has no rank-one extension at all (Theorem 2.1).  With --ndim k the
search is repeated for each k-dimensional subspace N of rad(P) (listed in a
fixed order; --nindex selects some of them) with all images in N^0; this is
the kernel-case split of Section 5 (im M = N^0 with N = ker tau).
Results are appended to run_exact.jsonl.
"""
import argparse, itertools, json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rankone_exact as rx
from survey_rankone import CLASSES


def subspaces(basis, k):
    """All k-dimensional subspaces of span(basis), each as a reduced basis."""
    seen, out = set(), []
    words = [w for w in rx.span_words(basis) if w]
    for combo in itertools.combinations(words, k):
        if rx.rank(list(combo)) < k:
            continue
        sp = frozenset(rx.span_words(list(combo)))
        if sp in seen:
            continue
        seen.add(sp)
        out.append(list(combo))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('classes', type=int, nargs='+')
    ap.add_argument('--ndim', type=int, default=0)
    ap.add_argument('--nindex', default=None)
    a = ap.parse_args()
    cls = json.load(open(CLASSES))['classes']
    for ci in a.classes:
        S = rx.Suffix(cls[ci]['representative']['generators'], d=6)
        Ns = subspaces(S.rad, a.ndim) if a.ndim else [[]]
        idxs = range(len(Ns)) if a.nindex is None else [int(x) for x in a.nindex.split(',')]
        for i in idxs:
            t0, st, found = time.time(), {}, []
            ztuples = 0
            for xi in rx.enumerate_Z(S, N=Ns[i], stats=st):
                ztuples += 1
                sols, _ = rx.solve_X(S, xi)
                for e in sols:
                    ok, why = rx.full_check(S, xi, e)
                    assert ok, why
                    found.append([list(xi), list(e)])
            rec = dict(cls=ci, hull=S.hull, ndim=a.ndim, n_subspaces=len(Ns), n_index=i, N=Ns[i],
                       nodes=st['nodes'], z_tuples=ztuples, extensions=found,
                       seconds=round(time.time() - t0, 1))
            print(json.dumps(rec), flush=True)
            with open(os.path.join(HERE, 'run_exact.jsonl'), 'a') as f:
                f.write(json.dumps(rec) + '\n')


if __name__ == '__main__':
    main()
