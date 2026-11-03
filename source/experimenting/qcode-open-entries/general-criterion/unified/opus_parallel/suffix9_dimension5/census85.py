"""Mass-checked lengthening census of additive (n,2^k,>=5)_4 codes feeding (8,2^5).

Reuses, read-only, the recurrence and graph encoding of
release16/census/mass_census.py (loaded by exec, its __main__ block stripped).
Every (n,2^k,>=5) code shortened at its last coordinate (rank r of that column)
is an (n-1,2^(k-r),>=5) code; extensions enumerate ALL r-dim subspaces of the
quotient whose nonzero cosets have coset-minimum weight >=4 (complement of the
weight<=3 error syndromes). Each cell is checked by the orbit-stabilizer
labelled-mass identity. Layout here is INTERLEAVED (symbol i = bits 2i,2i+1);
outputs are additionally converted to release16 X|Z layout (x | z<<n).

Checkpoint: each finished cell is written to cells_d5.json immediately.
Usage: python3.11 census85.py [maxn] [target_k]
"""
import json, sys, time
from fractions import Fraction
from math import factorial
from pathlib import Path

OUT = Path(__file__).resolve().parent
SRC = Path('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/mass_census.py')
ns = {'__file__': str(SRC)}
exec(SRC.read_text().split('if __name__ == "__main__":')[0], ns)
span, wt, graph, aut_order, extensions, certificate = (ns[k] for k in
    ['span', 'wt', 'graph', 'aut_order', 'extensions', 'certificate'])
assert ns['D'] == 5


def to_xz(v, n):
    out = 0
    for i in range(n):
        out |= ((v >> (2 * i)) & 1) << i
        out |= ((v >> (2 * i + 1)) & 1) << (n + i)
    return out


def main(maxn=8, target_k=5):
    ck = OUT / f'cells_d5_{maxn}_{target_k}.json'
    store = json.loads(ck.read_text()) if ck.exists() else {}
    cells = {(4, 0): [([], aut_order(graph([], 4)))]}
    for key, v in store.items():
        n, k = map(int, key.split(','))
        cells[(n, k)] = [(g, a) for g, a in v['interleaved_with_aut']]
    for n in range(5, maxn + 1):
        for k in range(max(0, target_k - 2 * (maxn - n)), target_k + 1):
            if f'{n},{k}' in store:
                continue
            t0 = time.time(); classes = {}; lhs = Fraction(0); cand = 0
            for r in (0, 1, 2):
                for gens, aut in cells.get((n - 1, k - r), []):
                    exts, ext = extensions(gens, n - 1)
                    lhs += Fraction(factorial(n - 1) * 6 ** (n - 1), aut) * ext[r]
                    for g in exts[r]:
                        cand += 1
                        G = graph(g, n); c = certificate(G)
                        if c not in classes:
                            classes[c] = (g, aut_order(G))
            found = list(classes.values())
            for g, a in found:
                assert (factorial(n) * 6 ** n) % a == 0
                sp = span(g); assert len(set(sp)) == 1 << k
                assert k == 0 or min(wt(v, n) for v in sp if v) >= 5
            rhs = sum(Fraction(factorial(n) * 6 ** n, a) for _, a in found)
            ok = lhs == rhs
            print(json.dumps(dict(cell=[n, k], classes=len(found), candidates=cand,
                                  lhs=str(lhs), rhs=str(rhs), mass_equal=ok,
                                  seconds=round(time.time() - t0, 1))), flush=True)
            assert ok
            cells[(n, k)] = found
            store[f'{n},{k}'] = dict(classes=len(found), labelled_mass=str(rhs), mass_equal=ok,
                                     interleaved_with_aut=[[g, a] for g, a in found],
                                     xz=[[to_xz(v, n) for v in g] for g, _ in found])
            ck.write_text(json.dumps(store) + '\n')
    print('DONE', {k: v['classes'] for k, v in store.items()})


if __name__ == '__main__':
    main(*map(int, sys.argv[1:3]))
