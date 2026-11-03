"""SAT version of the generic prefix-map test (independent of the C searches).  Same mathematics
as prefix_test.py; phi is a 2m-2r column matrix over F2 into W/C (codim bits).  Symmetry
breaking is NOT used.  usage: prefix_sat.py n k d m r cellsfile key ids [timeout_s]"""
import sys, pickle, time, threading
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general')
from prefix_test import to_xz, coset_dmin, prefix_input
from pysat.solvers import Cadical153


def sat(dmin, codim, inp_text):
    lines = inp_text.split('\n'); mm, N = map(int, lines[0].split())
    nb = 2 * mm - 2; nl = int(lines[2]); low = [tuple(map(int, l.split())) for l in lines[3:3 + nl]]
    nv = [0]
    def new():
        nv[0] += 1; return nv[0]
    M = [[new() for _ in range(nb)] for _ in range(codim)]
    cl = []; y = {}
    for c in range(1, 1 << nb):
        hb = c.bit_length() - 1; rest = c ^ (1 << hb)
        if rest == 0: y[c] = [M[kk][hb] for kk in range(codim)]; continue
        out = []
        for kk in range(codim):
            a, b = y[rest][kk], M[kk][hb]; z = new()
            cl += [[-z, a, b], [-z, -a, -b], [z, -a, b], [z, a, -b]]; out.append(z)
        y[c] = out
    for c, t in low:
        for s in range(N):
            if dmin[s] < t:
                cl.append([(-y[c][kk] if s >> kk & 1 else y[c][kk]) for kk in range(codim)])
    return cl


if __name__ == '__main__':
    n, k, d, m, r = map(int, sys.argv[1:6]); fn, key, ids = sys.argv[6:9]
    TO = float(sys.argv[9]) if len(sys.argv) > 9 else 0
    raw = pickle.load(open(fn, 'rb'))[key]; f = n - m
    for i in [int(x) for x in ids.split(',')]:
        g = [to_xz(w, f) for w in raw[i]]; t = time.time()
        dmin, codim = coset_dmin(g, f)
        cl = sat(dmin, codim, prefix_input(dmin, m, r, d))
        S = Cadical153(bootstrap_with=cl)
        if TO:
            tm = threading.Timer(TO, S.interrupt); tm.start()
            res = S.solve_limited(expect_interrupt=True); tm.cancel()
        else:
            res = S.solve()
        print(f'[[{n},{k},{d}]] m={m} r={r} {key} class {i}: {"SAT" if res else ("UNSAT" if res is False else "TIMEOUT")} '
              f'({len(cl)} clauses, {time.time()-t:.0f}s)', flush=True)
