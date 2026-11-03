"""For open entries whose only/likely enumerator profile is pure: which minimum stabilizer
weights m survive the split LP (S contains Z^m, A_1..A_{m-1}=0, pure up to d-1)?  For survivors,
report the shortened cell (f, 2^{dimC}, >=d), dimC = n+k-2m+1 (assuming dim S_T = 1)."""
import sys, json, time; sys.path.insert(0, 'split')
from system import build, farkas_exact
for n, k, d in [tuple(map(int, a.split(','))) for a in sys.argv[1:]]:
    for m in range(d, n + 1):
        t = time.time()
        W = {j: 0 for j in range(1, m)}
        Meq, b, G, *_ = build(n, k, d, [], weights=W)
        c0, _ = farkas_exact(Meq, b, G)
        if c0: print(f'[[{n},{k},{d}]] m>={m}: excluded (plain CERT) -> stop', flush=True); break
        try:
            Meq, b, G, *_ = build(n, k, d, [(m, ['Z' * m])], weights=W)
            c, _ = farkas_exact(Meq, b, G)
        except Exception as e: c = None; print('err', e)
        f = n - m; dc = n + k - 2 * m + 1
        print(f'[[{n},{k},{d}]] m={m}: {"CERT" if c else "open"}  cell ({f}, 2^{dc}, >={d})  {time.time()-t:.0f}s', flush=True)
