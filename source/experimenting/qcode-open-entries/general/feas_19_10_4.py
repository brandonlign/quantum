import sys, time
sys.path.insert(0, '/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/split')
from system import build, farkas_exact
n, k, d = 19, 10, 4
cases = [('w1', [(1, ['Z'])], None), ('w2', [(2, ['ZZ'])], None), ('w3', [(3, ['ZZZ'])], None),
         ('w2 A2=1', [(2, ['ZZ'])], {2: 1}), ('w3 A1=A2=0', [(3, ['ZZZ'])], {1: 0, 2: 0}),
         ('pure w4', [(4, ['ZZZZ'])], {1: 0, 2: 0, 3: 0}), ('pure w5', [(5, ['ZZZZZ'])], {1: 0, 2: 0, 3: 0, 4: 0})]
cases += [(f'minwt>={t}', [], {w: 0 for w in range(1, t)}) for t in range(6, 12)]
for name, blocks, W in cases:
    t = time.time()
    try:
        Meq, b, G, cols, le, ll = build(n, k, d, blocks, weights=W)
        cert, ok = farkas_exact(Meq, b, G)
        res = 'CERT' if cert else 'none'
    except Exception as e:
        res = 'ERR ' + str(e)[:80]
    print(f'{name:14s} {res}  {time.time()-t:.1f}s', flush=True)
