#!/usr/bin/env python3
"""Split-enumerator Farkas screen of the five weight-eight local frames of a putative
pure [[17,3,6]] (WEIGHT8_RECENTER.md). CERT = exactly verified infeasibility of a
necessary relaxation; 'none' proves nothing."""
import sys, json, time
sys.path.insert(0, 'split')
from system import build, farkas_exact
n, k, d = 17, 3, 6
frames = {'r1': ['ZZZZZZZZ'], 'r2_4_2': ['ZZZZZZZZ', 'XXXXZZII'], 'r2_6_0': ['ZZZZZZZZ', 'XXXXXXII'],
          'r2_6_1': ['ZZZZZZZZ', 'XXXXXXZI'], 'r2_8_0': ['ZZZZZZZZ', 'XXXXXXXX']}
out = {}
for tag, g in frames.items():
    for pure in (False, True):
        t = time.time()
        w = {i: 0 for i in range(1, 6)} if pure else None
        Meq, b, G, cols, *_ = build(n, k, d, [(8, g)], weights=w)
        cert, msg = farkas_exact(Meq, b, G)
        key = f'{tag}{"_pure" if pure else ""}'
        out[key] = dict(cert=bool(cert), msg=str(msg), cols=len(cols), sec=round(time.time()-t, 1))
        print(key, out[key], flush=True)
json.dump(out, open('general-criterion/unified/weight8_split/screen.json', 'w'), indent=1)
