"""Tiny fixture for NEXT.md: rigid local logical plane H=R+K for (q,t)=(3,2).

H=<Z^8,r,u,v> on 8 sites, polar radical <Z^8,r>, <u,v>=1.  NEXT.md Lemma 2
forces: non-Z^8 even words weight 6 (X-weight 4), odd words weight 7
(X-weight 6), odd words exactly the coset u+v+<Z^8,r> (four H words = two
K classes), every site full.  X weights on K: five 4s and two 6s (total 32).
Enumerates X-column multisets (up to site order) and all Z parts.
"""
import itertools, json
from collections import Counter

pc = lambda v: bin(v).count("1")
def wt(x, z): return pc(x | z)
def sp(a, b): return (pc(a[0] & b[1]) + pc(a[1] & b[0])) & 1
def add(a, b): return (a[0] ^ b[0], a[1] ^ b[1])
Z8 = (0, 255)

def xcodes():
    cols = range(1, 8)
    seen = []
    for comb in itertools.combinations_with_replacement(cols, 8):
        w = {y: sum((pc(y & c) & 1) for c in comb) for y in range(1, 8)}
        # functionals e1=r, e2=u, e3=v
        # integer coefficient order r,u,r+u,v,r+v,u+v,r+u+v
        if [w[y] for y in range(1, 8)] == [4, 4, 4, 4, 4, 6, 6]:
            seen.append(comb)
    return seen

def xvec(comb, bit):
    return sum(1 << i for i, c in enumerate(comb) if c >> bit & 1)

def words(gens):
    out = []
    for m in range(1, 16):
        w = (0, 0)
        for i, g in enumerate(gens):
            if m >> i & 1: w = add(w, g)
        out.append(w)
    return out

results = []
for comb in xcodes():
    xr, xu, xv = (xvec(comb, b) for b in range(3))
    sols = 0; example = None
    # local X->Y freedom on X-support of r: fix z_r = 0 there
    for zr in range(256):
        if zr & xr: continue
        r = (xr, zr)
        if wt(*r) != 6 or wt(*add(r, Z8)) != 6: continue
        for zu in range(256):
            u = (xu, zu)
            if sp(r, u): continue
            if any(wt(*w) != 6 for w in (u, add(u, Z8), add(u, r), add(add(u, r), Z8))): continue
            for zv in range(256):
                v = (xv, zv)
                if sp(r, v) or not sp(u, v): continue
                ws = words([Z8, r, u, v])
                if Counter(wt(*w) for w in ws) != Counter({6: 10, 7: 4, 8: 1}): continue
                sols += 1
                if example is None: example = [r, u, v]
    results.append({"columns": comb, "solutions": sols, "example": example})

print(json.dumps(results, indent=1))
json.dump(results, open(__file__.replace(".py", ".json"), "w"), indent=1)
