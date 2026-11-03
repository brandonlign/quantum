"""Cost profiles of E'=L^perp/L for every admissible isotropic prefix-logical subspace L.

R=<Z^8> on m=8 sites. E=R^perp/R (2^14 classes, each represented by the physical
word with X even and Z bit 7 cleared... we store min-cost over R). L ranges over
isotropic subspaces of E of dimension ell in {1,2,3} whose nonzero elements all have
cost>=6 (prefix-only logicals of a distance-6 code). Output: for each ell the set of
distinct cumulative tails T_ell[c]=#{nonzero e' in L^perp/L: cost'(e')>=c}.
Exhaustive over ALL such L (no symmetry assumption) for ell=1,2; ell=3 by extension of all ell=2.
"""
import json,itertools
import numpy as np
m=8;full=(1<<m)-1
# physical representatives of R^perp: x even weight; class = (x, z) modulo z->z^full
reps=[(x,z) for x in range(256) if bin(x).count('1')%2==0 for z in range(128)]  # z bit7 =0 normal form
idx={r:i for i,r in enumerate(reps)}
def norm(x,z):
    if z>>7&1:z^=full
    return idx[(x,z)]
def wt(x,z):return bin(x|z).count('1')
cost=np.array([min(wt(x,z),wt(x,z^full)) for x,z in reps])
X=np.array([x for x,z in reps]);Z=np.array([z for x,z in reps])
def sym(i,j):return (bin(X[i]&Z[j]).count('1')+bin(X[j]&Z[i]).count('1'))&1
# vectorised add table via (x,z) xor with normalisation
def add(i,j):return norm(X[i]^X[j],Z[i]^Z[j])
allx=X;allz=Z
def sym_vec(i):return (np.vectorize(lambda a:bin(a).count('1'))((allx&Z[i]))+np.vectorize(lambda a:bin(a).count('1'))((X[i]&allz)))&1
pc=np.array([bin(v).count('1') for v in range(256)])
def symv(i):return (pc[allx&Z[i]]+pc[X[i]&allz])&1
def addv(i):
    x=allx^X[i];z=allz^Z[i];fl=(z>>7)&1;z=np.where(fl==1,z^full,z)
    return x*128+z   # idx = position: reps ordered x-major over even x... build map
# position map from (x,z) key to index
key2=np.full(256*128,-1);
for i,(x,z) in enumerate(reps):key2[x*128+z]=i
def addvi(i):return key2[addv(i)]
heavy=[i for i in range(len(reps)) if cost[i]>=6]
print('E size',len(reps),'cost>=6',len(heavy))
def tails(Lspan):
    """Lspan: list of indices of span(L) incl 0. returns tail counts c=0..8 for L^perp/L."""
    perp=np.ones(len(reps),bool)
    for i in Lspan:perp&=symv(i)==0
    c=cost.copy()
    for i in Lspan:c=np.minimum(c,cost[addvi(i)])
    # every class of L^perp/L appears |L| times in perp; cost' constant on class
    zero=np.zeros(len(reps),bool);zero[Lspan]=True
    vals=c[perp&~zero]
    t=[int((vals>=k).sum())//len(Lspan) for k in range(9)]
    return tuple(t)
out={}
# ell=1
prof1=set();ok1=[]
for i in heavy:
    if i==0:continue
    prof1.add(tails([0,i]))
out[1]=sorted(prof1);print('ell1 profiles',out[1])
# ell=2: unordered pairs {i,j} commuting, i^j heavy
heavyset=np.zeros(len(reps),bool);heavyset[heavy]=True
prof2=set();planes=set()
for i in heavy:
    s=symv(i);a=addvi(i)
    for j in heavy:
        if j<=i or s[j] or not heavyset[a[j]]:continue
        k=int(a[j]);pl=tuple(sorted((i,j,k)))
        if pl in planes:continue
        planes.add(pl)
print('planes',len(planes),flush=True)
for pl in planes:prof2.add(tails([0,*pl]))
out[2]=sorted(prof2);print('ell2 profiles',len(prof2))
json.dump({str(k):v for k,v in out.items()},open('fixed_profiles_12.json','w'))
