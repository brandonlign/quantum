import sys; sys.path.insert(0, 'split')
from system import build, farkas_exact
n,k,d=17,3,6
for w in range(1,6):
    Meq,b,G,*_=build(n,k,d,[(w,['Z'*w])]); c,_=farkas_exact(Meq,b,G); print('S has weight',w,'word:','CERT' if c else 'open',flush=True)
for w in range(1,6):
    W={j:0 for j in range(1,6) if j!=w}
    Meq,b,G,*_=build(n,k,d,[],weights={**W}); c,_=farkas_exact(Meq,b,G)
