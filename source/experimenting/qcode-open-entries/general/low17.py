import sys; sys.path.insert(0, 'split')
from system import build, farkas_exact
n,k,d=17,6,5
for name,blocks,W in [('w1',[(1,['Z'])],{}),('w2',[(2,['ZZ'])],{}),('w3',[(3,['ZZZ'])],{}),('w4',[(4,['ZZZZ'])],{}),
                      ('A4>=1 plain',[], None)]:
    if W is None: continue
    Meq,b,G,*_=build(n,k,d,blocks,weights=W or None); c,_=farkas_exact(Meq,b,G); print(name,'CERT' if c else 'open',flush=True)
for a4 in range(0,4):
    Meq,b,G,*_=build(n,k,d,[],weights={1:0,2:0,3:0,4:a4}); c,_=farkas_exact(Meq,b,G); print('A1-3=0,A4=',a4,'CERT' if c else 'open',flush=True)
