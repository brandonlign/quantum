"""Exact split Farkas screen for [[17,6,5]], testing the already certified pure minimum-5 branch.
CERT means the necessary relaxation is impossible. This is not a complete proof."""
import sys,json;sys.path.insert(0,'split')
from system import build,farkas_exact
n,k,d=17,6,5
out={}
for m in range(5,10):
  w={i:0 for i in range(1,5)}
  A,b,G,*_=build(n,k,d,[(m,['Z'*m])],weights=w)
  cert,msg=farkas_exact(A,b,G);out[m]=dict(cert=bool(cert),msg=str(msg));print(m,out[m],flush=True)
json.dump(out,open('general-criterion/unified/weight8_split/minword1765.json','w'),indent=1)
