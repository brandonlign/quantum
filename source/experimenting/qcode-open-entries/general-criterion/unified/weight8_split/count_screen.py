"""Exact Farkas screen: pure [[17,3,6]] with fixed (A6,A7), with and without a frame block."""
import sys, json; sys.path.insert(0,'split')
import numpy as np
from system import build, farkas_exact
n,k,d=17,3,6
frames={'none':None,'r1':['ZZZZZZZZ'],'668':['ZZZZZZZZ','XXXXZZII'],'688':['ZZZZZZZZ','XXXXXXII'],
        '778':['ZZZZZZZZ','XXXXXXZI'],'888':['ZZZZZZZZ','XXXXXXXX']}
out={}
for a6 in range(3):
  for a7 in range(3):
    for tag,g in frames.items():
        w={i:0 for i in range(1,6)}; w[6]=a6; w[7]=a7
        blocks=[] if g is None else [(8,g)]
        Meq,b,G,*_=build(n,k,d,blocks,weights=w)
        cert,msg=farkas_exact(Meq,b,G)
        out[f'{a6}{a7}{tag}']=bool(cert)
    print(a6,a7,{t:out[f'{a6}{a7}{t}'] for t in frames},flush=True)
json.dump(out,open('general-criterion/unified/weight8_split/count_screen.json','w'),indent=1)
