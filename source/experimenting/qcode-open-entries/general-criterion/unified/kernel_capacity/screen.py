"""Exact necessary capacity counts over every C,N in the saved suffix cell."""
import json,sys
from pathlib import Path
from collections import Counter
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import coordinates,gram,lift_coordinate,nullspace,reduced_basis,span,weight
from mixed_frames.stage3.common import atomic_json,sha

def main():
 source=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json';classes=json.loads(source.read_text())['classes'];meta=json.loads((ROOT/'m8r1/delta_meta.json').read_text());assert meta['source_sha256']==sha(source)
 tables=np.load(ROOT/'m8r1/delta_tables.npy',mmap_mode='r');assert tables.shape==(len(classes),8192)
 profile=json.loads((HERE/'profiles.json').read_text());caps={}
 for r in profile['records']:caps.setdefault((r['q'],r['t']),set()).add(tuple(r['caps']))
 labels=np.arange(8192,dtype=np.int32);parity=np.array([v.bit_count()%2 for v in range(8192)],dtype=np.uint8)
 records=[];total=Counter();cases={1:[(1,0,0),(2,2,1)],3:[(2,0,1),(3,2,2)],5:[(3,0,2)]}
 for ci,c in enumerate(classes):
  C=c['representative']['generators'];h=c['representative']['hull'];rad=reduced_basis([lift_coordinate(v,C) for v in nullspace(gram(C,9),5)]);assert len(rad)==h
  words=sorted(w for w in span(rad) if w);planes=sorted({tuple(sorted((a,b,a^b))) for a in words for b in words if a!=b}) if h>1 else []
  rejected=[];counts=Counter()
  for q,t,ell in cases[h]:
   Ns=[()] if ell==0 else ([(w,) for w in words] if ell==1 else [p[:2] for p in planes])
   for ni,N in enumerate(Ns):
    mask=np.ones(8192,dtype=bool)
    for w in N:mask&=parity[labels&coordinates(w,meta['P_bases'][ci])]==0
    values=tables[ci][mask];assert len(values)==1<<(13-ell)
    actual=tuple(int(np.count_nonzero((values<=j)&(values>0))) for j in range(6))
    assert np.count_nonzero(values==0)==1
    ok=[p for p in caps[(q,t)] if all(actual[j]<=p[j] for j in range(6))]
    structural=None
    if ell==2:
     ws=list(N)+[N[0]^N[1]];wt=[weight(w,9) for w in ws];supports=[(w&511)|(w>>9) for w in ws]
     if sum(v<8 for v in wt)>2:structural='three_light'
     elif any(wt[a]==8 and any(b!=a and supports[b]&~supports[a]==0 for b in range(3)) for a in range(3)):structural='exclusive_rank_two_support'
    status=structural or ('capacity' if not ok else 'remaining');counts[(q,t,status)]+=1
    if status=='capacity':rejected.append(dict(q=q,t=t,N=list(N),actual=actual))
  total.update(counts);records.append(dict(cls=ci,hull=h,counts={str(k):v for k,v in counts.items()},capacity_rejections=rejected))
  if ci%5000==0:print(ci,dict(total),flush=True)
 report=dict(status='COMPLETE_RELATIVE_TO_SAVED_TABLES',records=records,counts={str(k):v for k,v in total.items()},hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'profiles.json',source,ROOT/'m8r1/delta_tables.npy',ROOT/'m8r1/delta_meta.json']},scope='Necessary C,N capacities conditional on saved classification and coset tables; no extension or target exclusion')
 atomic_json(HERE/'screen.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','hashes')}),flush=True)
if __name__=='__main__':main()
