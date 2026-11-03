"""Complete small isotropic H candidates for the 888 frame."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import dual,gram,lift_coordinate,reduced_basis,span,symp,weight
from suffix_mass_check.common import canonical
from mixed_frames.stage3.common import atomic_json,sha

def normalize(H):
 full=sorted(w for w in span(H) if weight(w,8)==8);assert len(full)==3 and full[0]^full[1]==full[2]
 p,q=full[:2]
 def apply(w):
  out=0
  for i in range(8):
   a=((p>>i)&1)|(((p>>(8+i))&1)<<1);b=((q>>i)&1)|(((q>>(8+i))&1)<<1)
   assert a and b and a!=b
   c=((w>>i)&1)|(((w>>(8+i))&1)<<1);v={0:0,a:2,b:1,a^b:3}[c];out|=(v&1)<<i;out|=(v>>1)<<(8+i)
  return out
 H=reduced_basis([apply(h) for h in H]);assert 255<<8 in span(H) and 255 in span(H)
 assert not any(gram(H,8)) and min(weight(w,8) for w in span(H) if w)>=6
 return H

def complement(space,sub):
 seen=list(sub);out=[]
 for v in space:
  if len(reduced_basis(seen+[v]))>len(seen):seen.append(v);out.append(v)
 return out

def main():
 start=time.monotonic();R=[255<<8,255];h3=R+[sum(1<<i for i in [2,3,6,7])|sum(1<<(8+i) for i in [4,5,6,7])];h3=normalize(h3)
 rows=json.loads((ROOT/'gauge_rankone_h5/all_critical_dim4.json').read_text())['classes'];h4=[]
 for i,row in enumerate(rows):
  H=row['generators'];full=sorted(w for w in span(H) if weight(w,8)==8)
  if len(full)==3 and full[0]^full[1]==full[2]:
   assert row['gram_rank']==0;h4.append(dict(source_index=i,H=normalize(H)))
 assert len(h4)==3 and len({canonical(r['H'],8) for r in h4})==3
 h5={};valid=0
 for row in h4:
  B=row['H'];L=complement(dual(B,8),B);assert len(L)==8
  for v in span(L):
   if min(weight(v^h,8) for h in span(B))<6:continue
   H=normalize(B+[v]);assert len(H)==5;key=canonical(H,8);valid+=1
   h5.setdefault(key,dict(H=H,parent_indices=[]))['parent_indices'].append(row['source_index'])
 assert len(h5)==1
 result=dict(status='PASS_RELATIVE_TO_H4_CENSUS',prefix={0:[dict(H=R)],1:[dict(H=h3)],2:h4,3:list(h5.values())},valid_H5_extensions=valid,seconds=time.monotonic()-start,scope='Local H candidates only, not target constructions or exclusions',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'PREFIX.md',ROOT/'gauge_rankone_h5/all_critical_dim4.json',ROOT/'model.py',ROOT/'suffix_mass_check/common.py']})
 atomic_json(HERE/'prefix.json',result)
 summary={k:v for k,v in result.items() if k not in ('prefix','hashes')};summary['prefix_counts']={q:len(v) for q,v in result['prefix'].items()}
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
