"""Independent complete labelled-mass verification of the (8,5,d5) cell."""
import itertools,json,sys,time
from pathlib import Path
from math import factorial
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import reduced_basis,span,weight
from suffix_mass_check.common import canonical,coset_minima,fwht
from mixed_frames.stage3.common import atomic_json,sha
DATA=ROOT/'opus_parallel/suffix9_dimension5'

def type_check(parents):
 types=[(0,0)]+[(a,0) for a in range(1,8)]+[(a,b) for a in range(1,8) for b in range(a+1,8) if b<a^b];assert len(types)==15
 masks=[sum(bool((a&c).bit_count()%2 or (b&c).bit_count()%2)<<c for c in range(1,8)) for a,b in types]
 hits=[[bool(v>>c&1) for c in range(1,8)] for v in masks];found=set();feasible=0;total=0
 for ms in itertools.combinations_with_replacement(range(15),7):
  total+=1
  if any(sum(hits[t][c] for t in ms)<5 for c in range(7)):continue
  feasible+=1;g=[]
  for bit in range(3):
   x=sum(((types[t][0]>>bit)&1)<<j for j,t in enumerate(ms));z=sum(((types[t][1]>>bit)&1)<<j for j,t in enumerate(ms));g.append(x|(z<<7))
  assert len(reduced_basis(g))==3 and min(weight(w,7) for w in span(g) if w)>=5
  found.add(canonical(g,7))
 assert found=={canonical(g,7) for g in parents} and total==116280 and feasible==2171
 return dict(multisets=total,feasible=feasible,classes=len(found),star_certificate_sets_equal=True)

def main():
 start=time.monotonic();cell=json.loads((DATA/'cell_8_5.json').read_text());parents={};orders={}
 for r,dim in [(0,5),(1,4),(2,3)]:
  c=json.loads((DATA/f'cells_d5_7_{dim}.json').read_text())[f'7,{dim}'];assert c['mass_equal'];parents[r]=c['xz'];orders[r]=[v[1] for v in c['interleaved_with_aut']]
 assert [len(parents[r]) for r in range(3)]==[1,43,32]
 independent3=type_check(parents[2]);print('dimension3 parent type enumeration passed',flush=True)
 release=json.loads((ROOT.parents[1]/'release16/census/cells.json').read_text())
 for r,dim in [(0,5),(1,4)]:assert {canonical(g,7) for g in parents[r]}=={canonical(g,7) for g in release[f'7,{dim}']}
 saved={(p['r'],p['parent']):p['labelled_ext'] for p in cell['per_parent']};records=[];lhs=0
 for r in range(3):
  for i,g in enumerate(parents[r]):
   assert len(reduced_basis(g))==5-r and min(weight(w,7) for w in span(g) if w)>=5
   key,order=canonical(g,7,True);assert order==orders[r][i]
   delta=coset_minima(g,7);D=(delta>=4).astype(np.int64);deep=np.flatnonzero(D)
   if r==0:count=1
   elif r==1:count=3*len(deep)
   else:
    spectrum=fwht(D);numerator=sum(int(v)**3 for v in spectrum);assert numerator%len(D)==0;count=numerator//len(D)
    # Independent direct ordered-pair count checks the Fourier cubic formula.
    direct=sum(int(np.count_nonzero(D[int(a)^deep])) for a in deep);assert direct==count and count%6==0
   assert count==saved[(r,i)]
   lhs+=(factorial(7)*6**7//order)*count;records.append(dict(r=r,parent=i,deep_cosets=len(deep),labelled_extensions=count,aut_order=order))
 rhs=0;seen=set()
 for i,(g,order) in enumerate(zip(cell['xz'],cell['aut_orders'])):
  assert len(reduced_basis(g))==5 and min(weight(w,8) for w in span(g) if w)>=5
  key,actual=canonical(g,8,True);assert actual==order and key not in seen;seen.add(key);rhs+=factorial(8)*6**8//actual
  if i%5000==0:print(json.dumps(dict(children_verified=i,total=len(cell['xz']))),flush=True)
 assert len(seen)==26299 and lhs==rhs==1593099450132480
 report=dict(status='PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS',parent_dimension3=independent3,parents=records,children=len(seen),labelled_mass_lhs=str(lhs),labelled_mass_rhs=str(rhs),seconds=time.monotonic()-start,hashes={str(p):sha(p) for p in [Path(__file__),HERE/'common.py',DATA/'cell_8_5.json',DATA/'cells_d5_7_3.json',DATA/'cells_d5_7_4.json',DATA/'cells_d5_7_5.json',ROOT.parents[1]/'release16/census/cells.json',ROOT/'model.py']},scope='Independent parent type enumeration, exhaustive physical coset minima, Fourier/direct labelled extension counts and a different star-graph certificate encoding; inherits certified release16 length7 dim4/5 parent completeness')
 atomic_json(HERE/'verification8.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('parents','hashes')}),flush=True)
if __name__=='__main__':main()
