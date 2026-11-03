"""Genuine pure-distance-three controls for every q=0..3 in a rank-two R."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,dual,gram,lift_coordinate,nullspace,reduced_basis,span,symp,weight
from rank2_kernel.common import local_R,kernel,contexts,quotient_input
from fixed_kernel.verify_large import paired_basis
from mixed_frames.suffix_tables import SOURCE
from mixed_frames.stage3.common import atomic_json,sha

def pauli(s):return sum((c in 'XY')<<i | (c in 'ZY')<<(len(s)+i) for i,c in enumerate(s))

def fixture_inputs():
 R4=local_R(4)
 yield 8,0,4,R4,R4,R4
 H=R4+[pauli('IXZY')];C=[pauli(s) for s in ['XZZXI','IXZZX','XIXZZ','ZXIXZ']]
 assert len(reduced_basis(C))==4 and not any(gram(C,5)) and min(weight(w,5) for w in span(C) if w)==4
 yield 9,1,4,H,C,C[:3]
 R6=local_R(6);chosen=None
 for v in dual(R6,6):
  B=reduced_basis(R6+[v])
  if len(B)!=3 or min(weight(w,6) for w in span(B) if w)<3:continue
  for u in dual(B,6):
   H=reduced_basis(B+[u])
   if len(H)==4 and min(weight(w,6) for w in span(H) if w)>=3:chosen=H;break
  if chosen is not None:break
 if chosen is None:
  # The normalizer's returned BASIS need not itself contain a good extension;
  # enumerate its small full space if the cheap first pass did not find one.
  for v in span(dual(R6,6)):
   B=reduced_basis(R6+[v])
   if len(B)!=3 or min(weight(w,6) for w in span(B) if w)<3:continue
   for u in span(dual(B,6)):
    H=reduced_basis(B+[u])
    if len(H)==4 and min(weight(w,6) for w in span(H) if w)>=3:chosen=H;break
   if chosen is not None:break
 assert chosen is not None
 yield 12,2,6,chosen,chosen,R6
 H=json.loads((HERE/'prefix.json').read_text())['prefix']['3'][0]['H'];row=next(r for r in json.loads(SOURCE.read_text())['classes'] if r['representative']['hull']==6);C=row['representative']['generators']
 yield 17,3,8,H,C,C[:3]

def make(n,k,m,H,C,N):
 H=reduced_basis(H);model=ExtensionModel(n,k,3,m,local_R(m),C);data=quotient_input(model,N);frame,cs=contexts(H,m,data['Gamma']);p=frame.p;mask=(1<<p)-1
 Ks=[coordinates(h,frame.R+frame.E)>>2 for h in H];Ks=reduced_basis(Ks)
 sp=lambda a,b:(((a&mask)&(b>>p)).bit_count()+((b&mask)&(a>>p)).bit_count())%2
 A=nullspace([(c>>p)|((c&mask)<<p) for c in Ks],2*p);BA,g=paired_basis(A,sp);G=data['Gamma'];b=len(G)
 BQ,gq=paired_basis([1<<j for j in range(b)],lambda a,c:(a&lift_coordinate(c,G)).bit_count()%2);assert g==gq==len(Ks)
 tau=[lift_coordinate(coordinates(1<<j,BQ),BA) for j in range(b)]
 assert gram([lift_coordinate(v,frame.E) for v in tau],m)==G
 small=[sum(sp(v,1<<i)<<j for j,v in enumerate(tau)) for i in range(2*p)]
 adapted=[lift_coordinate(v,data['units']) for v in small];columns=frame.model_columns(model,adapted);assert model.criterion(columns)
 images=[lift_coordinate(lift_coordinate(c,tau),frame.E) for c in cs[0]['radical_coeff']]
 Rspan=set(span(frame.R));ci=next(i for i,c in enumerate(cs) if all(a^b in Rspan for a,b in zip(c['radical_images'],images)))
 S,D=model.construct(columns);minimum=min(weight(w,n) for w in span(D) if w);assert minimum>=3
 from gauge_fix import supported_subspace
 from model import split
 assert reduced_basis([split(w,m,n-m)[0] for w in supported_subspace(D,n,range(m))])==H
 assert reduced_basis([split(w,m,n-m)[1] for w in supported_subspace(S,n,range(m,n))])==reduced_basis(N)
 return model,frame,data,cs,dict(parameters=[n,k,3],m=m,H=H,C=C,N=reduced_basis(N),context_index=ci,small_columns=small,model_columns=columns,S=S,physical_dual_minimum=minimum,words_checked=1<<len(D),scope='Genuine distance-three control only')

def main():
 start=time.monotonic();records=[]
 for fixture in fixture_inputs():records.append(make(*fixture)[-1])
 result=dict(status='PASS',records=records,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'prefix.json',ROOT/'model.py',ROOT/'fixed_kernel/verify_large.py',SOURCE]})
 atomic_json(HERE/'controls.json',result);print(json.dumps(dict(status='PASS',records=[{k:v for k,v in r.items() if k in ('parameters','m','physical_dual_minimum','words_checked','context_index')} for r in records],seconds=result['seconds'])),flush=True)
if __name__=='__main__':main()
