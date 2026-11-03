"""Real seventeen-site positives for both q3 radical types (distance THREE)."""
import json,random,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,dual,gram,lift_coordinate,nullspace,reduced_basis,span,split,symp,weight
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,standard,permute_word
from fixed_kernel.solve import solve
from fixed_kernel.verify_generic import native
from gauge_fix import supported_subspace
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

def paired_basis(basis,pair):
 remaining=list(basis);pairs=[]
 while True:
  ij=next(((i,j) for i in range(len(remaining)) for j in range(i+1,len(remaining)) if pair(remaining[i],remaining[j])),None)
  if ij is None:break
  i,j=ij;u,v=remaining[i],remaining[j];pairs.extend([u,v])
  remaining=[w^(u if pair(w,v) else 0)^(v if pair(w,u) else 0) for k,w in enumerate(remaining) if k not in (i,j)]
 out=remaining+pairs;assert len(reduced_basis(out))==len(basis)
 return out,len(remaining)

def main():
 start=time.monotonic();source=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json';classes=json.loads(source.read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];records=[]
 for hull,t in [(3,2),(5,0)]:
  ci,C=next((i,r['representative']['generators']) for i,r in enumerate(classes) if r['representative']['hull']==hull)
  M6=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),5)];N=rad[:2];data=quotient_input(M6,N)
  profile=next(r for r in profiles if (r['q'],r['t'])==(3,t));H=[normalize_word(h,profile['p'],8) for h in profile['H']];H,order,cs=kernel_contexts(H,8,data['Gamma'])
  model=ExtensionModel(17,3,3,8,[255<<8],C);frame=standard(8);p=7;mask=127
  # Work in independent standard quotient coordinates, not the new native
  # preimage construction. Directly build an isometry of presymplectic spaces.
  K=[coordinates(h,frame.R+frame.E)>>1 for h in H];K=reduced_basis(K)
  sp=lambda a,b:(((a&mask)&(b>>p)).bit_count()+((b&mask)&(a>>p)).bit_count())%2
  A=nullspace([(k>>p)|((k&mask)<<p) for k in K],14)
  BA,g=paired_basis(A,sp);Gamma=data['Gamma'];b=len(Gamma)
  BQ,gq=paired_basis([1<<j for j in range(b)],lambda a,c:(a&lift_coordinate(c,Gamma)).bit_count()%2)
  assert g==gq==3-t and len(BA)==len(BQ)==b
  tau=[lift_coordinate(coordinates(1<<j,BQ),BA) for j in range(b)]
  assert gram([lift_coordinate(v,frame.E) for v in tau],8)==Gamma
  rng=random.Random(1736+hull)
  for attempt in range(100):
   xsmall=[sum(sp(v,1<<i)<<j for j,v in enumerate(tau)) for i in range(14)]
   columns=[lift_coordinate(x,data['units']) for x in xsmall]
   mc=frame.model_columns(model,columns)
   if model.criterion(mc):break
   v=lift_coordinate(rng.randrange(1<<len(A)),A)
   tau=[w^(v if sp(w,v) else 0) for w in tau]
  else:raise AssertionError('No distance-three control after 100 exact isometries')
  radical_images=[lift_coordinate(lift_coordinate(r,tau),frame.E) for r in cs[0]['radical_coeff']]
  context=next(c for c in cs if all(a^b in (0,255<<8) for a,b in zip(c['radical_images'],radical_images)))
  st,L=native(model,data,context,fixed=xsmall[:7]);assert st['status']=='COMPLETE' and L==[tuple(xsmall[:7])]
  answer=solve(model,H,columns[:7]);assert answer['status']=='FOUND'
  for candidate in [mc,answer['solutions'][0]['model_columns']]:
   S,D=model.construct(candidate);assert len(reduced_basis(candidate))==b
   # Full 2^20 physical dual words, independently of all syndrome costs.
   minimum=min(weight(w,17) for w in span(D) if w);assert minimum>=3
   actual=[split(w,8,9)[0] for w in supported_subspace(D,17,range(8))];assert reduced_basis(actual)==H
  records.append(dict(parameters=[17,3,3],hull=hull,q=3,t=t,g=g,cls=ci,N=N,H=H,known_columns=mc,found_columns=answer['solutions'][0]['model_columns'],physical_generators=S,dual_minimum=minimum,physical_words_checked=1<<20,Z_leaves=len(L),Y_nodes=answer['nodes'],isometry_attempt=attempt+1))
 report=dict(status='PASS',records=records,seconds=time.monotonic()-start,scope='Distance-three positive controls ONLY, not [[17,3,6]] constructions or exclusions',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'solve.py',HERE/'verify_generic.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed',ROOT/'model.py',source,ROOT/'kernel_capacity/profiles.json']})
 atomic_json(HERE/'large_verification.json',report)
 summary={k:v for k,v in report.items() if k not in ('hashes','records')};summary['records']=[{k:v for k,v in r.items() if k not in ('H','known_columns','found_columns','physical_generators')} for r in records]
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
