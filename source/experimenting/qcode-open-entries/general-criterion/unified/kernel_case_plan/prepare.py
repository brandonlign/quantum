"""Complete residual q1/q2 C,N plan and distinct distinguished-H classes."""
import json,sys,time
from collections import Counter
from pathlib import Path
from pynauty import certificate
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import gram,lift_coordinate,nullspace,reduced_basis,span,weight
from suffix_mass_check.common import star_graph
from m8r1.suffix_symmetries import automorphisms
from mixed_frames.stage3.common import atomic_json,sha

def distinguished(H,p):
 words=span(H);W=len(words);i=words.index(p);g=star_graph(H,8)
 g.set_vertex_coloring([{i},set(range(W))-{i},set(range(W,W+24)),set(range(W+24,W+32))])
 return certificate(g)

def main():
 start=time.monotonic();source=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json';rows=json.loads(source.read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'];groups={}
 for i,r in enumerate(profiles):
  assert len(reduced_basis(r['H']))==r['q']+1 and min(weight(w,8) for w in span(r['H']) if w)>=6
  assert r['p'] in span(r['H']) and weight(r['p'],8)==8
  key=(r['q'],r['t'],distinguished(r['H'],r['p']));groups.setdefault(key,[]).append(i)
 prefix=[];prefixcounts=Counter()
 for (q,t,key),members in groups.items():
  i=members[0];r=profiles[i];prefix.append(dict(q=q,t=t,profile_index=i,equivalent_profile_indices=members,H=r['H'],p=r['p'],distinguished_certificate=key.hex()));prefixcounts[f'{q},{t}']+=1
 assert dict(prefixcounts)=={'1,0':4,'2,0':6,'2,2':2,'3,0':7,'3,2':5}
 records=[];counts=Counter()
 for ci,row in enumerate(rows):
  h=row['representative']['hull']
  if h==5:continue
  assert h in (1,3);C=row['representative']['generators'];rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),5)];assert len(rad)==h
  if h==1:
   branches=[dict(q=1,t=0,N_orbits=[dict(N=[],members=[[]])]),dict(q=2,t=2,N_orbits=[dict(N=rad,members=[rad])])]
  else:
   words=sorted(w for w in span(rad) if w);lookup={w:i for i,w in enumerate(words)};parent=list(range(7));maps,order=automorphisms(C,9)
   def find(i):
    while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
    return i
   for apply in maps:
    image=[lookup[apply(w)] for w in words];assert sorted(image)==list(range(7))
    for i,j in enumerate(image):
     a,b=find(i),find(j)
     if a!=b:parent[max(a,b)]=min(a,b)
   orbits={}
   for i,w in enumerate(words):orbits.setdefault(find(i),[]).append(w)
   branches=[dict(q=2,t=0,N_orbits=[dict(N=[words[i]],members=[[w] for w in members]) for i,members in sorted(orbits.items())])]
  for b in branches:
   key=f'{b["q"]},{b["t"]}';counts[key+'_C']+=1;counts[key+'_N_orbits']+=len(b['N_orbits']);counts[key+'_physical_N']+=sum(len(o['members']) for o in b['N_orbits'])
  records.append(dict(cls=ci,C=C,hull=h,branches=branches))
 assert counts['1,0_C']==counts['2,2_C']==26752 and counts['2,0_C']==5995 and counts['2,0_physical_N']==5995*7
 paths=[Path(__file__),source,ROOT/'kernel_capacity/profiles.json',ROOT/'kernel_capacity/prepare.py',ROOT/'suffix_mass_check/common.py',ROOT/'m8r1/suffix_symmetries.py',ROOT/'model.py',ROOT/'r1_column/Q4_KERNEL.md']
 report=dict(status='COMPLETE_PLAN_RELATIVE_TO_CLASSIFICATIONS',prefix=prefix,prefix_counts=dict(prefixcounts),records=records,counts=dict(counts),seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in paths},scope='Exhaustive residual q1/q2 C,N,H plan only; no extension search or target exclusion')
 atomic_json(HERE/'plan.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','prefix','hashes')}),flush=True)
if __name__=='__main__':main()
