"""Complete physical N-plane orbits for q3,t2 in the hull-three column."""
import json,sys,time
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import gram,lift_coordinate,nullspace,reduced_basis,span,weight
from m8r1.suffix_symmetries import automorphisms
from mixed_frames.stage3.common import atomic_json,sha

def screen(N):
 words=[w for w in span(N) if w];weights=[weight(w,9) for w in words];supports=[(w&511)|(w>>9) for w in words]
 if sum(w<8 for w in weights)>2:return 'three_light'
 if any(weights[i]==8 and any(j!=i and supports[j]&~supports[i]==0 for j in range(3)) for i in range(3)):return 'exclusive_rank_two_support'
 return None

def main():
 start=time.monotonic();source=ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json';rows=json.loads(source.read_text())['classes'];records=[];counts=Counter()
 for ci,row in enumerate(rows):
  if row['representative']['hull']!=3:continue
  C=row['representative']['generators'];rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),5)];assert len(rad)==3
  words=sorted(w for w in span(rad) if w);planes=sorted({tuple(sorted((a,b,a^b))) for a in words for b in words if a!=b});assert len(planes)==7
  maps,order=automorphisms(C,9);lookup={p:i for i,p in enumerate(planes)};parent=list(range(7))
  def find(i):
   while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
   return i
  for apply in maps:
   image=[lookup[tuple(sorted(apply(w) for w in plane))] for plane in planes];assert sorted(image)==list(range(7))
   for i,j in enumerate(image):
    a,b=find(i),find(j)
    if a!=b:parent[max(a,b)]=min(a,b)
  groups={}
  for i in range(7):groups.setdefault(find(i),[]).append(i)
  orbits=[]
  for i,members in sorted(groups.items()):
   reason=screen(planes[i][:2]);assert all(screen(planes[j][:2])==reason for j in members)
   orbits.append(dict(N=list(planes[i][:2]),members=[list(planes[j]) for j in members],screen=reason));counts[reason or 'remaining']+=1;counts['physical_'+str(reason or 'remaining')]+=len(members)
  records.append(dict(cls=ci,C=C,physical_group_order=order,orbits=orbits))
 assert len(records)==5995 and sum(len(o['members']) for r in records for o in r['orbits'])==5995*7
 assert counts['physical_remaining']==20664 and counts['physical_three_light']==14018 and counts['physical_exclusive_rank_two_support']==7283
 profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text());prefix=[dict(profile_index=i,**r) for i,r in enumerate(profiles['records']) if (r['q'],r['t'])==(3,2)]
 assert len(prefix)==5 and len({r['index'] for r in prefix})==5
 report=dict(status='PASS_RELATIVE_TO_CLASSIFICATIONS',records=records,prefix=prefix,counts=dict(counts),seconds=time.monotonic()-start,scope='Complete C,N,H plan for ONLY q3,t2,hull3 exclusive rank-one column',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),source,ROOT/'kernel_capacity/profiles.json',ROOT/'kernel_capacity/prepare.py',ROOT/'gauge_rankone_h5/all_critical_dim4.json',ROOT/'r1_column/Q4_KERNEL.md',ROOT/'rankone_recenter/screen.py',ROOT/'m8r1/suffix_symmetries.py',ROOT/'model.py']})
 atomic_json(HERE/'hull3_plan.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','prefix','hashes')}),flush=True)
if __name__=='__main__':main()
