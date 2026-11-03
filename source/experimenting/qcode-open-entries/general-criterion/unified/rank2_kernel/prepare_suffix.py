"""All exact N orbits in the 56-class 888 suffix census; no extensions run."""
import itertools,json,sys,time
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import gram,lift_coordinate,nullspace,reduced_basis,span,weight
from m8r1.suffix_symmetries import automorphisms
from mixed_frames.suffix_tables import SOURCE
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();rows=json.loads(SOURCE.read_text())['classes'];records=[];counts=Counter()
 for row in rows:
  ci=row['class_index'];C=row['representative']['generators'];h=row['representative']['hull'];assert h in (0,2,6)
  assert len(reduced_basis(C))==6 and min(weight(w,9) for w in span(C) if w)>=6
  rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),6)];assert len(rad)==h;q=h//2
  if q==0:spaces=[()]
  elif q==1:spaces=[(w,) for w in sorted(span(rad)) if w]
  else:
   words=sorted(w for w in span(rad) if w);spaces=sorted({tuple(sorted([a,b,c,a^b,a^c,b^c,a^b^c])) for a,b,c in itertools.combinations(words,3) if a^b!=c});assert len(spaces)==1395
  lookup={N:i for i,N in enumerate(spaces)};parent=list(range(len(spaces)));maps,order=automorphisms(C,9)
  def find(i):
   while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
   return i
  for apply in maps:
   image=[lookup[tuple(sorted(apply(w) for w in N))] for N in spaces];assert sorted(image)==list(range(len(spaces)))
   for i,j in enumerate(image):
    a,b=find(i),find(j)
    if a!=b:parent[max(a,b)]=min(a,b)
  groups={}
  for i in range(len(spaces)):groups.setdefault(find(i),[]).append(i)
  def screen(N):return 'more_than_two_light_stabilizers' if sum(weight(w,9) in (6,7) for w in N)>2 else None
  orbits=[]
  for i,members in sorted(groups.items()):
   reason=screen(spaces[i]);assert all(screen(spaces[j])==reason for j in members)
   N=reduced_basis(spaces[i]);assert len(N)==q
   orbits.append(dict(N=N,members=[list(spaces[j]) for j in members],screen=reason,spectrum=dict(Counter(weight(w,9) for w in spaces[i]))))
   counts[f'q{q}_'+('screened' if reason else 'remaining')]+=1;counts[f'q{q}_physical_'+('screened' if reason else 'remaining')]+=len(members)
  counts[f'q{q}_classes']+=1;records.append(dict(cls=ci,C=C,hull=h,q=q,physical_group_order=order,orbits=orbits))
 report=dict(status='PASS_RELATIVE_TO_SUFFIX_CLASSIFICATION',records=records,counts=dict(counts),seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'PREFIX.md',HERE/'prefix.json',SOURCE,ROOT/'model.py',ROOT/'m8r1/suffix_symmetries.py',ROOT/'weight8_split/frame_profile_reduction.md']},scope='Exact-N plan only, no extension or target exclusion')
 atomic_json(HERE/'suffix_plan.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','hashes')}),flush=True)
if __name__=='__main__':main()
