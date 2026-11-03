"""Exact full-H isometry coordinates for a full-support rank-two local R."""
import functools,itertools
from model import coordinates,dual,gram,lift_coordinate,nullspace,reduced_basis,span,symp,weight
from mixed_frames.frames import Frame,adapted_frame,xor_table
from m8r1.suffix_symmetries import automorphisms
from gauge_rankone_h5.symmetry import closure

def local_R(m):return reduced_basis([((1<<m)-1)<<m,(1<<m)-1])

def kernel(H,R,m):
 H=reduced_basis(H);R=reduced_basis(R);seen=list(R);K=[]
 assert all(r in span(H) for r in R) and not any(gram(H,m))
 for h in H:
  if len(reduced_basis(seen+[h]))>len(seen):seen.append(h);K.append(h)
 assert len(K)==len(H)-len(R) and len(K)<=3
 assert all((h&((1<<m)-1)) or h in span(R) for h in span(H))
 return K

def frame_for(H,m):
 R=local_R(m);base=adapted_frame(m,R,'rank2_fixed');K=kernel(H,R,m);q=len(K);tags=[sum(symp(k,z,m)<<j for j,k in enumerate(K)) for z in base.E[:base.p]]
 anchors=[];seen=[]
 for i,t in enumerate(tags):
  if len(reduced_basis(seen+[t]))>len(seen):anchors.append(i);seen.append(t)
  if len(anchors)==q:break
 assert len(anchors)==q
 order=anchors+[i for i in range(base.p) if i not in anchors];E=[base.E[i] for i in order]+[base.E[base.p+i] for i in order]
 costs=[min(weight(w^r,m) for r in span(R)) for w in xor_table(E)]
 f=Frame(m,R,E,costs,'rank2_fixed');assert f.p==m-2
 return f

def quotient_input(model,N):
 N=reduced_basis(N);ell=len(N);assert all(not symp(v,w,model.f) for v in N for w in model.P)
 full=[coordinates(v,model.P) for v in N];quotient=[]
 for i in range(model.b):
  v=1<<i
  if len(reduced_basis(full+[v]))>len(full):full.append(v);quotient.append(v)
 assert len(full)==model.b
 units=[sum(((coordinates(1<<i,full)>>(ell+j))&1)<<i for i in range(model.b)) for j in range(model.b-ell)]
 labels=xor_table(units);Q=[lift_coordinate(v,model.P) for v in quotient]
 assert all(all((label&coordinates(n,model.P)).bit_count()%2==0 for n in N) for label in labels)
 return dict(Gamma=gram(Q,model.f),units=units,labels=labels,distances=[model.distances[v] for v in labels])

def contexts(H,m,Gamma,with_group=False):
 H=reduced_basis(H);frame=frame_for(H,m);K=kernel(H,frame.R,m);q=len(K);p=frame.p;rad=nullspace(Gamma,len(Gamma));assert len(rad)==q and len(Gamma)==2*p-q
 maps,order=automorphisms(H,m,[((1<<m)-1)<<m]) if with_group else ([],None)
 # Small physical groups can be closed before stabilizing K. Filtering only
 # the original generators is safe but can miss products fixing the slice.
 if with_group and order<=4096:
  matrices=closure([[apply(1<<i) for i in range(2*m)] for apply in maps],2*m)
  assert len(matrices)==order
  maps=[lambda w,A=A:lift_coordinate(w,A) for A in matrices]
 out=[];Rspan=set(span(frame.R))
 for rows in itertools.product(range(1,1<<q),repeat=q):
  if len(reduced_basis(rows))!=q:continue
  images=[lift_coordinate(a,K) for a in rows];tags=[sum(symp(k,z,m)<<j for j,k in enumerate(images)) for z in frame.E[:p]]
  assert len(reduced_basis(tags[:q]))==q
  generators=[]
  for apply in maps:
   if {apply(r) for r in Rspan}!=Rspan:continue
   if any(apply(k)^k not in Rspan for k in images):continue
   A=[coordinates(apply(z),frame.R+frame.E)>>2 for z in frame.E[:p]];assert all(a<1<<p for a in A);generators.append(A)
  group=closure(generators,p) if with_group else [tuple(1<<i for i in range(p))]
  out.append(dict(g=q,t=0,radical_coeff=rad,radical_images=images,tags=tags+[functools.reduce(int.__xor__,tags,0)],Gram=[0]*p,group=group,physical_group_order=order))
 assert len(out)==(1,1,6,168)[q]
 return frame,out

def payload(model,frame,data,context,reps):
 arrays=[[frame.p,len(data['labels']),model.d],data['distances'],data['Gamma'],frame.costs[:1<<frame.p],[context['g'],context['t']],context['radical_coeff'],context['tags'],context['Gram'],[len(reps)],reps]
 return '\n'.join(' '.join(map(str,a)) for a in arrays)+'\n'

def identification_orbits(H,frame,contexts):
 """Outer FULL-code action may mix Z/X; inner Z lex subgroup may not."""
 key=lambda images:tuple(coordinates(k,frame.R+frame.E)>>2 for k in images)
 lookup={key(c['radical_images']):i for i,c in enumerate(contexts)};assert len(lookup)==len(contexts)
 maps,order=automorphisms(H,frame.m);Rspan=set(span(frame.R));parent=list(range(len(contexts)))
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 used=0
 for apply in maps:
  if {apply(r) for r in Rspan}!=Rspan:continue
  used+=1;image=[lookup[key([apply(k) for k in c['radical_images']])] for c in contexts];assert sorted(image)==list(range(len(contexts)))
  for i,j in enumerate(image):
   a,b=find(i),find(j)
   if a!=b:parent[max(a,b)]=min(a,b)
 groups={}
 for i in range(len(contexts)):groups.setdefault(find(i),[]).append(i)
 return dict(orbits=[dict(representative=i,members=members) for i,members in sorted(groups.items())],full_H_physical_order=order,used_generators=used)
