"""Fixed complete local kernel, with every radical identification retained."""
import functools,itertools
from model import coordinates,dual,gram,lift_coordinate,nullspace,reduced_basis,span,symp,weight
from mixed_frames.frames import Frame,xor_table
from m8r1.suffix_symmetries import automorphisms
from gauge_rankone_h5.symmetry import closure

def kernel(H,m):
 H=reduced_basis(H);p=((1<<m)-1)<<m;K=[];seen=[p]
 assert p in span(H) and 2<=len(H)<=4 and all(symp(p,h,m)==0 for h in H)
 for h in H:
  if len(reduced_basis(seen+[h]))>len(seen):K.append(h);seen.append(h)
 assert len(K)==len(H)-1
 return K

def project(e,K,m):return e^(K[0] if symp(e,K[1],m) else 0)^(K[1] if symp(e,K[0],m) else 0)

def standard(m):
 p=m-1;R=[((1<<m)-1)<<m];E=[1<<(m+i) for i in range(p)]+[(1<<i)|(1<<p) for i in range(p)]
 costs=[min(weight(w^r,m) for r in span(R)) for w in xor_table(E)]
 return Frame(m,R,E,costs,'fixed_kernel')

def prefix_data(H,m):
 frame=standard(m);K=kernel(H,m);low=[project(w,K,m) for w in frame.E[:m-1]];G=gram(low,m)
 assert len(reduced_basis(G))==2
 maps,order=automorphisms(H,m,[frame.R[0]]);gens=[]
 for apply in maps:
  A=[coordinates(apply(z),frame.R+frame.E)>>1 for z in frame.E[:m-1]];assert all(a<1<<(m-1) for a in A);gens.append(A)
 group=closure(gens,m-1)
 for U in group:
  transformed=[lift_coordinate(a,low) for a in U];assert gram(transformed,m)==G
 return dict(K=K,Gram=G,group=group,physical_order=order)

def payload(model,data,reps,prefix):
 arrays=[[model.m,len(data['labels']),model.d],data['distances'],data['Gamma']]
 if 'g' in prefix:arrays += [[prefix['g'],prefix['t']],prefix['radical_coeff'],prefix['tags']]
 arrays += [prefix['Gram'],[len(reps)],reps]
 return '\n'.join(' '.join(map(str,row)) for row in arrays)+'\n'

def permute_word(w,order,m):
 return sum(((w>>old)&1)<<new | ((w>>(m+old))&1)<<(m+new) for new,old in enumerate(order))

def kernel_contexts(H,m,Gamma,*,with_group=False):
 """Return a site-normalized H and ALL GL(g,2) radical identifications.

 No suffix symmetry is assumed. A context fixes the images of the exact
 basis nullspace(Gamma); prefix actions must fix those classes individually.
 """
 H=reduced_basis(H);K=kernel(H,m);p=m-1;R=((1<<m)-1)<<m
 assert all((w&((1<<m)-1)) or w in (0,R) for w in span(H))
 Kr=[lift_coordinate(c,K) for c in nullspace(gram(K,m),len(K))]
 g=len(Kr);t=len(K)-g;assert t in (0,2)
 assert len(Gamma)==2*p-len(K) and len(nullspace(Gamma,len(Gamma)))==g
 tags=[sum(((r>>i)&1)<<j for j,r in enumerate(Kr)) for i in range(m)]
 anchors=[];seen=[]
 for i in range(p):
  if len(reduced_basis(seen+[tags[i]]))>len(seen):anchors.append(i);seen.append(tags[i])
  if len(seen)==g:break
 assert len(anchors)==g
 order=anchors+[i for i in range(m) if i not in anchors]
 H=reduced_basis([permute_word(h,order,m) for h in H]);K=kernel(H,m)
 Kr=[lift_coordinate(c,K) for c in nullspace(gram(K,m),len(K))]
 pair=next(((u,v) for u in K for v in K if symp(u,v,m)),None)
 assert (pair is not None)==bool(t)
 frame=standard(m);rad=nullspace(Gamma,len(Gamma));actions=[];physical_order=None
 if with_group:
  maps,physical_order=automorphisms(H,m,[R]);actions=maps
 contexts=[]
 for rows in itertools.product(range(1,1<<g),repeat=g):
  if len(reduced_basis(rows))!=g:continue
  images=[lift_coordinate(c,Kr) for c in rows]
  tags=[sum(((r>>i)&1)<<j for j,r in enumerate(images)) for i in range(m)]
  assert len(reduced_basis(tags[:g]))==g and functools.reduce(int.__xor__,tags,0)==0
  low=[]
  for i,z in enumerate(frame.E[:p]):
   a=coordinates(tags[i],tags[:g]);v=z^lift_coordinate(a,frame.E[:g])
   assert all(not symp(v,r,m) for r in images)
   if pair:v=project(v,pair,m)
   low.append(v)
  G=gram(low,m);assert len(reduced_basis(G))==t and not any(G[:g])
  group=[[1<<i for i in range(p)]]
  if with_group:
   gens=[]
   for apply in actions:
    if any(apply(r)^r not in (0,R) for r in images):continue
    A=[coordinates(apply(z),frame.R+frame.E)>>1 for z in frame.E[:p]]
    assert all(a<1<<p for a in A);gens.append(A)
   group=closure(gens,p)
  contexts.append(dict(g=g,t=t,radical_coeff=rad,radical_images=images,tags=tags,Gram=G,group=group,physical_order=physical_order))
 assert len(contexts)==(1,1,6,168)[g]
 return H,order,contexts

def identification_orbits(H,m,contexts):
 """Joint prefix physical action on ALL radical identifications.

 This is separate from within-context prefix pruning. Its generators need
 not fix individual radical images; every orbit member is recorded.
 """
 R=((1<<m)-1)<<m;frame=standard(m)
 key=lambda images:tuple(coordinates(r,frame.R+frame.E)>>1 for r in images)
 lookup={key(c['radical_images']):i for i,c in enumerate(contexts)}
 assert len(lookup)==len(contexts)
 maps,order=automorphisms(H,m,[R]);parent=list(range(len(contexts)))
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 for apply in maps:
  image=[lookup[key([apply(r) for r in c['radical_images']])] for c in contexts]
  assert sorted(image)==list(range(len(contexts)))
  for i,j in enumerate(image):
   a,b=find(i),find(j)
   if a!=b:parent[max(a,b)]=min(a,b)
 groups={}
 for i in range(len(contexts)):groups.setdefault(find(i),[]).append(i)
 return dict(orbits=[dict(representative=i,members=members) for i,members in sorted(groups.items())],physical_group_order=order)
