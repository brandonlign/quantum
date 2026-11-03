"""Kernel-zero quotient inputs for hull-one, three and five gauge bases."""
from model import coordinates,gram,lift_coordinate,nullspace,reduced_basis,span,symp,weight
from mixed_frames.frames import xor_table

def quotient(model):
 N=reduced_basis([lift_coordinate(c,model.C) for c in nullspace(gram(model.C,model.f),len(model.C))])
 assert all(symp(n,p,model.f)==0 for n in N for p in model.P)
 full=[coordinates(n,model.P) for n in N];complement=[]
 for i in range(model.b):
  v=1<<i
  if len(reduced_basis(full+[v]))>len(full):full.append(v);complement.append(v)
 a=len(complement);h=len(N);assert a==model.a and len(full)==model.b
 Q=[lift_coordinate(v,model.P) for v in complement];G=gram(Q,model.f)
 assert len(reduced_basis(G))==a
 units=[sum(((coordinates(1<<i,full)>>(h+j))&1)<<i for i in range(model.b)) for j in range(a)]
 labels=xor_table(units);distances=[model.distances[v] for v in labels]
 weights=[min(weight(w^c,model.f) for c in span(model.C)) for w in xor_table(Q)]
 for i in range(model.b):
  pi=sum(((u>>i)&1)<<j for j,u in enumerate(units))
  for j in range(model.b):
   pj=sum(((u>>j)&1)<<k for k,u in enumerate(units))
   assert ((model.Gamma[i]>>j)&1)==((pi&lift_coordinate(pj,G)).bit_count()%2)
 assert all(w==distances[lift_coordinate(c,G)] for c,w in enumerate(weights))
 return dict(Gamma=G,units=units,labels=labels,distances=distances,weights=weights,Q=Q,N=N)
