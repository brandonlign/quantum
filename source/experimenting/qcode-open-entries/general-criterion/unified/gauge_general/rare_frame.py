"""Expose every independent singleton before the body difference circuit."""
from model import coordinates,gram,lift_coordinate,reduced_basis,span,weight
from mixed_frames.frames import Frame,adapted_frame,xor_table,FRAMES

def rare_frame(name):
 body=4 if name=='668' else 6;assert name in ('668','688','778')
 R=FRAMES[name];original=adapted_frame(8,R,name);p=6;r=len(original.R)
 Z=[1<<(8+i) for i in range(body,8)]+[(1|(1<<i))<<8 for i in range(1,body-1)]
 columns=[coordinates(z,original.R+original.E[:p])>>r for z in Z]
 assert len(reduced_basis(columns))==p
 inverse=[coordinates(1<<i,columns) for i in range(p)]
 paired=[sum(((inverse[i]>>j)&1)<<i for i in range(p)) for j in range(p)]
 E=Z+[lift_coordinate(c,original.E[p:]) for c in paired]
 assert gram(E,8)==[1<<(p+i) for i in range(p)]+[1<<i for i in range(p)]
 costs=[min(weight(w^r,8) for r in span(original.R)) for w in xor_table(E)]
 f=Frame(8,original.R,E,costs,name+'_rare')
 assert [f.costs[1<<i] for i in range(p)]==[1]*(8-body)+[2]*(body-2)
 return f,body
