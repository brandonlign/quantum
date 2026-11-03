"""Independent physical word/site-star graph and exhaustive coset minima."""
from math import factorial
import numpy as np
from pynauty import Graph,autgrp,certificate
from model import span,reduced_basis,dual,symp
from mixed_frames.frames import xor_table

def star_graph(g,n):
 words=span(g);W=len(words);adj={i:set() for i in range(W+4*n)}
 def edge(a,b):adj[a].add(b);adj[b].add(a)
 for j in range(n):
  for letter in range(3):edge(W+3*j+letter,W+3*n+j)
 for i,w in enumerate(words):
  for j in range(n):
   letter=((w>>j)&1)|(((w>>(n+j))&1)<<1)
   if letter:edge(i,W+3*j+letter-1)
 return Graph(W+4*n,adjacency_dict=adj,vertex_coloring=[set(range(W)),set(range(W,W+3*n)),set(range(W+3*n,W+4*n))])

def canonical(g,n,need_order=False):
 graph=star_graph(g,n);key=certificate(graph)
 if not need_order:return key
 result=autgrp(graph);order=round(result[1]*10**result[2]);group=factorial(n)*6**n
 assert group%order==0
 return key,order

def coset_minima(g,n):
 g=reduced_basis(g);P=dual(g,n);bits=len(P);size=1<<bits
 cols=[sum(symp(h,1<<i,n)<<j for j,h in enumerate(P)) for i in range(2*n)]
 labels=np.asarray(xor_table(cols),dtype=np.int32);words=np.arange(1<<(2*n),dtype=np.int32);mask=(1<<n)-1
 pop=np.array([v.bit_count() for v in range(1<<n)],dtype=np.int8);weights=pop[(words&mask)|(words>>n)]
 distances=np.full(size,n+1,dtype=np.int8);np.minimum.at(distances,labels,weights)
 assert np.all(distances<=n) and distances[0]==0
 assert set(map(int,words[labels==0]))==set(span(g))
 return distances

def fwht(a):
 out=np.asarray(a,dtype=np.int64).copy();h=1
 while h<len(out):
  blocks=out.reshape(-1,2*h);u=blocks[:,:h].copy();v=blocks[:,h:].copy();blocks[:,:h]=u+v;blocks[:,h:]=u-v;h*=2
 return out
