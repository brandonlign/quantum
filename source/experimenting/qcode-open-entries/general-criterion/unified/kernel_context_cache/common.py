"""Reuse prefix-only data; replace ONLY the actual suffix radical basis."""
import copy
from functools import lru_cache
from model import gram,nullspace,reduced_basis
from fixed_kernel.common import kernel_contexts as original_contexts
from fixed_kernel.common import identification_orbits as original_orbits

@lru_cache(maxsize=32)
def template(H,m,with_group):
 q=len(H)-1;t=len(reduced_basis(gram(H,m)));g=q-t;b=2*(m-1)-q
 assert t in (0,2) and 0<=g<=3
 Gamma=[0]*b
 for i in range(0,b-g,2):Gamma[i]=1<<(i+1);Gamma[i+1]=1<<i
 return original_contexts(list(H),m,Gamma,with_group=with_group)

def kernel_contexts(H,m,Gamma,*,with_group=False):
 H=tuple(reduced_basis(H));q=len(H)-1;t=len(reduced_basis(gram(H,m)));g=q-t
 assert len(Gamma)==2*(m-1)-q
 rad=nullspace(Gamma,len(Gamma));assert len(rad)==g
 normalized,order,contexts=template(H,m,with_group)
 # Prefix images, tags, Gram and physical groups depend on H, never Gamma.
 # The basis AND ORDER of rad(Gamma) are the original function's actual basis.
 return list(normalized),list(order),[dict(c,radical_coeff=list(rad)) for c in contexts]

@lru_cache(maxsize=64)
def orbit_template(H,m,images):
 return original_orbits(list(H),m,[dict(radical_images=list(v)) for v in images])

def identification_orbits(H,m,contexts):
 return copy.deepcopy(orbit_template(tuple(H),m,tuple(tuple(c['radical_images']) for c in contexts)))
