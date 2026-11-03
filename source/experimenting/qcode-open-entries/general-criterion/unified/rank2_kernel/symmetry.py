"""Verified suffix subgroup for a fixed radical identification and Z frame."""
from model import lift_coordinate,reduced_basis
from m8r1.suffix_symmetries import syndrome_orbits
from gauge_rankone_h5.common import first_orbits
from gauge_rankone_h5.symmetry import closure

def suffix_slice(model,N,frame,data,context):
 cache=getattr(model,'_rank2_suffix_groups',None)
 if cache is None:cache={};model._rank2_suffix_groups=cache
 key=(tuple(reduced_basis(N)),tuple(data['Gamma']),tuple(data['units']))
 if key not in cache:
  full=syndrome_orbits(model,data,N);generators=full['verified_syndrome_generators']
  # These are compositions of verified physical automorphisms. Closing the
  # image first gives the whole pointwise radical stabilizer when small.
  if full['automorphism_order']<=4096:generators=closure(generators,len(data['Gamma']))
  cache[key]=(full['automorphism_order'],generators)
 order,maps=cache[key];rad=context['radical_coeff']
 generators=[S for S in maps if all(sum(((r&s).bit_count()%2)<<j for j,s in enumerate(S))==r for r in rad)]
 eligible=[v for v,d in enumerate(data['distances']) if d>=model.d-frame.costs[1] and all((r&v).bit_count()%2==((context['tags'][0]>>j)&1) for j,r in enumerate(rad))]
 orbits=first_orbits(eligible,generators)
 return dict(eligible=eligible,first_representatives=[o['representative'] for o in orbits],orbits=orbits,generators=generators,physical_CN_order=order)
