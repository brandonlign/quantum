"""Physical action replay for identification orbits and compressed lex tests."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,dual,gram,lift_coordinate,reduced_basis,span,split
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,identification_orbits,standard
from fixed_kernel.verify_generic import native
from gauge_fix import supported_subspace
from m8r1.suffix_symmetries import automorphisms
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];D=dual(S,8)
 H=[split(w,4,4)[0] for w in supported_subspace(D,8,range(4))];C=[split(w,4,4)[1] for w in supported_subspace(D,8,range(4,8))];N=[split(w,4,4)[1] for w in supported_subspace(S,8,range(4,8))]
 model=ExtensionModel(8,2,3,4,[240],C);data=quotient_input(model,N);H,order,cs=kernel_contexts(H,4,data['Gamma'],with_group=True);f=standard(4)
 key=lambda images:tuple(coordinates(r,f.R+f.E)>>1 for r in images)
 lookup={key(c['radical_images']):i for i,c in enumerate(cs)};leaves=[set(native(model,data,c)[1]) for c in cs];maps,physical_order=automorphisms(H,4,[240]);edges=0
 for apply in maps:
  A=[coordinates(apply(z),f.R+f.E)>>1 for z in f.E[:3]]
  inv=[coordinates(1<<i,A) for i in range(3)]
  for i,c in enumerate(cs):
   j=lookup[key([apply(r) for r in c['radical_images']])]
   transformed={tuple(lift_coordinate(u,xs) for u in inv) for xs in leaves[i]}
   assert transformed==leaves[j];edges+=len(transformed)
 orbits=identification_orbits(H,4,cs)
 # Full S8 control, with all 40320 physical site permutations. Force the Zs
 # so this tests pruning without launching a huge artificial enumeration.
 class Toy:m=8;d=6
 G=[1<<(i^1) for i in range(12)]+[0]
 data=dict(Gamma=G,labels=list(range(8192)),distances=[0]+[6]*8191)
 H,order,contexts=kernel_contexts([255<<8,255],8,G,with_group=True);c=contexts[0];assert len(c['group'])==40320
 path=HERE/'s8_group.txt';path.write_text('40320 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
 xs=[4096]+[4096^(1<<(2*i)) for i in range(6)]
 assert xs==sorted(xs+[lift_coordinate(127,xs)])[:7]
 st,L=native(Toy(),data,c,group=path,fixed=xs);assert st['status']=='COMPLETE' and L==[tuple(xs)]
 rejected=xs[:];rejected[0],rejected[1]=rejected[1],rejected[0]
 st,L=native(Toy(),data,c,group=path,fixed=rejected);assert st['status']=='COMPLETE' and L==[]
 path.unlink()
 report=dict(status='PASS',physical_identification_arcs_checked=edges,identification_orbits=orbits,all_GL_contexts=len(cs),S8_group=40320,S8_canonical_retained=True,S8_noncanonical_rejected=True,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'verify_generic.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed',ROOT/'model.py']})
 atomic_json(HERE/'symmetry_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
