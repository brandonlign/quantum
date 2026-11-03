"""Independent joint-orbit coverage check for prefix and first-label pruning."""
import json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,dual,lift_coordinate,split
from gauge_fix import supported_subspace
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload
from fixed_kernel.verify_generic import native
from fixed_kernel.symmetry import suffix_slice
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();fixtures=[]
 a=json.loads((HERE/'control.json').read_text());fixtures.append((a['H'],ExtensionModel(8,2,3,4,[240],a['C']),a['N']))
 S=json.loads((ROOT/'m8r1/fixed_radical_verification.json').read_text())['nonzero_radical_positive_control']['generators'];H=[split(w,4,3)[0] for w in supported_subspace(dual(S,7),7,range(4))];fixtures.append((H,ExtensionModel(7,1,3,4,[240],[7]),[]))
 S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];D=dual(S,8);H=[split(w,4,4)[0] for w in supported_subspace(D,8,range(4))];C=[split(w,4,4)[1] for w in supported_subspace(D,8,range(4,8))];N=[split(w,4,4)[1] for w in supported_subspace(S,8,range(4,8))];fixtures.append((H,ExtensionModel(8,2,3,4,[240],C),N))
 records=[]
 for H,model,N in fixtures:
  data=quotient_input(model,N);H,order,cs=kernel_contexts(H,4,data['Gamma'],with_group=True)
  for ci,c in enumerate(cs):
   st,L=native(model,data,c);assert st['status']=='COMPLETE';lookup={x:i for i,x in enumerate(L)};parent=list(range(len(L)))
   def find(i):
    while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
    return i
   ss=suffix_slice(model,N,data,c)
   actions=[lambda xs,A=A:tuple(lift_coordinate(a,xs) for a in A) for A in c['group']]+[lambda xs,S=S:tuple(lift_coordinate(x,S) for x in xs) for S in ss['generators']]
   for action in actions:
    image=[lookup[action(xs)] for xs in L];assert sorted(image)==list(range(len(L)))
    for i,j in enumerate(image):
     a,b=find(i),find(j)
     if a!=b:parent[max(a,b)]=min(a,b)
   path=HERE/'joint_group.txt';path.write_text(f'{len(c["group"])} 3\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
   cmd=[str(HERE/'z_enum_fixed'),'--gamma','--avoid-radical','--first-orbits','--kernel-data','--prefix-group',str(path)]
   proc=subprocess.run(cmd,input=payload(model,data,ss['first_representatives'],c),text=True,capture_output=True,check=True);s=json.loads(proc.stderr);kept=[tuple(map(int,l.split()[1:])) for l in proc.stdout.splitlines()];assert s['status']=='COMPLETE' and len(kept)==s['leaves']
   assert {find(lookup[x]) for x in kept}=={find(i) for i in range(len(L))}
   records.append(dict(parameters=[model.n,model.k,model.d],context=ci,plain=len(L),joint_pruned=len(kept),joint_orbits=len({find(i) for i in range(len(L))}),eligible=len(ss['eligible']),first_orbits=len(ss['first_representatives'])))
 report=dict(status='PASS',records=records,seconds=time.monotonic()-start,scope='Complete small-tree combined group-orbit coverage, no target exclusions',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'symmetry.py',HERE/'common.py',HERE/'verify_generic.py',HERE/'z_enum_fixed.c',HERE/'z_enum_fixed',ROOT/'m8r1/suffix_symmetries.py',ROOT/'gauge_rankone_h5/common.py']})
 atomic_json(HERE/'joint_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
