"""Complete solution-set equality against frozen Y, plus full-size positives."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,dual,lift_coordinate,reduced_basis,span,split,weight
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,standard
from fixed_kernel.solve import solve as old
from kernel_y_parity.solve import solve as new
from fixed_kernel.verify_generic import native
from gauge_fix import supported_subspace
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();fixtures=[]
 a=json.loads((ROOT/'fixed_kernel/control.json').read_text());fixtures.append((a['H'],ExtensionModel(8,2,3,4,[240],a['C']),a['N']))
 S=json.loads((ROOT/'m8r1/fixed_radical_verification.json').read_text())['nonzero_radical_positive_control']['generators'];H=[split(w,4,3)[0] for w in supported_subspace(dual(S,7),7,range(4))];fixtures.append((H,ExtensionModel(7,1,3,4,[240],[7]),[]))
 S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];D=dual(S,8);H=[split(w,4,4)[0] for w in supported_subspace(D,8,range(4))];C=[split(w,4,4)[1] for w in supported_subspace(D,8,range(4,8))];N=[split(w,4,4)[1] for w in supported_subspace(S,8,range(4,8))];fixtures.append((H,ExtensionModel(8,2,3,4,[240],C),N))
 records=[]
 for H,model,N in fixtures:
  data=quotient_input(model,N);H,order,cs=kernel_contexts(H,4,data['Gamma']);count=maps=onodes=nnodes=0;otime=ntime=0
  for c in cs:
   stats,L=native(model,data,c);assert stats['status']=='COMPLETE'
   for xs in L:
    full=[lift_coordinate(x,data['units']) for x in xs];a=old(model,H,full,find_all=True);b=new(model,H,full,find_all=True)
    assert a['status']==b['status']=='COMPLETE'
    assert {tuple(s['standard_columns']) for s in a['solutions']}=={tuple(s['standard_columns']) for s in b['solutions']}
    count+=1;maps+=len(a['solutions']);onodes+=a['nodes'];nnodes+=b['nodes'];otime+=a['seconds'];ntime+=b['seconds']
  records.append(dict(parameters=[model.n,model.k,model.d],q=len(H)-1,X_cases=count,physical_maps=maps,old_nodes=onodes,new_nodes=nnodes,old_seconds=otime,new_seconds=ntime))
 large=[];source=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
 for control in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records']:
  C=source[control['cls']]['representative']['generators'];model=ExtensionModel(17,3,3,8,[255<<8],C);f=standard(8);columns=f.adapted_columns(model,control['known_columns']);r=new(model,control['H'],columns[:7]);assert r['status']=='FOUND'
  candidate=r['solutions'][0]['model_columns'];assert model.criterion(candidate);S,D=model.construct(candidate);minimum=min(weight(w,17) for w in span(D) if w);assert minimum>=3
  assert new(model,control['H'],columns[:7],node_limit=0)['status']=='LIMIT'
  large.append(dict(parameters=[17,3,3],g=control['g'],old_first_nodes=next(r['Y_nodes'] for r in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records'] if r['g']==control['g']),new_first_nodes=r['nodes'],new_seconds=r['seconds'],physical_dual_minimum=minimum))
 report=dict(status='PASS',complete_mapsets=records,large=large,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'solve.py',HERE/'PROOF.md',ROOT/'fixed_kernel/solve.py',ROOT/'fixed_kernel/common.py',ROOT/'model.py',ROOT/'two_site.py']})
 atomic_json(HERE/'verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
