"""Complete-output equality and repeated target-leaf timing for frame cache."""
import json,os,signal,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,lift_coordinate,reduced_basis,span,weight
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,identification_orbits,standard
from fixed_kernel.verify_generic import native
from fixed_kernel.symmetry import suffix_slice
from gauge_rankone_h5.common import normalize_word
from kernel_y_parity.solve import solve as old
from kernel_y_cache.solve import solve as new,standard as cached_frame
from mixed_frames.stage3.common import atomic_json,sha

def main():
 records=[];start=time.monotonic()
 def pair(name,model,H,xs,all=False):
  a=old(model,H,xs,find_all=all);b=new(model,H,xs,find_all=all)
  assert {k:v for k,v in a.items() if k!='seconds'}=={k:v for k,v in b.items() if k!='seconds'}
  assert a['status'] in ('COMPLETE','FOUND')
  if b['status']=='FOUND':
   columns=b['solutions'][0]['model_columns'];assert model.criterion(columns)
   if (model.n,model.k,model.d)==(17,3,6):
    S,D=model.construct(columns);candidate=dict(status='FOUND_REQUIRES_PHYSICAL_CHECK',C=model.C,H=H,X=xs,result=b,S=S,D=D)
    atomic_json(HERE/'target_candidate.json',candidate)
    assert len(reduced_basis(S))==14 and len(reduced_basis(D))==20
    candidate['physical_dual_minimum']=min(weight(w,17) for w in span(D) if w)
    assert candidate['physical_dual_minimum']>=6
    candidate['status']='FOUND_PHYSICALLY_VERIFIED';atomic_json(HERE/'target_candidate.json',candidate)
    raise RuntimeError('Target witness persisted; inspect kernel_y_cache/target_candidate.json')
  records.append(dict(name=name,old_seconds=a['seconds'],new_seconds=b['seconds'],status=b['status'],nodes=b['nodes'],solutions=len(b['solutions'])))
 classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
 control=json.loads((ROOT/'fixed_kernel/control.json').read_text());model=ExtensionModel(8,2,3,4,[240],control['C']);data=quotient_input(model,control['N']);H,_,cs=kernel_contexts(control['H'],4,data['Gamma'])
 for c in cs:
  stats,leaves=native(model,data,c);assert stats['status']=='COMPLETE'
  for i,xs in enumerate(leaves[:8]):pair(f'small_{i}',model,H,[lift_coordinate(x,data['units']) for x in xs],True)
 before=cached_frame(8).record()
 for control in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records']:
  model=ExtensionModel(17,3,3,8,[255<<8],classes[control['cls']]['representative']['generators']);xs=standard(8).adapted_columns(model,control['known_columns'])[:7]
  pair(f'positive17_g{control["g"]}',model,control['H'],xs)
 C=classes[0]['representative']['generators'];model=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(model,[]);row=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'][0]
 H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']];reps=suffix_slice(model,[],data,c)['first_representatives']
 pool=int((ROOT/'radical_projection/pool.pid').read_text());assert os.getpgid(pool)==pool
 assert '-m radical_projection.pool_run' in subprocess.check_output(['ps','-p',str(pool),'-o','command='],text=True)
 os.killpg(pool,signal.SIGSTOP)
 try:
  with tempfile.TemporaryDirectory() as directory:
   g=Path(directory)/'group';g.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
   proc=subprocess.run([str(ROOT/'kernel_z_joint/z_enum_joint'),'--gamma','--avoid-radical','--first-orbits','--kernel-data','--prefix-group',str(g),'--mitm4','65536','--seconds','2'],input=payload(model,data,reps,c),capture_output=True,text=True,check=True,timeout=15)
  for i,line in enumerate(proc.stdout.splitlines()[:24]):
   xs=[lift_coordinate(int(x),data['units']) for x in line.split()[1:]];pair(f'target_leaf_{i}',model,H,xs)
   assert records[-1]['status']=='COMPLETE','Inspect possible target candidate; never discard a FOUND result'
 finally:os.killpg(pool,signal.SIGCONT)
 assert len([r for r in records if r['name'].startswith('target')])>=1
 assert before==cached_frame(8).record()
 info=cached_frame.cache_info();assert info.hits>0
 report=dict(status='PASS',records=records,cache_info=info._asdict(),seconds=time.monotonic()-start,scope='Exact complete-output/node equality; bounded native generation provides individual target Z leaves, not a full target-tree exclusion',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'build.py',HERE/'solve.py',HERE/'PROOF.md',ROOT/'kernel_y_parity/solve.py',ROOT/'kernel_y_parity/verification.json',ROOT/'fixed_kernel/common.py',ROOT/'two_site.py',ROOT/'model.py',ROOT/'kernel_z_joint/z_enum_joint']})
 atomic_json(HERE/'verification.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','hashes')}),flush=True)

if __name__=='__main__':main()
