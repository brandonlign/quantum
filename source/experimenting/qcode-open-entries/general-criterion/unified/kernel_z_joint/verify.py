"""Paired full-label MITM vs compact MITM on active m8 controls."""
import json,os,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
import fixed_kernel_mitm4.gate1 as G
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,identification_orbits
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();G.HERE=HERE;G.NEW=HERE/'z_enum_joint';G.FROZEN=ROOT/'fixed_kernel_mitm4/z_enum_mitm4';G.CAP='65536';G.CHECK='2000000'
 original_run=G.run
 def run(binary,model,data,c,extra,reps=None,timeout=25):
  extra=list(extra)
  if '--mitm4' not in extra:extra+=['--mitm4',G.CAP,'--mitm4-check',G.CHECK]
  return original_run(binary,model,data,c,extra,reps,timeout)
 G.run=run;records=[]
 pool=int((ROOT/'radical_projection/pool.pid').read_text());assert os.getpgid(pool)==pool
 assert '-m radical_projection.pool_run' in subprocess.check_output(['ps','-p',str(pool),'-o','command='],text=True)
 os.killpg(pool,signal.SIGSTOP)
 try:
  G.large(records)
  classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records']
  for ri in [0,1,2,3,4,8]:
   row=profiles[ri];hull=1 if ri!=4 else 3;cls,C=next((i,r['representative']['generators']) for i,r in enumerate(classes) if r['representative']['hull']==hull)
   model=ExtensionModel(17,3,6,8,[255<<8],C)
   from model import gram,nullspace,lift_coordinate
   rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),len(C))];N=[] if row['q']==1 else rad[:row['q']-1];data=quotient_input(model,N)
   H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
   reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
   first=115 if 115 in reps else reps[0]
   G.compare(f'target_c{cls}_p{ri}_first{first}',model,data,c,['--prefix',str(first)],records,reps=reps,group='only')
  assert all(r['nodes']==r['frozen_nodes'] for r in records)
 finally:os.killpg(pool,signal.SIGCONT)
 report=dict(status='PASS',scope='Active m8 MITM paired complete positive and target-root controls, identical streams and nodes; not a target exclusion',records=records,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'PROOF.md',HERE/'build.py',HERE/'z_enum_joint.c',HERE/'z_enum_joint',ROOT/'kernel_z_compact/z_enum_compact.c',ROOT/'fixed_kernel_mitm4/patch.py',ROOT/'fixed_kernel_mitm4/z_enum_mitm4',ROOT/'fixed_kernel_mitm4/gate1.py']})
 atomic_json(HERE/'verification.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('records','hashes')}),flush=True)

if __name__=='__main__':main()
