"""Artificial d6 distance tables stress singleton5 versus pair4 boundaries."""
import hashlib,json,random,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload,standard
from mixed_frames.stage3.common import atomic_json,sha

def main():
 rng=random.Random(2026101408);records=[]
 classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
 for control in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records']:
  C=classes[control['cls']]['representative']['generators'];model=ExtensionModel(17,3,6,8,[65280],C);data=quotient_input(model,control['N'])
  frame=standard(8);known=frame.adapted_columns(model,control['known_columns']);xs=[coordinates(v,data['units']) for v in known[:7]]
  H,_,contexts=kernel_contexts(control['H'],8,data['Gamma'])
  c=next(c for c in contexts if all((r&xs[i]).bit_count()%2==(c['tags'][i]>>j)&1 for i in range(7) for j,r in enumerate(c['radical_coeff'])))
  last=0
  for x in xs:last^=x
  for trial in range(3):
   # The real Z tuple has valid tags/Gram. These distances are deliberately
   # artificial: only high singleton labels are compacted; pair sums at4 remain.
   high=set(xs+[last]);extra=rng.sample([v for v in range(1,len(data['labels'])) if v not in high],4+trial*2);high.update(extra)
   delta=[4]*len(data['labels']);delta[0]=0
   for v in high:delta[v]=6
   artificial=dict(data,distances=delta);reps=sorted(high)
   for pin in [3,4]:
    args=['--gamma','--avoid-radical','--first-orbits','--kernel-data','--prefix',','.join(map(str,xs[:pin]))];outputs=[];stats=[]
    for binary in [ROOT/'fixed_kernel/z_enum_fixed',ROOT/'fixed_kernel_mitm4/z_enum_mitm4',HERE/'z_enum_joint']:
     extraargs=[] if binary.name=='z_enum_fixed' else ['--mitm4','65536','--mitm4-check','2000000']
     proc=subprocess.run([str(binary)]+args+extraargs,input=payload(model,artificial,reps,c),text=True,capture_output=True,check=True,timeout=20);outputs.append(proc.stdout);stats.append(json.loads(proc.stderr))
    assert outputs[0]==outputs[1]==outputs[2] and all(s['status']=='COMPLETE' for s in stats)
    assert 'LEAF '+' '.join(map(str,xs)) in outputs[2].splitlines()
    assert stats[1]['nodes']==stats[2]['nodes'] and stats[2]['mitm_tested']>0
    records.append(dict(g=c['g'],trial=trial,pin=pin,columns=len(high),leaves=stats[2]['leaves'],nodes=stats[2]['nodes'],mitm_tested=stats[2]['mitm_tested'],mitm_pruned=stats[2]['mitm_pruned'],mitm_checked=stats[2]['mitm_checked'],stream_sha256=hashlib.sha256(outputs[2].encode()).hexdigest()));print(json.dumps(records[-1]),flush=True)
 report=dict(status='PASS',records=records,scope='Artificial distance tables, genuine known tags/Gram and retained Z tuples; these are not physical distance-six code constructions',hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'z_enum_joint.c',HERE/'z_enum_joint',ROOT/'fixed_kernel/z_enum_fixed',ROOT/'fixed_kernel_mitm4/z_enum_mitm4',ROOT/'fixed_kernel/large_verification.json']})
 atomic_json(HERE/'threshold_verification.json',report)

if __name__=='__main__':main()
