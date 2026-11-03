"""Gate 1: byte-identical COMPLETE LEAF streams, frozen backend vs lookahead.
Frozen files are only read and executed. Writes gate1.json here."""
import hashlib,json,random,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,dual,gram,lift_coordinate,nullspace,reduced_basis,split
from fixed_suffix_search import quotient_input
from gauge_fix import supported_subspace
from gauge_rankone_h5.common import normalize_word
from fixed_kernel.common import kernel_contexts,payload,standard
from fixed_kernel.verify_large import paired_basis
FROZEN=ROOT/'fixed_kernel/z_enum_fixed';NEW=HERE/'z_enum_lookahead'
BASE=['--gamma','--avoid-radical','--first-orbits','--kernel-data']

def run(binary,model,data,c,extra,reps=None):
 if reps is None:reps=[v for v,d in enumerate(data['distances']) if d>=model.d-1]
 p=subprocess.run([str(binary)]+BASE+extra,input=payload(model,data,reps,c),capture_output=True,text=True,check=True,timeout=60)
 return p.stdout,json.loads(p.stderr)

def compare(name,model,data,c,extra,rows,need_complete=True,expect=None):
 group=HERE/'group_tmp.txt'
 for mode in (['--prefix-group',str(group)],[]) if 'group' in c else ([],):
  if mode:group.write_text(f'{len(c["group"])} {model.m-1}\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
  out0,s0=run(FROZEN,model,data,c,extra+mode)
  for la in ('2','3'):
   out1,s1=run(NEW,model,data,c,extra+mode+['--lookahead',la])
   ok=out0==out1 and s0['status']==s1['status'] and (s0['status']=='COMPLETE' or not need_complete)
   if expect is not None:ok=ok and out1.split('\n')[0]==expect
   rows.append(dict(name=name,g=c['g'],t=c['t'],group=bool(mode),extra=extra,lookahead_free=int(la),status=s1['status'],leaves=s1['leaves'],frozen_nodes=s0['nodes'],nodes=s1['nodes'],look_pruned=s1['look_pruned'],sha256=hashlib.sha256(out1.encode()).hexdigest(),identical=ok))
   assert ok,rows[-1]

def small(rows):
 a=json.loads((ROOT/'fixed_kernel/control.json').read_text());cases=[('control_8_2_3',a['H'],ExtensionModel(8,2,3,4,[240],a['C']),a['N'])]
 S=json.loads((ROOT/'m8r1/fixed_radical_verification.json').read_text())['nonzero_radical_positive_control']['generators']
 cases.append(('radical_7_1_3',[split(w,4,3)[0] for w in supported_subspace(dual(S,7),7,range(4))],ExtensionModel(7,1,3,4,[240],[7]),[]))
 S=json.loads((ROOT/'radical_projection/positive_8.json').read_text())['result']['generators'];D=dual(S,8)
 cases.append(('positive_8_2_3',[split(w,4,4)[0] for w in supported_subspace(D,8,range(4))],ExtensionModel(8,2,3,4,[240],[split(w,4,4)[1] for w in supported_subspace(D,8,range(4,8))]),[split(w,4,4)[1] for w in supported_subspace(S,8,range(4,8))]))
 for name,H,model,N in cases:
  data=quotient_input(model,N);H,_,cs=kernel_contexts(H,model.m,data['Gamma'],with_group=True)
  for c in cs:
   compare(name,model,data,c,[],rows)
   leaves=[l for l in run(FROZEN,model,data,c,[])[0].splitlines()]
   if leaves:compare(name+'_fixed',model,data,c,['--prefix',','.join(leaves[0].split()[1:])],rows,expect=leaves[0])
   for lim in ('0','1'):compare(name+'_limit'+lim,model,data,c,['--max-nodes',lim],rows,need_complete=False)
 class Toy:m=4;d=3
 for H,G in [([240,9,10,12],[0,0,0]),([240,9^32,10,12],[2,1,0])]:
  data=dict(labels=list(range(8)),distances=[0]+[3]*7,Gamma=G);_,_,cs=kernel_contexts(H,4,G,with_group=False)
  for c in cs:compare('toy_q3',Toy(),data,c,[],rows)

def large(rows):
 classes=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records']
 ref={r['hull']:r for r in json.loads((ROOT/'fixed_kernel/large_verification.json').read_text())['records']}
 for hull,t in [(3,2),(5,0)]:
  # Same construction as fixed_kernel/verify_large.py, minus its 2^20 physical check.
  ci,C=next((i,r['representative']['generators']) for i,r in enumerate(classes) if r['representative']['hull']==hull)
  M6=ExtensionModel(17,3,6,8,[255<<8],C);rad=[lift_coordinate(c,C) for c in nullspace(gram(C,9),5)];data=quotient_input(M6,rad[:2])
  prof=next(r for r in profiles if (r['q'],r['t'])==(3,t));H=[normalize_word(h,prof['p'],8) for h in prof['H']];H,_,cs=kernel_contexts(H,8,data['Gamma'])
  model=ExtensionModel(17,3,3,8,[255<<8],C);frame=standard(8);p=7;mask=127
  K=reduced_basis([coordinates(h,frame.R+frame.E)>>1 for h in H])
  sp=lambda a,b:(((a&mask)&(b>>p)).bit_count()+((b&mask)&(a>>p)).bit_count())%2
  A=nullspace([(k>>p)|((k&mask)<<p) for k in K],14);BA,g=paired_basis(A,sp);Gamma=data['Gamma'];b=len(Gamma)
  BQ,_=paired_basis([1<<j for j in range(b)],lambda a,c:(a&lift_coordinate(c,Gamma)).bit_count()%2)
  tau=[lift_coordinate(coordinates(1<<j,BQ),BA) for j in range(b)];rng=random.Random(1736+hull)
  for attempt in range(100):
   xs=[sum(sp(v,1<<i)<<j for j,v in enumerate(tau)) for i in range(14)];mc=frame.model_columns(model,[lift_coordinate(x,data['units']) for x in xs])
   if model.criterion(mc):break
   v=lift_coordinate(rng.randrange(1<<len(A)),A);tau=[w^(v if sp(w,v) else 0) for w in tau]
  assert mc==ref[hull]['known_columns'] and H==ref[hull]['H']
  ri=[lift_coordinate(lift_coordinate(r,tau),frame.E) for r in cs[0]['radical_coeff']]
  c=next(c for c in cs if all(a^b in (0,255<<8) for a,b in zip(c['radical_images'],ri)))
  want='LEAF '+' '.join(map(str,xs[:7]))
  for k in (7,6,5):
   compare(f'large17_h{hull}_prefix{k}',model,data,c,['--prefix',','.join(map(str,xs[:k]))],rows,expect=want if k==7 else None)
   out,_=run(NEW,model,data,c,['--prefix',','.join(map(str,xs[:k])),'--lookahead','3'])
   assert want in out.splitlines(),(hull,k);rows[-1]['known_leaf_retained']=True

def main():
 start=time.monotonic();rows=[];small(rows);large(rows)
 # Positive controls with partial prefixes must still emit the known leaf.
 report=dict(status='PASS',comparisons=len(rows),rows=rows,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [FROZEN,ROOT/'fixed_kernel/z_enum_fixed.c',NEW,HERE/'z_enum_lookahead.c',Path(__file__)]})
 (HERE/'gate1.json').write_text(json.dumps(report,indent=1)+'\n');(HERE/'group_tmp.txt').unlink(missing_ok=True)
 print(json.dumps({k:v for k,v in report.items() if k!='rows'}))
 for r in rows:print(r['name'],r['g'],r['t'],r['group'],r['lookahead_free'],r['status'],r['leaves'],r['frozen_nodes'],r['nodes'],r['look_pruned'])
if __name__=='__main__':main()
