"""Compare complete native X,Y maps with the independent exact Python solver."""
import json,sys,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from model import ExtensionModel,reduced_basis,dual,span,weight,split,lift_coordinate,embed
from gauge_fix import fix_prefix_kernel,supported_subspace
from mixed_frames.fixed_z import solve
from gauge_rankone_h5.common import cheap_frame,input_text
from gauge_general.common import quotient
from mixed_frames.stage3.common import atomic_json,sha

def native(model,frame,data,complete=False,fast=False,**options):
 path=HERE/'validation_input.txt';path.write_text(input_text(model,frame,data))
 cmd=[str(HERE/'enumerate'),str(path),'--fast-lagrangian',str(int(fast))]
 if complete:cmd+=['--symplectic-domains','1','--complete-y','1']
 for k,v in options.items():cmd+=['--'+k,str(v)]
 p=subprocess.run(cmd,text=True,capture_output=True,check=True,timeout=30);rows=[json.loads(l) for l in p.stdout.splitlines()]
 return rows[-1],rows[:-1]

def main():
 saved=json.loads((ROOT/'rankone_recenter/verification.json').read_text());S7=saved['positive_controls'][1]['S'];base7=fix_prefix_kernel(S7,7,range(4))['S_fixed']
 letters={'I':0,'X':1,'Z':2,'Y':3}
 def word(s):return sum((letters[v]&1)<<j for j,v in enumerate(s))|(sum((letters[v]>>1)<<j for j,v in enumerate(s))<<len(s))
 S5=[word(s) for s in ('XZZXI','IXZZX','XIXZZ','ZXIXZ')]
 base12=[embed(w,0,7,5) for w in base7]+[embed(0,w,7,5) for w in S5]
 cases=[(5,0,3,3,[7<<3],[])]
 for S,n,m in [(base7,7,4),(base12,12,4)]:
  R=[split(w,m,n-m)[0] for w in supported_subspace(S,n,range(m))]
  C=[split(w,m,n-m)[1] for w in supported_subspace(dual(S,n),n,range(m,n))]
  cases.append((n,n-len(S),3,m,R,C))
 records=[]
 for n,k,d,m,R,C in cases:
  model=ExtensionModel(n,k,d,m,R,C);frame=cheap_frame(m,R,f'control_{n}');data=quotient(model)
  old,L=native(model,frame,data);fast,O=native(model,frame,data,fast=True)
  assert old['status']==fast['status']=='COMPLETE' and L==O
  expected=set()
  for row in L:
   xs=[lift_coordinate(v,data['units']) for v in row['X']];a=solve(model,frame,xs,find_all=True)
   assert a['status']=='COMPLETE';expected.update(tuple(s['adapted_columns']) for s in a['solutions'])
  st,maps=native(model,frame,data,complete=True,fast=True);actual=set()
  assert st['status']=='COMPLETE' and st['y_maps']==len(maps)
  for row in maps:
   assert row['type']=='map'
   columns=[lift_coordinate(v,data['units']) for v in row['X']+row['Y']]
   assert model.criterion(frame.model_columns(model,columns));actual.add(tuple(columns))
  assert actual==expected and actual and len(actual)==len(maps)
  arcstats,arcmaps=native(model,frame,data,complete=True,fast=True,**{'arc-depth':1})
  assert arcstats['status']=='COMPLETE'
  arcset={tuple(lift_coordinate(v,data['units']) for v in row['X']+row['Y']) for row in arcmaps}
  assert arcset==expected and len(arcmaps)==len(expected)
  records.append(dict(parameters=[n,k,d],compressed_bits=len(data['Gamma']),Z_leaves=len(L),maps=len(actual),fast_rank_prunes=st['rank_prunes']))
  for key in ('nodes','leaves','seconds'):
   st,_=native(model,frame,data,complete=True,fast=True,**{key:0});assert st['status']=='LIMIT'
 # At the largest bit width, fixed prefixes and proper depth give a small
 # complete traversal that compares the old rank test to parity masks.
 class Input:pass
 model=Input();model.d=6
 class Frame:pass
 frame=Frame();frame.p=6;frame.costs=[0]+[6]*4095
 data=dict(Gamma=[1<<(i^1) for i in range(12)],distances=[0]+[6]*4095,weights=[0]+[6]*4095)
 path=HERE/'large_fixture.txt';path.write_text(input_text(model,frame,data))
 xs=[]
 for fast in (0,1):
  p=subprocess.run([str(HERE/'enumerate'),str(path),'--depth','2','--prefix','1','--fast-lagrangian',str(fast)],text=True,capture_output=True,check=True)
  r=[json.loads(l) for l in p.stdout.splitlines()];assert r[-1]['status']=='COMPLETE';xs.append(r[:-1])
 assert xs[0]==xs[1] and len(xs[0])==2046
 report=dict(status='PASS',controls=records,large_width_complete_leaves=len(xs[0]),hashes={str(p.relative_to(ROOT)):sha(p) for p in [HERE/'enumerate',HERE/'enumerate.cpp',HERE/'common.py',Path(__file__),ROOT/'mixed_frames/fixed_z.py',ROOT/'model.py']})
 atomic_json(HERE/'backend_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
