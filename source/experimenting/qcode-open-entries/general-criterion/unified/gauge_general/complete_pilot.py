"""Exact bounded Y and inverse checks on retained gauge samples."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,reduced_basis,span,weight,lift_coordinate
from mixed_frames.frames import mixed_frame
from mixed_frames.fixed_z import solve
from gauge_rankone_h5.common import cheap_frame,normalize_word
from gauge_general.common import quotient
from gauge_general.inverse import select
from mixed_frames.stage3.common import atomic_json

def main():
 source=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];prefix=json.loads((ROOT/'gauge_prefix_candidates.json').read_text())['groups'];pilot=json.loads((HERE/'pilot.json').read_text())['records'];out=[]
 for row in pilot:
  hull=row['hull'];C=source[row['cls']]['representative']['generators'];r=(hull+3)//2;i=row['prefix'];kfix=row['kfix']
  if r==2:frame=mixed_frame(('668','688','778')[i])
  else:
   R=prefix[str(r)]['classes'][i]['generators'];distinguished=next(w for w in span(R) if weight(w,8)==8);R=reduced_basis([normalize_word(w,distinguished,8) for w in R]);frame=cheap_frame(8,R,f'h{hull}_{i}')
  model=ExtensionModel(17,kfix,6,8,frame.R,C);data=quotient(model)
  for x in row['sample'][:3]:
   xs=[lift_coordinate(v,data['units']) for v in x['X']]
   result=solve(model,frame,xs,find_all=True,seconds=1,node_limit=10000,solution_limit=5)
   inverses=[]
   for sol in result['solutions']:
    S,D=model.construct(frame.model_columns(model,sol['adapted_columns']))
    a=select(S,17,8,6,8,2,seconds=1,nodes=1000000);inverses.append({k:v for k,v in a.items() if k not in ('generators','normalizer')})
   out.append(dict(hull=hull,cls=row['cls'],prefix=i,X=x['X'],Y={k:v for k,v in result.items() if k!='solutions'},solutions=len(result['solutions']),inverses=inverses))
   atomic_json(HERE/'complete_pilot.json',dict(records=out,scope='Bounded sample only; each base inverse is exact only when terminal'))
   print(json.dumps({k:v for k,v in out[-1].items() if k not in ('X','Y')}),flush=True)
if __name__=='__main__':main()
