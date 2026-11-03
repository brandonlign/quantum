"""Bounded gauge first-root pilots; no incomplete search excludes a branch."""
import json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,reduced_basis
from mixed_frames.frames import mixed_frame
from gauge_rankone_h5.common import cheap_frame,input_text,normalize_word,suffix_generators,first_orbits
from gauge_general.common import quotient
from mixed_frames.stage3.common import atomic_json

def main():
 source=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];prefix=json.loads((ROOT/'gauge_prefix_candidates.json').read_text())['groups'];out=[]
 for hull in (1,3):
  row=next(r for r in source if r['representative']['hull']==hull);C=row['representative']['generators'];r=(hull+3)//2;kfix=(5-hull)//2
  cases=prefix[str(r)]['classes'][:3]
  for i,p in enumerate(cases):
   R=p['generators']
   if r==2:frame=mixed_frame(('668','688','778')[i])
   else:
    from model import span,weight
    distinguished=next(w for w in span(R) if weight(w,8)==8)
    R=reduced_basis([normalize_word(w,distinguished,8) for w in R]);frame=cheap_frame(8,R,f'h{hull}_{i}')
   model=ExtensionModel(17,kfix,6,8,frame.R,C);data=quotient(model);generators,order=suffix_generators(model,data)
   eligible=[v for v,w in enumerate(data['distances']) if w>=6-frame.costs[1]]
   reps=[o['representative'] for o in first_orbits(eligible,generators)]
   path=HERE/'pilot_input.txt';path.write_text(input_text(model,frame,data))
   binary=ROOT/'speedup_review/enumerate_688' if hull==1 and i==1 else ROOT/'mixed_frames/stage3/enumerate_z'
   cmd=[str(binary),str(path),'--prefix',str(reps[0]),'--seconds','5']
   if hull==1 and i==1:cmd+=['--canonical','1']
   result=subprocess.run(cmd,capture_output=True,text=True,check=True);rows=[json.loads(l) for l in result.stdout.splitlines()];stats=rows[-1]
   out.append(dict(hull=hull,cls=row['class_index'],prefix=i,kfix=kfix,quotient_bits=model.a,roots=len(reps),stats=stats,sample=rows[:min(10,len(rows)-1)]))
   atomic_json(HERE/'pilot.json',dict(records=out,scope='Bounded gauge first-root Z pilots only; no branch exclusions'))
   print(json.dumps({k:v for k,v in out[-1].items() if k!='sample'}),flush=True)
if __name__=='__main__':main()
