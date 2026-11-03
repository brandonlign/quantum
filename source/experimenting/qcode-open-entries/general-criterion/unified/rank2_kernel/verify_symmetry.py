"""Retain physical positives after joint prefix/suffix slice normalization."""
import json,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,coordinates,lift_coordinate,span,weight
from rank2_kernel.common import local_R,contexts,quotient_input
from rank2_kernel.symmetry import suffix_slice
from rank2_kernel.solve import solve
from rank2_kernel.verify_z import run
from gauge_rankone_h5.symmetry import closure
from mixed_frames.stage3.common import atomic_json,sha

def main():
 start=time.monotonic();records=[]
 for row in json.loads((HERE/'controls.json').read_text())['records']:
  n,k,d=row['parameters'];model=ExtensionModel(n,k,d,row['m'],local_R(row['m']),row['C']);data=quotient_input(model,row['N']);frame,cs=contexts(row['H'],row['m'],data['Gamma'],with_group=True);c=cs[row['context_index']];ss=suffix_slice(model,row['N'],frame,data,c)
  B=closure(ss['generators'],len(data['Gamma']));xs=row['small_columns'][:frame.p]
  normalized=min(tuple(lift_coordinate(lift_coordinate(a,xs),S) for a in U) for U in c['group'] for S in B)
  # The generic runner's first mask below is built independently, so verify
  # its extra labels only weaken the actual suffix mask; the retained point
  # itself MUST be one of the true minimum orbit representatives.
  assert normalized[0] in ss['first_representatives']
  stats,L=run(model,frame,data,c,fixed=normalized,group=True);assert stats['status']=='COMPLETE' and L==[normalized]
  answer=solve(model,row['H'],frame,[lift_coordinate(x,data['units']) for x in normalized]);assert answer['status']=='FOUND'
  columns=answer['solutions'][0]['model_columns'];assert model.criterion(columns)
  S,D=model.construct(columns);minimum=min(weight(w,n) for w in span(D) if w);assert minimum>=3
  records.append(dict(parameters=row['parameters'],q=c['g'],prefix_group=len(c['group']),suffix_slice_group=len(B),eligible=len(ss['eligible']),first_orbits=len(ss['first_representatives']),physical_dual_minimum=minimum,known_positive_transported=True))
 report=dict(status='PASS',records=records,seconds=time.monotonic()-start,hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'symmetry.py',HERE/'solve.py',HERE/'verify_z.py',HERE/'controls.json',ROOT/'rank2_kernel_backend/z_enum_cost.c',ROOT/'rank2_kernel_backend/z_enum_cost',ROOT/'m8r1/suffix_symmetries.py']})
 atomic_json(HERE/'symmetry_verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)
if __name__=='__main__':main()
