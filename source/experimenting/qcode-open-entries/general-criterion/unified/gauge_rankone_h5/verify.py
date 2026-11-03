"""Gauge-quotient and physical positive controls, compared with the full model."""
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from model import ExtensionModel,split,span,weight,lift_coordinate
from gauge_fix import fix_prefix_kernel
from mixed_frames.fixed_z import solve
from gauge_rankone_h5.common import (BINARY,PREFIX,cheap_frame,prefix_candidates,
    quotient,input_text,atomic_json,sha)
from gauge_rankone_h5.inverse import ordered_basis,select


def enumerate_z(model,frame,data,mu):
    path=HERE/'validation_input.txt';path.write_text(input_text(model,frame,data,mu))
    p=subprocess.run([str(BINARY),str(path)],capture_output=True,text=True,check=True,timeout=30)
    rows=[json.loads(l) for l in p.stdout.splitlines()]
    assert rows[-1]['status']=='COMPLETE'
    xs=[r['X'] for r in rows[:-1]];assert len(xs)==rows[-1]['leaf_count']
    return xs


def main():
    prefixes=prefix_candidates(json.loads(PREFIX.read_text())['groups']['4']['classes'])
    physical=0
    for row in prefixes:
        f=cheap_frame(8,row['R'],f'h5_gauge_{row["prefix_index"]}')
        physical+=f.verify_physical()
        assert f.costs[1]==min(f.cost[0][1:])
    control=json.loads((ROOT/'rankone_recenter/verification.json').read_text())['positive_controls'][1]
    S=control['S'];fixed=fix_prefix_kernel(S,7,list(range(4)))
    H=[split(w,4,3)[0] for w in fixed['R_fixed']]
    C=[split(w,4,3)[1] for w in fixed['shortened_centralizer']]
    frame=cheap_frame(4,H,'gauge_positive_7')
    model=ExtensionModel(7,0,3,4,frame.R,C);data=quotient(model)
    sets=[];bases=0;inverses=0
    for mu in (0,3):
        found=set()
        for xs in enumerate_z(model,frame,data,mu):
            full=[lift_coordinate(x,data['units']) for x in xs]
            result=solve(model,frame,full,find_all=True)
            assert result['status']=='COMPLETE' and result['complete']
            for solution in result['solutions']:
                columns=solution['adapted_columns']
                assert model.criterion(frame.model_columns(model,columns))
                B=ordered_basis(model,frame,columns)
                assert min(weight(w,7) for w in span(B) if w)>=3
                inverse=select(B,7,4,2,3,4,1000)
                inverses+=inverse['status']=='FOUND';bases+=1
                found.add(tuple(columns))
        sets.append(found)
    assert sets[0]==sets[1] and sets[0] and inverses
    report=dict(status='PASS',prefix_candidates=16,physical_prefix_words=physical,
                full_small_base_maps=len(sets[0]),evaluated_bases=bases,physical_inverse_positives=inverses,
                outside_mu_prune_preserves_all_full_small_maps=True,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'common.py',HERE/'inverse.py',HERE/'inverse.cpp',HERE/'inverse',BINARY,PREFIX,ROOT/'model.py']})
    atomic_json(HERE/'verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)


if __name__=='__main__':main()
