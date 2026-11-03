"""Complete solution-set comparisons for symplectic complementary domains."""
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from model import ExtensionModel,lift_coordinate
from mixed_frames.frames import adapted_frame
from mixed_frames.fixed_z import solve
from mixed_frames.stage3.common import input_text as full_input,sha,atomic_json
from gauge_rankone_h5.common import cheap_frame,quotient,input_text,BINARY


def native(text,active,binary=None):
    p=HERE/'domain_fixture.txt';p.write_text(text)
    cmd=[str(binary or HERE/'enumerate'),str(p)]
    if binary is None:cmd+=['--symplectic-domains',str(int(active))]
    result=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=30)
    rows=[json.loads(l) for l in result.stdout.splitlines()]
    assert rows[-1]['status']=='COMPLETE'
    return rows[-1],[r['X'] for r in rows[:-1]]


def all_maps(model,frame,xs,units=None):
    maps=set()
    for x in xs:
        x=x if units is None else [lift_coordinate(v,units) for v in x]
        result=solve(model,frame,x,find_all=True)
        assert result['status']=='COMPLETE'
        maps.update(tuple(s['adapted_columns']) for s in result['solutions'])
    return maps


def main():
    report=dict(status='PASS',fixtures=[])
    for model in [ExtensionModel(6,0,3,4,[15,240],[]),ExtensionModel(7,1,2,4,[240,15],[24,3])]:
        frame=adapted_frame(4,model.R)
        for mu in (0,model.d):
            text=full_input(model,frame,mu)
            old,L=native(text,False,binary=BINARY)
            off,D=native(text,False)
            assert L==D and old['nodes']==off['nodes']
            on,K=native(text,True)
            assert all_maps(model,frame,L)==all_maps(model,frame,K)
            report['fixtures'].append(dict(n=model.n,hull=model.hull,mu=mu,baseline_leaves=len(L),pruned_leaves=len(K),symplectic_domain_prunes=on['symplectic_domain_prunes']))
    saved=json.loads((ROOT/'rankone_recenter/verification.json').read_text())['positive_controls'][1]
    from gauge_fix import fix_prefix_kernel
    from model import split
    fixed=fix_prefix_kernel(saved['S'],7,list(range(4)))
    H=[split(w,4,3)[0] for w in fixed['R_fixed']]
    C=[split(w,4,3)[1] for w in fixed['shortened_centralizer']]
    frame=cheap_frame(4,H,'domain_gauge_7');model=ExtensionModel(7,0,3,4,frame.R,C);data=quotient(model)
    text=input_text(model,frame,data)
    old,L=native(text,False,binary=BINARY);on,K=native(text,True)
    maps=all_maps(model,frame,L,data['units']);kept=all_maps(model,frame,K,data['units'])
    assert maps==kept and len(maps)==480
    report['fixtures'].append(dict(n=7,quotient_bits=4,baseline_leaves=len(L),pruned_leaves=len(K),full_maps=len(maps),symplectic_domain_prunes=on['symplectic_domain_prunes']))
    report['hashes']={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'enumerate.cpp',HERE/'enumerate',HERE/'common.py',BINARY,ROOT/'mixed_frames/fixed_z.py',ROOT/'model.py']}
    atomic_json(HERE/'domain_verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)


if __name__=='__main__':main()
