"""Complete joint-prefix/suffix coverage test on a genuine gauge control."""
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel,lift_coordinate,coordinates,split
from gauge_fix import fix_prefix_kernel
from mixed_frames.fixed_z import solve
from gauge_rankone_h5.common import cheap_frame,quotient,input_text,suffix_generators,first_orbits,atomic_json,sha
from gauge_rankone_h5.symmetry import prefix_group,closure


def main():
    control=json.loads((ROOT/'rankone_recenter/verification.json').read_text())['positive_controls'][1]
    fixed=fix_prefix_kernel(control['S'],7,list(range(4)))
    H=[split(w,4,3)[0] for w in fixed['R_fixed']];C=[split(w,4,3)[1] for w in fixed['shortened_centralizer']]
    f=cheap_frame(4,H,'joint_control_7');m=ExtensionModel(7,0,3,4,f.R,C);data=quotient(m)
    pg=prefix_group(f);sg,_=suffix_generators(m,data);suffix=closure(sg,4)
    text=input_text(m,f,data);inp=HERE/'joint_input.txt';inp.write_text(text)
    mask=HERE/'joint_prefix.txt';mask.write_text(f'{len(pg["Z_group"])} 2\n'+'\n'.join(' '.join(map(str,g)) for g in pg['Z_group'])+'\n')
    eligible=[s for s,d in enumerate(data['distances']) if s and d+f.costs[1]>=3]
    roots=[o['representative'] for o in first_orbits(eligible,sg)]
    baseline=set();reduced=set()
    for first in eligible:
        p=subprocess.run([str(HERE/'enumerate'),str(inp),'--prefix',str(first)],text=True,capture_output=True,check=True)
        rows=[json.loads(s) for s in p.stdout.splitlines()];assert rows[-1]['status']=='COMPLETE'
        for row in rows[:-1]:
            xs=[lift_coordinate(x,data['units']) for x in row['X']]
            result=solve(m,f,xs,find_all=True);assert result['status']=='COMPLETE'
            baseline.update(tuple(s['adapted_columns']) for s in result['solutions'])
    for first in roots:
        cmd=[str(HERE/'enumerate'),str(inp),'--prefix',str(first),'--symplectic-domains','1','--prefix-group',str(mask)]
        p=subprocess.run(cmd,text=True,capture_output=True,check=True);rows=[json.loads(s) for s in p.stdout.splitlines()];assert rows[-1]['status']=='COMPLETE'
        for row in rows[:-1]:
            xs=[lift_coordinate(x,data['units']) for x in row['X']]
            result=solve(m,f,xs,find_all=True);assert result['status']=='COMPLETE'
            reduced.update(tuple(s['adapted_columns']) for s in result['solutions'])
    assert len(baseline)==480 and reduced
    for columns in baseline:
        small=[coordinates(c,data['units']) for c in columns];transports=[]
        for U in pg['full_lifts']:
            changed=[lift_coordinate(a,small) for a in U]
            for S in suffix:
                transformed=tuple(lift_coordinate(a,S) for a in changed)
                transports.append(transformed)
        least=min(transports,key=lambda c:c[:2])
        full=tuple(lift_coordinate(c,data['units']) for c in least)
        assert full in reduced
    # Recheck every target prefix's physical actions against its stored group.
    from gauge_rankone_h5.common import prefix_candidates,PREFIX
    stored=json.loads((HERE/'prefix_groups.json').read_text())
    checked=0
    for pr,row in zip(prefix_candidates(json.loads(PREFIX.read_text())['groups']['4']['classes']),stored):
        frame=cheap_frame(8,pr['R'],f'h5_gauge_{pr["prefix_index"]}');fresh=prefix_group(frame)
        assert fresh['Z_group']==row['Z_group'];checked+=1
    report=dict(status='PASS',baseline_maps=len(baseline),representative_maps=len(reduced),
                eligible_roots=len(eligible),suffix_roots=len(roots),prefix_Z_group_order=len(pg['Z_group']),
                target_prefix_groups_verified=checked,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'symmetry.py',HERE/'common.py',HERE/'enumerate',HERE/'enumerate.cpp',HERE/'prefix_groups.json']})
    atomic_json(HERE/'symmetry_verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='hashes'}),flush=True)


if __name__=='__main__':main()
