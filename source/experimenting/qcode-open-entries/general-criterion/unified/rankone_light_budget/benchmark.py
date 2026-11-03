"""Paired bounded first-root measurements; COMPLETE and LIMIT stay distinct."""
import json
from pathlib import Path
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel
from run_hull5_planes import quotient_input
from suffix_symmetries import syndrome_orbits
from rankone_light_budget.common import physical_data,payload,BINARY,OLD
from mixed_frames.stage3.common import atomic_json,sha


def measure(model,data,reps,light,mode,seconds=10):
    cmd=[str(BINARY if mode!='old' else OLD),'--gamma','--avoid-radical','--first-orbits','--seconds',str(seconds)]
    active=light if mode=='light' else None
    if mode=='light':cmd+=['--light-budget',str(2-light['baseline']),'--light-cut','8']
    if mode=='xdom':cmd+=['--xdom']
    start=time.monotonic()
    p=subprocess.run(cmd,input=payload(model,data,reps,active),text=True,capture_output=True,check=True)
    xs=[list(map(int,l.split()[1:])) for l in p.stdout.splitlines()]
    stats=json.loads(p.stderr);assert stats['leaves']==len(xs)
    return dict(mode=mode,wall_seconds=time.monotonic()-start,stats=stats,leaves=xs)


def main():
    plan=json.loads((ROOT/'rankone_recenter/plane_plan.json').read_text())
    screen=json.loads((ROOT/'rankone_recenter/column_screen.json').read_text())
    lookup={r['cls']:r for r in screen['records']}
    verified=json.loads((HERE/'verification.json').read_text())
    assert verified['status']=='PASS' and verified['hashes'][str(BINARY.relative_to(ROOT))]==sha(BINARY)
    results=[]
    for row in plan['records'][:3]:
        model=ExtensionModel(17,3,6,8,[255<<8],row['C'])
        i=next(i for i,o in enumerate(row['orbits']) if lookup[row['cls']]['orbits'][i]['reason'] is None)
        N=row['orbits'][i]['representative'];data=quotient_input(model,N)
        light=physical_data(model,N,data,8);assert light['baseline']<=2
        orbits=syndrome_orbits(model,data,N);reps=orbits['first_representatives'][:1]
        records=[]
        for mode in ('old','light','xdom'):
            r=measure(model,data,reps,light,mode);records.append(r)
            print(json.dumps(dict(cls=row['cls'],orbit=i,first=reps,mode=mode,stats=r['stats'])),flush=True)
        results.append(dict(cls=row['cls'],orbit=i,N=N,first=reps,baseline=light['baseline'],measurements=records))
        atomic_json(HERE/'benchmark.json',dict(records=results,scope='Bounded paired first-root measurements, no exclusion claim',backend_sha256=sha(BINARY)))


if __name__=='__main__':main()
