"""Full depth-six canonical circuit check in a prescribed isotropic space."""
from pathlib import Path
import sys,json,itertools,subprocess,time,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare import write_input
BASE=Path(__file__).resolve().parent

def main():
    start=time.monotonic();flags=[0]*4096
    for v in range(1,64):
        size=v.bit_count();flags[v]=2 | int(size in [1,2,5,6])
    form=[((v&63)<<6)|(v>>6) for v in range(4096)]
    data=dict(flags=flags,form_rows=form,needs=[1,1,2,2,1,1]);path=BASE/'depth6_input.txt'
    write_input(path,data,6,0,0,0)
    events=[json.loads(l) for l in subprocess.check_output([str(BASE/'enumerator'),str(path)],text=True).splitlines()]
    actual=[r['Z'] for r in events if r['status']=='LEAF'];assert events[-1]['status']=='COMPLETE'
    high=[v for v in range(1,64) if flags[v]&1];expected=[];choices=0
    for choice in itertools.combinations(high,6):
        choices+=1;dependent=0
        for v in choice:dependent ^= v
        if dependent<=choice[-1] or not flags[dependent]&1:continue
        if any(not flags[a^b]&1 for a,b in itertools.combinations(choice,2)):continue
        valid=True
        for index in range(1,64):
            value=0
            for j in range(6):
                if index>>j&1:value ^= choice[j]
            need=1 if index.bit_count() in [1,2,5,6] else 2
            if not flags[value]&need:valid=False;break
        if valid:expected.append(list(choice))
    assert actual==expected and len(actual)==8
    report=dict(status='PASS',depth=6,combinations_checked=choices,canonical_leaves=8,
                node_count=events[-1]['nodes'],seconds=time.monotonic()-start,
                scope='Exhaustive independent combination check of all depth-six subset and '
                      'canonical dependence conditions in a prescribed isotropic space; '
                      'this is a backend control, not a target17 code.',
                sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [Path(__file__),BASE/'enumerator.cpp',BASE/'enumerator']})
    (BASE/'depth6_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
