"""Small exhaustive backend checks with independent direct pair enumeration."""
from pathlib import Path
import json,random,subprocess,sys,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare import table,write_input
BASE=Path(__file__).resolve().parent


def verify_case(flags,form):
    expected=[]
    for a in range(1,16):
        for b in range(a+1,16):
            if a^b<=b:continue
            if not all(flags[v]&1 for v in [a,b,a^b]):continue
            if (form[a]&b).bit_count()%2:continue
            expected.append([a,b])
    data=dict(flags=flags,form_rows=form,needs=[1,1]);p=BASE/'test_input.txt';write_input(p,data,2,0,0,0)
    run=subprocess.run([str(BASE/'enumerator'),str(p)],capture_output=True,text=True,check=True)
    events=[json.loads(l) for l in run.stdout.splitlines()];actual=[r['Z'] for r in events if r['status']=='LEAF']
    assert actual==expected,(actual,expected)
    assert events[-1]['status']=='COMPLETE' and events[-1]['leaves']==len(expected)
    return len(expected)


def main():
    # Standard two-qubit symplectic form in x|z physical labels.
    form=[((v&3)<<2)|(v>>2) for v in range(16)]
    rng=random.Random(51736);checked=leaves=0
    for case in range(100):
        flags=[0]+[rng.randrange(4) for _ in range(15)]
        leaves+=verify_case(flags,form);checked+=1
    model,data=table([],4,6,0,3,3)
    positive_leaves=verify_case(data['flags'],data['form_rows']);assert positive_leaves==15
    # Every symplectic isometry is a genuine graph stabilizer; all nonzero
    # prefix costs>=2 and all nonzero suffix P weights>=1 here.
    assert min(model.cost[1:])==2 and model.C_distance is None
    p=BASE/'limit_input.txt';write_input(p,data,2,1,0,0)
    events=[json.loads(l) for l in subprocess.check_output([str(BASE/'enumerator'),str(p)],text=True).splitlines()]
    assert events[-1]['status']=='LIMIT' and not events[-1]['Z_excluded']
    report=dict(status='PASS',random_small_cases=checked,random_canonical_leaves=leaves,
                rank_two_positive_control=dict(n=6,k=0,d=3,m=4,Z_leaves=positive_leaves),
                strict_node_limit_test=True,
                scope='Exhaustive backend comparison on16 coefficient labels; same prepared model '
                      'for the positive control, no full target17 extension claim.',
                sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [Path(__file__),BASE/'enumerator.cpp',BASE/'enumerator',BASE/'prepare.py']})
    (BASE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
