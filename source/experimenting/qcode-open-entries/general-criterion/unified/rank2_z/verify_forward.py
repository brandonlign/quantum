"""Independent exhaustive checks of early X0-domain pruning."""
from pathlib import Path
import sys,subprocess,json,random,itertools,time,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent))
from prepare import table
BASE=Path(__file__).resolve().parent


def physical_cost(m,x,z):
    mask=(1<<m)-1;zphys=z|((z.bit_count()%2)<<(m-1))
    return min(((x^rx)|(zphys^rz)).bit_count() for rx in [0,mask] for rz in [0,mask])


def write(path,flags,form,costmask,p):
    needs=[1 if min(t+(t%2),p+2-t-(t%2))==2 else 2 for t in range(1,p+1)]
    rows=[f'{len(flags)} {p} 0 0 0',' '.join(map(str,needs))]
    rows += [f'{f} {g} {c}' for f,g,c in zip(flags,form,costmask)]
    path.write_text('\n'.join(rows)+'\n')


def reference(flags,form,mask,p):
    high=[v for v in range(1,len(flags)) if flags[v]&1];expected=[];m=p+2;checked=0
    costs=[physical_cost(m,1|(1<<(m-1)),z) for z in range(1<<p)]
    for choice in itertools.combinations(high,p):
        checked+=1;dep=0
        for v in choice:dep ^= v
        if dep<=choice[-1] or not flags[dep]&1:continue
        if any(not flags[a^b]&1 or (form[a]&b).bit_count()%2 for a,b in itertools.combinations(choice,2)):continue
        translates=[0]
        for v in choice:translates += [z^v for z in translates]
        if any(not flags[translates[z]]&(1 if min(z.bit_count()+z.bit_count()%2,m-z.bit_count()-z.bit_count()%2)==2 else 2) for z in range(1,1<<p)):continue
        beta_found=False
        for beta in high:
            if any((form[a]&beta).bit_count()%2 != int(j!=0) for j,a in enumerate(choice)):continue
            if all(mask[beta^az]&(1<<costs[z]) for z,az in enumerate(translates)):
                beta_found=True;break
        if beta_found:expected.append(list(choice))
    return expected,checked


def check(flags,form,mask,p):
    path=BASE/'forward_test.txt';write(path,flags,form,mask,p)
    events=[json.loads(l) for l in subprocess.check_output([str(BASE/'enumerator_unary'),str(path)],text=True).splitlines()]
    actual=[e['Z'] for e in events if e['status']=='LEAF'];assert events[-1]['status']=='COMPLETE'
    expected,checked=reference(flags,form,mask,p)
    assert actual==expected,(p,actual,expected)
    return len(actual),checked


def main():
    start=time.monotonic();rng=random.Random(91736);random_leaves=0;choices=0
    form=[((v&3)<<2)|(v>>2) for v in range(16)]
    for _ in range(100):
        mask=[rng.randrange(32) for v in range(16)]
        flags=[(int(bool(v&(1<<2)))|2*int(bool(v&(1<<4)))) for v in mask];flags[0]=0
        count,seen=check(flags,form,mask,2);random_leaves+=count;choices+=seen
    controls=[]
    for d in [3,4]:
        model,data=table([],4,6,0,d,d)
        masks=[sum((1<<c) for c in range(5) if delta>=d-c and w>=d-c) for delta,w in zip(data['coset_distances'],data['suffix_weights'])]
        count,seen=check(data['flags'],data['form_rows'],masks,2);assert count>0
        controls.append(dict(n=6,k=0,m=4,d=d,Z_leaves=count))
    # A depth-six fixture contains a prescribed isotropic Z space plus one
    # consistent X0 image. Full independent combination+beta enumeration follows.
    p=6;m=8;beta=62<<6;mask=[0]*4096
    for z in range(1,64):
        cost=min(z.bit_count()+z.bit_count()%2,8-z.bit_count()-z.bit_count()%2)
        mask[z]=sum(1<<c for c in range(cost,9))
    for z in range(64):
        cost=physical_cost(8,129,z);mask[beta^z]=sum(1<<c for c in range(cost,9))
    flags=[int(bool(v&4))|2*int(bool(v&16)) for v in mask];flags[0]=0
    form=[((v&63)<<6)|(v>>6) for v in range(4096)]
    depth6,combinations=check(flags,form,mask,6);assert depth6>0
    report=dict(status='PASS',random_small_cases=100,random_leaves=random_leaves,
        small_combinations_checked=choices,positive_controls=controls,
        depth6_canonical_leaves=depth6,depth6_combinations_checked=combinations,
        seconds=time.monotonic()-start,
        scope='Exhaustive independent Z-combination and X0-image enumeration; '
              'depth6 fixture is a backend control rather than a target17 code.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [Path(__file__),BASE/'enumerator_unary.cpp',BASE/'enumerator_unary']})
    (BASE/'forward_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
