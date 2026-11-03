"""Complete physical-oracle validation of singleton/pair-X forward domains."""
from pathlib import Path
import sys,json,random,subprocess,time,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'rank2_z'))
from verify_forward import write,physical_cost
from prepare import table
BASE=Path(__file__).resolve().parent


def reference(flags,form,mask,p):
    high=[a for a in range(1,len(flags)) if flags[a]&1];expected=[];pure_leaves=0;m=p+2
    sectors=[s for s in range(1,1<<p) if physical_cost(m,s|((s.bit_count()%2)<<(m-1)),0)==2]
    costs={s:[physical_cost(m,s|((s.bit_count()%2)<<(m-1)),z) for z in range(1<<p)] for s in sectors}
    def dfs(chosen,candidates):
        nonlocal pure_leaves
        if len(chosen)==p:
            dep=0
            for a in chosen:dep ^= a
            if dep<=chosen[-1] or not flags[dep]&1:return
            translates=[0]
            for a in chosen:translates += [z^a for z in translates]
            for z in range(1,1<<p):
                support=z.bit_count()+z.bit_count()%2;cost=min(support,m-support)
                if not flags[translates[z]]&(1 if cost==2 else 2):return
            pure_leaves+=1
            for sector in sectors:
                found=False
                for beta in high:
                    if any((form[a]&beta).bit_count()%2 != ((sector.bit_count()%2)^((sector>>j)&1)) for j,a in enumerate(chosen)):continue
                    if all(mask[beta^az]&(1<<costs[sector][z]) for z,az in enumerate(translates)):
                        found=True;break
                if not found:return
            expected.append(chosen[:]);return
        for pos,a in enumerate(candidates):
            remainder=[b for b in candidates[pos+1:] if flags[a^b]&1 and not (form[a]&b).bit_count()%2]
            dfs(chosen+[a],remainder)
    dfs([],high);return expected,pure_leaves


def check(flags,form,mask,p):
    path=BASE/'pair_domains_test.txt';write(path,flags,form,mask,p)
    events=[json.loads(l) for l in subprocess.check_output([str(BASE/'enumerator_physical_pairs'),str(path)],text=True).splitlines()]
    assert events[-1]['status']=='COMPLETE'
    actual=[e['Z'] for e in events if e['status']=='LEAF'];expected,pure=reference(flags,form,mask,p)
    assert actual==expected,(p,actual,expected)
    return len(actual),pure


def main():
    start=time.monotonic();rng=random.Random(917368);total=0
    form=[((a&3)<<2)|(a>>2) for a in range(16)]
    for _ in range(200):
        mask=[rng.randrange(32) if _<100 else (31 if rng.random()<.8 else rng.randrange(32)) for j in range(16)]
        flags=[int(bool(v&4))|2*int(bool(v&16)) for v in mask];flags[0]=0
        n,_=check(flags,form,mask,2);total+=n
    controls=[]
    for d in [3,4]:
        _,data=table([],4,6,0,d,d)
        masks=[sum(1<<c for c in range(5) if delta>=d-c and lam>=d-c) for delta,lam in zip(data['coset_distances'],data['suffix_weights'])]
        n,_=check(data['flags'],data['form_rows'],masks,2);assert n>0
        controls.append(dict(n=6,k=0,m=4,d=d,Z_leaves=n))
    # Prescribed full isometry fixture: P's form and admissible cost levels
    # match the physical prefix quotient itself. The identity therefore survives.
    p=6;m=8;form=[];mask=[0]*4096
    B=lambda a:a^(63 if a.bit_count()%2 else 0)
    for a in range(4096):
        z=a&63;x=a>>6;form.append(B(x)|(B(z)<<6))
        if not a:continue
        physical_x=x|((x.bit_count()%2)<<7);c=physical_cost(8,physical_x,z)
        mask[a]=sum(1<<level for level in range(c,9))
    flags=[int(bool(v&4))|2*int(bool(v&16)) for v in mask];flags[0]=0
    depth6,pure=check(flags,form,mask,6);assert depth6==pure==24
    report=dict(status='PASS',random_small_cases=200,random_surviving_Z_maps=total,
        positive_controls=controls,depth6_isometry_fixture_leaves=depth6,
        independent_depth6_pure_leaves=pure,seconds=time.monotonic()-start,
        scope='Independent clique enumeration and physical-cost/exhaustive X-image checks for singleton and pair sectors; '
              'depth6 fixture is an abstract isometry control, not a target17 code.',
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                [Path(__file__),BASE/'enumerator_physical_pairs.cpp',BASE/'enumerator_physical_pairs']})
    (BASE/'physical_pair_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
