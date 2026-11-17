"""Reconstruct S from one weight-eight word and two explicitly supplied maps."""
import json
from pathlib import Path
P=Path(__file__).parent
C=json.loads((P/'a5_construction.json').read_text())
N=34;I=tuple(1<<i for i in range(N));M=(1<<17)-1
def need(ok,message):
    if not ok:raise RuntimeError(message)
def apply(v,images):
    out=0
    while v:
        b=v&-v;v^=b;out^=images[b.bit_length()-1]
    return out
def compose(a,b):return tuple(apply(v,a) for v in b)
def convert(c):
    need(sorted(c['permutation'])==list(range(17)),'permutation')
    need(all(m[0]==0 and sorted(m)==list(range(4)) for m in c['local_label_maps']),'local bijections')
    out=[]
    for k in range(34):
        i=k%17;a=1 if k<17 else 2;j=c['permutation'][i];b=c['local_label_maps'][i][a]
        out.append((b&1)<<j|(b>>1)<<(17+j))
    return tuple(out)
a=convert(C['generator_a']);b=convert(C['generator_b'])
need(compose(a,a)==I and a!=I,'order two')
need(compose(b,compose(b,b))==I and b!=I,'order three')
ab=compose(a,b);v=I
for k in range(5):v=compose(ab,v)
need(v==I and ab!=I,'product order five')
group={I};queue=[I]
for h in queue:
    for generator in [a,b]:
        gh=compose(generator,h)
        if gh not in group:
            group.add(gh);queue.append(gh)
    need(len(group)<=60,'group exceeds 60')
need(len(group)==60,'group size')
five=C['five_site_orbit'];restrictions=set()
for h in group:
    perm=[((h[i]&M)|(h[i]>>17)).bit_length()-1 for i in range(17)]
    r=tuple(five.index(perm[i]) for i in five)
    need(sum(r[i]>r[j] for i in range(5) for j in range(i+1,5))%2==0,'restriction parity')
    restrictions.add(r)
need(len(restrictions)==60,'faithful natural A5 action')
orbit={apply(C['seed'],h) for h in group}
need(len(orbit)==30,'orbit size')
need(all(((v&M)|(v>>17)).bit_count()==8 for v in orbit),'orbit weights')
sp=lambda v,w:(((v&M)&(w>>17)).bit_count()+((v>>17)&(w&M)).bit_count())%2
need(all(sp(v,w)==0 for v in orbit for w in orbit),'orbit isotropy')
span={0}
for v in sorted(orbit):span|={w^v for w in span}
need(len(span)==4096,'orbit span rank 12')
rows=json.loads((P/'generators.json').read_text())['rows'];original={0}
for v in rows:original|={w^v for w in original}
need(span==original,'matches independently verified code')
out={'status':'PASS','group_order':60,'identification':'faithful action as all even permutations on five sites, hence A5','orbit_size':30,'orbit_weight':8,'span_rank':12,'same_code':True,'scope':'An explicit A5 subgroup and orbit construction; not a full automorphism-group classification.'}
(P/'a5_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
