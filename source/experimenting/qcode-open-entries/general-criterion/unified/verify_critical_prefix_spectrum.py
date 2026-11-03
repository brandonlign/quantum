"""Physical checks of the critical prefix bound and complete one-step extensions."""
from collections import Counter
import hashlib
import json
from pathlib import Path
from model import span,weight,symp,quotient_basis,lift_coordinate,dual,reduced_basis
from gauge_prefix_candidates import cert

HERE=Path(__file__).resolve().parent


def check(G):
    words=span(G);M=len(words)
    assert M>=16 and len(reduced_basis(G))==len(G)
    assert all(symp(a,b,8)==0 for a in G for b in G)
    assert min(weight(w,8) for w in words if w)>=6
    assert any(weight(w,8)==8 for w in words)
    assert all(weight(w,8)%2==0 for w in words)
    a=sum(len({((w>>j)&1)|(((w>>(8+j))&1)<<1) for w in words})==2 for j in range(8))
    spectrum=Counter(weight(w,8) for w in words)
    assert a<=16//M and spectrum[8]==3-M*a//8
    numerator=sum(count*{0:252,6:-12,8:28}[w] for w,count in spectrum.items())
    assert numerator%M==0
    coefficient=numerator//M
    assert coefficient>=0
    if M>=32:assert M==32 and a==0 and spectrum=={0:1,6:28,8:3}
    return dict(size=M,spectrum=dict(spectrum),rank_one_sites=a,dual_weight_two_coefficient=coefficient)


def main():
    source=HERE/'gauge_prefix_candidates.json'
    codes=json.loads(source.read_text())['groups']['4']['classes']
    assert len(codes)==6
    records=[];certificates=set();cosets=0;controls=[]
    for index,r in enumerate(codes):
        G=r['generators'];result=check(G);E=quotient_basis(G,8);words=span(G)
        assert len(E)==8
        extensions=[]
        for label in range(256):
            v=lift_coordinate(label,E);cosets+=1
            if min(weight(v^w,8) for w in words)<6:continue
            H=G+[v];physical=check(H);certificates.add(cert(H,8))
            normalizer=dual(H,8);S=set(span(H))
            distance=min(weight(w,8) for w in span(normalizer) if w not in S)
            assert distance==3
            extensions.append(dict(added_generator=v,physical=physical,quantum_distance=distance))
            controls.append(H)
        records.append(dict(prefix_class=index,physical=result,extensions=extensions))
    assert [len(r['extensions']) for r in records]==[2,0,0,0,0,2]
    result=dict(status='PASS',prefix_classes=6,normalizer_cosets=cosets,
                dimension_five_extension_classes=len(certificates),records=records,
                positive_control_generators=controls[0],
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [source,Path(__file__),HERE/'CRITICAL_PREFIX_SPECTRUM.md']},
                scope='All six existing prefix representatives and exhaustive one-step extensions; not an independent prefix census or target17 exclusion')
    (HERE/'critical_prefix_spectrum_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['records','sha256','positive_control_generators']})


if __name__=='__main__':main()
