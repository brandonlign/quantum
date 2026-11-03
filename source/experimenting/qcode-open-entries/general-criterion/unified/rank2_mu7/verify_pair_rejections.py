"""Direct physical pair-sector oracle on the four prior mu7 unary leaves."""
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from model import lift_coordinate,span,weight

HERE=Path(__file__).resolve().parent
BASE=HERE.parent


def cost(x,z):
    return min(((x^rx)|(z^rz)).bit_count() for rx in [0,255] for rz in [0,255])


def main():
    source=BASE/'rank2_mu7_orbits/leaf_solver.jsonl'
    leaves=[json.loads(line) for line in source.read_text().splitlines()]
    assert len(leaves)==4
    records=[]
    for leaf in leaves:
        index=leaf['class_index'];Z=leaf['Z']
        data=json.loads((BASE/f'rank2_z/class_{index}.json').read_text())
        physical=[lift_coordinate(v,data['P_generators']) for v in range(4096)]
        code=span(data['C_generators'])
        delta=[min(weight(w^c,9) for c in code) for w in physical]
        lam=[weight(w,9) for w in physical]
        assert delta==data['coset_distances'] and lam==data['suffix_weights']
        Zspan=[lift_coordinate(v,Z) for v in range(64)]
        results=[]
        for sector in range(1,64):
            if sector.bit_count()!=2:continue
            count=0
            for beta in range(4096):
                if any((data['form_rows'][a]&beta).bit_count()%2 != ((sector>>i)&1)
                       for i,a in enumerate(Z)):continue
                for z,az in enumerate(Zspan):
                    c=cost(sector,z|((z.bit_count()%2)<<7))
                    if delta[beta^az]+c<6 or lam[beta^az]+c<7:break
                else:count+=1
            results.append(dict(sector=sector,domain_size=count))
        empty=[r['sector'] for r in results if r['domain_size']==0]
        assert empty
        records.append(dict(class_index=index,Z=Z,empty_pair_sectors=empty,domains=results))
    report=dict(status='PASS',rejected_prior_leaves=4,records=records,
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [source,Path(__file__),BASE/'rank2_z/class_41.json',BASE/'rank2_z/class_50.json']},
                scope='All possible pair images physically checked for four prior leaves; no complete target Z search')
    (HERE/'pair_rejection_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(dict(status='PASS',rejected_prior_leaves=4,empty_pair_sector_counts=[len(r['empty_pair_sectors']) for r in records]))


if __name__=='__main__':main()
