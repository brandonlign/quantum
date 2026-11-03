"""Exact necessary screens for the exclusive rank-one weight-eight column."""
from collections import Counter
import json
from math import comb
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
from model import weight
from mixed_frames.stage3.common import sha,atomic_json


def kraw(n,w,j):
    return sum((-1)**a*3**(w-a)*comb(j,a)*comb(n-j,w-a)
               for a in range(max(0,w-n+j),min(w,j)+1))


def main():
    # Rebuild the A6+A7<=2 identity with integer arithmetic, independently
    # of LP or the saved screen booleans. Equality rows are unrestricted;
    # the negative shadow transform rows are <=0 and have positive multipliers.
    residual=[]
    for j in range(18):
        v=-1835008*int(j==0)+116
        v+=sum(m*(kraw(17,w,j)-16384*int(j==w)) for w,m in enumerate((60,25,11,3,1),1))
        v+=sum(m*int(j==w) for w,m in enumerate((-851968,-278528,-507904,-147456,-180224),1))
        v-=sum(m*((-1)**j)*kraw(17,w,j) for w,m in enumerate((88,32,18,4,2)))
        residual.append(v-32768*int(j in (6,7)))
    assert residual==[0]*18 and -1835008+116*16384==65536
    plan=json.loads((HERE/'plane_plan.json').read_text());records=[]
    for row in plan['records']:
        items=[]
        for i,orbit in enumerate(row['orbits']):
            reasons=[]
            for p in orbit['members']:
                assert len(set(p))==3 and p[0]^p[1]^p[2]==0
                ws=[weight(v,9) for v in p];supports=[(v&511)|(v>>9) for v in p]
                assert min(ws)>=6
                if sum(v<8 for v in ws)>2:reason='three_light'
                elif any(ws[a]==8 and any(b!=a and supports[b]&~supports[a]==0 for b in range(3)) for a in range(3)):
                    reason='rank_two_weight8_support'
                else:reason=None
                reasons.append(reason)
            assert len(set(reasons))==1
            items.append(dict(orbit_index=i,reason=reasons[0],covered_planes=len(orbit['members'])))
        records.append(dict(cls=row['cls'],orbits=items,remaining=sum(i['reason'] is None for i in items)))
    counts=Counter(i['reason'] or 'remaining' for r in records for i in r['orbits'])
    report=dict(status='PASS',records=records,screen_counts=dict(counts),
                identity_residual=residual,identity_rhs=65536,
                fully_screened_classes=[r['cls'] for r in records if r['remaining']==0],
                scope='Necessary screens only for the exclusive r=1 column: every weight-eight stabilizer support has local rank one; conditional on the classified suffix cell',
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'plane_plan.json',ROOT/'model.py']})
    atomic_json(HERE/'column_screen.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='records'}))


if __name__=='__main__':main()
