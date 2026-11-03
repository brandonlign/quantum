"""Short depth-four frontier comparisons and legacy full-Z corpus regressions."""
import json
from pathlib import Path
from model import coordinates,lift_coordinate
from .common import HERE,ROOT,BINARY,atomic_json,budget,profile_mu,sha,target
from .checks import cpp
from mixed_frames.stage2.common import BINARY as OLD_BINARY


def legacy_frontiers():
    rows=[]
    for index in [25,33,34,36,41,44,50]:
        m,f=target('888',index)
        zs=[((1<<i)|(1<<7))<<8 for i in range(6)]
        change=[coordinates(z,f.R+zs)>>2 for z in f.E[:6]]
        path=ROOT/f'rank2_mu6/class_{index}.jsonl'
        for event in [json.loads(l) for l in path.read_text().splitlines()]:
            if event['status']!='LEAF':continue
            old=[lift_coordinate(a,m.Gamma) for a in event['Z']]
            xs=[lift_coordinate(z,old) for z in change]
            rows.append(dict(class_index=index,X=xs,legacy_Z=event['Z'],source_sha256=sha(path)))
    return rows


def main():
    records=[]
    for name,index,xs in [('688',7,[1457,45,85,332]),('888',25,[3784,409,2105,3598])]:
        m,f=target(name,index)
        profiles=[(1,0),(1,1),(2,0)] if name=='688' else [(0,1),(0,2),(1,0),(1,1),(2,0)]
        for profile in profiles:
            for lane in ('stage2','stage3-profile'):
                mu=6 if lane=='stage2' else profile_mu(f,profile)
                counts=None if lane=='stage2' else budget(f,profile)
                r=cpp(m,f,mu=mu,counts=counts,binary=OLD_BINARY if lane=='stage2' else BINARY,
                      prefix=xs,nodes=50000,seconds=3,leaves=5)
                row=dict(frame=name,class_index=index,profile=list(profile),lane=lane,prefix=xs,**r)
                records.append(row);print(json.dumps(row),flush=True)
    atomic_json(HERE/'frontier_pilots.json',dict(records=records,scope='Specified depth-four frontier only; no class exclusion. Any LIMIT remains unresolved.'))

if __name__=='__main__':main()
