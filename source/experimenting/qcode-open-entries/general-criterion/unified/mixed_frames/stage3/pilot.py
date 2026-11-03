"""Bounded serial necessary-Z pilots; no completion or exhaustive batch."""
import argparse
import json
import time
from .common import HERE,ROOT,BINARY,atomic_json,budget,profile_mu,sha,target
from .checks import cpp
from mixed_frames.stage2.common import BINARY as OLD_BINARY

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--seconds',type=float,default=2)
    ap.add_argument('--nodes',type=int,default=20000);a=ap.parse_args()
    assert 0<a.seconds<=10 and 0<a.nodes<=100000
    report=[]
    for name,index,profile in [('688',7,(1,0)),('688',7,(1,1)),('688',7,(2,0)),
                               ('888',7,(0,1)),('888',7,(0,2)),('888',21,(1,1)),('888',0,(2,0)),('888',12,(2,0))]:
        m,f=target(name,index)
        for lane in ('stage2','stage3-infinity','stage3-profile'):
            mu=6 # Match the stage-2 688 baseline and supply the same 888 extension.
            counts=None;binary=OLD_BINARY if lane=='stage2' else BINARY
            if lane=='stage3-profile':mu=profile_mu(f,profile);counts=budget(f,profile)
            r=cpp(m,f,mu=mu,counts=counts,binary=binary,nodes=a.nodes,seconds=a.seconds,leaves=5)
            # Preserve logs as separate pilots, never as exhaustive coverage.
            file=HERE/f'pilot_{name}_{index}_{profile[0]}{profile[1]}_{lane}.json'
            atomic_json(file,r)
            row=dict(frame=name,class_index=index,hull=m.hull,profile=list(profile),lane=lane,
                     mu=mu,budget=list(counts) if counts is not None else None,
                     **{k:r[k] for k in ('status','reason','nodes','candidates','leaf_count','seconds','mu_prunes')},
                     nodes_per_second=r['nodes']/max(r['seconds'],1e-9),result_sha256=sha(file),
                     binary_sha256=sha(binary))
            for k in ('count_prunes','exact_prunes','proper_exact_prunes'):row[k]=r.get(k,0)
            report.append(row);print(json.dumps(row),flush=True)
    atomic_json(HERE/'pilots.json',dict(scope='Bounded necessary Z only; LIMIT never excludes',records=report))

if __name__=='__main__':main()
