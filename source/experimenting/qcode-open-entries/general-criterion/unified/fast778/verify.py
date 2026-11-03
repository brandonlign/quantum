"""Physical symmetry and positive control for the 778 S6 lane."""
import itertools
import json
import subprocess
from model import ExtensionModel,coordinates,lift_coordinate,span,weight
from mixed_frames.frames import mixed_frame,xor_table
from mixed_frames.fixed_z import validate_solution
from mixed_frames.stage3.common import graph_counts,input_text
from speedup_review.verify_688_symmetry import permute,normal
from .common import HERE,ROOT,BINARY,atomic_json,hashes


def main():
    f=mixed_frame('778')
    assert f.E[:6]==[0x4000,0x2100,0x2200,0x2400,0x2800,0x3000]
    checks=0
    for j in range(5):
        sites=list(range(8));sites[j],sites[j+1]=sites[j+1],sites[j]
        assert {permute(r,sites) for r in span(f.R)}==set(span(f.R))
        A=[coordinates(permute(w,sites),f.R+f.E)>>2 for w in f.E]
        assert A[0]==1 and all(v<64 for v in A[:6])
        for e,w in enumerate(xor_table(f.E)):
            ee=lift_coordinate(e,A)
            assert sorted(weight(w^r,8) for r in span(f.R))==sorted(weight(lift_coordinate(ee,f.E)^r,8) for r in span(f.R))
            checks+=1
    m=ExtensionModel(14,0,2,8,f.R,[])
    columns=[m.syndrome(1<<(6+i)) for i in range(6)]+[m.syndrome(1<<i) for i in range(6)]
    reverse=[5,4,3,2,1,0,6,7]
    A=[coordinates(permute(w,reverse),f.R+f.E)>>2 for w in f.E]
    original=[lift_coordinate(a,columns) for a in A]
    validate_solution(m,f,original)
    wanted=normal(original[:6]);assert wanted!=original[:6]
    for order in itertools.permutations(range(6)):
        A=[coordinates(permute(w,[*order,6,7]),f.R+f.E)>>2 for w in f.E]
        changed=[lift_coordinate(a,original) for a in A]
        if changed[:6]==wanted:break
    else:raise AssertionError('No transport to normal form')
    validate_solution(m,f,changed)
    measured=graph_counts(m,f,original,2,3)
    assert measured==graph_counts(m,f,changed,2,3)
    path=HERE/'positive_input.txt'
    path.write_text(input_text(m,f,2,counts=(measured['n6'],measured['n7']),w6=2,w7=3))
    proc=subprocess.run([str(BINARY),str(path),'--prefix',','.join(map(str,wanted)),'--canonical','1'],capture_output=True,text=True,check=True)
    rows=[json.loads(l) for l in proc.stdout.splitlines()]
    assert rows[-1]['status']=='COMPLETE' and [r['X'] for r in rows[:-1]]==[wanted]
    report=dict(status='PASS',hashes=hashes(),quotient_checks=checks,
                positive=dict(parameters=[14,0,2],original=original[:6],representative=wanted,exact_profile=True))
    atomic_json(HERE/'verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='hashes'}))


if __name__=='__main__':main()
