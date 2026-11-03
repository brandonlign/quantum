"""Direct physical verification of every suffix table used by the C++ search."""
from pathlib import Path
import json,hashlib,time
import numpy as np
BASE=Path(__file__).resolve().parent

def indexed(gens):
    out=np.array([0],dtype=np.uint32)
    for g in gens:out=np.concatenate([out,out^np.uint32(g)])
    return out


def main():
    t=time.monotonic();items=json.loads((BASE.parent/'rank2_z_inputs.json').read_text())['classes'];reports=[]
    counts=np.array([i.bit_count() for i in range(512)],dtype=np.uint8)
    def wt(words):return counts[(words&511)|(words>>9)]
    for item in items:
        data=json.loads(Path(item['table']).read_text());P=data['P_generators'];C=data['C_generators'];pw=indexed(P);cw=indexed(C)
        assert len(set(map(int,pw)))==4096 and len(set(map(int,cw)))==64
        assert int(wt(cw[1:]).min())==6
        lam=wt(pw);delta=wt(pw[:,None]^cw[None,:]).min(axis=1)
        assert list(map(int,lam))==data['suffix_weights'];assert list(map(int,delta))==data['coset_distances']
        def sp(a,b):return (((a&511)&(b>>9)).bit_count()+((b&511)&(a>>9)).bit_count())%2
        G=[sum(sp(a,b)<<j for j,b in enumerate(P)) for a in P];form=indexed(G)
        assert list(map(int,form))==data['form_rows']
        assert all(sp(a,b)==0 for a in C for b in P)
        assert set(map(int,pw))&set(map(int,cw))=={0}
        flags=(delta>=4)&(lam>=6);middle=(delta>=2)&(lam>=4)
        expected=flags.astype(int)+2*middle.astype(int);expected[0]=0
        assert list(map(int,expected))==data['flags']
        reports.append(dict(class_index=item['class_index'],cosets=4096,physical_words_per_coset=64,
                            table_sha256=hashlib.sha256(Path(item['table']).read_bytes()).hexdigest()))
    result=dict(status='PASS',classes=len(reports),total_physical_coset_word_checks=15*4096*64,
        reports=reports,seconds=time.monotonic()-t,
        scope='Direct suffix-word enumeration independent of syndrome BFS; exact P/C and flags checks.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (BASE/'table_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='reports'}))
if __name__=='__main__':main()
