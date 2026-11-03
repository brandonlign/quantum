#!/usr/bin/env python3
"""Check the matching constructions without a coding package.

Standard library only. Binary word layout: x | (z << n).
Enumerates every centralizer word with a Gray-code traversal.
"""
from pathlib import Path
import json,time
HERE=Path(__file__).resolve().parent

def nullspace(rows,width):
    rows=list(rows);piv=[];r=0
    for c in range(width):
        j=next((j for j in range(r,len(rows)) if (rows[j]>>c)&1),None)
        if j is None:continue
        rows[r],rows[j]=rows[j],rows[r]
        for j in range(len(rows)):
            if j!=r and (rows[j]>>c)&1:rows[j]^=rows[r]
        piv.append(c);r+=1
    basis=[]
    for c in range(width):
        if c in piv:continue
        v=1<<c
        for i,p in enumerate(piv):v|=((rows[i]>>c)&1)<<p
        basis.append(v)
    return r,basis

def main():
    (Path(__file__).resolve().parents[1]/'audit').mkdir(exist_ok=True)
    assert __debug__
    started=time.monotonic();records=[]
    for n,k in [(16,5),(17,3)]:
        data=json.loads((HERE/f'data/construction-{n}-{k}.json').read_text())
        g=data['generators'];mask=(1<<n)-1
        equations=[(v>>n)|((v&mask)<<n) for v in g]
        assert len(g)==nullspace(g,2*n)[0]==n-k
        assert all((a&b).bit_count()%2==0 for a in equations for b in g)
        S={0}
        for v in g:S|={w^v for w in S}
        rank,D=nullspace(equations,2*n);assert rank==n-k and len(D)==n+k
        minimum=n+1;witness=None;word=0
        for i in range(1,1<<len(D)):
            word^=D[(i&-i).bit_length()-1]
            if word in S:continue
            w=((word&mask)|(word>>n)).bit_count()
            if w<minimum:minimum=w;witness=word
        assert minimum==data['claimed_distance'],(n,k,minimum)
        records.append({'n':n,'k':k,'distance':minimum,'stabilizer_rank':n-k,'centralizer_rank':n+k,
            'centralizer_words_checked':1<<(n+k),'logical_word_witness':witness,'source':data['source']})
    out={'status':'PASS','scope':'Complete physical centralizer enumeration of both matching constructions.',
        'records':records,'seconds':time.monotonic()-started}
    (HERE.parent/'audit/review-constructions.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
