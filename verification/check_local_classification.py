#!/usr/bin/env python3
"""Reconstruct the full-plane local classifications from 43 seven-site parents.

Independent small enumerator: no production proof/search modules imported.
The complete parent classification is an explicit premise, proved by the
lengthening/incidence chain in the paper and the included release16 census.
Requires pynauty 2.8.8.1. Writes only the companion audit report.
"""
from collections import Counter
from pathlib import Path
import hashlib
import json
import time
import pynauty

HERE=Path(__file__).resolve().parent

def span(rows):
    words=[0]
    for row in rows:words += [w^row for w in words]
    assert len(set(words))==len(words)
    return sorted(words)

def weight(w,n):return ((w&((1<<n)-1))|(w>>n)).bit_count()

def symp(a,b,n):return (((a&((1<<n)-1))&(b>>n)).bit_count()+((b&((1<<n)-1))&(a>>n)).bit_count())%2

def canonical(rows,n):
    words=span(rows);W=len(words);letters=W+3*n
    adjacency={i:set() for i in range(letters+n)}
    def edge(a,b):adjacency[a].add(b);adjacency[b].add(a)
    for j,w in enumerate(words):
        for i in range(n):
            s=((w>>i)&1)|(((w>>(n+i))&1)<<1)
            if s:edge(j,W+3*i+s-1)
    for i in range(n):
        for letter in range(3):edge(letters+i,W+3*i+letter)
    graph=pynauty.Graph(letters+n,directed=False,adjacency_dict=adjacency,
        vertex_coloring=[set(range(W)),set(range(W,letters)),set(range(letters,letters+n))])
    return pynauty.certificate(graph)

def full_plane(rows):
    full=[w for w in span(rows) if weight(w,8)==8]
    return len(full)==3 and full[0]^full[1]==full[2]

def main():
    (Path(__file__).resolve().parents[1]/'audit').mkdir(exist_ok=True)
    assert __debug__
    start=time.monotonic();path=HERE/'data/parents-7-4-5.json'
    parents=json.loads(path.read_text());assert len(parents)==43
    assert len({canonical(g,7) for g in parents})==43
    found={};valid=0
    for rows in parents:
        words=span(rows);assert len(rows)==4 and min(weight(w,7) for w in words[1:])>=5
        for xform in range(16):
            for zform in range(16):
                lifted=[(v&127)|(((xform>>i)&1)<<7)|((v>>7)<<8)|(((zform>>i)&1)<<15) for i,v in enumerate(rows)]
                ws=span(lifted)
                if min(weight(w,8) for w in ws[1:])<6:continue
                valid+=1;found.setdefault(canonical(lifted,8),lifted)
    assert valid==861,(valid,len(found))
    critical=[g for g in found.values() if any(weight(w,8)==8 for w in span(g))]
    assert len(critical)==15,(len(found),len(critical))
    ranks=Counter()
    for rows in critical:
        # Gaussian elimination of the restricted Gram matrix.
        gram=[sum(symp(a,b,8)<<j for j,b in enumerate(rows)) for a in rows];piv={}
        for v in gram:
            while v:
                b=v.bit_length()-1
                if b in piv:v^=piv[b]
                else:piv[b]=v;break
        ranks[len(piv)]+=1
    assert ranks=={0:6,2:9},ranks
    bases=[g for g in found.values() if full_plane(g)];assert len(bases)==3
    assert all(g in critical for g in bases)
    h5={};accepted=0
    for rows in bases:
        words=span(rows);seen=set();cosets=0
        for w in range(1<<16):
            if w in seen or any(symp(w,h,8) for h in rows):continue
            coset=[w^h for h in words];seen.update(coset);cosets+=1
            if min(weight(v,8) for v in coset)<6:continue
            child=rows+[w];assert full_plane(child)
            accepted+=1;h5.setdefault(canonical(child,8),child)
        assert cosets==256
    assert accepted==2 and len(h5)==1,(accepted,len(h5))
    report={'status':'PASS','scope':'Fresh exhaustive coordinate lifts and physical cosets; completeness conditional on the 43 stated seven-site parents.',
        'parent_classes':43,'ordered_coordinate_maps':43*256,'minimum6_lifts':valid,
        'dimension4_classes':len(found),'dimension4_with_full_support_word':len(critical),
        'critical_dimension4_gram_rank_counts':dict(ranks),
        'dimension4_full_plane_classes':len(bases),'normalizer_cosets_per_parent':256,
        'minimum6_dimension5_extensions':accepted,'dimension5_full_plane_classes':len(h5),
        'dimension4_representatives':list(found.values()),'full_plane_dimension4_representatives':bases,
        'full_plane_dimension5_representatives':list(h5.values()),
        'parent_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'seconds':time.monotonic()-start}
    out=HERE.parent/'audit/review-local-classification.json';out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if 'representatives' not in k},indent=2))

if __name__=='__main__':main()
