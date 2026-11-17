#!/usr/bin/env python3
"""Exact portable certificate for all claims in the accompanying manuscript.
Python 3.10+, standard library only; run without -O. Does not edit inputs.
"""
from pathlib import Path
from collections import Counter
from math import comb
import json,runpy
O=Path(__file__).resolve().parent
if not __debug__: raise RuntimeError('Run without -O; assertions are proof checks.')
def span(rows):
    out=[0]
    for r in rows: out += [v^r for v in out]
    return out
def wt(v,n): return ((v&((1<<n)-1))|(v>>n)).bit_count()
def sp(a,b,n):
    m=(1<<n)-1
    return (((a&m)&(b>>n))^((b&m)&(a>>n))).bit_count()%2
def rank(rows):
    piv={}
    for v in rows:
        while v:
            i=v.bit_length()-1
            if i in piv: v^=piv[i]
            else: piv[i]=v;break
    return len(piv)
def solve(eq,vals,w):
    rows=[a|(b<<w) for a,b in zip(eq,vals)];piv=[];i=0
    for j in range(w):
        k=next((k for k in range(i,len(rows)) if rows[k]>>j&1),None)
        if k is None:continue
        rows[i],rows[k]=rows[k],rows[i]
        for k in range(len(rows)):
            if k!=i and rows[k]>>j&1:rows[k]^=rows[i]
        piv.append(j);i+=1
    assert all(r&((1<<w)-1) or not(r>>w) for r in rows)
    return sum((rows[k]>>w)<<j for k,j in enumerate(piv))
def gram(rows,n):return [sum(sp(a,b,n)<<j for j,b in enumerate(rows)) for a in rows]
def label_string(v,n=18):return ''.join('IXZY'[(v>>i&1)+2*(v>>(n+i)&1)] for i in range(n))
print('=== Independent full normalizer/stabilizer enumerations ===')
for c in 'ABC':
    result=runpy.run_path(str(O/('frozen_'+c)/'verify_standalone.py'),run_name='__main__')
    expected=json.loads((O/('frozen_'+c)/'STANDALONE_RESULT.json').read_text())
    for key,var in [('stabilizer','A'),('normalizer','B'),('logical','L')]:
        assert {str(w):count for w,count in result[var].items()}==expected[key+'_enumerator']
print('=== Local gluing, quotient dual forms and exact reconstruction ===')
runpy.run_path(str(O/'construction/reconstruct_from_gluing.py'),run_name='__main__')
ctx=json.loads((O/'construction/PHYSICAL_HIT_CONTEXT.json').read_text());F=json.loads((O/'construction/HIT.json').read_text())['map']
J=[2,1,8,4,32,16,128,64,0,0]
for name,dim,s,g in [('left',4,2,0),('right',7,3,2)]:
    z=ctx[name];n=z['n'];C=z['C'];Q=z['Q'];K=C[:s];m=(1<<n)-1
    assert rank(C)==dim and rank(K)==s and rank(gram(C,n))==g
    assert all(sp(a,b,n)==0 for a in K for b in C+Q)
    assert rank(C+Q)==2*n-s and len(Q)==10
    eq=[(r>>n)|((r&m)<<n) for r in C+Q]
    dual=[solve(eq,[0]*len(C)+[int(i==j) for i in range(10)],2*n) for j in range(10)]
    assert gram(dual,n)==J and rank(J)==8
    print(name,{'local_rank':dim,'marked_rank':s,'local_Gram_rank':g,'dual_Gram_rows':J,'dual_lifts':dual,'local_enumerator':dict(sorted(Counter(wt(v,n) for v in span(C)).items()))})
cols=[F[1<<i] for i in range(10)];assert cols==[240,31,248,100,144,202,242,97,389,678]
M=[[cols[j]>>i&1 for j in range(10)] for i in range(10)]
assert all(sum(M[i][a]*(J[a]>>b&1)*M[j][b] for a in range(10) for b in range(10))%2==(J[i]>>j&1) for i in range(10) for j in range(10))
counts=Counter(ctx['left']['cost'][t]+ctx['right']['cost'][F[t]] for t in range(1,1024))
assert counts=={6:261,7:333,8:337,9:73,10:19};print('nonzero_coset_minimum_sums',dict(sorted(counts.items())))
print('=== Exact binary transvections and Hermitian phases ===')
A=json.loads((O/'frozen_A/rows.json').read_text())
def phase(a,b,n):
    m=(1<<n)-1;c=a^b
    return (((a&m)&(a>>n)).bit_count()+((b&m)&(b>>n)).bit_count()-((c&m)&(c>>n)).bit_count()+2*((a>>n)&(b&m)).bit_count())%4
P=[((1,0),(0,1)),((0,1),(1,0)),((1,0),(0,-1)),((0,-1j),(1j,0))]
for a in range(4):
    for b in range(4):
        product=tuple(tuple(sum(P[a][i][k]*P[b][k][j] for k in range(2)) for j in range(2)) for i in range(2))
        assert product==tuple(tuple((1j**phase(a,b,1))*x for x in row) for row in P[a^b])
for c,v,w,neg in [('B',4328808577,17435,[9,12]),('C',8607273026,280636,[3,4,5])]:
    rows=json.loads((O/('frozen_'+c)/'rows.json').read_text());T=lambda r:r^(v if sp(r,v,18) else 0)
    assert rows==[T(r) for r in A] and wt(v,18)==4
    assert all(T(T(1<<i))==1<<i for i in range(36))
    assert all(sp(T(1<<i),T(1<<j),18)==sp(1<<i,1<<j,18) for i in range(36) for j in range(36))
    negative=[]
    for j,(a,b) in enumerate(zip(A,rows),1):
        e=(1+phase(v,a,18))%4 if sp(v,a,18) else 0
        assert e in [0,2]
        if e:negative.append(j)
        assert e//2==sp(w,b,18)
    assert negative==neg
    assert (O/('frozen_'+c)/'generators.txt').read_text().splitlines()==[label_string(r) for r in rows]
    print(c,{'v':v,'v_string':label_string(v),'correction':w,'correction_string':label_string(w),'bare_negative_rows_one_based':negative})
print('PASS: all three codes, gluing, complete enumerations, MacWilliams transforms, witnesses, transformations and phases.')
