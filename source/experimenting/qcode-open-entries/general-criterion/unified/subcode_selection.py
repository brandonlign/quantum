"""Exact bounded inverse of stabilizer augmentation via binary syndrome labels.

The theorem is general. This reference implementation bounds syndrome and
normalizer enumeration explicitly; it is not a large-code classification engine.
"""
from collections import deque
from model import reduced_basis,nullspace,dual,coordinates,lift_coordinate,symp,weight,span
from gauge_fix import supported_subspace


def ordered_span(gens):
    out=[0]
    for g in gens:out += [w^g for w in out]
    return out


def subcode(gens,rows):
    if len(reduced_basis(gens))!=len(gens):raise ValueError('basis must be independent')
    width=len(gens)
    if any(not isinstance(v,int) or not 0<=v<1<<width for v in rows):raise ValueError('invalid label row')
    if len(reduced_basis(rows))!=len(rows):raise ValueError('labels must be surjective')
    return [lift_coordinate(c,gens) for c in nullspace(rows,width)]


class SelectionContext:
    def __init__(self,gens,n,max_syndrome_bits=18,max_normalizer_bits=20):
        gens=list(gens)
        if not isinstance(n,int) or n<1:raise ValueError('positive length required')
        if any(not isinstance(v,int) or not 0<=v<1<<(2*n) for v in gens):raise ValueError('invalid word')
        if len(reduced_basis(gens))!=len(gens):raise ValueError('basis must be independent')
        if any(symp(a,b,n) for a in gens for b in gens):raise ValueError('stabilizer must be isotropic')
        if len(gens)>max_syndrome_bits or 2*n-len(gens)>max_normalizer_bits:
            raise ValueError('reference enumeration limit exceeded')
        self.gens,self.n=gens,n
        self.words=ordered_span(gens);self.weights=[weight(w,n) for w in self.words]
        st=set(self.words);normalizer=span(dual(gens,n))
        self.base_distance=min((weight(w,n) for w in normalizer if w not in st),default=None)
        N=1<<len(gens);self.distances=[-1]*N;self.distances[0]=0
        self.steps=[self.syndrome(v) for i in range(n) for v in [1<<i,1<<(n+i),(1<<i)|(1<<(n+i))]]
        queue=deque([0])
        while queue:
            a=queue.popleft()
            for t in self.steps:
                b=a^t
                if self.distances[b]<0:self.distances[b]=self.distances[a]+1;queue.append(b)
        assert all(v>=0 for v in self.distances)

    def syndrome(self,w):return sum(symp(g,w,self.n)<<i for i,g in enumerate(self.gens))

    def evaluate(self,rows,d,minimum_stabilizer_weight=1,T=None,R=None):
        if not isinstance(d,int) or d<1 or not isinstance(minimum_stabilizer_weight,int) or minimum_stabilizer_weight<1:
            raise ValueError('positive integer weights required')
        S=subcode(self.gens,rows);rowcode=span(rows)
        leader=min((self.distances[a] for a in rowcode if a),default=None)
        kept_low=True;avoids_low=True
        for a in range(1,len(self.words)):
            label=sum(((row&a).bit_count()%2)<<i for i,row in enumerate(rows));wt=self.weights[a]
            if wt<d and label:kept_low=False
            if wt<minimum_stabilizer_weight and not label:avoids_low=False
        distance_ok=(self.base_distance is None or self.base_distance>=d) and kept_low and (leader is None or leader>=d)
        support_ok=None
        if T is not None:
            if R is None:raise ValueError('desired local R required with T')
            Rp=supported_subspace(self.gens,self.n,T);R=reduced_basis(R)
            if any(w not in set(span(Rp)) for w in R):raise ValueError('R must be a subspace of the augmented local stabilizer')
            images=[sum(((row&coordinates(w,self.gens)).bit_count()%2)<<i for i,row in enumerate(rows)) for w in Rp]
            support_ok=(len(Rp)-len(R)==len(rows) and len(reduced_basis(images))==len(rows)
                        and all(not any((row&coordinates(w,self.gens)).bit_count()%2 for row in rows) for w in R))
            if support_ok:
                assert supported_subspace(S,self.n,T)==R
                W=[i for i in range(self.n) if i not in T]
                assert supported_subspace(dual(S,self.n),self.n,W)==supported_subspace(dual(self.gens,self.n),self.n,W)
        return dict(distance_ok=distance_ok,minimum_stabilizer_ok=avoids_low,
                    local_and_shortened_code_ok=support_ok,base_logical_distance=self.base_distance,
                    nonzero_syndrome_rowcode_distance=leader,low_base_stabilizers_retained=kept_low,
                    n=self.n,k=self.n-len(S),generators=S)
