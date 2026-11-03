"""Fixed radical/relative distance profiles for kernel-zero stabilizer extensions.

Works for impure codes as well as pure codes. Reference computation is bounded
explicitly. The unknown is a symplectic isometry, rather than an M-dependent
logical-distance profile.
"""
from model import reduced_basis,dual,gram,nullspace,lift_coordinate,coordinates,symp,weight,embed,split,span,quotient_basis
from gauge_fix import supported_subspace


def ordered_span(gens):
    out=[0]
    for g in gens:out += [w^g for w in out]
    return out


def quotient_in_space(P,N,f):
    seen=list(N);available=[]
    for w in P:
        if len(reduced_basis(seen+[w]))>len(seen):available.append(w);seen.append(w)
    zs=[];xs=[]
    while available:
        u=available.pop(0);j=next(i for i,v in enumerate(available) if symp(u,v,f))
        v=available.pop(j);zs.append(u);xs.append(v)
        available=[w^(u if symp(w,v,f) else 0)^(v if symp(w,u,f) else 0) for w in available]
    return zs+xs


def optional_min(values):return min((v for v in values if v is not None),default=None)


class ZeroKernelProfiles:
    def __init__(self,n,k,m,R,C,max_quotient_bits=14,max_code_bits=14):
        self.n,self.k,self.m=n,k,m;self.f=f=n-m
        if not all(isinstance(v,int) for v in [n,k,m]) or not 0<=k<=n or not 0<=m<=n:
            raise ValueError('invalid dimensions')
        if m>10 or f>10:raise ValueError('reference supports prefix and suffix lengths<=10')
        R=reduced_basis(R);C=reduced_basis(C)
        if any(not 0<=w<1<<(2*m) for w in R) or any(not 0<=w<1<<(2*f) for w in C):raise ValueError('invalid word')
        if any(symp(a,b,m) for a in R for b in R):raise ValueError('R must be isotropic')
        self.R,self.C=R,C;self.r=len(R)
        if len(C)!=n+k-2*m+self.r:raise ValueError('shortened dimension mismatch')
        self.P=dual(C,f);self.b=len(self.P)
        self.N=reduced_basis([lift_coordinate(v,C) for v in nullspace(gram(C,f),len(C))])
        self.E=quotient_basis(R,m);self.a=len(self.E);self.h=len(self.N)
        if self.b-self.h!=self.a:raise ValueError('requires kernel-zero Gram rank equality')
        if self.a>max_quotient_bits or len(C)>max_code_bits:raise ValueError('reference profile size limit exceeded')
        self.Q=quotient_in_space(self.P,self.N,f)
        assert len(self.Q)==self.a
        self.prefix_words=ordered_span(self.E);self.suffix_words=ordered_span(self.Q)
        Rs=ordered_span(R);Cs=ordered_span(C);Ns=ordered_span(self.N);Nset=set(Ns)
        logical=[c for c in Cs if c not in Nset]
        self.cost=[min(weight(w^r,m) for r in Rs) for w in self.prefix_words]
        self.stabilizer_profile=[min(weight(w^c,f) for c in Ns) for w in self.suffix_words]
        self.logical_profile=[min((weight(w^c,f) for c in logical),default=None) for w in self.suffix_words]
        self.absolute_profile=[min(weight(w^c,f) for c in Cs) for w in self.suffix_words]
        self.min_R=optional_min(weight(w,m) for w in Rs if w)
        self.min_N=optional_min(weight(w,f) for w in Ns if w)
        self.min_C=optional_min(weight(w,f) for w in Cs if w)

    def check_isometry(self,columns):
        if len(columns)!=self.a or len(reduced_basis(columns))!=self.a or any(not 0<=w<1<<self.a for w in columns):raise ValueError('isometry must be a bijection')
        images=[lift_coordinate(v,self.Q) for v in columns]
        if any(symp(self.E[i],self.E[j],self.m)!=symp(images[i],images[j],self.f) for i in range(self.a) for j in range(self.a)):
            raise ValueError('map does not preserve the symplectic form')
        return images

    def reconstruct(self,columns):
        images=self.check_isometry(columns)
        return reduced_basis([embed(r,0,self.m,self.f) for r in self.R]
             +[embed(0,w,self.m,self.f) for w in self.N]
             +[embed(e,w,self.m,self.f) for e,w in zip(self.E,images)])

    def metrics(self,columns):
        self.check_isometry(columns);values=ordered_span(columns)
        logical=optional_min(self.cost[e]+self.logical_profile[v] if self.logical_profile[v] is not None else None for e,v in enumerate(values))
        stab=optional_min([self.min_R,self.min_N]+[self.cost[e]+self.stabilizer_profile[v] for e,v in enumerate(values) if e])
        pure=optional_min([self.min_R,self.min_C]+[self.cost[e]+self.absolute_profile[v] for e,v in enumerate(values) if e])
        return dict(logical_distance=logical,minimum_stabilizer_weight=stab,pure_distance=pure)


def recover_profiles(S,n,T):
    T=sorted(T);W=[i for i in range(n) if i not in T];order=T+W;m=len(T);f=n-m
    normalized=[]
    for w in S:
        normalized.append(sum(((w>>i)&1)<<j for j,i in enumerate(order)) | (sum(((w>>(n+i))&1)<<j for j,i in enumerate(order))<<n))
    Rfull=supported_subspace(normalized,n,range(m));Cfull=supported_subspace(dual(normalized,n),n,range(m,n))
    R=[split(w,m,f)[0] for w in Rfull];C=[split(w,m,f)[1] for w in Cfull]
    model=ZeroKernelProfiles(n,n-len(reduced_basis(S)),m,R,C)
    qbasis=[];ebasis=[]
    for w in normalized:
        v,suffix=split(w,m,f);ecoef=coordinates(v,model.R+model.E)>>model.r
        qcoef=coordinates(suffix,model.N+model.Q)>>model.h
        if len(reduced_basis(qbasis+[qcoef]))>len(qbasis):qbasis.append(qcoef);ebasis.append(ecoef)
    assert len(qbasis)==model.a
    tau=[lift_coordinate(coordinates(1<<j,qbasis),ebasis) for j in range(model.a)]
    columns=[coordinates(1<<i,tau) for i in range(model.a)]
    assert model.reconstruct(columns)==reduced_basis(normalized)
    return model,columns
