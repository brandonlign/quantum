"""Eliminate a prefix-only logical kernel by fixing a maximal commuting part.

Exact binary algebra; preserves the shortened centralizer C. This does not
make the inverse distance test automatic or establish novelty.
"""
from model import reduced_basis,nullspace,dual,gram,symp,lift_coordinate


def supported_subspace(gens,n,sites):
    sites=set(sites);gens=reduced_basis(gens)
    equations=[sum(((g>>bit)&1)<<j for j,g in enumerate(gens))
               for i in range(n) if i not in sites for bit in [i,n+i]]
    return reduced_basis([lift_coordinate(c,gens) for c in nullspace(equations,len(gens))])


def maximal_isotropic(words,n):
    available=reduced_basis(words);out=[]
    while available:
        pair=next(((i,j) for i in range(len(available)) for j in range(i+1,len(available))
                   if symp(available[i],available[j],n)),None)
        if pair is None:out.extend(available);break
        i,j=pair;u,v=available[i],available[j]
        remaining=[w for k,w in enumerate(available) if k not in [i,j]]
        out.append(u)
        available=[w^(u if symp(w,v,n) else 0)^(v if symp(w,u,n) else 0) for w in remaining]
    assert all(symp(a,b,n)==0 for a in out for b in out)
    return reduced_basis(out)


def fix_prefix_kernel(S,n,T):
    if not isinstance(n,int) or n<1:raise ValueError('positive integer length required')
    S=list(S);T=list(T)
    if any(not isinstance(v,int) or not 0<=v<1<<(2*n) for v in S):raise ValueError('invalid Pauli word')
    if len(set(T))!=len(T) or any(not isinstance(i,int) or not 0<=i<n for i in T):raise ValueError('invalid support')
    S=reduced_basis(S)
    if any(symp(a,b,n) for a in S for b in S):raise ValueError('S is not isotropic')
    T=sorted(T);W=[i for i in range(n) if i not in T];D=dual(S,n)
    R=supported_subspace(S,n,T);H=supported_subspace(D,n,T)
    C=supported_subspace(D,n,W)
    seen=list(R);K=[]
    for v in H:
        if len(reduced_basis(seen+[v]))>len(seen):K.append(v);seen.append(v)
    q=len(K);t=len(reduced_basis(gram(K,n)));ell=q-t//2
    L=maximal_isotropic(K,n)
    assert len(L)==ell and ell<=n-len(S)
    Splus=reduced_basis(S+L);Dplus=dual(Splus,n)
    Rplus=supported_subspace(Splus,n,T);Hplus=supported_subspace(Dplus,n,T)
    Cplus=supported_subspace(Dplus,n,W)
    assert len(Splus)==len(S)+ell and len(Rplus)==len(R)+ell
    assert Rplus==reduced_basis(R+L) and Cplus==C
    assert Hplus==Rplus
    return dict(n=n,k=n-len(S),k_fixed=n-len(Splus),T=T,
                local_rank=len(R),kernel_dimension=q,kernel_form_rank=t,
                fixed_logicals=ell,local_rank_fixed=len(Rplus),
                S=S,S_fixed=Splus,R=R,R_fixed=Rplus,K_complement=K,L=L,
                shortened_centralizer=C,shortened_centralizer_fixed=Cplus,
                fixed_kernel_dimension=0)
