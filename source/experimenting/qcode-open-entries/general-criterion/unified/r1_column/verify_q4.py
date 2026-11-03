"""Checks for Q4_KERNEL.md: (1) wt mod 2 has polar form = symplectic product (all site pairs);
(2) Gauss-sum odd-count bound for every quadratic form on F2^5 (exhaustive over forms given by
upper-triangular matrices + linear part); (3) the arithmetic O<=4 and the Griesmer q<=4 input."""
import itertools, json
P=[(0,0),(1,0),(0,1),(1,1)]  # I X Z Y as (x,z)
for a in P:
  for b in P:
    s=(a[0]^b[0],a[1]^b[1]); wt=lambda u:int(u!=(0,0))
    assert (wt(s)-wt(a)-wt(b)-(a[0]*b[1]+a[1]*b[0]))%2==0
# every quadratic form on F2^5: q(v)=sum_{i<=j} c_ij v_i v_j
idx=[(i,j) for i in range(5) for j in range(i,5)]
vecs=list(itertools.product([0,1],repeat=5)); minodd_nonzero=99
for mask in range(1<<len(idx)):
  odd=0
  for v in vecs:
    odd+=sum(v[i]&v[j] for t,(i,j) in enumerate(idx) if mask>>t&1)&1
  if mask: minodd_nonzero=min(minodd_nonzero,odd)
assert minodd_nonzero==8
# arithmetic: |H|=32, weight sum <= 8*32*3/4, nonzero words >=6, one of weight 8, odd words >=7
upper=8*32*3//4; lower=30*6+8
assert upper-lower==4
griesmer=lambda m,D:max(q for q in range(m+1) if sum(-(-D//2**i) for i in range(q))<=m)
assert griesmer(8,4)==4
json.dump(dict(status='PASS',min_odd_nonzero_form=8,O_max=upper-lower,griesmer_q_bound=4),
          open('verify_q4.json','w'),indent=1); print('PASS')
