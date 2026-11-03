"""Exact necessary-enumerator witnesses; no code construction or exclusion claim."""
from math import comb
from pathlib import Path
import hashlib,json

def kraw(n,j,w):
    return sum((-1)**i * 3**(j-i)*comb(w,i)*comb(n-w,j-i)
               for i in range(max(0,j-(n-w)),min(j,w)+1))

def main():
    n=9; dimension=6; size=1<<dimension
    A=[1,0,0,0,0,0,35,15,0,13]
    nums=[sum(count*kraw(n,j,w) for w,count in enumerate(A)) for j in range(n+1)]
    assert sum(A)==size and all(v>=0 and v%size==0 for v in nums)
    B=[v//size for v in nums]
    assert B==[1,0,3,14,231,420,917,1326,864,320]
    assert sum(B)==1<<(2*n-dimension)
    inverse=[sum(count*kraw(n,j,w) for w,count in enumerate(B)) for j in range(n+1)]
    assert [v//sum(B) for v in inverse]==A and all(v%sum(B)==0 for v in inverse)
    odd=sum(A[1::2]); hull=0
    assert size-2*odd==1<<((dimension+hull)//2)
    report=dict(status='PASS',candidate='[[17,3,6]]',local_branch='m=8,r=2,hull=0',
                suffix_length=n,suffix_dimension=dimension,weight_enumerator=A,
                dual_weight_enumerator=B,odd_weight_words=odd,
                nondegenerate_quadratic_gauss_sum=size-2*odd,
                scope='An exact witness satisfying distance, integral nonnegative MacWilliams '
                      'coefficients and the hull-zero parity Gauss sum. This is not an actual '
                      'additive code, stabilizer code, extension, or proof of feasibility.',
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_name('target17_enumerator_witness.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
if __name__=='__main__':main()
