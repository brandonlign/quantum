from itertools import product
from pathlib import Path
import hashlib,json,time

def image(mask,columns):
    out=0
    for i,v in enumerate(columns):
        if mask>>i&1:out^=v
    return out

def symp(a,b,p):
    mask=(1<<p)-1
    return (((a&mask)&(b>>p)).bit_count()+((b&mask)&(a>>p)).bit_count())&1

def main():
    (Path(__file__).resolve().parents[1]/'audit').mkdir(exist_ok=True)
    assert __debug__, "Run without Python -O: assertions are required."
    start=time.monotonic();total=0;cases=[]
    for p,b in [(1,1),(1,2),(1,3),(1,4),(2,2),(2,3),(2,4)]:
        count=0;limit=1<<b
        xs=((x,) for x in range(1,limit)) if p==1 else ((x,y) for x in range(1,limit) for y in range(1,limit) if y!=x)
        for X in xs:
            for Y in product(range(limit),repeat=p):
                M=X+Y
                kernel={c for c in range(1<<(2*p)) if image(c,M)==0}
                n_actual=set();n_recovered=set();rad_images=set()
                for c in range(limit):
                    xx=sum(((x&c).bit_count()&1)<<j for j,x in enumerate(X))
                    yy=sum(((y&c).bit_count()&1)<<j for j,y in enumerate(Y))
                    gamma_c=image(yy,X)^image(xx,Y)
                    if xx==0 and yy==0:n_actual.add(c)
                    if gamma_c==0:
                        rad_images.add(yy|(xx<<p))
                        if xx==0:n_recovered.add(c)
                assert n_actual==n_recovered
                rad_kernel={k for k in kernel if all(symp(k,l,p)==0 for l in kernel)}
                assert rad_images==rad_kernel
                assert len(kernel)*(1<<b)==(1<<(2*p))*len(n_actual)
                count+=1
        total+=count;cases.append({'half_dimension':p,'syndrome_dimension':b,'matrices':count})
    # Local coset used in the worked example: ZIII and IZZZ.
    coset=[0x10,0x10^0xf0]
    weights=[((v&15)|(v>>4)).bit_count() for v in coset]
    assert weights==[1,3] and min(weights)+3==4
    report={'status':'PASS','scope':'Exhaustive finite sanity check of the expanded exact-kernel argument, not a proof of the large code exclusions.','matrices':total,'cases':cases,'identities':['ker(M^T)=ker(Gamma) intersect ker(X^T) for independent columns of X','tau(rad(P))=rad(ker M)','rank-nullity for both kernels'],'worked_local_coset_minimum_weight':1,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'seconds':time.monotonic()-start}
    (Path(__file__).resolve().parents[1]/'audit/human-kernel-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
