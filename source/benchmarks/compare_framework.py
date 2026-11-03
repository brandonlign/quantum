"""Complete matched extension problems: full-map baseline versus exact pruning.

This measures small, fully exhausted instances of the paper's general model.
It does not extrapolate their speedups to the 16/17-qubit proof searches.
"""
from pathlib import Path
import argparse,collections,hashlib,itertools,json,platform,statistics,sys,time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'experimenting/qcode-open-entries/general-criterion/unified'
sys.path.insert(0,str(BASE))
from model import ExtensionModel,pullback_gram,reduced_basis,nullspace,span

def image(mask,cols):
    ans=0
    while mask:
        bit=mask&-mask;ans^=cols[bit.bit_length()-1];mask^=bit
    return ans

def digest(solutions):
    return hashlib.sha256(json.dumps(sorted(solutions),separators=(',',':')).encode()).hexdigest()

def run(model,variant):
    start=time.perf_counter();N=1<<model.b;a=model.a;p=a//2
    stats=collections.Counter(nodes=0,complete_maps_tested=0)
    solutions=[]
    if variant=='full_map_baseline':
        for columns in itertools.product(range(N),repeat=a):
            stats['nodes']+=1;stats['complete_maps_tested']+=1
            if model.criterion(columns):solutions.append(columns)
    else:
        if variant=='all_structural_tests' and max(model.cost)<model.d:
            obstruction=model.small_prefix_obstruction()
            if obstruction['excluded']:
                stats['root_'+obstruction['reason']+'_rejections']+=1
                return dict(variant=variant,seconds=time.perf_counter()-start,**stats,solutions=0,solution_sha256=digest([])),[]
        def dfs(columns,values):
            stats['nodes']+=1;depth=len(columns)
            if depth==a:
                stats['complete_maps_tested']+=1
                if model.criterion(columns):solutions.append(tuple(columns))
                return
            for v in range(N):
                newvalues=[x^v for x in values]
                offset=1<<depth
                if any(x==0 and model.cost[offset+i]<model.d for i,x in enumerate(newvalues)):
                    stats['kernel_rejections']+=1;continue
                if variant in ('kernel_and_distance','all_structural_tests'):
                    if any(model.cost[offset+i]+model.distances[x]<model.d for i,x in enumerate(newvalues)):
                        stats['distance_rejections']+=1;continue
                cols=columns+[v]
                if variant=='all_structural_tests':
                    j=len(cols)
                    if j<=p:
                        ann=nullspace(cols,model.b)
                        products=[image(c,model.Gamma) for c in ann]
                        restricted=[sum(((x&y).bit_count()&1)<<k for k,y in enumerate(products)) for x in ann]
                        if len(reduced_basis(restricted))>2*(p-j):
                            stats['partial_Gram_rank_rejections']+=1;continue
                    else:
                        pairs=j-p
                        partial=pullback_gram(cols[:pairs]+cols[p:],model.b)
                        residual=[x^y for x,y in zip(model.Gamma,partial)]
                        if len(reduced_basis(residual))>2*(p-pairs):
                            stats['remaining_Gram_rank_rejections']+=1;continue
                dfs(cols,values+newvalues)
        dfs([], [0])
    return dict(variant=variant,seconds=time.perf_counter()-start,**stats,solutions=len(solutions),solution_sha256=digest(solutions)),solutions

def physical_check(model,columns):
    """Check every accepted graph as physical Pauli words, not coset distances."""
    def words(gens):
        result={0}
        for g in gens:result|={x^g for x in result}
        return result
    def embed(pre,suf):
        x=(pre&((1<<model.m)-1))|((suf&((1<<model.f)-1))<<model.m)
        z=(pre>>model.m)|((suf>>model.f)<<model.m)
        return x|(z<<model.n)
    def wt(v):return ((v&((1<<model.n)-1))|(v>>model.n)).bit_count()
    D=[embed(r,0) for r in model.R]+[embed(0,c) for c in model.C]+[embed(e,model.representatives[s]) for e,s in zip(model.E,columns)]
    Dwords=words(D);assert len(Dwords)==1<<(model.n+model.k)
    Rwords=words([embed(r,0) for r in model.R])
    assert all(wt(w)>=model.d for w in Dwords-Rwords)
    # An ordinary dot-product nullspace of symplectically swapped graph rows.
    rows=[(w>>model.n)|((w&((1<<model.n)-1))<<model.n) for w in D]
    pivot={}
    for row in rows:
        while row:
            bit=row.bit_length()-1
            if bit in pivot:row^=pivot[bit]
            else:pivot[bit]=row;break
    free=[j for j in range(2*model.n) if j not in pivot]
    S=[]
    for j in free:
        v=1<<j
        for k in sorted(pivot):
            if (pivot[k]&v).bit_count()&1:v|=1<<k
        S.append(v)
    Swords=words(S);assert len(Swords)==1<<(model.n-model.k)
    assert Swords<=Dwords
    assert Rwords<=Swords
    local_words={w for w in Swords if ((w&((1<<model.n)-1))>>model.m)==0 and (w>>(model.n+model.m))==0}
    assert local_words==Rwords
    def commute(u,v):
        mask=(1<<model.n)-1
        left=((u&mask)&(v>>model.n)).bit_count()
        right=((v&mask)&(u>>model.n)).bit_count()
        return (left+right)%2==0
    assert all(commute(u,v) for u in S for v in S)

def cases():
    R=[15,15<<4]
    return [
        ('five_qubits_d2',ExtensionModel(5,1,2,4,R,[])),
        ('five_qubits_d3',ExtensionModel(5,1,3,4,R,[])),
        ('five_qubits_d4_impossible',ExtensionModel(5,1,4,4,R,[])),
        ('six_qubits_d3',ExtensionModel(6,0,3,4,R,[])),
        ('six_qubits_d4',ExtensionModel(6,0,4,4,R,[])),
        ('six_qubits_d5_impossible',ExtensionModel(6,0,5,4,R,[])),
        ('seven_qubits_d3',ExtensionModel(7,1,3,4,R,[7,7<<3])),
        ('six_qubits_one_suffix_generator_d2',ExtensionModel(6,1,2,4,R,[3<<2])),
    ]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repeats',type=int,default=3);args=p.parse_args()
    assert args.repeats>0 and __debug__
    output=ROOT/'benchmarks/results';output.mkdir(exist_ok=True)
    variants=['full_map_baseline','kernel_only','kernel_and_distance','all_structural_tests']
    report={'status':'RUNNING','scope':__doc__,'python':sys.version,'platform':platform.platform(),'repeats':args.repeats,'cases':[]}
    for name,model in cases():
        assert model.C_distance is None or model.C_distance>=model.d
        case={'name':name,'parameters':{'n':model.n,'k':model.k,'target_distance':model.d,'local_sites':model.m,'R':model.R,'C':model.C,'domain_dimension':model.a,'syndrome_dimension':model.b},'all_possible_maps':1<<(model.a*model.b),'prefix_cost_histogram':dict(collections.Counter(model.cost)),'variants':[]}
        reference=None
        for variant in variants:
            timings=[];deterministic=None
            for repeat in range(args.repeats):
                result,solutions=run(model,variant);timings.append(result['seconds'])
                stable={k:v for k,v in result.items() if k!='seconds'}
                if deterministic is None:deterministic=stable
                assert stable==deterministic
                if reference is None:reference=set(solutions)
                assert set(solutions)==reference,(name,variant)
            result['seconds_repetitions']=timings;result['median_seconds']=statistics.median(timings)
            case['variants'].append(result)
            print(name,variant,'tested',result['complete_maps_tested'],'solutions',result['solutions'],'median',round(result['median_seconds'],4),flush=True)
        for columns in reference:physical_check(model,columns)
        case['accepted_maps_checked_as_physical_Pauli_spaces']=len(reference)
        case['all_variants_identical_solution_sets']=True
        report['cases'].append(case)
        (output/'framework_comparison.json').write_text(json.dumps(report,indent=2))
    report['status']='PASS';report['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output/'framework_comparison.json').write_text(json.dumps(report,indent=2))
    print('COMPLETE MATCHED COMPARISON PASS',flush=True)
if __name__=='__main__':main()
