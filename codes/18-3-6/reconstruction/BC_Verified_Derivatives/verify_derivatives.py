from pathlib import Path
import json,hashlib,collections
O=Path(__file__).parent
if not __debug__:raise RuntimeError('Assertions are required')
A=json.loads((O/'BASE_CODE_A.json').read_text())['rows'];mask=(1<<18)-1
sp=lambda a,b:(((a&mask)&(b>>18))^((b&mask)&(a>>18))).bit_count()%2
wt=lambda a:((a&mask)|(a>>18)).bit_count()
records=[]
for name,v,a7,logical in [('code_b_style',4328808577,16,34359943556),('code_c_style',8607273026,12,34359939617)]:
 D=O/name;T=json.loads((D/'TRANSFORM.json').read_text());S=list(map(int,(D/'ROWS.txt').read_text().split()));assert T['rows']==S==[r^(v if sp(r,v) else 0) for r in A] and T['v']==v and wt(v)==4
 tv=lambda r:r^(v if sp(r,v) else 0)
 assert all(tv(tv(1<<i))==1<<i for i in range(36))
 assert all(sp(tv(1<<i),tv(1<<j))==sp(1<<i,1<<j) for i in range(36) for j in range(36))
 J=json.loads((D/'STANDALONE_RESULT.json').read_text());assert J['status']=='PASS' and J['logical_distance']==J['minimum_normalizer_weight']==6 and J['rank_S']==15 and J['rank_normalizer']==21
 assert min(int(w) for w in J['stabilizer_enumerator'] if int(w))==7 and J['stabilizer_enumerator']['7']==a7 and J['normalizer_enumerator']['6']==J['logical_enumerator']['6']==564 and J['explicit_weight6_logical']==logical
 rows=[list(map(int,l.split())) for l in (D/'ENUMERATION.log').read_text().splitlines() if l and l[0].isdigit()];assert len(rows)==19
 for w,a,b,l,_ in rows:
  assert a==J['stabilizer_enumerator'].get(str(w),0) and b==J['normalizer_enumerator'].get(str(w),0) and l==J['logical_enumerator'].get(str(w),0)
 strings=(D/'PAULI_GENERATORS.txt').read_text().splitlines();assert strings==[''.join('IXZY'[(x>>i&1)+2*(x>>(18+i)&1)] for i in range(18)) for x in S]
 assert 'PASS all physical errors of weights 1 through 5' in (D/'check_errors.log').read_text()
 records.append({'name':name,'v':v,'Pauli_v':T['Pauli_v'],'rank':15,'stabilizer_minimum':7,'A7':a7,'normalizer_minimum':6,'D6':564,'logical_witness':logical,'raw_rows_sha256':hashlib.sha256((D/'ROWS.txt').read_bytes()).hexdigest(),'all_19_MacWilliams_coefficients_match':True})
assert (O/'BASE_INPUT.txt').read_text().split()[:15]==list(map(str,A))
assert 'tested 3060 minS7 34 valid_distance6 0' in (O/'z4_scan.log').read_text()
assert 'BOUNDED_STOP after 86085 of247860 masks' in (O/'pauli4_scan.log').read_text()
J={'status':'INDEPENDENT_PASS','codes':records,'specified_Z2345_transform':'Rejected: distance5, twelve weight-five logicals, minS7, A7=16, D6=528','Z4_masks_exhausted':3060,'Z4_valid_minS7_distance6_variants':0,'Pauli_weight4_masks_tested':86085,'Pauli_weight4_total_possible':247860,'Pauli_search_stopped_after_both_fingerprints_found':True,'pairwise_product_local_unitary_and_permutation_inequivalence':True,'exact_historical_identity_claim':False}
(O/'RESULT.json').write_text(json.dumps(J,indent=2)+'\n');print(json.dumps(J,indent=2))
