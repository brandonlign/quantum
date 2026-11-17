"""Standalone exact audit of the FINAL 1999 IEEE Example 1, printed p.2496.
Source: https://biblio.telecom-paristech.fr/cgi-bin/download.cgi?id=9892
The displayed rows were independently inspected in a rendered final PDF.
Optional argument: the downloaded final PDF; requires pdftotext for extraction.
The default audit needs only Python 3 standard library and no input files.
"""
from pathlib import Path
from collections import Counter
import hashlib,json,re,subprocess,sys
if not __debug__: raise RuntimeError('Run without -O; assertions perform verification')
FINAL_ROWS='''111111000000000000
001011100000101000
001001011000110000
011100000100110000
011010000000011100
011001001010001000
001000001100101010
001000001100010101
011000000111000100
110000101000000000
000000001000110001
000000010100001001'''.splitlines()
PDF_SHA256='9e9f47adf47e625563ca7f2c3f279d6ceb47a1cdf55a05c22ce36acc14f67395'
def rref(rows,width):
    a=rows[:];p=[];i=0
    for j in range(width):
        k=next((k for k in range(i,len(a)) if (a[k]>>j)&1),None)
        if k is None:continue
        a[i],a[k]=a[k],a[i]
        for k in range(len(a)):
            if k!=i and (a[k]>>j)&1:a[k]^=a[i]
        p.append(j);i+=1
    return a[:i],p

def inspect(rowstrings):
    rows=[int(s,2) for s in rowstrings]
    rr,piv=rref(rows,18)
    assert len(rows)==len(rr)
    gram=[sum((((a&b).bit_count()&1)<<j) for j,b in enumerate(rows)) for a in rows]
    gg,gp=rref(gram,len(rows))
    words=[0]
    for r in rows:words.extend([w^r for w in words])
    assert len(set(words))==1<<len(rows)
    counts=Counter(w.bit_count() for w in words)
    assert all(w%2==0 for w in counts)
    rank=len(gg)
    assert rank%2==0
    rad=[]
    for free in range(len(rows)):
        if free in gp:continue
        x=(1<<free)|sum(((r>>free)&1)<<j for r,j in zip(gg,gp))
        assert all((r&x).bit_count()%2==0 for r in gram)
        rad.append(x)
    return {'dimension':len(rr),'minimum_nonzero_distance':min(w for w in counts if w),'weight_enumerator':dict(sorted(counts.items())),'Euclidean_Gram_rows':gram,'Euclidean_Gram_rank':rank,'Euclidean_Gram_pivot_columns_zero_based':gp,'radical_coefficient_basis':rad,'radical_dimension':len(rad),'maximum_totally_isotropic_dimension':len(rows)-rank//2}

report={'status':'PASS','publication':{'title':'On Binary Constructions of Quantum Codes','authors':['Gerard Cohen','Sylvia Encheva','Simon Litsyn'],'journal':'IEEE Transactions on Information Theory','volume':45,'issue':7,'pages':'2495–2498','year':1999,'DOI':'10.1109/18.796389','public_final_PDF_url':'https://biblio.telecom-paristech.fr/cgi-bin/download.cgi?id=9892','PDF_bytes':128437,'PDF_SHA256':PDF_SHA256,'inspected_page':2496,'PDF_page_index_zero_based':1},'final_printed_rows':FINAL_ROWS,'first_nine':inspect(FINAL_ROWS[:9]),'full_twelve':inspect(FINAL_ROWS),'first_two_overlap_positions_one_based':[i+1 for i,(a,b) in enumerate(zip(FINAL_ROWS[0],FINAL_ROWS[1])) if a==b=='1'],'theorem_condition':'Theorem 1 requires C^perp subseteq C and C subseteq Cprime. At length 18 and dimension 9, the first condition forces self-duality.','conclusion':'The FINAL published extension matrix has no 9-dimensional totally isotropic subcode. Therefore its displayed data do not satisfy the hypothesis of the theorem invoked for Example 1. The same displayed-data obstruction persists from the preprint into the final article.','limits':'This does not establish the absence of a later erratum, exclude other constructions, invalidate all results in the article, or settle priority.'}
assert report['first_nine']['dimension']==9 and report['first_nine']['minimum_nonzero_distance']==6
assert report['first_nine']['Euclidean_Gram_rank']==8
assert report['full_twelve']['dimension']==12 and report['full_twelve']['minimum_nonzero_distance']==4
assert report['full_twelve']['Euclidean_Gram_rank']==10
assert report['full_twelve']['radical_dimension']==2
assert report['full_twelve']['maximum_totally_isotropic_dimension']==7
# Optional direct verification against the identified final published PDF.
if len(sys.argv)>1:
    pdf=Path(sys.argv[1]);data=pdf.read_bytes()
    assert len(data)==128437 and hashlib.sha256(data).hexdigest()==PDF_SHA256
    txt=subprocess.check_output(['pdftotext','-f','2','-l','2','-layout',str(pdf),'-'],text=True)
    extracted=re.findall(r'(?<![01])[01]{18}(?![01])',txt)
    assert extracted==FINAL_ROWS,(extracted,FINAL_ROWS)
    report['optional_direct_final_PDF_check']={'SHA256_match':True,'all_twelve_extracted_rows_match':True,'page':2496}
# Compare with independently frozen preprint audit when present.
old=Path(__file__).with_name('CEL_PREPRINT_MATRIX_AUDIT.json')
if old.exists():
    assert json.loads(old.read_text())['rows']==FINAL_ROWS
    report['all_twelve_rows_identical_to_preprint_audit']=True
out=Path(__file__).with_name('CEL_PUBLISHED_MATRIX_AUDIT.json')
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
