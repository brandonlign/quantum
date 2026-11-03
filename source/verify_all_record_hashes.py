#!/usr/bin/env python3
"""Check every archived payload and every hash in the frozen residual manifest.

This is a complete archival integrity/status check. It is distinct from the
mathematical per-record reconstruction performed by verify.py --deep.
"""
from pathlib import Path
import collections,hashlib,itertools,json,re,sqlite3,time,zlib

ROOT=Path(__file__).resolve().parent
BASE=ROOT/'experimenting/qcode-open-entries/general-criterion/unified'
def main():
    start=time.monotonic();saved_path=BASE/'oracle_residual/aggregate.json'
    original_sha=hashlib.sha256(saved_path.read_bytes()).hexdigest()
    expected_original='9aeae1ae01235e42ee92d1f854346c0ea7d8c2c281c3a264e23ed26df89083d0'
    assert original_sha==expected_original,'Frozen aggregate differs from original evidence audit'
    saved=json.loads(saved_path.read_text());assert saved['counts']['EXCLUDED_SLICE_CHUNK']==2419796
    expected={}
    marker='/experimenting/qcode-open-entries/general-criterion/unified/'
    for name,h in saved['record_hashes'].items():
        assert marker in name
        rel=name.split(marker,1)[1];assert rel not in expected;expected[rel]=h
    assert len(expected)==2419796
    assignment=json.loads((BASE/'oracle_residual/assignment.json').read_text())
    print('Assignment keys',list(assignment),flush=True)
    remote=set(assignment['remote']);local=set(assignment['local'])
    assert not local&remote and len(local|remote)==32747
    store=ROOT/'evidence/search_records.sqlite'
    db=sqlite3.connect('file:'+str(store)+'?mode=ro',uri=True)
    counts=collections.Counter();classes=set();bytes_checked=0;matched=0;statuses=collections.Counter()
    legacy_names=[name for name in expected if name.startswith('rankone_residual/runs/')]
    assert legacy_names and all((BASE/name).is_file() for name in legacy_names)
    print('Additional materialized legacy manifest records',len(legacy_names),flush=True)
    def legacy_rows():
        for name in legacy_names:
            data=(BASE/name).read_bytes()
            yield name,0,zlib.crc32(data)&0xffffffff,len(data),data
    for name,method,crc,size,packed in itertools.chain(db.execute('SELECT path,method,crc,nbytes,packed FROM files'),legacy_rows()):
        data=zlib.decompress(packed,-15) if method==8 else packed
        assert len(data)==size and zlib.crc32(data)&0xffffffff==crc,name
        counts['all_stored_files']+=1;bytes_checked+=size
        counts['materialized_legacy' if name.startswith('rankone_residual/runs/') else 'indexed_payloads']+=1
        if name in expected:
            assert hashlib.sha256(data).hexdigest()==expected.pop(name),name
            r=json.loads(data)
            assert r['status']=='EXCLUDED_SLICE_CHUNK',name
            assert r['z']['status']=='COMPLETE' and r['z'].get('reason') is None,name
            assert r['z']['leaves']==len(r['leaves']),name
            assert r.get('found') is None and 'candidate_file' not in r and 'physical' not in r,name
            cls=int(name.rsplit('/',1)[1].split('_',1)[0]);assert cls==r['cls']
            assert cls in (remote if name.startswith('oracle_residual/runs/') else local),name
            for leaf in r['leaves']:
                result=leaf['result']
                assert result['status']=='COMPLETE' and result['complete'] and result['fixed_z_excluded'] and not result['solutions'],name
                assert result['reason']!='search_limit',name
            classes.add(cls);statuses[r['status']]+=1;matched+=1
        if counts['all_stored_files']%100000==0:
            print('Verified',counts['all_stored_files'],'payloads;',matched,'residual manifest hashes;',round(time.monotonic()-start),'s',flush=True)
    assert not expected,list(expected)[:5]
    assert matched==2419796 and len(classes)==32747
    assert classes==local|remote
    result={'status':'PASS','stored_payloads_size_and_CRC_checked':counts['indexed_payloads'],
            'materialized_legacy_records_SHA256_checked':counts['materialized_legacy'],
            'residual_manifest_SHA256_matches':matched,'residual_classes':len(classes),
            'uncompressed_bytes_checked':bytes_checked,'source_aggregate_sha256':original_sha,
            'completed_record_statuses':dict(statuses),'elapsed_seconds':time.monotonic()-start,
            'scope':__doc__}
    (ROOT/'provenance').mkdir(parents=True,exist_ok=True)
    (ROOT/'provenance/all_record_integrity.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
