#!/usr/bin/env python3
"""Run the unchanged quantum evidence checker using portable record storage."""
from pathlib import Path
import argparse,json,runpy,sys,time
from record_store import RecordStore

ROOT=Path(__file__).resolve().parent
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--deep',action='store_true',help='Revalidate every saved residual record, not only the default sample.')
    parser.add_argument('--script',default='quantum-extension/verification/verify_evidence.py')
    parser.add_argument('--filesystem',action='store_true',help='Use materialized records; only relocate historical file paths.')
    args=parser.parse_args()
    if not __debug__:raise RuntimeError('Run without Python -O')
    script=(ROOT/args.script).resolve();script.relative_to(ROOT)
    store=RecordStore(ROOT,False if args.filesystem else None)
    sys.path.insert(0,str(script.parent))
    sys.argv=[str(script)]+(['--deep'] if args.deep else [])
    started=time.monotonic()
    with store.overlay():runpy.run_path(str(script),run_name='__main__')
    report={'status':'PASS','script':args.script,'deep':args.deep,'record_payloads_read':store.reads,'elapsed_seconds':time.monotonic()-started,
            'scope':'Original mathematical verifier and source bytes unchanged. Historical absolute file reads relocated; original search-record bytes supplied from an indexed store with size/CRC checks. This revalidates saved computations; it does not rerun every original search.'}
    (ROOT/'provenance').mkdir(parents=True,exist_ok=True)
    (ROOT/'provenance'/('portable_deep_replay.json' if args.deep else 'portable_replay.json')).write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
