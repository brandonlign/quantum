#!/usr/bin/env python3
"""Extract byte-preserved records for use without the storage adapter."""
from pathlib import Path
import argparse,json,time
from record_store import RecordStore

ROOT=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix',default='',help='Optional stored relative prefix, e.g. oracle_residual/runs/')
    args=p.parse_args();store=RecordStore(ROOT);count=0;start=time.monotonic()
    if '..' in Path(args.prefix).parts:raise ValueError('Invalid prefix')
    cursor=store.db.execute('SELECT path FROM files WHERE path>=? AND path<? ORDER BY path',(args.prefix,args.prefix+'\uffff'))
    for (key,) in cursor:
        target=store.base/key;target.resolve().relative_to(store.base)
        data=store.read(key)
        if target.exists():
            if target.read_bytes()!=data:raise ValueError('Conflicting existing file: '+key)
        else:
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        count+=1
        if count%100000==0:print('Extracted',count,'files',flush=True)
    print(json.dumps({'files':count,'elapsed_seconds':time.monotonic()-start,'prefix':args.prefix}),flush=True)
if __name__=='__main__':main()
