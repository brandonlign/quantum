"""Acceptance/rejection checks of the storage adapter, separate from mathematics."""
from pathlib import Path
import hashlib,json,sqlite3,tempfile,zlib
import os
from record_store import RecordStore

def main():
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);database=root/'fixture.sqlite'
        paths={'rankone_reuse/runs/1_ok.json':b'{"status":"COMPLETE","solutions":[]}',
               'oracle_residual/runs/2_ok.json':b'{"value":17}',
               'mixed_frames/stage4/668_1/data.txt':b'leaf 1 2 3\n'}
        with sqlite3.connect(database) as db:
            db.execute('CREATE TABLE files(path TEXT PRIMARY KEY,method INT,crc INT,nbytes INT,packed BLOB)')
            for name,data in paths.items():
                c=zlib.compressobj(wbits=-15);packed=c.compress(data)+c.flush()
                db.execute('INSERT INTO files VALUES(?,?,?,?,?)',(name,8,zlib.crc32(data)&0xffffffff,len(data),packed))
            db.execute('INSERT INTO files VALUES(?,?,?,?,?)',('oracle_residual/runs/3_corrupt.json',0,123,3,b'bad'))
        store=RecordStore(root,database);base=store.base
        original=root/'unchanged.txt';original.write_text('ordinary file')
        # Exercise the actual nested Desktop/research/quantum layout, where
        # the canonical root is also below a historical relocation prefix.
        nested=RecordStore(Path('/Users/brandonli/Desktop/research/quantum'),database)
        current_nested=str(nested.root/'quantum-extension/audit/result.json')
        old_nested='/Users/Brandon/Desktop/research/experimenting/result.json'
        mapped_nested=nested.relocate(old_nested)
        assert nested.relocate(current_nested)==current_nested
        assert mapped_nested==str(nested.root)+'/experimenting/result.json'
        assert nested.relocate(mapped_nested)==mapped_nested
        assert nested.relocate('/var/tmp/record-store-unrelated.json')=='/var/tmp/record-store-unrelated.json'
        assert nested.relocate('relative/unrelated.json')=='relative/unrelated.json'
        nested.db.close()
        with store.overlay():
            for name,data in paths.items():
                p=base/name
                assert p.is_file() and p.exists() and not p.is_dir()
                assert p.read_bytes()==data and p.read_text()==data.decode()
                assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(data).digest()
                assert p.stat().st_size==len(data)
                with open(p,'rb') as f:assert f.read()==data
            names=sorted(p.name for p in (base/'oracle_residual/runs').glob('*.json'))
            assert names==['2_ok.json','3_corrupt.json']
            assert not (base/'oracle_residual/runs/missing.json').exists()
            assert original.read_text()=='ordinary file'
            descriptor=os.open(root,os.O_RDONLY)
            try:
                with os.scandir(descriptor) as entries:assert 'unchanged.txt' in {e.name for e in entries}
            finally:os.close(descriptor)
            try:(base/'oracle_residual/runs/3_corrupt.json').read_bytes()
            except ValueError:pass
            else:raise AssertionError('Corrupt CRC accepted')
            # Historical absolute paths and stored JSON keys stay compatible.
            assert Path('/Users/Brandon/Desktop/research/unchanged.txt').read_text()=='ordinary file'
            encoded='{"/Users/Brandon/Desktop/research/input.py":"unchanged-hash"}'
            decoded=json.loads(encoded)
            assert decoded[str(store.root/'input.py')]=='unchanged-hash'
            assert list(decoded.items())==[('/Users/Brandon/Desktop/research/input.py','unchanged-hash')]
        assert original.read_text()=='ordinary file'
    report={'status':'PASS','checks':['original byte/text reads','builtins open','file existence/type/size','glob','missing-file rejection','CRC-corruption rejection','historical path relocation','new-path idempotence','old-to-new relocation idempotence','unrelated path preservation','frozen historical JSON-key lookup','ordinary-file fallback','hooks restored on exit']}
    root=Path(__file__).resolve().parent
    (root/'provenance').mkdir(parents=True,exist_ok=True)
    (root/'provenance/storage_adapter_checks.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
