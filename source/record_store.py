"""Read-only access to byte-preserved search files in the compact evidence store.

The mathematical checkers remain unchanged. This module only supplies their
original file bytes and relocates historical absolute paths to this folder.
Run materialize_records.py to use an ordinary filesystem instead.
"""
from __future__ import annotations
import builtins
import contextlib
import functools
import io
import json
import os
from pathlib import Path
import sqlite3
import stat
import zlib

class RecordStore:
    def __init__(self, root, database=None):
        self.root=Path(root).resolve()
        self.base=self.root/'experimenting/qcode-open-entries/general-criterion/unified'
        self.database=None if database is False else Path(database or self.root/'evidence/search_records.sqlite')
        self.db=None if self.database is None else sqlite3.connect('file:'+str(self.database)+'?mode=ro',uri=True)
        if self.db is not None:self.db.execute('PRAGMA cache_size=-131072')
        self.original={n:getattr(os,n) for n in ['stat','lstat','scandir']}
        self.original['open']=builtins.open
        self.original['io_open']=io.open
        self.original['json_loads']=json.loads
        self.reads=0

    def relocate(self,path):
        if isinstance(path,int):return path
        p=os.fsdecode(os.fspath(path))
        # The current workspace can itself live below a historical workspace
        # root (for example Desktop/research/quantum). Preserve paths already
        # inside this root before translating legacy absolute paths.
        current=str(self.root)
        if p==current or p.startswith(current+'/'):return p
        for old in ['/Users/Brandon/Desktop/research','/Users/brandonli/Desktop/research']:
            if p==old or p.startswith(old+'/'):
                return current+p[len(old):]
        return p

    def key(self,path):
        if isinstance(path,int):return None
        p=os.path.abspath(self.relocate(path))
        base=str(self.base)+'/'
        return p[len(base):] if p.startswith(base) else None

    @functools.lru_cache(maxsize=8192)
    def info(self,key):
        if key is None or self.db is None:return None
        return self.db.execute('SELECT method,crc,nbytes FROM files WHERE path=?',(key,)).fetchone()

    @functools.lru_cache(maxsize=8192)
    def directory(self,key):
        if key is None or self.db is None:return False
        prefix=key.rstrip('/')+'/'
        return self.db.execute('SELECT 1 FROM files WHERE path>=? AND path<? LIMIT 1',(prefix,prefix+'\uffff')).fetchone() is not None

    @functools.lru_cache(maxsize=128)
    def read(self,key):
        if self.db is None:raise FileNotFoundError(key)
        row=self.db.execute('SELECT method,crc,nbytes,packed FROM files WHERE path=?',(key,)).fetchone()
        if row is None:raise FileNotFoundError(key)
        method,crc,size,packed=row
        data=zlib.decompress(packed,-15) if method==8 else packed
        if len(data)!=size or zlib.crc32(data)&0xffffffff!=crc:
            raise ValueError('Corrupt stored record: '+key)
        self.reads+=1
        return data

    def open(self,file,mode='r',buffering=-1,encoding=None,errors=None,newline=None,closefd=True,opener=None):
        if isinstance(file,int):
            return self.original['io_open'](file,mode,buffering,encoding,errors,newline,closefd,opener)
        mapped=self.relocate(file)
        try:return self.original['io_open'](mapped,mode,buffering,encoding,errors,newline,closefd,opener)
        except FileNotFoundError:
            key=self.key(mapped)
            if self.info(key) is None:raise
            if any(c in mode for c in 'wax+'):
                raise PermissionError('Frozen evidence is read-only: '+str(file))
            stream=io.BytesIO(self.read(key))
            return stream if 'b' in mode else io.TextIOWrapper(stream,encoding=encoding or 'utf-8',errors=errors,newline=newline)

    def stat(self,path,*args,**kwargs):
        mapped=self.relocate(path)
        try:return self.original['stat'](mapped,*args,**kwargs)
        except FileNotFoundError:
            key=self.key(mapped);info=self.info(key)
            if info is not None:return os.stat_result((stat.S_IFREG|0o444,0,0,1,0,0,info[2],0,0,0))
            if self.directory(key):return os.stat_result((stat.S_IFDIR|0o555,0,0,1,0,0,0,0,0,0))
            raise

    def lstat(self,path,*args,**kwargs):
        kwargs.pop('follow_symlinks',None)
        return self.stat(path,*args,follow_symlinks=False,**kwargs)

    def scandir(self,path='.'):
        if isinstance(path,int):return self.original['scandir'](path)
        return _Scan(self,path)

    @contextlib.contextmanager
    def overlay(self):
        store=self
        class HistoricalPathKeys(dict):
            def __missing__(self,key):
                if isinstance(key,str) and key.startswith(str(store.root)+'/'):
                    suffix=key[len(str(store.root)):]
                    matches=[dict.__getitem__(self,old+suffix) for old in
                             ['/Users/Brandon/Desktop/research','/Users/brandonli/Desktop/research']
                             if dict.__contains__(self,old+suffix)]
                    if matches:
                        if any(x!=matches[0] for x in matches):raise ValueError('Conflicting historical path hashes: '+key)
                        return matches[0]
                raise KeyError(key)
        def load_json(*args,**kwargs):
            # Preserve every stored key/value and raw byte hash. Only missing
            # lookups under the relocated root may use the historical key.
            if kwargs.get('object_pairs_hook') is None and kwargs.get('object_hook') is None:
                kwargs['object_pairs_hook']=HistoricalPathKeys
            return self.original['json_loads'](*args,**kwargs)
        builtins.open=self.open;io.open=self.open
        json.loads=load_json
        os.stat=self.stat;os.lstat=self.lstat;os.scandir=self.scandir
        try:yield self
        finally:
            builtins.open=self.original['open'];io.open=self.original['io_open']
            json.loads=self.original['json_loads']
            for n in ['stat','lstat','scandir']:setattr(os,n,self.original[n])

class _Entry:
    def __init__(self,store,parent,name,is_dir):
        self.store=store;self.name=name;self.path=os.path.join(parent,name);self.isdir=is_dir
    def is_dir(self,follow_symlinks=True):return self.isdir
    def is_file(self,follow_symlinks=True):return not self.isdir
    def is_symlink(self):return False
    def stat(self,follow_symlinks=True):return self.store.stat(self.path)
    def inode(self):return 0
    def __fspath__(self):return self.path

class _Scan:
    def __init__(self,store,path):
        self.store=store;self.path=store.relocate(path);self.it=self.generate()
    def generate(self):
        seen=set()
        try:
            with self.store.original['scandir'](self.path) as entries:
                for entry in entries:
                    seen.add(entry.name);yield entry
        except FileNotFoundError:
            if not self.store.directory(self.store.key(self.path)):raise
        key=self.store.key(self.path)
        if key is None or self.store.db is None:return
        prefix=key.rstrip('/')+'/'
        cursor=self.store.db.execute('SELECT path FROM files WHERE path>=? AND path<?',(prefix,prefix+'\uffff'))
        for (stored,) in cursor:
            tail=stored[len(prefix):];name,sep,_=tail.partition('/')
            if name not in seen:
                seen.add(name);yield _Entry(self.store,self.path,name,bool(sep))
    def __iter__(self):return self
    def __next__(self):return next(self.it)
    def close(self):self.it.close()
    def __enter__(self):return self
    def __exit__(self,*args):self.close()
