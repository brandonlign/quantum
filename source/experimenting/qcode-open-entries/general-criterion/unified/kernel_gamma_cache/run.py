"""Bounded exact-stream benchmark of memoized fixed Gamma multiplication."""
import hashlib
import json
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from mixed_frames.stage3.common import atomic_json,sha

OLD='static int gapply(int c) { int o = 0; for (int i = 0; i < bdim; ++i) o |= par(gam[i] & c) << i; return o; }'
NEW='''/* Gamma is read once and never changes within this process. N<=8192,
 * all coefficients are b-bit words, and every image is in [0,N). Zero
 * denotes an uncached entry; value+1 fits uint16_t including N=8192. */
static uint16_t gamma_image_cache[MAX_N];
static int gapply(int c) {
    unsigned cached=gamma_image_cache[c];
    if(cached)return (int)cached-1;
    int o=0;
    for(int i=0;i<bdim;++i)o|=par(gam[i]&c)<<i;
    gamma_image_cache[c]=(uint16_t)(o+1);
    return o;
}'''


def call(binary,f):
    with tempfile.TemporaryDirectory() as td:
        gp=Path(td)/'group.txt';gp.write_text(f['group'])
        start=time.monotonic()
        p=subprocess.run([str(binary),*f['argv'],'--prefix-group',str(gp)],
            input=f['stdin'],text=True,capture_output=True,check=True,timeout=180)
        elapsed=time.monotonic()-start
    st=json.loads(p.stderr)
    return dict(seconds=elapsed,native=st,stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest())


def main():
    fixtures=json.loads((ROOT/'oracle_residual/fixtures.json').read_text())
    paths={'compact':ROOT/'kernel_z_compact/z_enum_compact',
           'joint':ROOT/'kernel_z_joint/z_enum_joint'}
    frozen={str(p.relative_to(ROOT)):sha(p) for p in
        [*paths.values(),*(p.with_suffix('.c') for p in paths.values())]}
    for name,old in paths.items():
        source=old.with_suffix('.c').read_text();assert source.count(OLD)==1
        new=HERE/f'{name}.c';new.write_text(source.replace(OLD,NEW))
        subprocess.run(['cc','-O3','-std=c11',str(new),'-o',str(HERE/name)],check=True)
    records=[]
    for round_number in (0,1):
        for f in fixtures['cases']:
            results={}
            for mode in (['baseline','cached'] if round_number==0 else ['cached','baseline']):
                binary=paths[f['backend']] if mode=='baseline' else HERE/f['backend']
                results[mode]=call(binary,f)
                assert results[mode]['stdout_sha256']==f['expected_stdout_sha256']
                assert {k:v for k,v in results[mode]['native'].items() if k!='seconds'}==f['expected_stats']
            r=dict(round=round_number,cls=f['cls'],profile=f['profile'],backend=f['backend'],
                   results=results)
            records.append(r)
            print(json.dumps(dict(round=round_number,cls=f['cls'],profile=f['profile'],
                ratio=results['cached']['seconds']/results['baseline']['seconds'],
                seconds={k:v['seconds'] for k,v in results.items()},status='COMPLETE_STREAM_AND_COUNTER_MATCH')),flush=True)
            atomic_json(HERE/'benchmark.json',dict(status='BOUNDED_CACHE_BENCHMARK_IN_PROGRESS',records=records))
    assert frozen=={p:sha(ROOT/p) for p in frozen}
    ratios=[r['results']['cached']['seconds']/r['results']['baseline']['seconds'] for r in records]
    report=dict(status='BOUNDED_CACHE_BENCHMARK_COMPLETE',records=records,
        summary=dict(median_ratio=statistics.median(ratios),total_ratio=
            sum(r['results']['cached']['seconds'] for r in records)/sum(r['results']['baseline']['seconds'] for r in records)),
        frozen_hashes=frozen,hashes={str(p.relative_to(ROOT)):sha(p) for p in
            [Path(__file__),HERE/'compact.c',HERE/'joint.c',HERE/'compact',HERE/'joint']},
        scope='Exact same finite search with a pure-function lookup cache. Bounded fixture controls only; no production deployment or new exclusion.')
    atomic_json(HERE/'benchmark.json',report)
    print(json.dumps(dict(status=report['status'],summary=report['summary'])),flush=True)


if __name__=='__main__':main()
