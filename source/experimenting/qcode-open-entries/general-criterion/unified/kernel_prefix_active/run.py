"""Cache which immutable prefix witnesses apply at each DFS depth."""
import json
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from kernel_gamma_cache import run as harness
from mixed_frames.stage3.common import atomic_json,sha

source=(ROOT/'kernel_z_compact/z_enum_compact.c').read_text()
start=source.index('static int prefix_canonical(int depth){')
end=source.index('\nstatic int par(',start)
harness.HERE=Path(__file__).resolve().parent
harness.OLD=source[start:end]
harness.NEW='''/* The prefix action list is immutable after input admission. Cache
 * exactly the original applicability predicate for each depth. */
static unsigned short active_witness[MAX_M][MAX_M*(1<<MAX_M)];
static int active_count[MAX_M],active_ready[MAX_M];
static int prefix_canonical(int depth){
    if(!active_ready[depth]){
        int assigned=(1<<depth)-1;
        for(int k=0;k<witness_count;k++){
            int i=witness_index[k],c=witness_mask[k];
            if(i>=depth || (c&~assigned))continue;
            active_witness[depth][active_count[depth]++]=(unsigned short)k;
        }
        active_ready[depth]=1;
    }
    int last=depth-1;
    for(int z=0;z<active_count[depth];z++){
        int k=active_witness[depth][z];
        int i=witness_index[k],c=witness_mask[k];
        int value=sums[last][c&((1<<last)-1)]^((c>>last&1)?xs[last]:0);
        if(value<xs[i])return 0;
    }
    return 1;
}'''


if __name__=='__main__':
    harness.main()
    path=harness.HERE/'benchmark.json'
    d=json.loads(path.read_text())
    d['hashes'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    d['scope']='Bounded complete stream/counter comparisons for caching the exact depth-only witness predicate and using the already maintained subset-sum table. No production deployment or new exclusion.'
    atomic_json(path,d)
