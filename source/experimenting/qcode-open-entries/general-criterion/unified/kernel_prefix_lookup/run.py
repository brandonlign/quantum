"""Use the DFS's already computed subset sums for prefix witnesses."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from kernel_gamma_cache import run as harness
from mixed_frames.stage3.common import atomic_json,sha
import json

harness.HERE=Path(__file__).resolve().parent
harness.OLD='int value=0;for(int j=0;j<depth;j++)if(c>>j&1)value^=xs[j];'
harness.NEW='''/* At this call xs[depth-1] has just been assigned; sums[depth-1]
         * contains every subset sum of the earlier columns, including zero. */
        int last=depth-1;
        int value=sums[last][c&((1<<last)-1)]^((c>>last&1)?xs[last]:0);'''


if __name__=='__main__':
    harness.main()
    path=harness.HERE/'benchmark.json'
    d=json.loads(path.read_text())
    d['hashes'][str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
    d['scope']='Bounded complete stream/counter comparisons for replacing repeated prefix XOR loops with the already maintained subset-sum table. No production deployment or new exclusion.'
    atomic_json(path,d)
