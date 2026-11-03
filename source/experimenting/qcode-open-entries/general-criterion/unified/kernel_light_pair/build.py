"""Generate separate compact/joint prototypes; never modify frozen sources."""
import json
import subprocess
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
from mixed_frames.stage3.common import atomic_json, sha

FUNCTION = r'''
/* Every entry represents a distinct physical suffix p, not just a quotient. */
static int pair_on,pair_budget,pair_cut,pair_purity,pair_count,pair_min_weight=999;
static int pair_a[MAX_N],pair_weight[MAX_N];
static uint64_t pair_pruned;
static int pair_ok(int depth){
    int minimum_x=m-depth;minimum_x-=minimum_x%2;
    int lower=m+minimum_x+2*pair_min_weight;
    if(lower>=2*pair_cut && lower>=2*pair_purity)return 1;
    int forced=0;
    for(int k=0;k<pair_count;k++){
        int zeros=0;
        for(int j=0;j<depth;j++)zeros+=!par(pair_a[k]&xs[j]);
        int max_x=m-zeros;max_x-=max_x%2;
        int sum=m+max_x+2*pair_weight[k];
        if(sum<2*pair_purity)return 0;
        if(sum<pair_cut+pair_purity)forced+=2;
        else if(sum<2*pair_cut)forced++;
        if(forced>pair_budget)return 0;
    }
    return 1;
}
'''


def once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)


def main():
    records = []
    for backend, source in [('compact', ROOT / 'kernel_z_compact/z_enum_compact.c'),
                            ('joint', ROOT / 'kernel_z_joint/z_enum_joint.c')]:
        text = source.read_text()
        text = once(text, 'static int gapply(int c)', FUNCTION + '\nstatic int gapply(int c)')
        text = once(text, '        } else if (!strcmp(argv[i], "--kernel-data")) {',
                    '''        } else if (!strcmp(argv[i], "--pair-budget") && i+1<argc) {
            char *end;long value=strtol(argv[++i],&end,10);
            if(!*argv[i]||*end||value<0||value>1000000)return error("bad pair budget");
            pair_on=1;pair_budget=(int)value;
        } else if (!strcmp(argv[i], "--pair-cut") && i+1<argc) {
            char *end;long value=strtol(argv[++i],&end,10);
            if(!*argv[i]||*end||value<1||value>100)return error("bad pair cutoff");
            pair_cut=(int)value;
        } else if (!strcmp(argv[i], "--pair-purity") && i+1<argc) {
            char *end;long value=strtol(argv[++i],&end,10);
            if(!*argv[i]||*end||value<1||value>100)return error("bad pair purity");
            pair_purity=(int)value;
        } else if (!strcmp(argv[i], "--kernel-data")) {''')
        text = once(text, '    word_count = (ZN + 63) / 64;',
                    '''    if(pair_on){
        if(pair_cut<1||pair_purity<1||pair_purity>pair_cut)return error("pair thresholds required");
        if(scanf("%d",&pair_count)!=1||pair_count<0||pair_count>MAX_N)return error("bad pair count");
        for(int k=0;k<pair_count;k++){
            if(scanf("%d %d",pair_a+k,pair_weight+k)!=2||pair_a[k]<=0||pair_a[k]>=N||pair_weight[k]<0)return error("bad pair entry");
            if(pair_weight[k]<pair_min_weight)pair_min_weight=pair_weight[k];
        }
    }
    word_count = (ZN + 63) / 64;''')
        text = once(text, '            if(use_projection && !projection_step(depth,x)){++projection_pruned;continue;}',
                    '            if(pair_on && !pair_ok(depth+1)){++pair_pruned;continue;}\n'
                    '            if(use_projection && !projection_step(depth,x)){++projection_pruned;continue;}')
        final = '    else fputs("null", stderr);'
        text = once(text, final, final + '\n    fprintf(stderr,' + json.dumps(',"pair_pruned":%') +
                    ' PRIu64,pair_pruned);')
        generated = HERE / f'z_{backend}.c'
        binary = HERE / f'z_{backend}'
        generated.write_text(text)
        cmd = ['cc', '-O3', '-std=c11', '-Wall', '-Wextra', str(generated), '-o', str(binary)]
        subprocess.run(cmd, check=True)
        records.append(dict(backend=backend, command=cmd,
                            hashes={str(p.relative_to(ROOT)): sha(p) for p in [source, generated, binary]}))
    atomic_json(HERE / 'build.json', dict(status='BUILT_UNVALIDATED_PROTOTYPES', records=records,
                build_source_sha256=sha(Path(__file__))))


if __name__ == '__main__':
    main()
