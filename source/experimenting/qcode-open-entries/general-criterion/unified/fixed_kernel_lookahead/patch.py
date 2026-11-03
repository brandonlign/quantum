"""Regenerate z_enum_lookahead.c from the frozen native source."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'fixed_kernel/z_enum_fixed.c').read_text()
def rep(a,b):
    global s
    assert s.count(a)==1,a[:60]; s=s.replace(a,b)
rep("/* Exact streaming canonical Z enumeration; reference: z_enumerator.py.",
"/* fixed_kernel_lookahead: copy of frozen fixed_kernel/z_enum_fixed.c plus an\n * optional sum-feasibility lookahead (--lookahead F, F in {2,3}).\n * Exact streaming canonical Z enumeration; reference: z_enumerator.py.")
rep("static uint64_t prefix_pruned;","static uint64_t prefix_pruned;\nstatic uint64_t depth_nodes[MAX_M+1];\nstatic int lookahead_free;\nstatic uint64_t look_tested,look_pruned;\nstatic uint64_t ldom[3][MAX_WORDS];")
helper=r'''
/* Lookahead (only with --kernel-data). At depth d the free sites are d..m-1;
 * x_{m-1}=xs[0]^...^xs[m-2], so the free values XOR to s=sums[d][last].
 * D_f = native tests on site f involving only assigned sites: subset-distance
 * rows (cost symmetric under complement, hence valid for site m-1 too), tag
 * parities, independence from the assigned span, fixed-form parity vs
 * assigned non-anchor sites (f<m-1, as native), --prefix pins. Any completion
 * picks f_k in D_k with XOR s, so an infeasible node has no completion. */
static int free_domain(int depth,int span_size,int site,uint64_t *out){
    for(int w=0;w<word_count;w++)out[w]=UINT64_MAX;
    if(N<64)out[0]=(UINT64_C(1)<<N)-1;
    for(int j=0;j<span_size;j++){
        int a=(int)sizes[depth][j]+1,t=a<m-a?a:m-a;
        if(distance-t<=0)continue;
        const uint64_t *row=translated_row(t,sums[depth][j]);
        for(int w=0;w<word_count;w++)out[w]&=row[w];
    }
    int tag=site_tags[site];
    for(int j=0;j<kernel_g;j++){
        const uint64_t *mask=parity_allowed+(size_t)radical_coeff[j]*word_count;int wanted=(tag>>j)&1;
        for(int w=0;w<word_count;w++)out[w]&=wanted?~mask[w]:mask[w];
    }
    for(int j=0;j<span_size;j++)out[sums[depth][j]/64]&=~(UINT64_C(1)<<(sums[depth][j]%64));
    if(site<m-1&&depth>=kernel_g){
        int shift=0;for(int j=0;j<kernel_g;j++)if(anchor_coeff[tag]>>j&1)shift^=xs[j];
        for(int i=kernel_g;i<depth;i++){
            int prev=xs[i];for(int j=0;j<kernel_g;j++)if(anchor_coeff[site_tags[i]]>>j&1)prev^=xs[j];
            int a=inverse_coeff[prev],wanted=((fixed_form[site]>>i)&1)^par(a&shift);
            const uint64_t *mask=parity_allowed+(size_t)a*word_count;
            for(int w=0;w<word_count;w++)out[w]&=wanted?~mask[w]:mask[w];
        }
    }
    if(site<fixed_prefix_len){int label=fixed_prefix[site];uint64_t kept=out[label/64]&(UINT64_C(1)<<(label%64));
        for(int w=0;w<word_count;w++)out[w]=0;out[label/64]=kept;}
    uint64_t any=0;for(int w=0;w<word_count;w++)any|=out[w];
    return any!=0;
}
static int has(const uint64_t *b,int v){return (int)((b[v>>6]>>(v&63))&1);}
static int count_bits(const uint64_t *b){int c=0;for(int w=0;w<word_count;w++)c+=__builtin_popcountll(b[w]);return c;}
static int sum_feasible(int depth,int span_size){
    int s=sums[depth][span_size-1],free=m-depth,o[3]={0,1,2};
    for(int k=0;k<free;k++)if(!free_domain(depth,span_size,depth+k,ldom[k]))return 0;
    if(free==2){
        for(int w=0;w<word_count;w++){uint64_t b=ldom[0][w];
            while(b){int v=64*w+__builtin_ctzll(b);b&=b-1;if(has(ldom[1],v^s))return 1;}}
        return 0;
    }
    int c[3];for(int k=0;k<3;k++)c[k]=count_bits(ldom[k]);
    for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(c[o[j]]<c[o[i]]){int t=o[i];o[i]=o[j];o[j]=t;}
    const uint64_t *A=ldom[o[0]],*B=ldom[o[1]],*C=ldom[o[2]];
    for(int w=0;w<word_count;w++){uint64_t a=A[w];
        while(a){int v=64*w+__builtin_ctzll(a);a&=a-1;
            for(int w2=0;w2<word_count;w2++){uint64_t b=B[w2];
                while(b){int u=64*w2+__builtin_ctzll(b);b&=b-1;if(has(C,v^u^s))return 1;}}}}
    return 0;
}
'''
rep("static void recurse(int depth, int span_size, int lo) {",helper+"static void recurse(int depth, int span_size, int lo) {")
rep("    ++nodes;\n","    ++nodes; ++depth_nodes[depth];\n")
rep("    uint64_t allowed[MAX_WORDS];\n","    if(lookahead_free&&m-depth>=2&&m-depth<=lookahead_free){\n        ++look_tested;if(!sum_feasible(depth,span_size)){++look_pruned;return;}\n    }\n    uint64_t allowed[MAX_WORDS];\n")
rep('''        } else if (!strcmp(argv[i], "--fixed-form")) {''','''        } else if (!strcmp(argv[i], "--lookahead") && i+1<argc) {
            lookahead_free=atoi(argv[++i]);
            if(lookahead_free!=2&&lookahead_free!=3)return error("lookahead must be 2 or 3");
        } else if (!strcmp(argv[i], "--fixed-form")) {''')
rep('''    if (avoid_radical && !use_gamma) return error("avoid-radical requires --gamma");''','''    if (lookahead_free && !use_kernel_data) return error("lookahead requires --kernel-data");
    if (avoid_radical && !use_gamma) return error("avoid-radical requires --gamma");''')
rep('''    else fputs("null", stderr);
    fputs("}\\n", stderr);''','''    else fputs("null", stderr);
    fprintf(stderr, ",\\"lookahead_free\\":%d,\\"look_tested\\":%" PRIu64 ",\\"look_pruned\\":%" PRIu64 ",\\"depth_nodes\\":[", lookahead_free, look_tested, look_pruned);
    for(int i=0;i<=m;i++)fprintf(stderr,"%s%" PRIu64,i?",":"",depth_nodes[i]);
    fputs("]}\\n", stderr);''')
(HERE/'z_enum_lookahead.c').write_text(s)
