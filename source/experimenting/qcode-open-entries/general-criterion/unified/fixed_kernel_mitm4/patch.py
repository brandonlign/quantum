"""Regenerate z_enum_mitm4.c from the frozen fixed_kernel/z_enum_fixed.c."""
from pathlib import Path
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'fixed_kernel/z_enum_fixed.c').read_text()
def rep(a,b):
    global s
    assert s.count(a)==1,a[:60]; s=s.replace(a,b)
rep("/* Exact streaming canonical Z enumeration; reference: z_enumerator.py.",
"/* fixed_kernel_mitm4: frozen fixed_kernel/z_enum_fixed.c plus an optional\n * four-free-site meet-in-middle feasibility test at depth 4 (m=8, --kernel-data).\n * Exact streaming canonical Z enumeration; reference: z_enumerator.py.")
rep("static uint64_t prefix_pruned;","""static uint64_t prefix_pruned;
static uint64_t depth_nodes[MAX_M+1];
static int mitm_cap; static uint64_t mitm_check_bound, mitm_work_limit=UINT64_C(1)<<22;
static uint64_t mitm_tested,mitm_pruned,mitm_skipped,mitm_overflow,mitm_checked,mitm_matches;""")
helper=r'''
/* ---- four-free MITM (README/REPORT in fixed_kernel_mitm4/) ----
 * Depth 4, m=8: assigned x0..x3 (all anchors, g<=3), W=span, s=x0^..^x3.
 * Free a=x4,b=x5,c=x6,d=x7 (x7 implied) satisfy a^b^c^d=s. With
 * T_f=x_f^shift(tag_f) and B(u,v)=par(inverse_coeff[u]&v), every native leaf
 * has B(T_i,T_j)=F[i][j] for all i,j<8, F[7]=XOR_k<7 F[k]; tags hold for all
 * eight sites; x_f notin W; every subset-distance row holds. */
#define MITM_MAX_CAP 65536
static uint64_t mdom[4][MAX_WORDS], mbase[MAX_WORDS], mH2[MAX_WORDS], minW[MAX_WORDS];
static int mlist[4][MAX_N], mT[4][MAX_N], mI[4][MAX_N], mcount[4];
static int bucket_head[MAX_N]; static unsigned bucket_epoch[MAX_N], mepoch;
static uint16_t pu[MITM_MAX_CAP], pv[MITM_MAX_CAP]; static int pnext[MITM_MAX_CAP];
static int ptu[MITM_MAX_CAP], piu[MITM_MAX_CAP];
static int F7, Fbit[8][8];
static int mbit(const uint64_t *b,int v){return (int)((b[v>>6]>>(v&63))&1);}
static int site_shift(int site){int sh=0,tag=site_tags[site];for(int j=0;j<kernel_g;j++)if(anchor_coeff[tag]>>j&1)sh^=xs[j];return sh;}
static int bpair(int inv_u,int tv){return par(inv_u&tv);}
static int mitm_quad_ok(int a,int b,int c,int d,int Ia,int Ib,int Tc,int Td){
    int ac=a^c,ad=a^d;
    if(mbit(minW,ac)||mbit(minW,ad)||!mbit(mH2,ac)||!mbit(mH2,ad))return 0;
    if(bpair(Ia,Tc)!=Fbit[4][6]||bpair(Ia,Td)!=Fbit[4][7]||bpair(Ib,Tc)!=Fbit[5][6]||bpair(Ib,Td)!=Fbit[5][7])return 0;
    return 1;
}
/* returns 1 feasible/skipped, 0 = provably no completion */
static int mitm_core(int s,int *skipped){
    int n4=mcount[0],n5=mcount[1],n6=mcount[2],n7=mcount[3];
    uint64_t p1=(uint64_t)n4*n5,p2=(uint64_t)n6*n7;
    *skipped=0;
    if(p1>(uint64_t)mitm_cap||p2>(uint64_t)mitm_cap){*skipped=1;return 1;}
    /* hash side 0: (x4,x5) key a^b ; side 1: (x6,x7) key c^d^s. both equal h. */
    int hs=(p1<=p2)?0:1,u0=hs?2:0,ps=hs?0:2,np=0;
    int Fu=hs?Fbit[6][7]:Fbit[4][5],Fp=hs?Fbit[4][5]:Fbit[6][7],ksh=hs?s:0,kpsh=hs?0:s;
    if(++mepoch==0){memset(bucket_epoch,0,sizeof bucket_epoch);mepoch=1;}
    for(int i=0;i<mcount[u0];i++)for(int j=0;j<mcount[u0+1];j++){
        int u=mlist[u0][i],v=mlist[u0+1][j],h=u^v;
        if(mbit(minW,h)||!mbit(mH2,h)||bpair(mI[u0][i],mT[u0+1][j])!=Fu)continue;
        int key=h^ksh;
        if(bucket_epoch[key]!=mepoch){bucket_epoch[key]=mepoch;bucket_head[key]=-1;}
        pu[np]=(uint16_t)u;pv[np]=(uint16_t)v;ptu[np]=mT[u0][i];piu[np]=(i<<16)|j;pnext[np]=bucket_head[key];bucket_head[key]=np;np++;
    }
    uint64_t work=0;
    for(int i=0;i<mcount[ps];i++)for(int j=0;j<mcount[ps+1];j++){
        int u=mlist[ps][i],v=mlist[ps+1][j],h=u^v;
        if(mbit(minW,h)||!mbit(mH2,h)||bpair(mI[ps][i],mT[ps+1][j])!=Fp)continue;
        int key=h^kpsh;
        if(bucket_epoch[key]!=mepoch)continue;
        for(int k=bucket_head[key];k>=0;k=pnext[k]){
            if(++work>mitm_work_limit){++mitm_overflow;mitm_matches+=work;return 1;}
            int hi=piu[k]>>16,hj=piu[k]&0xffff,ok;
            if(hs==0) ok=mitm_quad_ok(pu[k],pv[k],u,v,mI[0][hi],mI[1][hj],mT[2][i],mT[3][j]);
            else      ok=mitm_quad_ok(u,v,pu[k],pv[k],mI[0][i],mI[1][j],mT[2][hi],mT[3][hj]);
            if(ok){mitm_matches+=work;return 1;}
        }
    }
    mitm_matches+=work;return 0;
}
static int mitm_brute(int s){
    for(int i=0;i<mcount[0];i++)for(int j=0;j<mcount[1];j++){
        int a=mlist[0][i],b=mlist[1][j],h=a^b;
        if(mbit(minW,h)||!mbit(mH2,h)||bpair(mI[0][i],mT[1][j])!=Fbit[4][5])continue;
        for(int k=0;k<mcount[2];k++){
            int c=mlist[2][k],d=s^h^c;
            if(!mbit(mdom[3],d))continue;
            int Td=d^site_shift(7);
            if(bpair(inverse_coeff[mT[2][k]],Td)!=Fbit[6][7])continue;
            if(mitm_quad_ok(a,b,c,d,mI[0][i],mI[1][j],mT[2][k],Td))return 1;
        }
    }
    return 0;
}
static int mitm4(int span_size){
    int depth=4,s=sums[depth][span_size-1];
    /* shared base: subset-distance rows (cost symmetric in complements) and x notin W */
    for(int w=0;w<word_count;w++){mbase[w]=UINT64_MAX;mH2[w]=UINT64_MAX;minW[w]=0;}
    if(N<64){mbase[0]=mH2[0]=(UINT64_C(1)<<N)-1;}
    for(int j=0;j<span_size;j++){
        int k=(int)sizes[depth][j],t1=k+1<m-k-1?k+1:m-k-1,t2=k+2<m-k-2?k+2:m-k-2;
        if(distance-t1>0){const uint64_t *row=translated_row(t1,sums[depth][j]);for(int w=0;w<word_count;w++)mbase[w]&=row[w];}
        if(distance-t2>0){const uint64_t *row=translated_row(t2,sums[depth][j]);for(int w=0;w<word_count;w++)mH2[w]&=row[w];}
        minW[sums[depth][j]>>6]|=UINT64_C(1)<<(sums[depth][j]&63);
    }
    for(int w=0;w<word_count;w++)mbase[w]&=~minW[w];
    int Tas[8],Ias[8];
    for(int i=kernel_g;i<depth;i++){Tas[i]=xs[i]^site_shift(i);Ias[i]=inverse_coeff[Tas[i]];}
    for(int f=4;f<8;f++){
        uint64_t *out=mdom[f-4];memcpy(out,mbase,(size_t)word_count*sizeof(uint64_t));
        int tag=site_tags[f],sh=site_shift(f);
        for(int j=0;j<kernel_g;j++){const uint64_t *mask=parity_allowed+(size_t)radical_coeff[j]*word_count;int wanted=(tag>>j)&1;
            for(int w=0;w<word_count;w++)out[w]&=wanted?~mask[w]:mask[w];}
        for(int i=kernel_g;i<depth;i++){
            int wanted=Fbit[f][i]^par(Ias[i]&sh);const uint64_t *mask=parity_allowed+(size_t)Ias[i]*word_count;
            for(int w=0;w<word_count;w++)out[w]&=wanted?~mask[w]:mask[w];
        }
        if(f<fixed_prefix_len){int label=fixed_prefix[f];uint64_t kept=out[label/64]&(UINT64_C(1)<<(label%64));
            for(int w=0;w<word_count;w++)out[w]=0;out[label/64]=kept;}
        int n=0;
        for(int w=0;w<word_count;w++){uint64_t b=out[w];while(b){int v=64*w+__builtin_ctzll(b);b&=b-1;
            mlist[f-4][n]=v;mT[f-4][n]=v^sh;mI[f-4][n]=inverse_coeff[v^sh];n++;}}
        mcount[f-4]=n;
        if(!n)return 0;
    }
    int skipped,ans=mitm_core(s,&skipped);
    if(skipped){++mitm_skipped;return 1;}
    if(mitm_check_bound&&(uint64_t)mcount[0]*mcount[1]*mcount[2]<=mitm_check_bound){
        ++mitm_checked;int brute=mitm_brute(s);
        if(brute!=ans&&!(ans==1&&mitm_overflow)){fprintf(stderr,"{\"status\":\"ERROR\",\"error\":\"mitm4 brute mismatch %d %d\"}\n",ans,brute);exit(3);}
    }
    return ans;
}
'''
rep("static void recurse(int depth, int span_size, int lo) {",helper+"static void recurse(int depth, int span_size, int lo) {")
rep("    ++nodes;\n","    ++nodes; ++depth_nodes[depth];\n")
rep("    uint64_t allowed[MAX_WORDS];\n","    if(mitm_cap&&depth==4){++mitm_tested;if(!mitm4(span_size)){++mitm_pruned;return;}}\n    uint64_t allowed[MAX_WORDS];\n")
rep('''        } else if (!strcmp(argv[i], "--fixed-form")) {''','''        } else if (!strcmp(argv[i], "--mitm4") && i+1<argc) {
            mitm_cap=atoi(argv[++i]); if(mitm_cap<1||mitm_cap>MITM_MAX_CAP)return error("mitm4 cap must be 1..65536");
        } else if (!strcmp(argv[i], "--mitm4-check") && i+1<argc) {
            mitm_check_bound=strtoull(argv[++i],NULL,10);
        } else if (!strcmp(argv[i], "--mitm4-work") && i+1<argc) {
            mitm_work_limit=strtoull(argv[++i],NULL,10);
        } else if (!strcmp(argv[i], "--fixed-form")) {''')
rep('''    if (avoid_radical && !use_gamma) return error("avoid-radical requires --gamma");''','''    if (mitm_cap && (!use_kernel_data || m!=8 || kernel_g>3 || N>65536)) mitm_cap=0; /* inapplicable: native path */
    if (mitm_cap){
        for(int i=0;i<m-1;i++)for(int j=0;j<m-1;j++)Fbit[i][j]=(fixed_form[i]>>j)&1;
        F7=0;for(int k=0;k<m-1;k++)F7^=fixed_form[k];
        for(int i=0;i<m-1;i++){Fbit[7][i]=(F7>>i)&1;Fbit[i][7]=Fbit[7][i];}
        Fbit[7][7]=0;
    }
    if (avoid_radical && !use_gamma) return error("avoid-radical requires --gamma");''')
rep('''    else fputs("null", stderr);
    fputs("}\\n", stderr);''','''    else fputs("null", stderr);
    fprintf(stderr, ",\\"mitm_cap\\":%d,\\"mitm_tested\\":%" PRIu64 ",\\"mitm_pruned\\":%" PRIu64 ",\\"mitm_skipped\\":%" PRIu64 ",\\"mitm_overflow\\":%" PRIu64 ",\\"mitm_checked\\":%" PRIu64 ",\\"mitm_matches\\":%" PRIu64 ",\\"depth_nodes\\":[", mitm_cap, mitm_tested, mitm_pruned, mitm_skipped, mitm_overflow, mitm_checked, mitm_matches);
    for(int i=0;i<=m;i++)fprintf(stderr,"%s%" PRIu64,i?",":"",depth_nodes[i]);
    fputs("]}\\n", stderr);''')
(HERE/'z_enum_mitm4.c').write_text(s)
