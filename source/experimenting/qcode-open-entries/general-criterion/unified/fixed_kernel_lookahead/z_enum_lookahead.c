/* fixed_kernel_lookahead: copy of frozen fixed_kernel/z_enum_fixed.c plus an
 * optional sum-feasibility lookahead (--lookahead F, F in {2,3}).
 * Exact streaming canonical Z enumeration; reference: z_enumerator.py.
 * stdin: m N d, then N distances. stdout: flushed LEAF x0 ... x(m-2).
 * stderr: final JSON statistics; interrupted enumeration is always LIMIT.
 */
#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_M 10
#define MAX_N 8192
#define MAX_WORDS (MAX_N / 64)
#define MAX_SPAN (1 << (MAX_M - 2))

static int m, N, distance, word_count, max_t;
static int *distances;
static uint64_t *translated;
static int sums[MAX_M][MAX_SPAN];
static unsigned char sizes[MAX_M][MAX_SPAN];
static int xs[MAX_M];
static uint64_t nodes, leaf_count, max_nodes;
static int has_max_nodes, has_seconds, limited, output_error;
static double seconds_limit, start_time;
static const char *limit_reason;
/* Gamma pruning: ann[depth] = basis of {c: <x_i,c>=0, i<depth}; restricted form
 * rank must be <= 2*(m-1-depth) (final annihilator, dim b-(m-1), is totally
 * isotropic). At the last image the condition is applied exactly via the radical. */
static int use_gamma, bdim;
static int avoid_radical, global_rad[16], global_radn;
static int use_first_orbits;
static uint64_t first_mask[MAX_WORDS];
static int gam[16];
static int ann[MAX_M][16], ann_n[MAX_M];
static uint64_t gamma_pruned, xdom_pruned;
static uint64_t *parity_allowed;
/* Unary X-domain forward check (option --xdom): for X sector X_i X_{m-1}
 * (i<m-1) and assigned Z subset A, the image y_i + zsyn(A) must have coset
 * distance >= d - (2 + min(t, m-2-t)), t = |A minus {i}|. dom[depth][i] is the
 * bitset of y_i values meeting every constraint over assigned subsets. */
static int use_xdom;
static uint64_t dom[MAX_M + 1][MAX_M][MAX_N / 64];
static int use_fixed_form, use_kernel_data, fixed_form[16], inverse_coeff[MAX_N];
static int kernel_g,kernel_t,radical_coeff[3],site_tags[MAX_M],anchor_coeff[8];
static int fixed_prefix[MAX_M],fixed_prefix_len;
static int prefix_count, prefix_actions[40320][MAX_M];
static int witness_count,witness_index[MAX_M*(1<<MAX_M)],witness_mask[MAX_M*(1<<MAX_M)];
static const char *prefix_path;
static uint64_t prefix_pruned;
static uint64_t depth_nodes[MAX_M+1];
static int lookahead_free;
static uint64_t look_tested,look_pruned;
static uint64_t ldom[3][MAX_WORDS];
static int prefix_canonical(int depth){
    int assigned=(1<<depth)-1;
    /* Independent X images make distinct coefficient masks unequal. Thus
     * only the first nonidentity row of each action can decide its lex test.
     * Deduplicating these witnesses preserves every group inequality. */
    for(int k=0;k<witness_count;k++){
        int i=witness_index[k],c=witness_mask[k];
        if(i>=depth || (c&~assigned))continue;
        int value=0;for(int j=0;j<depth;j++)if(c>>j&1)value^=xs[j];
        if(value<xs[i])return 0;
    }
    return 1;
}
static int par(int v) { return __builtin_parity((unsigned)v); }
static int gapply(int c) { int o = 0; for (int i = 0; i < bdim; ++i) o |= par(gam[i] & c) << i; return o; }
static int rank_rows(int *rows, int n) {
    int r = 0;
    for (int bit = 15; bit >= 0; --bit) {
        int p = -1;
        for (int i = r; i < n; ++i) if (rows[i] >> bit & 1) { p = i; break; }
        if (p < 0) continue;
        int t = rows[p]; rows[p] = rows[r]; rows[r] = t;
        for (int i = 0; i < n; ++i) if (i != r && (rows[i] >> bit & 1)) rows[i] ^= rows[r];
        ++r;
    }
    return r;
}
/* restricted Gram rows in ann coordinates */
static int form_rank(int depth, int *G) {
    int n = ann_n[depth], g[16], rows[16];
    for (int j = 0; j < n; ++j) g[j] = gapply(ann[depth][j]);
    for (int i = 0; i < n; ++i) { G[i] = 0; for (int j = 0; j < n; ++j) G[i] |= par(ann[depth][i] & g[j]) << j; rows[i] = G[i]; }
    return rank_rows(rows, n);
}
static void restrict_ann(int depth, int x) {
    int n = ann_n[depth], piv = -1, k = 0;
    for (int i = 0; i < n; ++i) if (par(ann[depth][i] & x)) { piv = i; break; }
    for (int i = 0; i < n; ++i) {
        if (i == piv) continue;
        int c = ann[depth][i];
        if (piv >= 0 && par(c & x)) c ^= ann[depth][piv];
        ann[depth + 1][k++] = c;
    }
    ann_n[depth + 1] = k;
}

/* Radical prefix words and their translates by Z^m both have weight
 * >=d. Hence their binary X supports have even weight >=2*d-m. */
#define MAX_PATTERNS 2048
#define PATTERN_WORDS (MAX_PATTERNS/64)
static int use_projection, projection_min, projection_tags, pattern_count, pattern_words;
static unsigned char patterns[MAX_PATTERNS][8], composition[8];
static uint64_t capacity[8][MAX_M+1][PATTERN_WORDS];
static uint64_t viable[MAX_M+1][PATTERN_WORDS];
static unsigned char tag_counts[MAX_M+1][8];
static unsigned char syndrome_tag[MAX_N];
static uint64_t projection_pruned;
static void pattern_rec(int v,int left) {
    if(v==projection_tags-1){
        composition[v]=(unsigned char)left;
        for(int u=1;u<projection_tags;++u){int w=0;for(int t=0;t<projection_tags;++t)if(par(u&t))w+=composition[t];
            if(w<projection_min || w%2)return;}
        if(pattern_count>=MAX_PATTERNS){fprintf(stderr,"projection pattern bound exceeded\n");exit(2);}
        memcpy(patterns[pattern_count++],composition,8);return;
    }
    for(int c=0;c<=left;++c){composition[v]=(unsigned char)c;pattern_rec(v+1,left-c);}
}
static void initialize_projection(void) {
    projection_min=2*distance-m;projection_tags=1<<global_radn;
    pattern_rec(0,m);pattern_words=(pattern_count+63)/64;
    for(int s=0;s<N;++s)for(int i=0;i<global_radn;++i)syndrome_tag[s]|=par(s&global_rad[i])<<i;
    for(int p=0;p<pattern_count;++p){
        viable[0][p/64]|=UINT64_C(1)<<(p%64);
        for(int t=0;t<projection_tags;++t)for(int c=0;c<=patterns[p][t];++c)
            capacity[t][c][p/64]|=UINT64_C(1)<<(p%64);
    }
}
static int projection_step(int depth,int x) {
    int t=syndrome_tag[x],c=tag_counts[depth][t]+1;uint64_t any=0;
    memcpy(tag_counts[depth+1],tag_counts[depth],8);tag_counts[depth+1][t]=(unsigned char)c;
    for(int w=0;w<pattern_words;++w){viable[depth+1][w]=viable[depth][w]&capacity[t][c][w];any|=viable[depth+1][w];}
    return any!=0;
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.0;
}

static const uint64_t *translated_row(int t, int syndrome) {
    return translated + (((size_t)(t - 1) * N + syndrome) * word_count);
}

static int xdom_step(int depth, int span_size, int x) {
    for (int i = 0; i < m - 1; ++i) {
        uint64_t *nd = dom[depth + 1][i];
        memcpy(nd, dom[depth][i], (size_t)word_count * sizeof(uint64_t));
        for (int j = 0; j < span_size; ++j) {
            int mask = j | (1 << depth);
            int t = __builtin_popcount((unsigned)(mask & ~(1 << i)));
            int cost = 2 + (t < m - 2 - t ? t : m - 2 - t);
            if (cost >= distance) continue;
            const uint64_t *row = translated_row(cost, sums[depth][j] ^ x);
            uint64_t any = 0;
            for (int w = 0; w < word_count; ++w) { nd[w] &= row[w]; any |= nd[w]; }
            if (!any) return 0;
        }
    }
    return 1;
}

static int radical_evaluation_rank(int depth) {
    int rows[16];
    for (int i=0;i<global_radn;++i) {
        rows[i]=0;
        for (int j=0;j<depth;++j) rows[i] |= par(global_rad[i]&xs[j])<<j;
    }
    return rank_rows(rows,global_radn);
}
static void initialize_global_radical(void) {
    int piv[16]={0};
    for (int v=1;v<N;++v) if (!gapply(v)) {
        int x=v;
        for (int bit=bdim-1;bit>=0;--bit) if (x>>bit&1) {
            if (piv[bit]) x ^= piv[bit];
            else { piv[bit]=x;global_rad[global_radn++]=v;break; }
        }
    }
}


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
static void recurse(int depth, int span_size, int lo) {
    if (limited || output_error) return;
    ++nodes; ++depth_nodes[depth];
    if (has_max_nodes && nodes > max_nodes) {
        limited = 1; limit_reason = "node_limit"; return;
    }
    if (has_seconds && nodes % 1024 == 0 && now() - start_time >= seconds_limit) {
        limited = 1; limit_reason = "time_limit"; return;
    }
    int radn = 0, rad[16];
    if (use_gamma && !use_fixed_form && !use_kernel_data) {
        int G[16];
        int rho = form_rank(depth, G);
        int bound=2*(m-1-depth);
        if (avoid_radical) bound=2*(m-1-depth+radical_evaluation_rank(depth)-global_radn);
        if (rho > bound) { ++gamma_pruned; return; }
        if (rho == 2 * (m - 1 - depth)) {
            /* radical basis: nullspace of symmetric G in ann coordinates */
            int n = ann_n[depth], rows[16], pc[16], r = 0;
            for (int i = 0; i < n; ++i) rows[i] = G[i];
            for (int bit = 0; bit < n && r < n; ++bit) {
                int p = -1;
                for (int i = r; i < n; ++i) if (rows[i] >> bit & 1) { p = i; break; }
                if (p < 0) continue;
                int t = rows[p]; rows[p] = rows[r]; rows[r] = t;
                for (int i = 0; i < n; ++i) if (i != r && (rows[i] >> bit & 1)) rows[i] ^= rows[r];
                pc[r++] = bit;
            }
            for (int fcol = 0; fcol < n; ++fcol) {
                int isp = 0; for (int i = 0; i < r; ++i) if (pc[i] == fcol) isp = 1;
                if (isp) continue;
                int u = 1 << fcol;
                for (int i = 0; i < r; ++i) if (rows[i] >> fcol & 1) u |= 1 << pc[i];
                int v = 0; for (int i = 0; i < n; ++i) if (u >> i & 1) v ^= ann[depth][i];
                rad[radn++] = v;
            }
        }
    }
    if(lookahead_free&&m-depth>=2&&m-depth<=lookahead_free){
        ++look_tested;if(!sum_feasible(depth,span_size)){++look_pruned;return;}
    }
    uint64_t allowed[MAX_WORDS];
    for (int w = 0; w < word_count; ++w) allowed[w] = UINT64_MAX;
    if (N < 64) allowed[0] = (UINT64_C(1) << N) - 1;
    /* Remove labels <= lo; lo=-1 deliberately removes none. */
    for (int w = 0; w < word_count && 64 * w <= lo; ++w) {
        int last = lo - 64 * w;
        if (last >= 63) allowed[w] = 0;
        else allowed[w] &= UINT64_MAX << (last + 1);
    }
    for (int j = 0; j < span_size; ++j) {
        int a = (int)sizes[depth][j] + 1;
        int t = a < m - a ? a : m - a;
        if (distance - t <= 0) continue;
        const uint64_t *row = translated_row(t, sums[depth][j]);
        uint64_t any = 0;
        for (int w = 0; w < word_count; ++w) {
            allowed[w] &= row[w];
            any |= allowed[w];
        }
        if (!any) return;
    }
    /* Independent pruning: use the radical of the full restricted Gram
     * matrix, not the small rank-budget identity. At saturation each next
     * hyperplane must drop the restricted rank by two, hence contain the
     * radical. Independence is guaranteed by d>floor(m/2) and delta(0)=0.
     */
    if (use_first_orbits && depth==0)
        for (int w=0;w<word_count;++w) allowed[w] &= first_mask[w];
    for (int j = 0; j < radn; ++j) {
        const uint64_t *row = parity_allowed + (size_t)rad[j] * word_count;
        for (int w = 0; w < word_count; ++w) allowed[w] &= row[w];
    }
    if(use_kernel_data){
        int tag=site_tags[depth];
        for(int j=0;j<kernel_g;j++){
            const uint64_t *mask=parity_allowed+(size_t)radical_coeff[j]*word_count;
            int wanted=(tag>>j)&1;
            for(int w=0;w<word_count;w++)allowed[w]&=wanted?~mask[w]:mask[w];
        }
        if(depth>=kernel_g){
            int shift=0;for(int j=0;j<kernel_g;j++)if(anchor_coeff[tag]>>j&1)shift^=xs[j];
            for(int i=kernel_g;i<depth;i++){
                int prev=xs[i];for(int j=0;j<kernel_g;j++)if(anchor_coeff[site_tags[i]]>>j&1)prev^=xs[j];
                int a=inverse_coeff[prev],wanted=((fixed_form[depth]>>i)&1)^par(a&shift);
                const uint64_t *mask=parity_allowed+(size_t)a*word_count;
                for(int w=0;w<word_count;w++)allowed[w]&=wanted?~mask[w]:mask[w];
            }
        }
        // Exact fixed H has no nonzero pure-Z class: enforce independent X
        // even for lower-distance physical controls.
        for(int j=0;j<span_size;j++)allowed[sums[depth][j]/64]&=~(UINT64_C(1)<<(sums[depth][j]%64));
    }
    if(depth<fixed_prefix_len){
        int label=fixed_prefix[depth];uint64_t kept=allowed[label/64]&(UINT64_C(1)<<(label%64));
        for(int w=0;w<word_count;w++)allowed[w]=0;allowed[label/64]=kept;
    }
    if(use_fixed_form)for(int i=0;i<depth;i++){
        int a=inverse_coeff[xs[i]],wanted=(fixed_form[depth]>>i)&1;
        const uint64_t *mask=parity_allowed+(size_t)a*word_count;
        for(int w=0;w<word_count;w++)allowed[w]&=wanted?~mask[w]:mask[w];
    }
    for (int w = 0; w < word_count; ++w) {
        while (allowed[w]) {
            int bit = __builtin_ctzll(allowed[w]);
            int x = 64 * w + bit;
            allowed[w] &= allowed[w] - 1;
            xs[depth] = x;
            if((use_fixed_form||use_kernel_data) && !prefix_canonical(depth+1)){++prefix_pruned;continue;}
            if(use_projection && !projection_step(depth,x)){++projection_pruned;continue;}
            if (depth == m - 2) {
                int last = x ^ sums[depth][span_size - 1];
                if (!use_fixed_form && !use_kernel_data && last <= x) continue;
                if(use_projection && !projection_step(depth+1,last)){++projection_pruned;continue;}
                if (use_xdom && !xdom_step(depth, span_size, x)) { ++xdom_pruned; continue; }
                if (avoid_radical && radical_evaluation_rank(depth+1)!=global_radn) {
                    ++gamma_pruned; continue;
                }
                ++leaf_count;
                if (fputs("LEAF", stdout) == EOF) { output_error = 1; return; }
                for (int i = 0; i < m - 1; ++i) {
                    if (fprintf(stdout, " %d", xs[i]) < 0) { output_error = 1; return; }
                }
                if (fputc('\n', stdout) == EOF || fflush(stdout) == EOF) {
                    output_error = 1; return;
                }
            } else {
                for (int j = 0; j < span_size; ++j) {
                    sums[depth + 1][j] = sums[depth][j];
                    sizes[depth + 1][j] = sizes[depth][j];
                    sums[depth + 1][span_size + j] = sums[depth][j] ^ x;
                    sizes[depth + 1][span_size + j] = sizes[depth][j] + 1;
                }
                if (use_xdom && !xdom_step(depth, span_size, x)) { ++xdom_pruned; continue; }
                if (use_gamma && !use_fixed_form && !use_kernel_data) restrict_ann(depth, x);
                recurse(depth + 1, span_size * 2, (use_fixed_form||use_kernel_data)?-1:x);
                if (limited || output_error) return;
            }
        }
    }
}

static int error(const char *message) {
    fprintf(stderr, "{\"status\":\"ERROR\",\"error\":\"%s\"}\n", message);
    return 2;
}

int main(int argc, char **argv) {
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--max-nodes") && i + 1 < argc) {
            char *end; errno = 0;
            if (argv[i + 1][0] == '-') return error("invalid node limit");
            max_nodes = strtoull(argv[++i], &end, 10);
            if (errno || !*argv[i] || *end) return error("invalid node limit");
            has_max_nodes = 1;
        } else if (!strcmp(argv[i], "--seconds") && i + 1 < argc) {
            char *end; errno = 0;
            seconds_limit = strtod(argv[++i], &end);
            if (errno || !*argv[i] || *end || !(seconds_limit >= 0.0))
                return error("invalid time limit");
            has_seconds = 1;
        } else if (!strcmp(argv[i], "--kernel-data")) {
            use_kernel_data=1;
        } else if (!strcmp(argv[i], "--prefix") && i+1<argc) {
            const char *v=argv[++i];
            while(*v){char *end;long x=strtol(v,&end,10);if(end==v||x<0||x>=MAX_N||fixed_prefix_len>=MAX_M)return 2;fixed_prefix[fixed_prefix_len++]=(int)x;if(!*end)break;if(*end!=',')return 2;v=end+1;}
        } else if (!strcmp(argv[i], "--lookahead") && i+1<argc) {
            lookahead_free=atoi(argv[++i]);
            if(lookahead_free!=2&&lookahead_free!=3)return error("lookahead must be 2 or 3");
        } else if (!strcmp(argv[i], "--fixed-form")) {
            use_fixed_form=1;
        } else if (!strcmp(argv[i], "--prefix-group") && i+1<argc) {
            prefix_path=argv[++i];
        } else if (!strcmp(argv[i], "--projection")) {
            use_projection=1;
        } else if (!strcmp(argv[i], "--first-orbits")) {
            use_first_orbits=1;
        } else if (!strcmp(argv[i], "--avoid-radical")) {
            avoid_radical=1;
        } else if (!strcmp(argv[i], "--xdom")) {
            use_xdom = 1;
        } else if (!strcmp(argv[i], "--gamma")) {
            use_gamma = 1;
        } else {
            return error("usage: z_enumerator [--seconds S] [--max-nodes K]");
        }
    }
    if (scanf("%d %d %d", &m, &N, &distance) != 3)
        return error("expected m N d and N distances on stdin");
    if (m < 3 || m > MAX_M || N < 1 || N > MAX_N || (N & (N - 1)) || distance < 3)
        return error("require m=3..10, power-of-two N<=8192, d>=3");
    if (use_gamma && !use_kernel_data && distance <= m / 2)
        return error("independent backend requires d>floor(m/2)");
    distances = malloc((size_t)N * sizeof(*distances));
    if (!distances) return error("distance allocation failed");
    int deep_count = 0;
    for (int i = 0; i < N; ++i) {
        if (scanf("%d", distances + i) != 1 || distances[i] < 0)
            return error("expected N nonnegative integer distances");
        if (distances[i] >= distance - 1) ++deep_count;
    }
    if (use_gamma && distances[0] != 0)
        return error("independent backend requires delta(0)=0");
    /* Match both reference global early-return checks, before rec(). */
    if (deep_count < m) {
        fprintf(stderr, "{\"status\":\"COMPLETE\",\"nodes\":0,\"leaves\":0,\"seconds\":0.0}\n");
        free(distances); return 0;
    }
    if (use_gamma) {
        bdim = __builtin_ctz((unsigned)N);
        for (int i = 0; i < bdim; ++i) if (scanf("%d", gam + i) != 1) return error("expected b Gamma rows");
        for (int i = 0; i < bdim; ++i) {
            if (gam[i]<0 || gam[i]>=N || (gam[i]>>i&1)) return error("Gamma must be alternating");
            for (int j=0;j<bdim;++j) if ((gam[i]>>j&1)!=(gam[j]>>i&1)) return error("Gamma must be symmetric");
        }
        for (int i = 0; i < bdim; ++i) ann[0][i] = 1 << i;
        ann_n[0] = bdim;
    }
    if(use_fixed_form){
        if(!use_gamma||bdim!=2*(m-2))return error("fixed form needs nonsingular quotient of dimension 2m-4");
        int rows[16];for(int i=0;i<bdim;i++)rows[i]=gam[i];if(rank_rows(rows,bdim)!=bdim)return error("fixed form needs nonsingular Gamma");
        for(int i=0;i<m-1;i++)if(scanf("%d",fixed_form+i)!=1||fixed_form[i]<0||fixed_form[i]>=(1<<(m-1)))return error("bad prefix Gram");
        for(int i=0;i<m-1;i++){if(fixed_form[i]>>i&1)return error("prefix Gram must alternate");for(int j=0;j<m-1;j++)if((fixed_form[i]>>j&1)!=(fixed_form[j]>>i&1))return error("prefix Gram must be symmetric");rows[i]=fixed_form[i];}
        if(rank_rows(rows,m-1)!=2)return error("prefix Gram must have rank two");
        for(int c=0;c<N;c++)inverse_coeff[gapply(c)]=c;
    }
    if(use_kernel_data){
        if(!use_gamma||use_fixed_form||scanf("%d %d",&kernel_g,&kernel_t)!=2||kernel_g<0||kernel_g>3||(kernel_t!=0&&kernel_t!=2)||kernel_g+kernel_t<1||kernel_g+kernel_t>3||bdim!=2*(m-1)-kernel_g-kernel_t)return error("bad fixed kernel dimensions");
        int rows[16];for(int i=0;i<bdim;i++)rows[i]=gam[i];if(rank_rows(rows,bdim)!=bdim-kernel_g)return error("wrong fixed-kernel radical dimension");
        for(int i=0;i<kernel_g;i++){if(scanf("%d",radical_coeff+i)!=1||radical_coeff[i]<=0||radical_coeff[i]>=N||gapply(radical_coeff[i]))return error("bad radical coefficient");rows[i]=radical_coeff[i];}
        if(rank_rows(rows,kernel_g)!=kernel_g)return error("dependent radical basis");
        int total=0;for(int i=0;i<m;i++){if(scanf("%d",site_tags+i)!=1||site_tags[i]<0||site_tags[i]>=(1<<kernel_g))return error("bad kernel tags");total^=site_tags[i];}
        if(total)return error("kernel tags must sum to zero");
        for(int i=0;i<kernel_g;i++)rows[i]=site_tags[i];if(rank_rows(rows,kernel_g)!=kernel_g)return error("first sites must give independent tags");
        for(int c=0;c<(1<<kernel_g);c++){int value=0;for(int j=0;j<kernel_g;j++)if(c>>j&1)value^=site_tags[j];anchor_coeff[value]=c;}
        for(int i=0;i<m-1;i++){if(scanf("%d",fixed_form+i)!=1||fixed_form[i]<0||fixed_form[i]>=(1<<(m-1)))return error("bad null-tag Gram");rows[i]=fixed_form[i];}
        for(int i=0;i<m-1;i++){if(fixed_form[i]>>i&1)return error("null-tag Gram must alternate");if(i<kernel_g&&fixed_form[i])return error("anchor Gram must vanish");for(int j=0;j<m-1;j++)if((fixed_form[i]>>j&1)!=(fixed_form[j]>>i&1))return error("null-tag Gram must be symmetric");}
        if(rank_rows(rows,m-1)!=kernel_t)return error("wrong null-tag Gram rank");
        for(int c=0;c<N;c++)inverse_coeff[gapply(c)]=c;
    }
    if(fixed_prefix_len>m-1)return error("too many fixed prefix images");
    for(int i=0;i<fixed_prefix_len;i++)if(fixed_prefix[i]>=N)return error("fixed prefix image out of range");
    if(prefix_path){
        FILE *f=fopen(prefix_path,"r");int width;
        if((!use_fixed_form&&!use_kernel_data)||!f||fscanf(f,"%d %d",&prefix_count,&width)!=2||prefix_count<1||prefix_count>40320||width!=m-1)return error("bad prefix group");
        for(int g=0;g<prefix_count;g++)for(int j=0;j<m-1;j++)if(fscanf(f,"%d",&prefix_actions[g][j])!=1||prefix_actions[g][j]<=0||prefix_actions[g][j]>=(1<<(m-1)))return error("bad prefix action");
        for(int g=0;g<prefix_count;g++){
            int rows[16];for(int j=0;j<m-1;j++)rows[j]=prefix_actions[g][j];
            if(rank_rows(rows,m-1)!=m-1)return error("singular prefix action");
            for(int i=0;i<m-1;i++)if(prefix_actions[g][i]!=(1<<i)){
                int c=prefix_actions[g][i],seen=0;
                for(int k=0;k<witness_count;k++)if(witness_index[k]==i&&witness_mask[k]==c){seen=1;break;}
                if(!seen){witness_index[witness_count]=i;witness_mask[witness_count++]=c;}
                break;
            }
        }
        fclose(f);
    }
    if (use_first_orbits) {
        int count;
        if (scanf("%d",&count)!=1 || count<0 || count>N) return error("expected valid orbit representative count");
        for (int i=0;i<count;++i) {
            int v;
            if (scanf("%d",&v)!=1 || v<0 || v>=N || distances[v]<distance-1) return error("invalid first orbit representative");
            first_mask[v>>6] |= UINT64_C(1)<<(v&63);
        }
    }
    word_count = (N + 63) / 64;
    max_t = m / 2;
    if (use_xdom) max_t = distance - 1;
    if (max_t >= distance) max_t = distance - 1;
    translated = calloc((size_t)max_t * N * word_count, sizeof(*translated));
    if (!translated) return error("translated bitset allocation failed");
    for (int t = 1; t <= max_t; ++t) {
        int need = distance - t;
        for (int s = 0; s < N; ++s) {
            uint64_t *row = translated + (((size_t)(t - 1) * N + s) * word_count);
            for (int x = 0; x < N; ++x) {
                if (distances[x ^ s] >= need) row[x >> 6] |= UINT64_C(1) << (x & 63);
            }
        }
    }
    if (use_xdom) {
        for (int i = 0; i < m - 1; ++i) {
            const uint64_t *row = translated_row(2, 0);
            memcpy(dom[0][i], row, (size_t)word_count * sizeof(uint64_t));
        }
    }
    if (lookahead_free && !use_kernel_data) return error("lookahead requires --kernel-data");
    if (avoid_radical && !use_gamma) return error("avoid-radical requires --gamma");
    if (avoid_radical) initialize_global_radical();
    if (use_gamma) {
        parity_allowed=calloc((size_t)N*word_count,sizeof(uint64_t));
        if (!parity_allowed) return error("parity allocation failed");
        for (int row=0;row<N;++row) for (int x=0;x<N;++x) if (!par(row&x))
            parity_allowed[(size_t)row*word_count+(x>>6)] |= UINT64_C(1)<<(x&63);
    }
    if(use_projection){
        if(!use_gamma||!avoid_radical||global_radn>3)return error("projection requires gamma, exact radical and radical bits<=3");
        initialize_projection();
    }
    start_time = now();
    sums[0][0] = sizes[0][0] = 0;
    recurse(0, 1, -1);
    double elapsed = now() - start_time;
    if (output_error) {
        free(translated); free(distances); return error("stdout write failed");
    }
    fprintf(stderr, "{\"status\":\"%s\",\"nodes\":%" PRIu64
            ",\"leaves\":%" PRIu64 ",\"gamma_pruned\":%" PRIu64 ",\"xdom_pruned\":%" PRIu64 ",\"prefix_pruned\":%" PRIu64 ",\"projection_pruned\":%" PRIu64 ",\"seconds\":%.9f,\"reason\":",
            limited ? "LIMIT" : "COMPLETE", nodes, leaf_count, gamma_pruned, xdom_pruned, prefix_pruned, projection_pruned, elapsed);
    if (limited) fprintf(stderr, "\"%s\"", limit_reason);
    else fputs("null", stderr);
    fprintf(stderr, ",\"lookahead_free\":%d,\"look_tested\":%" PRIu64 ",\"look_pruned\":%" PRIu64 ",\"depth_nodes\":[", lookahead_free, look_tested, look_pruned);
    for(int i=0;i<=m;i++)fprintf(stderr,"%s%" PRIu64,i?",":"",depth_nodes[i]);
    fputs("]}\n", stderr);
    free(parity_allowed); free(translated); free(distances);
    return 0;
}
