/* Exact streaming canonical Z enumeration; reference: z_enumerator.py.
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

static void recurse(int depth, int span_size, int lo) {
    if (limited || output_error) return;
    ++nodes;
    if (has_max_nodes && nodes > max_nodes) {
        limited = 1; limit_reason = "node_limit"; return;
    }
    if (has_seconds && nodes % 1024 == 0 && now() - start_time >= seconds_limit) {
        limited = 1; limit_reason = "time_limit"; return;
    }
    int radn = 0, rad[16];
    if (use_gamma) {
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
    for (int j = 0; j < radn; ++j) {
        const uint64_t *row = parity_allowed + (size_t)rad[j] * word_count;
        for (int w = 0; w < word_count; ++w) allowed[w] &= row[w];
    }
    for (int w = 0; w < word_count; ++w) {
        while (allowed[w]) {
            int bit = __builtin_ctzll(allowed[w]);
            int x = 64 * w + bit;
            allowed[w] &= allowed[w] - 1;
            xs[depth] = x;
            if (depth == m - 2) {
                int last = x ^ sums[depth][span_size - 1];
                if (last <= x) continue;
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
                if (use_gamma) restrict_ann(depth, x);
                recurse(depth + 1, span_size * 2, x);
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
    if (use_gamma && distance <= m / 2)
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
    if (avoid_radical && !use_gamma) return error("avoid-radical requires --gamma");
    if (avoid_radical) initialize_global_radical();
    if (use_gamma) {
        parity_allowed=calloc((size_t)N*word_count,sizeof(uint64_t));
        if (!parity_allowed) return error("parity allocation failed");
        for (int row=0;row<N;++row) for (int x=0;x<N;++x) if (!par(row&x))
            parity_allowed[(size_t)row*word_count+(x>>6)] |= UINT64_C(1)<<(x&63);
    }
    start_time = now();
    sums[0][0] = sizes[0][0] = 0;
    recurse(0, 1, -1);
    double elapsed = now() - start_time;
    if (output_error) {
        free(translated); free(distances); return error("stdout write failed");
    }
    fprintf(stderr, "{\"status\":\"%s\",\"nodes\":%" PRIu64
            ",\"leaves\":%" PRIu64 ",\"gamma_pruned\":%" PRIu64 ",\"xdom_pruned\":%" PRIu64 ",\"seconds\":%.9f,\"reason\":",
            limited ? "LIMIT" : "COMPLETE", nodes, leaf_count, gamma_pruned, xdom_pruned, elapsed);
    if (limited) fprintf(stderr, "\"%s\"", limit_reason);
    else fputs("null", stderr);
    fputs("}\n", stderr);
    free(parity_allowed); free(translated); free(distances);
    return 0;
}
