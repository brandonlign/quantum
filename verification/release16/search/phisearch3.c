/* phisearch3: exhaustive prefix-map search = phisearch2's forward checking + phisearch's symmetry
   breaking (proof: Lemma 10 of the paper):  (a) xi_0 < ... < xi_{m-2} < xi_{m-1} := xor of the others,
   (b) phi(X_iX_{m-1}) < phi(X_iX_{m-1}) ^ xi_i.
   Order: Z basis vectors 0..m-2 first (fixed order, needed for (a)), then the X-type vectors
   most-constrained-first ((b) is per-vector so any order is fine).
   Input format identical to phisearch/phisearch2.  Optional argv[1] = node limit. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int m, N, nb;
static int *dmin; static unsigned char *T;
static int *bad[8], nbad[8];
static long long nodes = 0, nodelimit = 0;
static int spanC[1 << 16]; static unsigned short spanS[1 << 16];
static int y[32];
static int dfs(int nspan, unsigned rem, unsigned char *allow) {
    nodes++;
    if (nodelimit && nodes > nodelimit) { printf("LIMIT nodes=%lld\n", nodes); exit(2); }
    if (!rem) return 1;
    int g = -1;
    for (int z = 0; z < m - 1; z++) if (rem >> z & 1) { g = z; break; }
    if (g < 0) {
        int bc = N + 1;
        for (int h = 0; h < nb; h++) if (rem >> h & 1) {
            int c = 0; unsigned char *a = allow + h * N;
            for (int v = 0; v < N; v++) c += a[v];
            if (c == 0) return 0;
            if (c < bc) { bc = c; g = h; }
        }
    }
    unsigned rem2 = rem & ~(1u << g);
    unsigned char *na = malloc(nb * N);
    for (int v = 0; v < N; v++) {
        if (!allow[g * N + v]) continue;
        if (g < m - 1) {
            if (g > 0 && v <= y[g - 1]) continue;
            if (g == m - 2) { int x = v; for (int i = 0; i < g; i++) x ^= y[i]; if (x <= v) continue; }
        } else {
            int i = g - (m - 1);
            if (v > (v ^ y[i])) continue;
        }
        y[g] = v;
        memcpy(na, allow, nb * N);
        for (int i = 0; i < nspan; i++) { spanC[nspan + i] = spanC[i] | (1 << g); spanS[nspan + i] = spanS[i] ^ v; }
        int dead = 0;
        for (int h = 0; h < nb && !dead; h++) if (rem2 >> h & 1) {
            unsigned char *a = na + h * N;
            for (int i = nspan; i < 2 * nspan; i++) {
                int t = T[spanC[i] | (1 << h)];
                if (!t) continue;
                int s = spanS[i];
                for (int b = 0; b < nbad[t]; b++) a[s ^ bad[t][b]] = 0;
            }
            int c = 0; for (int w = 0; w < N && !c; w++) c = a[w];
            if (!c) dead = 1;
        }
        if (!dead && dfs(2 * nspan, rem2, na)) { free(na); return 1; }
    }
    free(na);
    return 0;
}
int main(int argc, char **argv) {
    if (argc > 1) nodelimit = atoll(argv[1]);
    if (scanf("%d %d", &m, &N) != 2) return 1;
    nb = 2 * m - 2;
    dmin = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++) if (scanf("%d", &dmin[i]) != 1) return 1;
    for (int t = 0; t < 8; t++) { bad[t] = malloc(sizeof(int) * N); nbad[t] = 0;
        for (int i = 0; i < N; i++) if (dmin[i] < t) bad[t][nbad[t]++] = i; }
    T = calloc(1 << nb, 1);
    int nl; if (scanf("%d", &nl) != 1) return 1;
    for (int k = 0; k < nl; k++) { int c, t; if (scanf("%d %d", &c, &t) != 2) return 1; T[c] = t; }
    unsigned char *allow = malloc(nb * N);
    for (int g = 0; g < nb; g++) for (int v = 0; v < N; v++) allow[g * N + v] = dmin[v] >= T[1 << g];
    spanC[0] = 0; spanS[0] = 0;
    int r = dfs(1, (1u << nb) - 1, allow);
    printf("%s nodes=%lld\n", r ? "FOUND" : "NONE", nodes);
    return 0;
}
