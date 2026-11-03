/* Exact exhaustive search for the prefix map phi : E -> W/C (see full.py for the mathematics).
   Input (stdin): m N  then N values dmin[0..N-1]  then  nlow, then nlow lines "c t" (coordinate of a
   nonzero element of E in the basis Z_0..Z_{m-2}, X_0X_{m-1}..X_{m-2}X_{m-1}; t = 5 - cost > 0).
   Symmetry breaking (sound: Lemma 10 of the paper): xi_0 < xi_1 < ... < xi_{m-1} where xi_i = phi(Z_i),
   xi_{m-1} = xor of the others;  and phi(X_iX_{m-1}) < phi(X_iX_{m-1}) ^ xi_i  for i < m-1.
   Output: number of nodes and FOUND/NONE (stops at first solution unless argv[1] == "all"). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int m, N, nb, all = 0;
static int *dmin;
static int *lowc[64], *lowt[64], nlow[64];
static int *bad[8], nbad[8];
static unsigned short *phi;
static int y[64];
static long long nodes = 0, sols = 0;
static int dfs(int j) {
    nodes++;
    if (j == nb) { sols++; return !all; }
    unsigned char *allowed = malloc(N);
    memset(allowed, 1, N);
    for (int k = 0; k < nlow[j]; k++) {
        int s = phi[lowc[j][k]], t = lowt[j][k];
        for (int b = 0; b < nbad[t]; b++) allowed[s ^ bad[t][b]] = 0;
    }
    int res = 0;
    for (int v = 0; v < N && !res; v++) {
        if (!allowed[v]) continue;
        if (j < m - 1) {
            if (j > 0 && v <= y[j - 1]) continue;
            if (j == m - 2) { int x = v; for (int i = 0; i < j; i++) x ^= y[i]; if (x <= v) continue; }
        } else {
            int i = j - (m - 1);
            if (v > (v ^ y[i])) continue;
        }
        y[j] = v;
        int L = 1 << j;
        for (int c = 0; c < L; c++) phi[c | L] = phi[c] ^ v;
        res = dfs(j + 1);
    }
    free(allowed);
    return res;
}
int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "all")) all = 1;
    if (scanf("%d %d", &m, &N) != 2) return 1;
    nb = 2 * m - 2;
    dmin = malloc(sizeof(int) * N);
    for (int i = 0; i < N; i++) if (scanf("%d", &dmin[i]) != 1) return 1;
    for (int t = 0; t < 8; t++) { bad[t] = malloc(sizeof(int) * N); nbad[t] = 0;
        for (int i = 0; i < N; i++) if (dmin[i] < t) bad[t][nbad[t]++] = i; }
    int nl; if (scanf("%d", &nl) != 1) return 1;
    for (int j = 0; j < nb; j++) { lowc[j] = malloc(sizeof(int) * nl); lowt[j] = malloc(sizeof(int) * nl); nlow[j] = 0; }
    for (int k = 0; k < nl; k++) {
        int c, t; if (scanf("%d %d", &c, &t) != 2) return 1;
        int hb = 31 - __builtin_clz(c);
        lowc[hb][nlow[hb]] = c ^ (1 << hb); lowt[hb][nlow[hb]] = t; nlow[hb]++;
    }
    phi = calloc(1 << nb, sizeof(unsigned short));
    dfs(0);
    printf("%s nodes=%lld sols=%lld\n", sols ? "FOUND" : "NONE", nodes, sols);
    return 0;
}
