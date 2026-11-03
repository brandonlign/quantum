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
#define MAX_N 4096
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

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec / 1000000000.0;
}

static const uint64_t *translated_row(int t, int syndrome) {
    return translated + (((size_t)(t - 1) * N + syndrome) * word_count);
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
    for (int w = 0; w < word_count; ++w) {
        while (allowed[w]) {
            int bit = __builtin_ctzll(allowed[w]);
            int x = 64 * w + bit;
            allowed[w] &= allowed[w] - 1;
            xs[depth] = x;
            if (depth == m - 2) {
                int last = x ^ sums[depth][span_size - 1];
                if (last <= x) continue;
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
        } else {
            return error("usage: z_enumerator [--seconds S] [--max-nodes K]");
        }
    }
    if (scanf("%d %d %d", &m, &N, &distance) != 3)
        return error("expected m N d and N distances on stdin");
    if (m < 3 || m > MAX_M || N < 1 || N > MAX_N || (N & (N - 1)) || distance < 3)
        return error("require m=3..10, power-of-two N<=4096, d>=3");
    distances = malloc((size_t)N * sizeof(*distances));
    if (!distances) return error("distance allocation failed");
    int deep_count = 0;
    for (int i = 0; i < N; ++i) {
        if (scanf("%d", distances + i) != 1 || distances[i] < 0)
            return error("expected N nonnegative integer distances");
        if (distances[i] >= distance - 1) ++deep_count;
    }
    /* Match both reference global early-return checks, before rec(). */
    if (deep_count < m) {
        fprintf(stderr, "{\"status\":\"COMPLETE\",\"nodes\":0,\"leaves\":0,\"seconds\":0.0}\n");
        free(distances); return 0;
    }
    word_count = (N + 63) / 64;
    max_t = m / 2;
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
    start_time = now();
    sums[0][0] = sizes[0][0] = 0;
    recurse(0, 1, -1);
    double elapsed = now() - start_time;
    if (output_error) {
        free(translated); free(distances); return error("stdout write failed");
    }
    fprintf(stderr, "{\"status\":\"%s\",\"nodes\":%" PRIu64
            ",\"leaves\":%" PRIu64 ",\"seconds\":%.9f,\"reason\":",
            limited ? "LIMIT" : "COMPLETE", nodes, leaf_count, elapsed);
    if (limited) fprintf(stderr, "\"%s\"", limit_reason);
    else fputs("null", stderr);
    fputs("}\n", stderr);
    free(translated); free(distances);
    return 0;
}
