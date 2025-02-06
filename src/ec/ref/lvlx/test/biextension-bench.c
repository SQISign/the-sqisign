#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <tools.h>
#include <mp.h>
#include "biextension.h"
#include <rng.h>
#include "bench.h"

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

void
biextension_bench(uint64_t bench)
{
    uint64_t t0, t1;
    uint32_t e = TORSION_EVEN_POWER;

    fp2_t r1;
    ec_curve_t curve;
    ec_point_t tmp;

    digit_t scal_r1[NWORDS_ORDER];
    digit_t scal_r2[NWORDS_ORDER];
    digit_t scal_s1[NWORDS_ORDER];
    digit_t scal_s2[NWORDS_ORDER];

    ec_basis_t BPQ, BRS;

    // Get constants form curve E6 : y^2 = x^3 + 6*x^2 + x
    ec_curve_init(&curve);
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Compute 2^e torsion on curve and copy to a second basis
    (void)ec_curve_to_basis_2f_to_hint(&BPQ, &curve, e);
    copy_basis(&BRS, &BPQ);

    // Benchmark doubling on the curve
    printf("\n\nBenchmarking doublings\n");
    t0 = cpucycles();
    for (uint64_t i = 0; i < bench; ++i) {
        ec_dbl_iter(&tmp, e, &BPQ.P, &curve);
    }
    t1 = cpucycles();
    printf("\x1b[34mAvg doubling: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

    printf("\n\nBenchmarking (Weil) pairings\n");
    t0 = cpucycles();
    for (uint64_t i = 0; i < bench; ++i) {
        weil(&r1, e, &BPQ.P, &BPQ.Q, &BPQ.PmQ, &curve);
    }
    t1 = cpucycles();
    printf("\x1b[34mAvg pairing: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

    printf("\n\nBenchmarking (Weil) dlogs\n");
    t0 = cpucycles();
    for (uint64_t i = 0; i < bench; ++i) {
        ec_dlog_2_weil(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e);
    }
    t1 = cpucycles();
    printf("\x1b[34mAvg pairing dlog: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

    printf("\n\nBenchmarking (Tate) dlogs\n");
    t0 = cpucycles();
    for (uint64_t i = 0; i < bench; ++i) {
        ec_dlog_2_tate(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e);
    }
    t1 = cpucycles();
    printf("\x1b[34mAvg Tate dlog: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);
}

int
main(int argc, char *argv[])
{
    int iterations = 1000 * SQISIGN_TEST_REPS;
    int help = 0;

#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        return 1;
    }

    cpucycles_init();

    printf("Running biextension benchmarks for " STRINGIFY(SQISIGN_VARIANT) ":\n\n");

    biextension_bench(iterations);

    return 0;
}
