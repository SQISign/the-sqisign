#include <bench.h>
#include <bench_test_arguments.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "test_utils.h"
#include <rng.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

bool
fp2_run(int iterations)
{
    bool OK = true;
    int n, i;
    uint64_t cycles1, cycles2;
    fp2_t a, b;
    uint8_t tmp[2 * FP_ENCODED_BYTES];

    fp2_random_test(&a);
    fp2_random_test(&b);

    printf("\n-------------------------------------------------------------------------------------"
           "-------------------\n\n");
    printf("Benchmarking GF(p^2) field arithmetic for " STRINGIFY(SQISIGN_VARIANT) ": \n\n");

    // GF(p^2) addition
    uint64_t cycle_runs[20];

    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_add(&a, &a, &b);
            fp2_add(&b, &b, &a);
            fp2_add(&a, &a, &b);
            fp2_add(&b, &b, &a);
            fp2_add(&a, &a, &b);
            fp2_add(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) addition runs in .......................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p^2) subtraction
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_sub(&a, &a, &b);
            fp2_sub(&b, &b, &a);
            fp2_sub(&a, &a, &b);
            fp2_sub(&b, &b, &a);
            fp2_sub(&a, &a, &b);
            fp2_sub(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) subtraction runs in ....................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p^2) multiplication
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_mul(&a, &a, &b);
            fp2_mul(&b, &b, &a);
            fp2_mul(&a, &a, &b);
            fp2_mul(&b, &b, &a);
            fp2_mul(&a, &a, &b);
            fp2_mul(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) multiplication runs in .................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p^2) squaring
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_sqr(&a, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) squaring runs in .......................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (iterations),
           tmp[0]);

    // GF(p^2) small multiplication
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        uint32_t val = cycles1;
        for (n = 0; n < iterations; n++) {
            fp2_mul_small(&a, &a, val);
            fp2_mul_small(&a, &a, val);
            fp2_mul_small(&a, &a, val);
            fp2_mul_small(&a, &a, val);
            fp2_mul_small(&a, &a, val);
            fp2_mul_small(&a, &a, val);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_add(&a.re, &a.re, &a.im);
    fp2_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) small multiplication runs in .............................. %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p^2) inversion
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_inv(&a);
            fp2_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) inversion runs in ......................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / iterations,
           tmp[0]);

    // GF(p^2) sqrt
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_sqrt(&a);
            fp2_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p^2) sqrt runs in .............................................. %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / iterations,
           tmp[0]);

    // GF(p^2) is_square
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp2_is_square(&a);
            fp2_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp2_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  Square checking runs in ........................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / iterations,
           tmp[0]);

    return OK;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 1000 * SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;

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

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);
    cpucycles_init();

    return !fp2_run(iterations);
}
