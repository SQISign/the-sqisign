#include <bench.h>
#include <bench_test_arguments.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "test_utils.h"
#include <rng.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

bool
fp_run(int iterations)
{
    bool OK = true;
    int n, i;
    uint64_t cycles1, cycles2;
    fp_t a, b;
    uint8_t tmp[32];

    fp_random_test(&a);
    fp_random_test(&b);

    printf("\n-------------------------------------------------------------------------------------"
           "-------------------\n\n");
    printf("Benchmarking GF(p) field arithmetic for " STRINGIFY(SQISIGN_VARIANT) ": \n\n");

    // GF(p) addition
    uint64_t cycle_runs[20];

    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_add(&a, &a, &b);
            fp_add(&b, &b, &a);
            fp_add(&a, &a, &b);
            fp_add(&b, &b, &a);
            fp_add(&a, &a, &b);
            fp_add(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) addition runs in .......................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p) subtraction
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_sub(&a, &a, &b);
            fp_sub(&b, &b, &a);
            fp_sub(&a, &a, &b);
            fp_sub(&b, &b, &a);
            fp_sub(&a, &a, &b);
            fp_sub(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) subtraction runs in ....................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p) division by 3
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_div3(&a, &a);
            fp_div3(&b, &b);
            fp_div3(&a, &a);
            fp_div3(&b, &b);
            fp_div3(&a, &a);
            fp_div3(&b, &b);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) division by 3 runs in ..................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p) multiplication
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_mul(&a, &a, &b);
            fp_mul(&b, &b, &a);
            fp_mul(&a, &a, &b);
            fp_mul(&b, &b, &a);
            fp_mul(&a, &a, &b);
            fp_mul(&b, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &b);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) multiplication runs in .................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p) small multiplication
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        uint32_t val = cycles1;
        for (n = 0; n < iterations; n++) {
            fp_mul_small(&a, &a, val);
            fp_mul_small(&a, &a, val);
            fp_mul_small(&a, &a, val);
            fp_mul_small(&a, &a, val);
            fp_mul_small(&a, &a, val);
            fp_mul_small(&a, &a, val);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) small multiplication runs in .............................. %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (6 * iterations),
           tmp[0]);

    // GF(p) squaring
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_sqr(&a, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) squaring runs in .......................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / (iterations),
           tmp[0]);

    // GF(p) inversion
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_inv(&a);
            fp_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) inversion runs in ......................................... %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / iterations,
           tmp[0]);

    // GF(p) sqrt
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_sqrt(&a);
            fp_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  GF(p) sqrt runs in .............................................. %" PRIu64 " cycles, (%u ignore me)\n",
           cycle_runs[4] / iterations,
           tmp[0]);

    // GF(p) is_square
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (n = 0; n < iterations; n++) {
            fp_is_square(&a);
            fp_add(&a, &b, &a);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    fp_encode(tmp, &a);
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  Square checking runs in ......................................... %" PRIu64 " cycles, (%u ignore me)\n",
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

    return !fp_run(iterations);
}
