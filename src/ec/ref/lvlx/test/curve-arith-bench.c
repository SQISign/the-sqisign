#include <bench.h>
#include <bench_test_arguments.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>

#include "test_extras.h"
#include <ec.h>
#include <isog.h>
#include <rng.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

uint64_t
bench_xDBL(unsigned int Nbench)
{
    uint64_t cycles0, cycles1;
    unsigned int i;
    ec_point_t P[Nbench], A24[Nbench];
    for (i = 0; i < Nbench; i++) {
        fp2_random_test(&(P[i].x));
        fp2_random_test(&(P[i].z));
        fp2_random_test(&(A24[i].x));
        fp2_random_test(&(A24[i].z));
    }
    cycles0 = cpucycles();
    for (i = 0; i < Nbench; i++) {
        xDBL(&P[i], &P[i], &A24[i]);
    }
    cycles1 = cpucycles();
    return cycles1 - cycles0;
}

uint64_t
bench_xEVAL4(unsigned int Nbench)
{
    uint64_t cycles0, cycles1;
    unsigned int i;
    ec_point_t P[Nbench];
    ec_kps4_t KPS[Nbench];
    for (i = 0; i < Nbench; i++) {
        fp2_random_test(&(P[i].x));
        fp2_random_test(&(P[i].z));
        for (int j = 0; j < 3; j++) {
            fp2_random_test(&(KPS[i].K[j].x));
            fp2_random_test(&(KPS[i].K[j].z));
        }
    }
    cycles0 = cpucycles();
    for (i = 0; i < Nbench; i++) {
        xeval_4(&P[i], &P[i], 1, &KPS[i]);
    }
    cycles1 = cpucycles();
    return cycles1 - cycles0;
}

uint64_t
bench_isog_strategy(unsigned int Nbench)
{
    uint64_t cycles0, cycles1;
    unsigned int i;
    ec_curve_t E0;
    ec_isog_even_t phi[Nbench];
    ec_basis_t basis2;
    ec_curve_init(&E0);
    fp2_set_small(&(E0.A), 6);
    fp2_set_one(&(E0.C));
    (void)ec_curve_to_basis_2f_to_hint(&basis2, &E0, TORSION_EVEN_POWER);
    for (i = 0; i < Nbench; i++) {
        copy_curve(&phi[i].curve, &E0);
        phi[i].length = TORSION_EVEN_POWER;
        if (i == 0) {
            xADD(&phi[i].kernel, &basis2.P, &basis2.Q, &basis2.PmQ);
        }
        if (i == 1) {
            xADD(&phi[i].kernel, &phi[i - 1].kernel, &basis2.Q, &basis2.P);
        }
        if (i > 1) {
            xADD(&phi[i].kernel, &phi[i - 1].kernel, &basis2.Q, &phi[i - 2].kernel);
        }
    }
    cycles0 = cpucycles();
    for (i = 2; i < Nbench; i++) {
        if (ec_eval_even(&phi[i].curve, &phi[i], NULL, 0)) {
            printf("Failed isogeny strategy\n");
            return 0;
        }
    }
    cycles1 = cpucycles();
    return cycles1 - cycles0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 100 * SQISIGN_TEST_REPS;
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

    printf("Benchmarking elliptic curve arithmetic for " STRINGIFY(SQISIGN_VARIANT) ":\n\n");

    uint64_t cycles;

    cycles = bench_xDBL(10 * iterations);
    printf("Bench xDBL_A24:\t%" PRIu64 " cycles\n", cycles / (10 * iterations));

    cycles = bench_xEVAL4(iterations);
    printf("Bench xEVAL4:\t%" PRIu64 " cycles\n", cycles / iterations);

    cycles = bench_isog_strategy(iterations);
    printf("Bench isog strategy:\t%" PRIu64 " cycles\n", cycles / iterations);

    return 0;
}
