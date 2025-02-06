#include "lll_internals.h"
#include "quaternion_tests.h"
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>

// norms must be either a vector of length iterations of norms or NULL
int
quat_bench_lll_benchmark_lll_core(int bitsize,
                                  int iterations,
                                  int warmup_loops, //(int)(100000 / (bitsize*bitsize)) + 1;
                                  quat_lattice_t *lattices,
                                  const ibz_t *norms,
                                  const quat_alg_t *alg,
                                  int add_tests)
{
    int res = 0;
    uint64_t start, end, time;
    quat_lattice_t test;
    ibz_t num, denom;
    ibq_t eta, delta;
    ibz_mat_4x4_t *reds;
    ibz_mat_4x4_t gram;
    ibz_init(&num);
    ibz_init(&denom);
    ibq_init(&eta);
    ibq_init(&delta);
    quat_lattice_init(&test);
    ibz_mat_4x4_init(&gram);
    quat_lll_set_ibq_parameters(&delta, &eta);

    reds = (ibz_mat_4x4_t *)malloc(iterations * sizeof(ibz_mat_4x4_t));
    for (int i = 0; i < iterations; i++)
        ibz_mat_4x4_init(&(reds[i]));

    // warmup setup
    quat_lattice_t *wu;
    ibz_mat_4x4_t *wur;
    wu = (quat_lattice_t *)malloc(warmup_loops * sizeof(quat_lattice_t));
    wur = (ibz_mat_4x4_t *)malloc(warmup_loops * sizeof(ibz_mat_4x4_t));
    for (int i = 0; i < warmup_loops; i++) {
        quat_lattice_init(&(wu[i]));
        ibz_mat_4x4_init(&(wur[i]));
    }
    quat_test_input_random_lattice_generation(wu, bitsize, warmup_loops, 1);
    // Case for unknown norms
    if (norms == NULL) {
        printf("Start warmup and measures for bitsize %d iterations %d with non-ideals\n", bitsize, iterations);
        // warmup loop
        for (int iter = 0; iter < warmup_loops; iter++) {
            quat_lattice_gram(&gram, &(lattices[iter]), alg);
            ibz_mat_4x4_copy(&(reds[iter]), &(lattices[iter].basis));
            quat_lll_core(&gram, &(reds[iter]));
        }
        // benchmark loop
        start = cpucycles();
        for (int iter = 0; iter < iterations; iter++) {
            quat_lattice_gram(&gram, &(lattices[iter]), alg);
            ibz_mat_4x4_copy(&(reds[iter]), &(lattices[iter].basis));
            quat_lll_core(&gram, &(reds[iter]));
        }
        end = cpucycles();
    } else {
        // Using division by norm
        quat_lattice_t O0;
        quat_left_ideal_t *ideals;
        quat_lattice_init(&O0);
        quat_lattice_O0_set(&O0);
        ideals = (quat_left_ideal_t *)malloc(iterations * sizeof(quat_left_ideal_t));
        for (int i = 0; i < iterations; i++) {
            quat_left_ideal_init(&(ideals[i]));
            ibz_mat_4x4_copy(&(ideals[i].lattice.basis), &(lattices[i].basis));
            ibz_copy(&(ideals[i].lattice.denom), &(lattices[i].denom));
            ibz_copy(&(ideals[i].norm), &(norms[i]));
            (ideals[i].parent_order) = &O0;
        }

        printf("Start warmup and measures for bitsize %d iterations %d with ideals\n", bitsize, iterations);
        // warmup loop
        for (int iter = 0; iter < warmup_loops; iter++) {
            quat_lideal_class_gram(&gram, &(ideals[iter % iterations]), alg);
            ibz_mat_4x4_copy(&(reds[iter % iterations]), &(ideals[iter % iterations].lattice.basis));
            quat_lll_core(&gram, &(reds[iter % iterations]));
        }
        // benchmark loop
        start = cpucycles();
        for (int iter = 0; iter < iterations; iter++) {
            quat_lideal_class_gram(&gram, &(ideals[iter]), alg);
            ibz_mat_4x4_copy(&(reds[iter]), &(ideals[iter].lattice.basis));
            quat_lll_core(&gram, &(reds[iter]));
        }
        end = cpucycles();
        for (int i = 0; i < iterations; i++)
            quat_left_ideal_finalize(&(ideals[i]));
        free(ideals);
        quat_lattice_finalize(&O0);
    }
    // results output
    time = (end - start);
    printf("%" PRIu64 " cycles per lattice\n%d Lattices %" PRIu64 " cycles total\n",
           (uint64_t)(time / iterations),
           iterations,
           time);

    // test loop
    if (add_tests) {
        for (int iter = 0; iter < iterations; iter++) {
            // test lll reduced
            res = res || !quat_lll_verify(&(reds[iter]), &delta, &eta, alg);
            // test lattice equality
            ibz_copy(&(test.denom), &((lattices[iter]).denom));
            ibz_mat_4x4_copy(&(test.basis), &(reds[iter]));
            quat_lattice_hnf(&test);
            res = res || !quat_lattice_equal(&test, &(lattices[iter]));
        }

        if (res != 0) {
            printf("Quaternion benchmark tests for lll_core failed\n");
        }
    }
    printf("\n");

    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibq_finalize(&delta);
    ibq_finalize(&eta);
    quat_lattice_finalize(&test);
    ibz_mat_4x4_finalize(&gram);
    for (int i = 0; i < warmup_loops; i++) {
        quat_lattice_finalize(&(wu[i]));
        ibz_mat_4x4_finalize(&(wur[i]));
    }
    for (int i = 0; i < iterations; i++)
        ibz_mat_4x4_finalize(&(reds[i]));
    free(reds);
    free(wu);
    free(wur);
    return (res);
}

// ideals determines if generated as ideals: 0 means ideals, 1 lattices without HNF, 2 and any other
// number lattices in HNF tests if tests are run
int
quat_bench_lll_one_level(const ibz_t *prime,
                         const int *norm_bitsizes,
                         int nb_bitsizes,
                         int iterations,
                         int ideals,
                         int tests)
{
    // initializations
    quat_alg_t alg;
    quat_p_extremal_maximal_order_t order;
    quat_lattice_t *lattices;
    ibz_t *norms;
    const ibz_t *used_norms;
    quat_represent_integer_params_t params;
    quat_alg_init_set(&alg, prime);
    lattices = malloc(iterations * sizeof(quat_lattice_t));
    norms = malloc(iterations * sizeof(ibz_t));
    // initialize params:
    quat_lattice_init(&(order.order));
    quat_alg_elem_init(&(order.t));
    quat_alg_elem_init(&(order.z));
    quat_lattice_O0_set_extremal(&order);
    params.algebra = &alg;
    params.order = &order;
    params.primality_test_iterations = 30;
    int res = 0;
    int randret = 0;
    if (ideals == 0) {
        used_norms = norms;
    } else {
        used_norms = NULL;
    }

    // run benchmarks
    for (int i = 0; i < iterations; i++) {
        quat_lattice_init(&(lattices[i]));
        ibz_init(&(norms[i]));
    }
    for (int i = 0; i < nb_bitsizes; i++) {
        int warmups = (int)(100000 / (norm_bitsizes[i] * norm_bitsizes[i])) + 1;
        if (ideals == 0) {
            randret = randret || quat_test_input_random_ideal_lattice_generation(
                                     lattices, norms, norm_bitsizes[i], iterations, &params);
        } else {
            randret = randret || quat_test_input_random_lattice_generation(
                                     lattices, norm_bitsizes[i], iterations, ideals - 1); // in HNF
        }
        if (!randret) {
            res = res || quat_bench_lll_benchmark_lll_core(
                             norm_bitsizes[i], iterations, warmups, lattices, used_norms, &alg, tests);
        }
    }
    quat_alg_finalize(&alg);
    for (int i = 0; i < iterations; i++) {
        ibz_finalize(&(norms[i]));
        quat_lattice_finalize(&(lattices[i]));
    }
    quat_lattice_finalize(&(order.order));
    quat_alg_elem_finalize(&(order.t));
    quat_alg_elem_finalize(&(order.z));
    free(norms);
    free(lattices);
    return (res + 2 * (randret));
}

// this function must be adapted if algorithms or parameters change
int
quat_bench_lll_level(int lvl, int iterations, int ideals, int tests)
{
    if (!((lvl == 1) || (lvl == 3) || (lvl == 5))) {
        printf("Invalid input level to quat_bench_lll_level: %d\n", lvl);
        return (1);
    }
    ibz_t prime;
    int norm_bitsizes[4] = { 0 };
    int nb_bitsizes = 4;
    ibz_init(&prime);
    if (lvl == 5) {
        ibz_set_from_str(&prime,
                         "1afffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
                         "ffffffffffffffffffffffffffffffffffffffffffffffffffffff",
                         16);
        norm_bitsizes[0] = 254;
        norm_bitsizes[1] = 254;
        norm_bitsizes[2] = 747;
        norm_bitsizes[3] = 1006;
    }
    if (lvl == 3) {
        ibz_set_from_str(&prime,
                         "40fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff"
                         "fffffffffffffffffffffff",
                         16);
        norm_bitsizes[0] = 193;
        norm_bitsizes[1] = 193;
        norm_bitsizes[2] = 571;
        norm_bitsizes[3] = 761;
    }
    if (lvl == 1) {
        ibz_set_from_str(&prime, "4ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff", 16);
        norm_bitsizes[0] = 127;
        norm_bitsizes[1] = 127;
        norm_bitsizes[2] = 377;
        norm_bitsizes[3] = 501;
    }
    printf("Running benchmarks for lvl %d with bitsizes of [%d, %d, %d, %d] \n",
           lvl,
           norm_bitsizes[0],
           norm_bitsizes[1],
           norm_bitsizes[2],
           norm_bitsizes[3]);
    printf("Using prime ");
    ibz_print(&prime, 10);
    printf("\n");
    printf("With %d iterations, ", iterations);
    if (tests == 0)
        printf("no ");
    printf("tests, and lattices generated as ");
    if (ideals == 0)
        printf("ideals");
    else {
        printf("lattices");
        if (ideals != 1)
            printf("not in hnf");
    }
    printf("\n\n");
    int res = quat_bench_lll_one_level(&prime, norm_bitsizes, nb_bitsizes, iterations, ideals, tests);
    ibz_finalize(&prime);
    return (res);
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12];
    int iterations = SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    int tests = 0;
    int ideals = 0;
    int res = 0;
    int level = 1;
    int level_set = 0;

#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!tests && strcmp(argv[i], "--tests") == 0) {
            tests = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (level_set == 0 && sscanf(argv[i], "--level=%d", &level) == 1) {
            level_set = 1;
            continue;
        }

        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!help && strcmp(argv[i], "--not_ideals") == 0) {
            ideals = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || (level != 1 && level != 3 && level != 5)) {
        printf("Usage: %s [--level=<level>] [--iterations=<iterations>] [--tests] [--seed=<seed>]\n", argv[0]);
        printf("Where <level> is either 1 or 3 or 5; if not passed, runs on all levels\n");
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        printf("Additional verifications are run on each output if --tests is passed\n");
        printf("If --not_ideals is passed, input lattices are not generated as random O0 ideals of "
               "norm of a level-specific bitsize, but as  random lattices with entries of the same "
               "bitsize\n");
        printf("Output has last bit set if tests failed, second-to-last if randomness failed\n");
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

    if (level_set != 0) {
        res = quat_bench_lll_level(level, iterations, ideals, tests);
    } else {
        for (level = 1; level <= 5; level += 2) {
            res |= quat_bench_lll_level(level, iterations, ideals, tests);
        }
    }

    return (res);
}

// Notes on realistic benchmarks:
// Denominator is always 2
// Lattice numerators have following maximal bitsizes (when in HNF):
// lvl5: 254, 254, 755-757, ni 1006
// lvl3: 193, 193, 570-571, ni 760-761
// lvl1: 127, 127, 372-377, ni 500-501
// where only the "ni" ones are not calls through the ideal reduction function.
// The largest coeff of all ideals above is about 2 times its norm.
// The quat_bench_lll_level and main functions takes the above into account.
// In particular lattices are always generated as ideals of O0 in HNF.

// Measures obtained by prints in the L2 function and its applications
