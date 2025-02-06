#include <bench.h>
#include <bench_test_arguments.h>
#include <rng.h>
#include <id2iso.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

// return 0 if ok, anything else if error
// sample array of odd numbers
int
id2iso_represent_integer_benchmarks_input_generation(ibz_t *n_gammas, int bitsize, int iterations)
{
    int randret = 1;
    ibz_t max;
    ibz_init(&max);
    ibz_pow(&max, &ibz_const_two, bitsize - 1);
    for (int i = 0; i < iterations; i++) {
        randret = randret && ibz_rand_interval(&(n_gammas[i]), &ibz_const_one, &max);
        ibz_mul(&(n_gammas[i]), &ibz_const_two, &(n_gammas[i]));
        ibz_add(&(n_gammas[i]), &ibz_const_one, &(n_gammas[i]));
    }
    ibz_finalize(&max);
    return (!randret);
}

// 0 if ok, 1 otherwise
int
id2iso_represent_integer_benchmarks_test(const quat_alg_elem_t *gamma,
                                         const ibz_t *n_gamma,
                                         int non_diag,
                                         const quat_represent_integer_params_t *params)
{
    ibz_vec_4_t coord;
    ibz_t norm_d, norm_n;
    ibz_init(&norm_d);
    ibz_init(&norm_n);
    ibz_vec_4_init(&coord);
    int res = 0;
    quat_alg_norm(&norm_n, &norm_d, gamma, &QUATALG_PINFTY);
    res = res || !ibz_is_one(&norm_d);
    res = res || !(ibz_cmp(&norm_n, n_gamma) == 0);
    res = res || !quat_lattice_contains(NULL, &((*params).order->order), gamma);
    if (non_diag) {
        if ((*params).order->q == 1) {
            ibz_mul(&norm_n, &ibz_const_two, &ibz_const_two);
            ibz_add(&norm_d, &((*gamma).coord[0]), &((*gamma).coord[3]));
            ibz_mod(&norm_d, &norm_d, &norm_n);
            res = res || (0 != ibz_cmp(&ibz_const_two, &norm_d));
            ibz_sub(&norm_d, &((*gamma).coord[1]), &((*gamma).coord[2]));
            ibz_mod(&norm_d, &norm_d, &norm_n);
            res = res || (0 != ibz_cmp(&ibz_const_two, &norm_d));
        } else {
            quat_lattice_contains(&coord, &((*params).order->order), gamma);
            ibz_gcd(&norm_d, &(coord[1]), &(coord[2]));
            ibz_gcd(&norm_d, &norm_d, &(coord[3]));
            res = res || ibz_is_even(&norm_d);
        }
    }
    ibz_finalize(&norm_d);
    ibz_finalize(&norm_n);
    ibz_vec_4_finalize(&coord);
    return (res);
}

int
id2iso_represent_integer_benchmarks(int bitsize, int iterations, int non_diag, int test, int orders)
{
    int res = 0;
    int randret = 0;
    uint64_t start, end, sum;
    ibz_t *n_gammas;
    ibz_t aux;
    quat_alg_elem_t *gammas;
    n_gammas = malloc(iterations * sizeof(ibz_t));
    gammas = malloc(iterations * sizeof(quat_alg_elem_t));
    for (int i = 0; i < iterations; i++) {
        ibz_init(&(n_gammas[i]));
        quat_alg_elem_init(&(gammas[i]));
    }
    ibz_init(&aux);

    quat_represent_integer_params_t ri_params;
    quat_p_extremal_maximal_order_t order_hnf;
    quat_alg_elem_init(&order_hnf.z);
    quat_alg_elem_init(&order_hnf.t);
    quat_lattice_init(&order_hnf.order);

    printf("Running represent_integer benchmarks for " STRINGIFY(SQISIGN_VARIANT) " with %d iterations, %d orders, "
           "non_diag set to %d, inputs of bitsize %d and the %d-bit prime\n",
           iterations,
           orders,
           non_diag,
           bitsize,
           ibz_bitsize(&(QUATALG_PINFTY.p)));
    if (test) {
        printf("Tests are run on the outputs\n");
    }
    for (int index_alternate_order = 0; index_alternate_order < orders; index_alternate_order++) {

        ri_params.primality_test_iterations = QUAT_represent_integer_params.primality_test_iterations;

        quat_alg_elem_copy(&order_hnf.z, &EXTREMAL_ORDERS[index_alternate_order].z);
        quat_alg_elem_copy(&order_hnf.t, &EXTREMAL_ORDERS[index_alternate_order].t);
        ibz_copy(&order_hnf.order.denom, &EXTREMAL_ORDERS[index_alternate_order].order.denom);
        ibz_mat_4x4_copy(&order_hnf.order.basis, &EXTREMAL_ORDERS[index_alternate_order].order.basis);
        order_hnf.q = EXTREMAL_ORDERS[index_alternate_order].q;
        ri_params.order = &order_hnf;
        ri_params.algebra = &QUATALG_PINFTY;

        // adjust bitsize
        ibz_set(&aux, ri_params.order->q);
        int bitsize_adjustment = ibz_bitsize(&aux);

#ifndef NDEBUG
        assert(quat_lattice_contains(NULL, &ri_params.order->order, &ri_params.order->z));
        assert(quat_lattice_contains(NULL, &ri_params.order->order, &ri_params.order->t));
#endif

        randret = randret || id2iso_represent_integer_benchmarks_input_generation(
                                 n_gammas, bitsize + bitsize_adjustment, iterations);

        if (randret)
            goto fin;

        sum = 0;
        for (int iter = 0; iter < iterations; iter++) {
            start = cpucycles();
            res = res | !quat_represent_integer(&(gammas[iter]), &(n_gammas[iter]), non_diag, &ri_params);
            end = cpucycles();
            sum = sum + end - start;
        }

        printf("For order with q=%d: %" PRIu64 " cycles on average per iteration (adjusted bitsize %d)\n",
               (ri_params.order->q),
               sum / iterations,
               bitsize + bitsize_adjustment);
        if (test) {
            for (int iter = 0; iter < iterations; iter++) {
                res = res || id2iso_represent_integer_benchmarks_test(
                                 &(gammas[iter]), &(n_gammas[iter]), non_diag, &ri_params);
            }
        }
    }

fin:;
    if (randret) {
        printf("Randomness failure in id2iso_represent_integer_benchmarks\n");
    }
    if (res) {
        printf("Tests in id2iso_represent_integer_benchmarks failed\n");
    }

    quat_alg_elem_finalize(&order_hnf.z);
    quat_alg_elem_finalize(&order_hnf.t);
    quat_lattice_finalize(&order_hnf.order);
    ibz_finalize(&aux);
    for (int i = 0; i < iterations; i++) {
        ibz_finalize(&(n_gammas[i]));
        quat_alg_elem_finalize(&(gammas[i]));
    }
    free(n_gammas);
    free(gammas);

    return (res | 2 * randret);
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12];
    int iterations = 10 * SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    int max_orders = NUM_ALTERNATE_EXTREMAL_ORDERS + 1;
    int orders = max_orders;
    int tests = 0;
    int non_diag = 1;
    int bitsize = ibz_bitsize(&(QUATALG_PINFTY.p)) + 25;
    int invalid = 0;

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

        if (non_diag && (strcmp(argv[i], "--disable-non-diag") == 0)) {
            non_diag = 0;
            continue;
        }

        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }

        if (sscanf(argv[i], "--bitsize=%d", &bitsize) == 1) {
            continue;
        }

        if (sscanf(argv[i], "--orders=%d", &orders) == 1) {
            continue;
        }
    }

    invalid = invalid || (argc > 6);
    invalid = invalid || ((non_diag != 1) && (non_diag != 0)) || (iterations < 0);
    invalid = invalid || (bitsize < ibz_bitsize(&(QUATALG_PINFTY.p)) + 20);
    invalid = invalid || (orders > max_orders);
    invalid = invalid || (orders < 0);

    if (help || invalid) {
        if (invalid) {
            printf("Invalid input\n");
        }
        printf("Usage: %s [--bitsize=<bitsize>] [--orders=<orders>] [--iterations=<iterations>] "
               "[--tests] [--disable-non-diag] [--seed=<seed>]\n",
               argv[0]);
        printf("Where <bitsize> is least 20 higher than the bitsize of the prime used in the "
               "algebra(which depends on the level); if no present, uses the default %d; it will be "
               "adjusted for larger maximal orders by this program\n",
               bitsize);
        printf("Where <orders> is the number of maximal orders on which benchmarks are run. It "
               "must be lower than %d;"
               "if not present, uses the default: %d)\n",
               max_orders,
               orders);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        printf("Additional verifications are run on each output if --tests is passed\n");
        printf("If --disable-non-diag is passed, the non_diag flag in the calls is set to 0, otherwise it "
               "defaults to 1. If non_diag=1, benchmarks are only run for one order where q=1.\n");
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

    return id2iso_represent_integer_benchmarks(bitsize, iterations, non_diag, tests, orders);
}
