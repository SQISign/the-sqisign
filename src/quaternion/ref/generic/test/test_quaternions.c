#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "quaternion_tests.h"
#include <rng.h>
#include <bench_test_arguments.h>

// run all tests in module
int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;
    int res = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }
    }

    if (help) {
        printf("Usage: %s [--seed=<seed>]\n", argv[0]);
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

    printf("Running quaternion module unit tests\n");

    res = res | ibz_test_intbig();
    res = res | mini_gmp_test();
    res = res | quat_test_finit();
    res = res | quat_test_dim4();
    res = res | quat_test_dim2();
    res = res | quat_test_integers();
    res = res | quat_test_hnf();
    res = res | quat_test_algebra();
    res = res | quat_test_lattice();
    res = res | quat_test_lll();
    res = res | quat_test_lideal();
    res = res | quat_test_normeq();
    res = res | quat_test_lat_ball();
    res = res | quat_test_with_randomization();
    if (res != 0) {
        printf("\nSome tests failed!\n");
    }
    return (res);
}
