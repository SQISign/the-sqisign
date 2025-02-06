#include <inttypes.h>
#include <stdio.h>
#include <time.h>

#include "id2iso_tests.h"
#include "dim2id2iso_tests.h"
#include <hd.h>
#include <rng.h>
#include <bench_test_arguments.h>

// run all tests in module
int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;
    int res = 1;

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

    randombytes_init((unsigned char *)seed, NULL, 256);

    printf("Running id2iso module unit tests\n");

    res = res & id2iso_test_ker2id();

    printf("\nRunning dim2id2iso module unit tests\n");

    printf("\nRunning find uv tests \n");
    int number_test_find_uv = 5;
    for (int i = 0; i < number_test_find_uv; i++) {
        res = res & dim2id2iso_test_find_uv();
    }

    printf("\nRunning fixed degree tests\n");

    res = res & dim2id2iso_test_fixed_degree_isogeny();

    printf("\nRunning id2iso_clapotis tests\n");
    res = res & dim2id2iso_test_dimid2iso();

    if (!res) {
        printf("\nSome tests failed!\n");
    } else {
        printf("\nAll tests passed!\n");
    }
    return (!res);
}
