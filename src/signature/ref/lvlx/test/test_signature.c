#include <stdio.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include <verification.h>
#include <signature.h>

#include <tools.h>
#include <rng.h>
#include <bench_test_arguments.h>

int
test_sqisign(int repeat)
{
    int res = 1;

    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    const unsigned char msg[32] = { 0 };

    public_key_init(&pk);
    secret_key_init(&sk);

    printf("\n\nTesting signatures\n");
    for (int i = 0; i < repeat; ++i) {
        printf("#%d \n", i);

        protocols_keygen(&pk, &sk);
        protocols_sign(&sig, &pk, &sk, msg, 32);
        int check = protocols_verify(&sig, &pk, msg, 32);
        if (!check) {
            printf("verif failed ! \n");
        }
    }

    public_key_finalize(&pk);
    secret_key_finalize(&sk);

    return res;
}

// run all tests in module
int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    int res;

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
        printf("Where <iterations> is the number of iterations used for testing; if not "
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

    res = test_sqisign(iterations);

    if (!res) {
        printf("\nSome tests failed!\n");
    } else {
        printf("All tests passed!\n");
    }
    return (!res);
}
