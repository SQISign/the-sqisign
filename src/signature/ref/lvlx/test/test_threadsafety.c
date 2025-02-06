#include <inttypes.h>
#include <stdio.h>
#include <assert.h>

#include <verification.h>
#include <signature.h>
#include <rng.h>
#include <bench_test_arguments.h>

#include <pthread.h>

int threads = 4;
int iterations = SQISIGN_TEST_REPS;

char dummy = 0;

void *
test_sqisign(void *_)
{
    (void)_;

    bool res = 1;

    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    unsigned char msg[32] = { 0 };

    public_key_init(&pk);
    secret_key_init(&sk);

    for (int i = 0; i < iterations; ++i) {
        protocols_keygen(&pk, &sk);
        int scheck = protocols_sign(&sig, &pk, &sk, msg, sizeof(msg) / sizeof(*msg));
        int check = protocols_verify(&sig, &pk, msg, sizeof(msg) / sizeof(*msg));
        assert(scheck);
        assert(check);
        res = res && scheck && check;
    }

    public_key_finalize(&pk);
    secret_key_finalize(&sk);

    return res ? &dummy : NULL;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;

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

        if (sscanf(argv[i], "--threads=%d", &threads) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0 || threads <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--threads=<threads>] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for testing; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <threads> is the number of threads used for testing; if not "
               "present, uses the default: %d)\n",
               threads);
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

    pthread_t thread_handles[threads];

    {
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setstacksize(&attr, 8 << 20); // 8 MB
        for (int i = 0; i < threads; ++i)
            pthread_create(&thread_handles[i], &attr, &test_sqisign, NULL);
        pthread_attr_destroy(&attr);
    }

    bool ok = true;
    for (int i = 0; i < threads; ++i) {
        void *res;
        pthread_join(thread_handles[i], &res);
        ok = ok && res;
    }

    if (!ok) {
        printf("\nSome tests failed!\n");
    } else {
        printf("All tests passed!\n");
    }
    return !ok;
}
