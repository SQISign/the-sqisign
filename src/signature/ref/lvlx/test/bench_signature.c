#include <stdio.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include <verification.h>
#include <signature.h>

#include <tools.h>
#include <rng.h>
#include <bench.h>
#include <bench_test_arguments.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

static __inline__ uint64_t
rdtsc(void)
{
    return (uint64_t)cpucycles();
}

void
bench_sqisign(uint64_t bench)
{
    setlocale(LC_NUMERIC, "");
    uint64_t t0, t1;
    clock_t t;
    float ms;

    public_key_t pks[bench];
    secret_key_t sks[bench];
    signature_t sigs[bench];
    unsigned char msg[32] = { 0 };

    for (uint64_t i = 0; i < bench; i++) {
        public_key_init(&pks[i]);
        secret_key_init(&sks[i]);
    }

    printf("\n\nBenchmarking signatures for " STRINGIFY(SQISIGN_VARIANT) ":\n\n");
    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        protocols_keygen(&pks[i], &sks[i]);
    }
    t1 = rdtsc();
    ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    printf("Average keygen time [%.2f ms]\n", (float)(ms / bench));
    printf("\x1b[34mAvg keygen: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        protocols_sign(&(sigs[i]), &(pks[i]), &(sks[i]), msg, sizeof(msg) / sizeof(*msg));
    }
    t1 = rdtsc();
    ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    printf("Average signature time [%.2f ms]\n", (float)(ms / bench));
    printf("\x1b[34mAvg signature: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

    t = tic();
    t0 = rdtsc();
    for (uint64_t i = 0; i < bench; ++i) {
        int check = protocols_verify(&(sigs[i]), &(pks[i]), msg, sizeof(msg) / sizeof(*msg));
        if (!check) {
            printf("verif failed ! \n");
        }
    }
    t1 = rdtsc();
    ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
    printf("Average verification time [%.2f ms]\n", (float)(ms / bench));
    printf("\x1b[34mAvg verification: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

    for (uint64_t i = 0; i < bench; i++) {
        public_key_finalize(&pks[i]);
        secret_key_finalize(&sks[i]);
    }
}

// run all tests in module
int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = SQISIGN_TEST_REPS;
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

    bench_sqisign(iterations);

    return (0);
}
