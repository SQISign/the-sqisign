// SPDX-License-Identifier: Apache-2.0

#include <sig.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>
#include <api.h>


#if defined(TARGET_OS_UNIX) && (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_OTHER))
#include <time.h>
#endif
#if (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_S390X) || defined(TARGET_OTHER))
#define print_unit printf("nsec\n");
#else
#define print_unit printf("cycles\n");
#endif

static int bench_sig(int runs, int csv);
static inline int64_t cpucycles(void);

int main(int argc, char *argv[]) {
    int rc = 0;

#ifndef NDEBUG
fprintf(stderr, "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
                "This will severely impact performance measurements.\x1b[0m\n");
#endif

    if (argc < 2) {
        printf("One argument needed\n");
        rc = 1;
        goto end;
    }
    int runs = atoi(argv[1]);
    rc = bench_sig(runs, 0);
end:
    return rc;
}

#if (defined(TARGET_ARM) || defined(TARGET_ARM64) || defined(TARGET_S390X))
#define BENCH_UNITS "nsec"
#else
#define BENCH_UNITS "cycles"
#endif

int cmpfunc (const void *a, const void *b) {
    return ( *(uint64_t *)a - * (uint64_t *)b );
}

#define BENCH_CODE_1(r) \
    cycles = 0; \
    for (i = 0; i < (r); ++i) { \
        cycles1 = cpucycles();

#define BENCH_CODE_2(name, csv) \
        cycles2 = cpucycles(); \
        if(i < LIST_SIZE) \
          cycles_list[i] = (cycles2 - cycles1);\
        cycles = cycles + (cycles2 - cycles1); \
    } \
    qsort(cycles_list, (runs < LIST_SIZE)? runs : LIST_SIZE, sizeof(uint64_t), cmpfunc);\
    if (csv) \
      printf("%2" PRId64 ",", cycles_list[(runs < LIST_SIZE)? runs/2 : LIST_SIZE/2]); \
    else { \
      printf("  %-20s-> median: %2" PRId64 ", average: %2" PRId64 " ", name, \
      cycles_list[(runs < LIST_SIZE)? runs/2 : LIST_SIZE/2], (cycles / runs)); \
      printf("%s\n", BENCH_UNITS); \
    }

#define LIST_SIZE 10000

static int bench_sig(int runs, int csv) {

    int rc = 0;
    int i;

    int64_t cycles, cycles1, cycles2;
    int64_t cycles_list[10000];

    const int m_len = 32;

    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);
    unsigned char *sig = calloc(CRYPTO_BYTES + m_len, 1);
    unsigned char *m   = calloc(m_len, 1);
    unsigned long long len = CRYPTO_BYTES;

    if (csv) {
        printf("%s,", CRYPTO_ALGNAME);
    } else {
        printf("Benchmarking %s\n", CRYPTO_ALGNAME);
    }

    BENCH_CODE_1(runs);
    sqisign_keypair(pk, sk);
    BENCH_CODE_2("sqisign_keypair", csv);

    BENCH_CODE_1(runs);
    sqisign_sign(sig, &len, m, m_len, sk);
    BENCH_CODE_2("sqisign_sign", csv);

    len = 32;
    BENCH_CODE_1(runs);
    sqisign_open(m, &len, sig, CRYPTO_BYTES, pk);
    BENCH_CODE_2("sqisign_verify", csv);

    if (csv) {
        printf("\n");
    }

    free(pk);
    free(sk);
    free(sig);
    free(m);
    return rc;
}

static inline int64_t cpucycles(void) {
#if (defined(TARGET_AMD64) || defined(TARGET_X86))
    unsigned int hi, lo;

    asm volatile ("rdtsc" : "=a" (lo), "=d"(hi));
    return ((int64_t) lo) | (((int64_t) hi) << 32);
#elif (defined(TARGET_S390X))
    uint64_t tod;
    asm volatile("stckf %0\n" : "=Q" (tod) : : "cc");
    return (tod * 1000 / 4096);
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}
