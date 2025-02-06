#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "bench.h"

#define RANDOMBYTES_MAX_LENGTH 131072
#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

void
randombytes_init_nist(unsigned char *entropy_input,
                      unsigned char *personalization_string,
                      int security_strength);

int
randombytes_nist(unsigned char *x, size_t xlen);

void
RANDOMBYTES_INIT_PLATFORM(unsigned char *entropy_input,
                          unsigned char *personalization_string,
                          int security_strength);

int
RANDOMBYTES_PLATFORM(unsigned char *x, size_t xlen);

int
randombytes_select(void *buf, size_t n);

// run all tests in module
int main(int argc, char *argv[]) {
#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

  printf("Running AES-CTR-DRBG benchmarks\n");

  unsigned char x[RANDOMBYTES_MAX_LENGTH];

  cpucycles_init();

  BENCH_CODE_1(1000 * SQISIGN_TEST_REPS);
  RANDOMBYTES_PLATFORM(x, RANDOMBYTES_MAX_LENGTH);
  BENCH_CODE_2(STRINGIFY(RANDOMBYTES_PLATFORM));

  BENCH_CODE_1(SQISIGN_TEST_REPS);
  randombytes_nist(x, RANDOMBYTES_MAX_LENGTH);
  BENCH_CODE_2("randombytes_nist");

  BENCH_CODE_1(1000 * SQISIGN_TEST_REPS);
  randombytes_select(x, RANDOMBYTES_MAX_LENGTH);
  BENCH_CODE_2("randombytes_system");

  return 0;
}
