#include <stddef.h>
#include <stdio.h>
#include <string.h>

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
  int res = 1;

  printf("Running AES-CTR-DRBG unit tests\n");

  unsigned char seed[48];
  unsigned char x_nist[RANDOMBYTES_MAX_LENGTH], x_platform[RANDOMBYTES_MAX_LENGTH];

  for (int i = 0; i < 8; i++) {
    for (unsigned j = 0; j < sizeof(seed); j++) {
      seed[j] = 1 << i;
    }

    RANDOMBYTES_INIT_PLATFORM(seed, NULL, 256);
    randombytes_init_nist(seed, NULL, 256);

    for (int j = RANDOMBYTES_MAX_LENGTH; j <= RANDOMBYTES_MAX_LENGTH; j *= 2) {
      RANDOMBYTES_PLATFORM(x_platform, j);
      randombytes_nist(x_nist, j);

      if (memcmp(x_platform, x_nist, j) != 0) {
        for (int k = 0; k < j; k++) {
          if (x_platform[k] != x_nist[k]) {
            printf("Test failed for seed = %d, length = %d bytes: mismatch at index %d: %d != %d\n", i, j, k, x_platform[k], x_nist[k]);
            break;
          }
        }
        res = 0;
      }
    }
  }

  if (!res) {
    printf("\nSome tests failed!\n");
  } else {
    printf("\nAll tests passed!\n");
  }

  return (!res);
}
