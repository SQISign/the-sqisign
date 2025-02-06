/*
 * A custom SHA-3 / SHAKE implementation is used for pseudorandom (but
 * reproducible) generation of test values.
 */

#include "test_utils.h"
#include "rng.h"

// Make n random-ish field elements (for tests only!).
void
fp_random_test(fp_t *a)
{
    uint8_t tmp[FP_ENCODED_BYTES];

    randombytes(tmp, sizeof(tmp));

    fp_decode_reduce(a, tmp, sizeof(tmp));
}

void
fp2_random_test(fp2_t *a)
{
    fp_random_test(&(a->re));
    fp_random_test(&(a->im));
}

int
cmp_u64(const void *v1, const void *v2)
{
    uint64_t x1 = *(const uint64_t *)v1;
    uint64_t x2 = *(const uint64_t *)v2;
    if (x1 < x2) {
        return -1;
    } else if (x1 == x2) {
        return 0;
    } else {
        return 1;
    }
}
