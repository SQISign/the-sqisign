#ifndef test_utils_h__
#define test_utils_h__

#include "fp.h"
#include "fp2.h"
#include <encoded_sizes.h>
#include <stdlib.h>

#define PASSED 0
#define FAILED 1

// Random elements of fp and fp2, only suitable for testing
void fp_random_test(fp_t *a);
void fp2_random_test(fp2_t *a);

// Comparison of u64 for qsort
int cmp_u64(const void *v1, const void *v2);

#endif
