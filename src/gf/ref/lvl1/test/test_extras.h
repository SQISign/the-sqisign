
#ifndef TEST_EXTRAS_H
#define TEST_EXTRAS_H

#include <time.h>
#include <stdlib.h>
#include "../include/fp.h"
#include "../include/fp2.h"

#define PASSED    0
#define FAILED    1
    
// Access system counter for benchmarking
//int64_t cpucycles(void);

// Comparing "nword" elements, a=b? : (1) a!=b, (0) a=b
int compare_words(digit_t* a, digit_t* b, unsigned int nwords);

// Generating a pseudo-random field element in [0, p-1] 
void fprandom_test(digit_t* a);

// Generating a pseudo-random element in GF(p^2)
void fp2random_test(fp2_t* a);

#endif