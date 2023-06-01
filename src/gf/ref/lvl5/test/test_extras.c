#include "test_extras.h"
#include <bench.h>

// Global constants
extern const digit_t p[NWORDS_FIELD];
extern const digit_t R2[NWORDS_FIELD];

#if 0
int64_t cpucycles(void)
{ // Access system counter for benchmarking
    unsigned int hi, lo;

    asm volatile ("rdtsc\n\t" : "=a" (lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
}
#endif


int compare_words(digit_t* a, digit_t* b, unsigned int nwords)
{ // Comparing "nword" elements, a=b? : (1) a>b, (0) a=b, (-1) a<b
  // SECURITY NOTE: this function does not have constant-time execution. TO BE USED FOR TESTING ONLY.
    int i;

    for (i = nwords-1; i >= 0; i--)
    {
        if (a[i] > b[i]) return 1;
        else if (a[i] < b[i]) return -1;
    }

    return 0; 
}


static void sub_test(digit_t* out, digit_t* a, digit_t* b, unsigned int nwords)
{ // Subtraction without borrow, out = a-b where a>b
  // SECURITY NOTE: this function does not have constant-time execution. It is for TESTING ONLY.     
    unsigned int i;
    digit_t res, carry, borrow = 0;
  
    for (i = 0; i < nwords; i++)
    {
        res = a[i] - b[i];
        carry = (a[i] < b[i]);
        out[i] = res - borrow;
        borrow = carry || (res < borrow);
    } 
}


void fprandom_test(digit_t* a)
{ // Generating a pseudo-random field element in [0, p-1] 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.
    unsigned int i, diff = 256-254, nwords = NWORDS_FIELD;
    unsigned char* string = NULL;

    string = (unsigned char*)a;
    for (i = 0; i < sizeof(digit_t)*nwords; i++) {
        *(string + i) = (unsigned char)rand();              // Obtain 256-bit number
    }
    a[nwords-1] &= (((digit_t)(-1) << diff) >> diff);

    while (compare_words((digit_t*)p, a, nwords) < 1) {  // Force it to [0, modulus-1]
        sub_test(a, a, (digit_t*)p, nwords);
    }
}


void fp2random_test(fp2_t* a)
{ // Generating a pseudo-random element in GF(p^2) 
  // SECURITY NOTE: distribution is not fully uniform. TO BE USED FOR TESTING ONLY.

    fprandom_test(a->re);
    fprandom_test(a->im);
}