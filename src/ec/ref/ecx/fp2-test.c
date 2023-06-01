#include <assert.h>
#include <time.h>
#include <stdio.h>
#include "../generic/include/fp2_tmp.h"

int main()
{
	fp2_t fp2_0, fp2_1;
	// ------------
	fp2_set0(fp2_0);
	fp2_set1(fp2_1);
	// ------------

	int i;
	fp2_t a, b, c, d;
	fp_t e;

	for (i = 0; i < 1024; i++)
	{
		printf("[%3d%%] Testing fp2_t arithmetic", 100 * i / (int)1024);
		fflush(stdout);
		printf("\r\x1b[K");
                
		// Random elements of fp
		fp2_random(a);
		fp2_random(b);
		fp2_copy(c, a);
		c.re[0] += 1;
		fp2_copy(d, b);
		d.re[0] -= 1;

		assert(fp2_isequal(a,b) == 0);		// different values check --> (a != b)
		assert(fp2_isequal(c,c) == 1);		// equal values check --> 1 (c == c)

		// Testing neg
		fp2_set0(b);
		fp2_copy(c, a);
		fp2_neg(a, a);
		fp2_sub(c, b, c);
		assert(fp2_isequal(a,c) == 1);

		fp2_set1(a);	// Now a == 1
		fp2_set0(b);	// Now b == 0

		assert(fp2_is_zero(a) == 0);
		assert(fp2_is_zero(b) == 1);

		// testing c - c
		fp2_sub(d, c, c);
		assert(fp2_is_zero(d) == 1);

		// tetsing c * 0
		fp2_mul(d, c, b);
		assert(fp2_is_zero(d) == 1);

		// tetsing c * 1 ... recall, in Montgomery domain R mod p plays the role of the 1
		fp2_set1(a);
		fp2_mul(d, c, a);
		assert(fp2_isequal(d, c) == 1);

		// fp_set(e, 1);	// Now e == 1
		// fp2_pow(d, e, c);
		// assert(fp2_isequal(d, c) == 1);
		
		// fp_set(e, 0);	// Now e == 0
		// fp2_pow(d, e, c);
		// assert(fp2_isone(d) == 1);

		// fp2_set(a, 1);	// Now e == R mod p
		// fp_random(e);
		// fp2_pow(d, e, a);
		// assert(fp2_isone(d) == 1);

		// Testing 1/a by computing (1/a) x a
		fp2_random(a);
		fp2_copy(b, a);
		fp2_inv(a);
		fp2_mul(c, a, b);
		assert(fp2_isone(c) == 1);

		fp2_random(a);
		fp2_sqr(b, a);
		assert( fp2_issquare(b) );

	};

	printf("[%2d%%] Tested fp2_t arithmetic:\tNo errors!\n", 100 * i / (int)1024);
	printf("-- All tests passed.\n");
	return 0;
}
