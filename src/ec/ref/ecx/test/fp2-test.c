#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <fp2.h>
#include <inttypes.h>

static int BENCH_LOOPS = 1000;       // Number of iterations per bench
static int TEST_LOOPS  = 512;       // Number of iterations per test

bool fp2_isequal(fp2_t a, fp2_t b){
    return fp_is_equal(a.re, b.re) && fp_is_equal(a.im, b.im);
}

bool fp2_isone(fp2_t a){
    fp_t one;
    bool res = 1;
    fp_mont_setone(one);
    for(int i = 0; i < NWORDS_FIELD; i++){
        res = res && (a.re[i] == one[i]);
        res = res && (a.im[i] == 0);
    }
    return res;
}

void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.im[i]);
    printf("\n");
}

// VERY NOT SECURE (testing only)
void fp2_random(fp2_t *a){
    for(int i = 0; i < NWORDS_FIELD; i++){
        a->re[i] = rand();
        a->im[i] = rand();
    }
    // Normalize
    fp2_t one;
    fp_mont_setone(one.re);fp_set(one.im,0);
    fp2_mul(&*a, &*a, &one);
    // Update seed
    srand((unsigned) a->re[0]);
}

int main(int argc, char* argv[])
{
	if (argc > 1) {
		TEST_LOOPS = atoi(argv[1]);
	}

	fp2_t fp2_0, fp2_1;
	// ------------
	fp2_set(&fp2_0, 0);
	fp_mont_setone(fp2_1.re);fp_set(fp2_1.im,0);
	// ------------

	int i;
	fp2_t a, b, c, d;
	fp_t e;

	for (i = 0; i < TEST_LOOPS; i++)
	{
		printf("[%3d%%] Testing fp2_t arithmetic", 100 * i / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");
                
		// Random elements of fp
		fp2_random(&a);
		fp2_random(&b);
		fp2_copy(&c, &a);
		c.re[0] += 1;
		fp2_copy(&d, &b);
		d.re[0] -= 1;

		assert(fp2_isequal(a,b) == 0);		// different values check --> (a != b)
		assert(fp2_isequal(c,c) == 1);		// equal values check --> 1 (c == c)

		// Testing neg
		fp2_set(&b, 0);
		fp2_copy(&c, &a);
		fp2_neg(&a, &a);
		fp2_sub(&c, &b, &c);
		assert(fp2_isequal(a,c) == 1);

		fp_mont_setone(a.re);fp_set(a.im,0);	// Now a == 1
		fp2_set(&b, 0);	// Now b == 0

		assert(fp2_is_zero(&a) == 0);
		assert(fp2_is_zero(&b) == 1);

		// testing c - c
		fp2_sub(&d, &c, &c);
		assert(fp2_is_zero(&d) == 1);

		// tetsing c * 0
		fp2_mul(&d, &c, &b);
		assert(fp2_is_zero(&d) == 1);

		// tetsing c * 1 ... recall, in Montgomery domain R mod p plays the role of the 1
		fp_mont_setone(a.re);fp_set(a.im,0);
		fp2_mul(&d, &c, &a);
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
		fp2_random(&a);
		fp2_copy(&b, &a);
		fp2_inv(&a);
		fp2_mul(&c, &a, &b);
		assert(fp2_isone(c) == 1);

		fp2_random(&a);
		fp2_sqr(&b, &a);
		assert( fp2_is_square(&b) );

	};

	if(TEST_LOOPS){
		printf("[%2d%%] Tested fp2_t arithmetic:\tNo errors!\n", 100 * i /TEST_LOOPS);
	}
	printf("-- All tests passed.\n");
	return 0;
}
