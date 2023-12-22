#include<time.h>
#include <stdio.h>
#include <assert.h>
#include <inttypes.h>

#include "isog.h"
#include "sdacs.h"
#include "ec.h"
#include "test-basis.h"

void random_scalar(fp_t k, const uint8_t j)
{
    for(int i = 0; i < NWORDS_FIELD; i++)
        k[i] = rand();
}

// Affine Montgomery coefficient computation (A + 2C : 4C) --> A/C
void coeff(fp2_t *B, ec_point_t const A)
{
	fp2_t t;
	fp2_add(&t, &A.x, &A.x);	// (2 * A24)
	fp2_sub(&t, &t, &A.z);	// (2 * A24) - C24

	fp2_copy(&*B, &A.z);
	fp2_inv(&*B);		// 1 / (C24)
	fp2_add(&t, &t, &t);	// 4*A = 2[(2 * A24) - C24]
	fp2_mul(&*B, &t, &*B);	// A/C = 2[(2 * A24) - C24] / C24
}

// Determines if point is fp2-rational (if not, then it must be a zero trace point)
uint8_t isrational(ec_point_t const T, fp2_t const a)
{
	fp2_t XT, tmp, aux, YT_squared;

	fp2_copy(&XT, &T.z);
	fp2_inv(&XT);

	fp2_mul(&XT, &XT, &T.x);

	fp2_sqr(&tmp, &XT);
	fp2_mul(&aux, &tmp, &XT);
	fp2_mul(&tmp, &tmp, &a);
	fp2_add(&YT_squared, &tmp, &aux);
	fp2_add(&YT_squared, &YT_squared, &XT);

	return fp2_is_square(&YT_squared);
}

// ladder3pt computes x(P + [m]Q)
void ladder3pt(ec_point_t *R, fp_t const m, ec_point_t const *P, ec_point_t const *Q, ec_point_t const *PQ, ec_point_t const *A)
{
	ec_point_t X0, X1, X2;
	copy_point(&X0, Q);
	copy_point(&X1, P);
	copy_point(&X2, PQ);

	int i,j;
	digit_t t;
	for (i = 0; i < NWORDS_FIELD; i++)
	{
		t = 1;
		for (j = 0 ; j < RADIX; j++)
		{
			swap_points(&X1, &X2, -((t & m[i]) == 0));
			xDBLADD(&X0, &X1, &X0, &X1, &X2, A);
			swap_points(&X1, &X2, -((t & m[i]) == 0));
			t <<= 1;
		};
	};
	copy_point(R, &X1);
}

// The projective x-coordinate point (X : Z) at infinity is such that Z == 0
static inline int isinfinity(ec_point_t const P)
{
	return fp2_is_zero(&P.z);
}

int main()
{
	
	fp2_t fp2_0, fp2_1;
	fp2_set(&fp2_0, 0);
	fp_mont_setone(fp2_1.re);fp_set(fp2_1.im,0);

	int i, j;

	ec_point_t A, B, T;
	fp2_set(&A.x, 0);
	fp_mont_setone(A.z.re);fp_set(A.z.im,0);
	
	// fp2_add(&A.x, &A.z, &A.x);	// 1
	// fp2_add(&A.x, &A.x, &A.x);	// 2
	// fp2_add(&A.x, &A.z, &A.x);	// 3
	// fp2_add(&A.x, &A.x, &A.x);	// 6

	fp2_add(&A.z, &A.z, &A.z);	// 2C
	fp2_add(&A.x, &A.x, &A.z);	// A' + 2C
	fp2_add(&A.z, &A.z, &A.z);	// 4C

	// Just to ensure the projective curve coeffientes are different from zero
	assert( !fp2_is_zero(&A.x) & !fp2_is_zero(&A.x) );

	fp2_t a;
	coeff(&a, A);

	ec_point_t PA, QA, PQA, PB, QB, PQB, RA, RB;

	// Writing the public projective x-coordinate points into Montogmery domain
	fp2_tomont(&(PA.x), &(xPA));
	fp_mont_setone(PA.z.re);fp_set(PA.z.im,0);
	fp2_tomont(&(QA.x), &(xQA));
	fp_mont_setone(QA.z.re);fp_set(QA.z.im,0);
	fp2_tomont(&(PQA.x), &(xPQA));
	fp_mont_setone(PQA.z.re);fp_set(PQA.z.im,0);

	assert( isrational(PA, a) );
	assert( isrational(QA, a) );
	assert( isrational(PQA, a) );

	fp2_tomont(&(PB.x), &(xPB));
	fp_mont_setone(PB.z.re);fp_set(PB.z.im,0);
	fp2_tomont(&(QB.x), &(xQB));
	fp_mont_setone(QB.z.re);fp_set(QB.z.im,0);
	fp2_tomont(&(PQB.x), &(xPQB));
	fp_mont_setone(PQB.z.re);fp_set(PQB.z.im,0);

	assert( !isrational(PB, a) );
	assert( !isrational(QB, a) );
	assert( !isrational(PQB, a) );
	// ======================================================================================================
	// Recall, PA, QA, and PQA are expeted to be N-order points, but we require to ensure they are of order N
	for (j = 0; j < P_LEN; j++)
	{
		for (i = 1; i < TORSION_ODD_POWERS[j]; i++)
		{
			xMULv2(&PA, &PA, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
			xMULv2(&QA, &QA, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
			xMULv2(&PQA, &PQA, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);

			assert( isrational(PA, a) );
			assert( isrational(QA, a) );
			assert( isrational(PQA, a) );
		};
	};

	assert( !isinfinity(PA) );
	assert( !isinfinity(QA) );
	assert( !isinfinity(PQA) );

	// --------------------------------------------------------------
	fp_t m;
	random_scalar(m, 0);
	ladder3pt(&RA, m, &PA, &QA, &PQA, &A);
	for (i = 0; i < P_LEN; i++)
	{
		printf("// Processing the %d-th prime:\t", i + 1);
		printf("%2d%%", 100 * i / (int)P_LEN);
		fflush(stdout);
		printf("\r\x1b[K");

		copy_point(&T, &RA);
		for (j = (i+1); j < P_LEN; j++)
			xMULv2(&T, &T, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);

		assert( !isinfinity(T) );

		kps(i, T, A);
		if (TORSION_ODD_PRIMES[i] > gap)
			printf("[\033[0;31m%7" PRId "\033[0m] (#I: %3d, #J: %3d, #K: %3d) \n", TORSION_ODD_PRIMES[i], sI, sJ, sK);
		else
			printf("[\033[0;31m%7" PRId "\033[0m] --------------------------- \n", TORSION_ODD_PRIMES[i]);

		xisog(&B, i, A);

		xeval(&PB, i, PB, A);
		coeff(&a, B);
		assert( !isinfinity(PB) );
		assert( !isrational(PB, a) );

		xeval(&RA, i, RA, A);
		assert( (!isinfinity(RA) && (i < (P_LEN - 1))) || (isinfinity(RA) && (i == (P_LEN - 1))) );
		assert( (isrational(RA, a) && (i < (P_LEN - 1))) || (isinfinity(RA) && (i == (P_LEN - 1))) );
		
		copy_point(&A, &B);
		// Verifying the order of the image point of  PA has been reduced 
		copy_point(&T, &RA);
		for (j = (i+1); j < P_LEN; j++)
			xMULv2(&T, &T, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);

		assert( isinfinity(T) );
		kps_clear(i);
	};

	fp2_set(&A.x, 0);
	fp_mont_setone(A.z.re);fp_set(A.z.im,0);
	
	// fp2_add(&A.x, &A.z, &A.x);	// 1
	// fp2_add(&A.x, &A.x, &A.x);	// 2
	// fp2_add(&A.x, &A.z, &A.x);	// 3
	// fp2_add(&A.x, &A.x, &A.x);	// 6

	fp2_add(&A.z, &A.z, &A.z);	// 2C
	fp2_add(&A.x, &A.x, &A.z);	// A' + 2C
	fp2_add(&A.z, &A.z, &A.z);	// 4C

	// Just to ensure the projective curve coeffientes are different from zero
	assert( !fp2_is_zero(&A.x) & !fp2_is_zero(&A.x) );

	coeff(&a, A);
	// Writing the public projective x-coordinate points into Montogmery domain
	fp2_tomont(&(PA.x), &(xPA));
	fp_mont_setone(PA.z.re);fp_set(PA.z.im,0);
	fp2_tomont(&(QA.x), &(xQA));
	fp_mont_setone(QA.z.re);fp_set(QA.z.im,0);
	fp2_tomont(&(PQA.x), &(xPQA));
	fp_mont_setone(PQA.z.re);fp_set(PQA.z.im,0);

	assert( isrational(PA, a) );
	assert( isrational(QA, a) );
	assert( isrational(PQA, a) );

	fp2_tomont(&(PB.x), &(xPB));
	fp_mont_setone(PB.z.re);fp_set(PB.z.im,0);
	fp2_tomont(&(QB.x), &(xQB));
	fp_mont_setone(QB.z.re);fp_set(QB.z.im,0);
	fp2_tomont(&(PQB.x), &(xPQB));
	fp_mont_setone(PQB.z.re);fp_set(PQB.z.im,0);

	assert( !isrational(PB, a) );
	assert( !isrational(QB, a) );
	assert( !isrational(PQB, a) );

	// ======================================================================================================
	// Recall, PA, QA, and PQA are expeted to be N-order points, but we require to ensure they are of order N
	for (j = P_LEN; j < (P_LEN+M_LEN); j++)
	{
		for (i = 1; i < TORSION_ODD_POWERS[j]; i++)
		{
			xMULv2(&PB, &PB, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
			xMULv2(&QB, &QB, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
			xMULv2(&PQB, &PQB, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);

			assert( !isrational(PB, a) );
			assert( !isrational(QB, a) );
			assert( !isrational(PQB, a) );
		};
	};

	assert( !isinfinity(PB) );
	assert( !isinfinity(QB) );
	assert( !isinfinity(PQB) );

	random_scalar(m, 1);
	ladder3pt(&RB, m, &PB, &QB, &PQB, &A);
	for (i = P_LEN; i < (P_LEN+M_LEN); i++)
	{
		printf("// Processing the %d-th prime:\t", i + 1);
		printf("%2d%%", 100 * i / (int)(P_LEN+M_LEN));
		fflush(stdout);
		printf("\r\x1b[K");

		copy_point(&T, &RB);
		for (j = (i+1); j < (P_LEN+M_LEN); j++)
			xMULv2(&T, &T, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);

		assert( !isinfinity(T) );

		kps(i, T, A);
		if (TORSION_ODD_PRIMES[i] > gap)
			printf("[\033[0;31m%7" PRId "\033[0m] (#I: %3d, #J: %3d, #K: %3d) \n", TORSION_ODD_PRIMES[i], sI, sJ, sK);
		else
			printf("[\033[0;31m%7" PRId "\033[0m] --------------------------- \n", TORSION_ODD_PRIMES[i]);

		xisog(&B, i, A);

		xeval(&PA, i, PA, A);
		coeff(&a, B);
		assert( !isinfinity(PA) );
		assert( isrational(PA, a) );

		xeval(&RB, i, RB, A);
		assert( (!isinfinity(RB) && (i < (P_LEN + M_LEN - 1))) || (isinfinity(RB) && (i == (P_LEN + M_LEN - 1))) );
		assert( (!isrational(RB, a) && (i < (P_LEN + M_LEN - 1))) || (isinfinity(RB) && (i == (P_LEN + M_LEN - 1))) );
	
		copy_point(&A, &B);
		// Verifying the order of the image point of  PB has been reduced 
		copy_point(&T, &RB);
		for (j = (i+1); j < (P_LEN+M_LEN); j++)
			xMULv2(&T, &T, &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);

		assert( isinfinity(T) );
		kps_clear(i);
	};

	printf("-- All tests passed!\n");
	return 0;
}
