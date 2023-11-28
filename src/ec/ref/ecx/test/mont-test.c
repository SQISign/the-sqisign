#include <time.h>
#include <assert.h>
#include <stdio.h>

#include "ec.h"
#include "isog.h"
#include "test-basis.h"
#include <bench.h> 

static int BENCH_LOOPS = 1000;       // Number of iterations per bench
static int TEST_LOOPS  = 128;       // Number of iterations per test

// void random_scalar(fp_t k, const uint8_t j)
// {
//         // To implement a better random function (We must use some of the SHAKE family functions)
//         do
//         {
//                 randombytes((void *)k, keyspace_bytes[j]);
//         } while (fp_issmaller((uint64_t *)k, keyspace_size[j]));
// }

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
void ladder3pt(ec_point_t* R, fp_t const m, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A)
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

// For computing [(p + 1) / l_i]P, i:=0, ..., (N - 1)
void cofactor_multiples(ec_point_t P[], ec_point_t const* A, size_t lower, size_t upper)
{
	assert(lower < upper);
	if (upper - lower == 1)
		return ;

	int i;
	size_t mid = lower + (upper - lower + 1) / 2;
	copy_point(&(P[mid]), &(P[lower]));
	for (i = lower; i < (int)mid; i++)
		xMULv2(&(P[mid]), &(P[mid]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], A);
	for (i = (int)mid; i < (int)upper; i++)
		xMULv2(&(P[lower]), &(P[lower]), &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], A);

	cofactor_multiples(P, A, lower, mid);
	cofactor_multiples(P, A, mid, upper);
}

// The projective x-coordinate point (X : Z) at infinity is such that Z == 0
static inline int isinfinity(ec_point_t const P)
{
	return fp2_is_zero(&P.z);
}

int main(int argc, char* argv[])
{
	if (argc > 1) {
		TEST_LOOPS = atoi(argv[1]);
	}

	fp2_t fp2_0, fp2_1;
	fp2_set(&fp2_0, 0);
	fp_mont_setone(fp2_1.re);fp_set(fp2_1.im,0);

	int i, j;

	ec_point_t A;
	fp2_set(&A.x, 0);
	fp_mont_setone(A.z.re);fp_set(A.z.im,0);

	fp2_add(&A.z, &A.z, &A.z);	// 2C
	fp2_add(&A.x, &A.x, &A.z);	// A' + 2C
	fp2_add(&A.z, &A.z, &A.z);	// 4C

	// Just to ensure the projective curve coeffientes are different from zero
	assert( !fp2_is_zero(&A.x) & !fp2_is_zero(&A.x) );

	fp2_t a;
	coeff(&a, A);

	ec_point_t PA, QA, PQA, PB, QB, PQB;

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
	
	ec_point_t P[P_LEN + M_LEN], Q[P_LEN + M_LEN], PQ[P_LEN + M_LEN];
	copy_point(&(P[0]), &PA);
	cofactor_multiples(P, &A, 0, P_LEN);
	copy_point(&(Q[0]), &QA);
	cofactor_multiples(Q, &A, 0, P_LEN);
	copy_point(&(PQ[0]), &PQA);
	cofactor_multiples(PQ, &A, 0, P_LEN);
	for (j = 0; j < P_LEN; j++)
	{
		// x(PA)
		assert( !isinfinity(P[j]) );	// It must be different from the point at infinity
		assert( isrational(P[j], a) );
		xMULv2(&P[j], &P[j], &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
		assert( isinfinity(P[j]) );		// It must be now the point at infinity
		// x(QA)
		assert( !isinfinity(Q[j]) );	// It must be different from the point at infinity
		assert( isrational(Q[j], a) );
		xMULv2(&Q[j], &Q[j], &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
		assert( isinfinity(Q[j]) );		// It must be now the point at infinity
		// x(PQA)
		assert( !isinfinity(PQ[j]) );	// It must be different from the point at infinity
		assert( isrational(PQ[j], a) );
		xMULv2(&PQ[j], &PQ[j], &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
		assert( isinfinity(PQ[j]) );	// It must be now the point at infinity
	};
	// Writing the public projective x-coordinate points into Montogmery domain
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
	// Recall, PB, QB, and PQB are expeted to be M-order points, but we require to ensure they are of order M
	for (j = P_LEN; j < (P_LEN + M_LEN); j++)
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

	copy_point(&(P[P_LEN]), &PB);
	cofactor_multiples(P, &A, P_LEN, P_LEN + M_LEN);
	copy_point(&(Q[P_LEN]), &QB);
	cofactor_multiples(Q, &A, P_LEN, P_LEN + M_LEN);
	copy_point(&(PQ[P_LEN]), &PQB);
	cofactor_multiples(PQ, &A, P_LEN, P_LEN + M_LEN);
	for (j = P_LEN; j < (P_LEN+M_LEN); j++)
	{
		// x(PB)
		assert( !isinfinity(P[j]) );	// It must be different from the point at infinity
		assert( !isrational(P[j], a) );
		xMULv2(&P[j], &P[j], &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
		assert( isinfinity(P[j]) );		// It must be now the point at infinity
		// x(QB)
		assert( !isinfinity(Q[j]) );	// It must be different from the point at infinity
		assert( !isrational(Q[j], a) );
		xMULv2(&Q[j], &Q[j], &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
		assert( isinfinity(Q[j]) );		// It must be now the point at infinity
		// x(PQB)
		assert( !isinfinity(PQ[j]) );	// It must be different from the point at infinity
		assert( !isrational(PQ[j], a) );
		xMULv2(&PQ[j], &PQ[j], &(TORSION_ODD_PRIMES[j]), p_plus_minus_bitlength[j], &A);
		assert( isinfinity(PQ[j]) );	// It must be now the point at infinity
	};

	fp2_t m;

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

	ec_point_t R[P_LEN + M_LEN];
	int k;
	for (j = 0; j < TEST_LOOPS; j++)
	{
		printf("[%3d%%] Testing EC differential arithmetic", 100 * j / TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");
		fp2_random(&m);
		ladder3pt(&(R[0]), m.re, &PA, &QA, &PQA, &A);
		assert( isrational(R[0], a) );
		for (k = 0; k < P_LEN; k++)
		{
			for (i = 1; i < TORSION_ODD_POWERS[k]; i++)
			{
				xMULv2(&R[0], &R[0], &(TORSION_ODD_PRIMES[k]), p_plus_minus_bitlength[k], &A);
				assert( isrational(R[0], a) );
			};
		};
		cofactor_multiples(R, &A, 0, P_LEN);
		for (i = 0; i < P_LEN; i++)
		{
			assert( !isinfinity(R[i]) );	// It must be different from the point at infinity
			assert( isrational(R[i], a) );
			xMULv2(&R[i], &R[i], &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &A);
			assert( isinfinity(R[i]) );		// It must be now the point at infinity
		};

		fp2_random(&m);
		ladder3pt(&(R[P_LEN]), m.re, &PB, &QB, &PQB, &A);
		assert( !isrational(R[P_LEN], a) );
		for (k = P_LEN; k < (P_LEN+M_LEN); k++)
		{
			for (i = 1; i < TORSION_ODD_POWERS[k]; i++)
			{
				xMULv2(&R[P_LEN], &R[P_LEN], &(TORSION_ODD_PRIMES[k]), p_plus_minus_bitlength[k], &A);
				assert( !isrational(R[P_LEN], a) );
			};
		};
		cofactor_multiples(R, &A, P_LEN, P_LEN + M_LEN);
		for (i = P_LEN; i < (P_LEN+M_LEN); i++)
		{
			assert( !isinfinity(R[i]) );	// It must be different from the point at infinity
			assert( !isrational(R[i], a) );
			xMULv2(&R[i], &R[i], &(TORSION_ODD_PRIMES[i]), p_plus_minus_bitlength[i], &A);
			assert( isinfinity(R[i]) );		// It must be now the point at infinity
		};
	};

	if(TEST_LOOPS)
		printf("[%3d%%] Tested EC differential arithmetic:\tNo errors!\n", 100 * j / TEST_LOOPS);
	printf("-- All tests passed.\n");

	// BENCHMARK xDBLv2
    unsigned long long cycles, cycles1, cycles2;
    cycles = 0;
	ec_point_t PP[TEST_LOOPS], EE[TEST_LOOPS];
	for(int i = 0; i < TEST_LOOPS; i++){
		fp2_random(&PP[i].x);
		fp2_random(&PP[i].z);
		fp2_random(&EE[i].x);
		fp2_random(&EE[i].z);
	}
    cycles1 = cpucycles(); 
	for(int i = 0; i < TEST_LOOPS; i++){
		xDBLv2(&PP[i], &PP[i], &EE[i]);
	}
    cycles2 = cpucycles();
    cycles = cycles+(cycles2-cycles1);
	
	printf("xDBLv2 bench: %7lld cycles\n", cycles/TEST_LOOPS);

	// BENCHMARK xIsog4
    cycles = 0;
	ec_point_t KK0[TEST_LOOPS], KK1[TEST_LOOPS], KK2[TEST_LOOPS];
	for(int i = 0; i < TEST_LOOPS; i++){
		fp2_random(&KK0[i].x);
		fp2_random(&KK0[i].z);
		fp2_random(&KK1[i].x);
		fp2_random(&KK1[i].z);
		fp2_random(&KK2[i].x);
		fp2_random(&KK2[i].z);
	}
    cycles1 = cpucycles(); 
	for(int i = 0; i < TEST_LOOPS; i++){
	fp2_t t0, t1;
	fp2_add(&t0, &PP[i].x, &PP[i].z);
	fp2_sub(&t1, &PP[i].x, &PP[i].z);
	fp2_mul(&(EE[i].x), &t0, &KK1[i].x);
	fp2_mul(&(EE[i].z), &t1, &KK2[i].x);
	fp2_mul(&t0, &t0, &t1);
	fp2_mul(&t0, &t0, &KK0[i].x); 
	fp2_add(&t1, &(EE[i].x), &(EE[i].z));
	fp2_sub(&(EE[i].z), &(EE[i].x), &(EE[i].z));
	fp2_sqr(&t1, &t1);
	fp2_sqr(&(EE[i].z), &(EE[i].z));
	fp2_add(&(EE[i].x), &t0, &t1);
	fp2_sub(&t0, &(EE[i].z), &t0);
	fp2_mul(&(EE[i].x), &(EE[i].x), &t1);
	fp2_mul(&(EE[i].z), &(EE[i].z), &t0);
	}
    cycles2 = cpucycles();
    cycles = cycles+(cycles2-cycles1);
	printf("xeval_4 bench: %7lld cycles\n", cycles/TEST_LOOPS);

	return 0;
}
