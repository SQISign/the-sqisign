#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <inttypes.h>

#include "isog.h"
#include "test-basis.h"
#include <bench.h> 

static int BENCH_LOOPS = 1000;       // Number of iterations per bench
static int TEST_LOOPS  = 128;       // Number of iterations per test

bool curve_equal(ec_curve_t* E1, ec_curve_t* E2){
	fp2_t a, b;
	fp2_mul(&a, &E1->A, &E2->C);
	fp2_mul(&b, &E2->A, &E1->C);
	return fp2_is_equal(&a, &b);
}

void random_scalar(fp_t k)
{
    for(int i = 0; i < NWORDS_FIELD; i++)
        k[i] = rand();
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

void point_print(char *name, ec_point_t P){
	fp2_t a;
	if(fp2_is_zero(&P.z)){
		printf("%s = INF\n", name);
	}
	else{
	fp2_copy(&a, &P.z);
	fp2_inv(&a);
	fp2_mul(&a, &a, &P.x);
	fp2_print(name, a);
	}
}

void curve_print(char *name, ec_curve_t E){
	fp2_t a;
	fp2_copy(&a, &E.C);
	fp2_inv(&a);
	fp2_mul(&a, &a, &E.A);
	fp2_print(name, a);
}

// Affine Montgomery coefficient computation (A + 2C : 4C) --> A/C
void coeff(fp2_t *B, ec_curve_t const *E)
{
	fp2_t t;
	fp2_add(&t, &E->A, &E->A);	// (2 * A24)
	fp2_sub(&t, &t, &E->C);	// (2 * A24) - C24

	fp2_copy(&*B, &E->C);
	fp2_inv(&*B);		// 1 / (C24)
	fp2_add(&t, &t, &t);	// 4*A = 2[(2 * A24) - C24]
	fp2_mul(&*B, &t, &*B);	// A/C = 2[(2 * A24) - C24] / C24
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

static void xTPL(ec_point_t* Q, const ec_point_t* P, const ec_point_t* A3)
{
    /* ----------------------------------------------------------------------------- *
     * Differential point tripling given the montgomery coefficient A3 = (A+2C:A-2C)
     * ----------------------------------------------------------------------------- */
     
    fp2_t t0, t1, t2, t3, t4;
    fp2_sub(&t0, &P->x, &P->z);
    fp2_sqr(&t2, &t0);
    fp2_add(&t1, &P->x, &P->z);
    fp2_sqr(&t3, &t1);
    fp2_add(&t4, &t1, &t0);
    fp2_sub(&t0, &t1, &t0);
    fp2_sqr(&t1, &t4);
    fp2_sub(&t1, &t1, &t3);
    fp2_sub(&t1, &t1, &t2);
    fp2_mul(&Q->x, &t3, &A3->x);
    fp2_mul(&t3, &Q->x, &t3);
    fp2_mul(&Q->z, &t2, &A3->z);
    fp2_mul(&t2, &t2, &Q->z);
    fp2_sub(&t3, &t2, &t3);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_mul(&t1, &t2, &t1);
    fp2_add(&t2, &t3, &t1);
    fp2_sqr(&t2, &t2);
    fp2_mul(&Q->x, &t2, &t4);
    fp2_sub(&t1, &t3, &t1);
    fp2_sqr(&t1, &t1);
    fp2_mul(&Q->z, &t1, &t0);
}

static void xisog_2_singular(ec_point_t* B24, ec_point_t A24){
	fp2_t t0, four;
	fp_mont_setone(four.re);
	fp_set(four.im, 0);
	fp2_add(&four, &four, &four);
	fp2_add(&four, &four, &four);
	fp2_add(&t0, &A24.x, &A24.x);
	fp2_sub(&t0, &t0, &A24.z);
	fp2_add(&t0, &t0, &t0);
	fp2_inv(&A24.z);
	fp2_mul(&t0, &t0, &A24.z);
	fp2_copy(&K[0].x, &t0);
	fp2_add(&B24->x, &t0, &t0);
	fp2_sqr(&t0, &t0);
	fp2_sub(&t0, &t0, &four);
	fp2_sqrt(&t0);
	fp2_neg(&K[0].z, &t0);
	fp2_add(&B24->z, &t0, &t0);
	fp2_add(&B24->x, &B24->x, &B24->z);
	fp2_add(&B24->z, &B24->z, &B24->z);
}

static void xeval_2_singular(ec_point_t* R, const ec_point_t* Q, const int lenQ){
	fp2_t t0, t1;
	for(int i = 0; i < lenQ; i++){
		fp2_mul(&t0, &Q[i].x, &Q[i].z);
		fp2_mul(&t1, &K[0].x, &Q[i].z);
		fp2_add(&t1, &t1, &Q[i].x);
		fp2_mul(&t1, &t1, &Q[i].x);
		fp2_sqr(&R[i].x, &Q[i].z);
		fp2_add(&R[i].x, &R[i].x, &t1);
		fp2_mul(&R[i].z, &t0, &K[0].z);
	}
}

int main(int argc, char* argv[])
{
    unsigned long long cycles, cycles1, cycles2;

	// ----------------- TEST FOR BASIS GENERATION ----------------- //
	if (argc > 1) {
		TEST_LOOPS = atoi(argv[1]);
	}

	// Initial curve with A = 0
	ec_curve_t E;
	fp2_set(&E.A, 0);
	fp_mont_setone(E.C.re);
	fp_set(E.C.im, 0);

	
	for(int iter = 0; iter < TEST_LOOPS; iter++){
		
		printf("[%3d%%] Testing basis generation", 100 * iter / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");

		// Curve coefficient A24=(A+2C:4C)
		ec_point_t A24;
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		fp2_t j;
		ec_j_inv(&j, &E);

		// Construct basis for 2^f torsion
		ec_basis_t B2;
		ec_curve_to_basis_2(&B2, &E);

		// Check that basis is rational
		assert(ec_is_on_curve(&E, &B2.P));
		assert(ec_is_on_curve(&E, &B2.Q));
		assert(ec_is_on_curve(&E, &B2.PmQ));

		// Compute P+Q
		ec_point_t PpQ;
		xADD(&PpQ, &B2.P, &B2.Q, &B2.PmQ);

		// Check the order
		ec_point_t P2, Q2, PmQ2, PpQ2, R;
		copy_point(&P2, &B2.P);
		copy_point(&Q2, &B2.Q);
		copy_point(&PmQ2, &B2.PmQ);
		copy_point(&PpQ2, &PpQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			assert(!ec_is_zero(&P2));
			xDBLv2(&Q2, &Q2, &A24);
			assert(!ec_is_zero(&Q2));
			xDBLv2(&PmQ2, &PmQ2, &A24);
			assert(!ec_is_zero(&PmQ2));
			xDBLv2(&PpQ2, &PpQ2, &A24);
			assert(!ec_is_zero(&PpQ2));
		}
		assert(is_point_equal(&PmQ2, &PpQ2));
		xDBLv2(&P2, &P2, &A24);
		assert(ec_is_zero(&P2));
		xDBLv2(&Q2, &Q2, &A24);
		assert(ec_is_zero(&Q2));
		xDBLv2(&PmQ2, &PmQ2, &A24);
		assert(ec_is_zero(&PmQ2));
		xDBLv2(&PpQ2, &PpQ2, &A24);
		assert(ec_is_zero(&PpQ2));

		// Check the complete_basis function
		fp_t k;
		random_scalar(k);
		ladder3pt(&R, k, &B2.P, &B2.Q, &B2.PmQ, &A24);
		ec_complete_basis_2(&B2, &E, &R);
		assert(is_point_equal(&R, &B2.P));

		// Check that basis is rational
		assert(ec_is_on_curve(&E, &B2.P));
		assert(ec_is_on_curve(&E, &B2.Q));
		assert(ec_is_on_curve(&E, &B2.PmQ));

		// Compute P+Q
		xADD(&PpQ, &B2.P, &B2.Q, &B2.PmQ);

		// Check the order
		copy_point(&P2, &B2.P);
		copy_point(&Q2, &B2.Q);
		copy_point(&PmQ2, &B2.PmQ);
		copy_point(&PpQ2, &PpQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			assert(!ec_is_zero(&P2));
			xDBLv2(&Q2, &Q2, &A24);
			assert(!ec_is_zero(&Q2));
			xDBLv2(&PmQ2, &PmQ2, &A24);
			assert(!ec_is_zero(&PmQ2));
			xDBLv2(&PpQ2, &PpQ2, &A24);
			assert(!ec_is_zero(&PpQ2));
		}
		assert(is_point_equal(&PmQ2, &PpQ2));
		xDBLv2(&P2, &P2, &A24);
		assert(ec_is_zero(&P2));
		xDBLv2(&Q2, &Q2, &A24);
		assert(ec_is_zero(&Q2));
		xDBLv2(&PmQ2, &PmQ2, &A24);
		assert(ec_is_zero(&PmQ2));
		xDBLv2(&PpQ2, &PpQ2, &A24);
		assert(ec_is_zero(&PpQ2));

		// Curve coefficient A3=(A+2C:A-2C)
		ec_point_t A3;
		fp2_copy(&A3.x, &A24.x);
		fp2_sub(&A3.z, &A3.x, &A24.z);

		// Construct basis for 3^g torsion
		ec_basis_t B3;
		ec_curve_to_basis_3(&B3, &E);

		// Check that basis is rational
		assert(ec_is_on_curve(&E, &B3.P));
		assert(ec_is_on_curve(&E, &B3.Q));
		assert(ec_is_on_curve(&E, &B3.PmQ));

		// Compute P+Q
		xADD(&PpQ, &B3.P, &B3.Q, &B3.PmQ);

		// Check the order
		ec_point_t P3, Q3, PmQ3, PpQ3;
		copy_point(&P3, &B3.P);
		copy_point(&Q3, &B3.Q);
		copy_point(&PmQ3, &B3.PmQ);
		copy_point(&PpQ3, &PpQ);
		for(int i = 0; i < POWER_OF_3 - 1; i++){
			xTPL(&P3, &P3, &A3);
			assert(!ec_is_zero(&P3));
			xTPL(&Q3, &Q3, &A3);
			assert(!ec_is_zero(&Q3));
			xTPL(&PmQ3, &PmQ3, &A3);
			assert(!ec_is_zero(&PmQ3));
			xTPL(&PpQ3, &PpQ3, &A3);
			assert(!ec_is_zero(&PpQ3));
		}
		xADD(&P2, &PpQ3, &Q3, &P3);
		assert(is_point_equal(&PmQ3, &P2));
		xTPL(&P3, &P3, &A3);
		assert(ec_is_zero(&P3));
		xTPL(&Q3, &Q3, &A3);
		assert(ec_is_zero(&Q3));
		xTPL(&PmQ3, &PmQ3, &A3);
		assert(ec_is_zero(&PmQ3));
		xTPL(&PpQ3, &PpQ3, &A3);
		assert(ec_is_zero(&PpQ3));

        // Construct basis for 2^f*3^g torsion
        ec_basis_t B6;
        ec_curve_to_basis_6(&B6, &E);

        // Check that basis is rational
        assert(ec_is_on_curve(&E, &B6.P));
        assert(ec_is_on_curve(&E, &B6.Q));
        assert(ec_is_on_curve(&E, &B6.PmQ));

        // Compute P+Q
        xADD(&PpQ, &B6.P, &B6.Q, &B6.PmQ);

        // Check the order
        ec_point_t P6, Q6, PmQ6, PpQ6;
        copy_point(&P6, &B6.P);
        copy_point(&Q6, &B6.Q);
        copy_point(&PmQ6, &B6.PmQ);
        copy_point(&PpQ6, &PpQ);

        for(int i = 0; i < POWER_OF_2 - 1; i++){
            xDBLv2(&P6, &P6, &A24);
            assert(!ec_is_zero(&P6));
            xDBLv2(&Q6, &Q6, &A24);
            assert(!ec_is_zero(&Q6));
            xDBLv2(&PmQ6, &PmQ6, &A24);
            assert(!ec_is_zero(&PmQ6));
            xDBLv2(&PpQ6, &PpQ6, &A24);
            assert(!ec_is_zero(&PpQ6));
        }

        for(int i = 0; i < POWER_OF_3 - 1; i++){
            xTPL(&P6, &P6, &A3);
            assert(!ec_is_zero(&P6));
            xTPL(&Q6, &Q6, &A3);
            assert(!ec_is_zero(&Q6));
            xTPL(&PmQ6, &PmQ6, &A3);
            assert(!ec_is_zero(&PmQ6));
            xTPL(&PpQ6, &PpQ6, &A3);
            assert(!ec_is_zero(&PpQ6));
        }
        copy_point(&R, &P6);
        xDBLv2(&P6, &P6, &A24);
        assert(!ec_is_zero(&P6));
        xTPL(&P6, &P6, &A3);
        assert(ec_is_zero(&P6));
        xTPL(&R, &R, &A3);
        assert(!ec_is_zero(&R));
        xDBLv2(&R, &R, &A24);
        assert(ec_is_zero(&R));


        copy_point(&R, &Q6);
        xDBLv2(&Q6, &Q6, &A24);
        assert(!ec_is_zero(&Q6));
        xTPL(&Q6, &Q6, &A3);
        assert(ec_is_zero(&Q6));
        xTPL(&R, &R, &A3);
        assert(!ec_is_zero(&R));
        xDBLv2(&R, &R, &A24);
        assert(ec_is_zero(&R));


        copy_point(&R, &PmQ6);
        xDBLv2(&PmQ6, &PmQ6, &A24);
        assert(!ec_is_zero(&PmQ6));
        xTPL(&PmQ6, &PmQ6, &A3);
        assert(ec_is_zero(&PmQ6));
        xTPL(&R, &R, &A3);
        assert(!ec_is_zero(&R));
        xDBLv2(&R, &R, &A24);
        assert(ec_is_zero(&R));


        copy_point(&R, &PpQ6);
        xDBLv2(&PpQ6, &PpQ6, &A24);
        assert(!ec_is_zero(&PpQ6));
        xTPL(&PpQ6, &PpQ6, &A3);
        assert(ec_is_zero(&PpQ6));
        xTPL(&R, &R, &A3);
        assert(!ec_is_zero(&R));
        xDBLv2(&R, &R, &A24);
        assert(ec_is_zero(&R));


		// Compute a 2^e-basis with 2^(e-1)*Q=(0,0)
		copy_point(&P2, &B2.P);
		copy_point(&Q2, &B2.Q);
		copy_point(&PmQ2, &B2.PmQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			xDBLv2(&Q2, &Q2, &A24);
			xDBLv2(&PmQ2, &PmQ2, &A24);
		}
		if(fp2_is_zero(&P2.x)){
			copy_point(&Q2, &B2.Q);
			copy_point(&B2.Q, &B2.P);
			copy_point(&B2.P, &Q2);
		}
		else if(fp2_is_zero(&PmQ2.x)){
			copy_point(&Q2, &B2.Q);
			copy_point(&B2.Q, &B2.PmQ);
			copy_point(&B2.PmQ, &Q2);
		}

		// Compute a 2^e-isogeny
		ec_isog_even_t isog;
		random_scalar(k);
		ladder3pt(&R, k, &B2.P, &B2.Q, &B2.PmQ, &A24);
		fp2_copy(&isog.curve.A, &E.A);
		fp2_copy(&isog.curve.C, &E.C);
		copy_point(&isog.kernel, &R);
		isog.length = POWER_OF_2;
		ec_eval_even_nonzero(&E, &isog, &B3.P, 1);
	}
	printf("[%2d%%] Tested basis generation:\t\tNo errors!\n", 100);
	


	// ----------------- TEST FOR NONZERO 2-ISOGENIES VS 4-ISOGENIES----------------- //

	// Initial curve with A = 0
	fp2_set(&E.A, 0);
	fp_mont_setone(E.C.re);
	fp_set(E.C.im, 0);

	for(int iter = 0; iter < TEST_LOOPS; iter++){

		printf("[%3d%%] Testing 2-isog vs 4-isog", 100 * iter / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");

		// Curve coefficient A24=(A+2C:4C)
		ec_point_t A24, B24, C24;
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Compute a 2^e-basis with 2^(e-1)*Q=(0,0)
		ec_basis_t B0, B1, B2;
		ec_point_t P4, Q4, PmQ4, P2, Q2, PmQ2, tmp;
		ec_curve_to_basis_2(&B2, &E);
		copy_point(&P4, &B2.P);
		copy_point(&Q4, &B2.Q);
		copy_point(&PmQ4, &B2.PmQ);
		for(int i = 0; i < POWER_OF_2 - 2; i++){
			xDBLv2(&P4, &P4, &A24);
			xDBLv2(&Q4, &Q4, &A24);
			xDBLv2(&PmQ4, &PmQ4, &A24);
		}
		xDBLv2(&P2, &P4, &A24);
		xDBLv2(&Q2, &Q4, &A24);
		xDBLv2(&PmQ2, &PmQ4, &A24);
		if(fp2_is_zero(&P2.x)){
			copy_point(&tmp, &Q4);
			copy_point(&Q4, &P4);
			copy_point(&P4, &tmp);
		}
		else if(fp2_is_zero(&PmQ2.x)){
			copy_point(&tmp, &Q4);
			copy_point(&Q4, &PmQ4);
			copy_point(&PmQ4, &tmp);
		}

		// Non-singular 2-isogenies
		xDBLv2(&P2, &P4, &A24);
		xisog_2(&B24, P2);
		xeval_2((ec_point_t*)&B0, (ec_point_t*)&B2, 3);
		xeval_2(&P2, &P4, 1);
		xisog_2(&B24, P2);
		xeval_2((ec_point_t*)&B0, (ec_point_t*)&B0, 3);

		// Non-singular 4-isogeny
		xisog_4(&C24, P4);
		xeval_4((ec_point_t*)&B1, (ec_point_t*)&B2, 3);

		// Compare results
		assert(ec_is_equal(&B24, &C24));
		assert(ec_is_equal(&B0.P, &B1.P));
		assert(ec_is_equal(&B0.Q, &B1.Q));
		assert(ec_is_equal(&B0.PmQ, &B1.PmQ));

		// Singular 2-isogenies case 1
		xisog_2_singular(&B24, A24);
		xeval_2_singular((ec_point_t*)&B0, (ec_point_t*)&B2, 3);
		xeval_2_singular(&Q2, &Q4, 1);
		xisog_2(&B24, Q2);
		xeval_2((ec_point_t*)&B0, (ec_point_t*)&B0, 3);

		// Singular 4-isogeny case 1
		xisog_4_singular(&C24, Q4, A24);
		xeval_4_singular((ec_point_t*)&B1, (ec_point_t*)&B2, 3, Q4);

		// Compare results
		assert(ec_is_equal(&B24, &C24));
		assert(ec_is_equal(&B0.P, &B1.P));
		assert(ec_is_equal(&B0.Q, &B1.Q));
		assert(ec_is_equal(&B0.PmQ, &B1.PmQ));

		// Singular 2-isogenies case 2
		fp2_sub(&A24.x, &A24.z, &A24.x);
		fp2_neg(&Q4.x, &Q4.x);
		fp2_neg(&P4.x, &P4.x);
		fp2_neg(&B2.P.x, &B2.P.x);
		fp2_neg(&B2.Q.x, &B2.Q.x);
		fp2_neg(&B2.PmQ.x, &B2.PmQ.x);
		xisog_2_singular(&B24, A24);
		xeval_2_singular((ec_point_t*)&B0, (ec_point_t*)&B2, 3);
		xeval_2_singular(&Q2, &Q4, 1);
		xisog_2(&B24, Q2);
		xeval_2((ec_point_t*)&B0, (ec_point_t*)&B0, 3);

		// Singular 4-isogeny case 2
		xisog_4_singular(&C24, Q4, A24);
		xeval_4_singular((ec_point_t*)&B1, (ec_point_t*)&B2, 3, Q4);

		// Compare results
		assert(ec_is_equal(&B24, &C24));
		assert(ec_is_equal(&B0.P, &B1.P));
		assert(ec_is_equal(&B0.Q, &B1.Q));
		assert(ec_is_equal(&B0.PmQ, &B1.PmQ));

		// Move to next curve
		xisog_4(&A24, P4);
		fp2_add(&E.A, &A24.x, &A24.x);
		fp2_sub(&E.A, &E.A, &A24.z);
		fp2_sub(&E.A, &E.A, &E.A);
		fp2_copy(&E.C, &A24.z);
	}
	printf("[%2d%%] Tested 2-isog vs 4-isog:\t\tNo errors!\n", 100);
	


	// ----------------- TEST FOR NONZERO 2^f ISOGENIES ----------------- //

	// Initial curve with A = 0
	fp2_set(&E.A, 0);
	fp_mont_setone(E.C.re);
	fp_set(E.C.im, 0);

	for(int iter = 0; iter < TEST_LOOPS; iter++){

		printf("[%3d%%] Testing 2^f-isogenies (nonzero)", 100 * iter / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");

		// Curve coefficient A24=(A+2C:4C)
		ec_point_t A24;
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Compute a 2^e-basis with 2^(e-1)*Q=(0,0)
		ec_basis_t B2;
		ec_point_t P2, Q2, PmQ2;
		ec_curve_to_basis_2(&B2, &E);
		copy_point(&P2, &B2.P);
		copy_point(&Q2, &B2.Q);
		copy_point(&PmQ2, &B2.PmQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			xDBLv2(&Q2, &Q2, &A24);
			xDBLv2(&PmQ2, &PmQ2, &A24);
		}
		if(fp2_is_zero(&P2.x)){
			copy_point(&Q2, &B2.Q);
			copy_point(&B2.Q, &B2.P);
			copy_point(&B2.P, &Q2);
		}
		else if(fp2_is_zero(&PmQ2.x)){
			copy_point(&Q2, &B2.Q);
			copy_point(&B2.Q, &B2.PmQ);
			copy_point(&B2.PmQ, &Q2);
		}

		// Generate 3^g-basis and a 2^f-kernel point
		ec_basis_t B3;
		fp_t k;
		ec_point_t R;
		ec_curve_to_basis_3(&B3, &E);
		random_scalar(k);
		ladder3pt(&R, k, &B2.P, &B2.Q, &B2.PmQ, &A24);

		// Evaluate 2^f-isogeny
		ec_isog_even_t isog;
		fp2_copy(&isog.curve.A, &E.A);
		fp2_copy(&isog.curve.C, &E.C);
		isog.length = POWER_OF_2;
		for(int i = isog.length; i < POWER_OF_2; i++)
			xDBLv2(&R, &R, &A24);
		copy_point(&isog.kernel, &R);
		ec_eval_even_nonzero(&E, &isog, (ec_point_t*)&B3, 3);

		// Curve coefficient A3 = (A+2C:A-2C)
		ec_point_t A3;
		fp2_add(&A3.z, &E.C, &E.C);
		fp2_add(&A3.x, &E.A, &A3.z);
		fp2_sub(&A3.z, &E.A, &A3.z);

		// Check order of the pushed 3^g-basis
		ec_point_t P3, Q3, PmQ3, PpQ, PpQ3;
		xADD(&PpQ, &B3.P, &B3.Q, &B3.PmQ);
		copy_point(&P3, &B3.P);
		copy_point(&Q3, &B3.Q);
		copy_point(&PmQ3, &B3.PmQ);
		copy_point(&PpQ3, &PpQ);
		for(int i = 0; i < POWER_OF_3 - 1; i++){
			xTPL(&P3, &P3, &A3);
			assert(!ec_is_zero(&P3));
			xTPL(&Q3, &Q3, &A3);
			assert(!ec_is_zero(&Q3));
			xTPL(&PmQ3, &PmQ3, &A3);
			assert(!ec_is_zero(&PmQ3));
			xTPL(&PpQ3, &PpQ3, &A3);
			assert(!ec_is_zero(&PpQ3));
		}
		xADD(&R, &PpQ3, &Q3, &P3);
		assert(is_point_equal(&PmQ3, &R));
		xTPL(&P3, &P3, &A3);
		assert(ec_is_zero(&P3));
		xTPL(&Q3, &Q3, &A3);
		assert(ec_is_zero(&Q3));
		xTPL(&PmQ3, &PmQ3, &A3);
		assert(ec_is_zero(&PmQ3));
		xTPL(&PpQ3, &PpQ3, &A3);
		assert(ec_is_zero(&PpQ3));
	}
	printf("[%2d%%] Tested 2^f-isogenies (nonzero):\tNo errors!\n", 100);
	


	// ----------------- TEST FOR 2^f ISOGENIES ----------------- //

	// Initial curve with A = 0
	fp2_set(&E.A, 0);
	fp_mont_setone(E.C.re);
	fp_set(E.C.im, 0);

    cycles = 0;
	for(int iter = 0; iter < TEST_LOOPS; iter++){

		printf("[%3d%%] Testing 2^f-isogenies", 100 * iter / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");

		// Curve coefficient A24=(A+2C:4C)
		ec_point_t A24;
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Compute a 2^e-basis
		ec_basis_t B2;
		ec_curve_to_basis_2(&B2, &E);

		// Generate 3^g-basis and a 2^f-kernel point
		ec_basis_t B3;
		fp_t k;
		ec_point_t R;
		ec_curve_to_basis_3(&B3, &E);
		random_scalar(k);
		ladder3pt(&R, k, &B2.P, &B2.Q, &B2.PmQ, &A24);

		// Evaluate 2^f-isogeny
		ec_isog_even_t isog;
		fp2_copy(&isog.curve.A, &E.A);
		fp2_copy(&isog.curve.C, &E.C);
		isog.length = POWER_OF_2;
		for(int i = isog.length; i < POWER_OF_2; i++)
			xDBLv2(&R, &R, &A24);
		copy_point(&isog.kernel, &R);
        cycles1 = cpucycles(); 
		ec_eval_even_basis(&E, &isog, &B3, 1);
        cycles2 = cpucycles();
        cycles = cycles+(cycles2-cycles1);

		// Curve coefficient A3 = (A+2C:A-2C)
		ec_point_t A3;
		fp2_add(&A3.z, &E.C, &E.C);
		fp2_add(&A3.x, &E.A, &A3.z);
		fp2_sub(&A3.z, &E.A, &A3.z);

		// Check order of the pushed 3^g-basis
		ec_point_t P3, Q3, PmQ3, PpQ, PpQ3;
		xADD(&PpQ, &B3.P, &B3.Q, &B3.PmQ);
		copy_point(&P3, &B3.P);
		copy_point(&Q3, &B3.Q);
		copy_point(&PmQ3, &B3.PmQ);
		copy_point(&PpQ3, &PpQ);
		for(int i = 0; i < POWER_OF_3 - 1; i++){
			xTPL(&P3, &P3, &A3);
			assert(!ec_is_zero(&P3));
			xTPL(&Q3, &Q3, &A3);
			assert(!ec_is_zero(&Q3));
			xTPL(&PmQ3, &PmQ3, &A3);
			assert(!ec_is_zero(&PmQ3));
			xTPL(&PpQ3, &PpQ3, &A3);
			assert(!ec_is_zero(&PpQ3));
		}
		xADD(&R, &PpQ3, &Q3, &P3);
		assert(is_point_equal(&PmQ3, &R));
		xTPL(&P3, &P3, &A3);
		assert(ec_is_zero(&P3));
		xTPL(&Q3, &Q3, &A3);
		assert(ec_is_zero(&Q3));
		xTPL(&PmQ3, &PmQ3, &A3);
		assert(ec_is_zero(&PmQ3));
		xTPL(&PpQ3, &PpQ3, &A3);
		assert(ec_is_zero(&PpQ3));
	}
	printf("[%2d%%] Tested 2^f-isogenies:\t\tNo errors! (%7lld cycles)\n", 100, cycles/TEST_LOOPS);\

	// ----------------- TEST FOR 3^g ISOGENIES ----------------- //

	// Initial curve with A = 0
	fp2_set(&E.A, 0);
	fp_mont_setone(E.C.re);
	fp_set(E.C.im, 0);

	for(int iter = 0; iter < TEST_LOOPS; iter++){

		printf("[%3d%%] Testing 3^g-isogenies", 100 * iter / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");

		// Curve coefficient A24=(A+2C:4C)
		ec_point_t A24;
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Curve coefficient A3=(A+2C:A-2C)
		ec_point_t A3;
		fp2_copy(&A3.x, &A24.x);
		fp2_sub(&A3.z, &A24.x, &A24.z);

		// Compute bases and a kernel point
		ec_basis_t B2, B3;
		fp_t k;
		ec_point_t R;
		ec_curve_to_basis_2(&B2, &E);
		ec_curve_to_basis_3(&B3, &E);
		random_scalar(k);
		ladder3pt(&R, k, &B3.P, &B3.Q, &B3.PmQ, &A24);

		// Evaluate 3^g-isogeny
		ec_isog_odd_t isog;
		fp2_copy(&isog.curve.A, &E.A);
		fp2_copy(&isog.curve.C, &E.C);
		for(int i = 0; i < P_LEN+M_LEN; i++)
			isog.degree[i] = 0;
		isog.degree[0] = rand() % POWER_OF_3 + 1;
		for(int i = isog.degree[0]; i < POWER_OF_3; i++)
			xTPL(&R, &R, &A3);
		copy_point(&isog.ker_plus, &R);
		ec_eval_odd_basis(&E, &isog, &B2, 1);

		// Update A24
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Check order of the pushed 2^f-basis
		ec_point_t P2, Q2, PmQ2, PpQ, PpQ2;
		xADD(&PpQ, &B2.P, &B2.Q, &B2.PmQ);
		copy_point(&P2, &B2.P);
		copy_point(&Q2, &B2.Q);
		copy_point(&PmQ2, &B2.PmQ);
		copy_point(&PpQ2, &PpQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			assert(!ec_is_zero(&P2));
			xDBLv2(&Q2, &Q2, &A24);
			assert(!ec_is_zero(&Q2));
			xDBLv2(&PmQ2, &PmQ2, &A24);
			assert(!ec_is_zero(&PmQ2));
			xDBLv2(&PpQ2, &PpQ2, &A24);
			assert(!ec_is_zero(&PpQ2));
		}
		assert(is_point_equal(&PmQ2, &PpQ2));
		xDBLv2(&P2, &P2, &A24);
		assert(ec_is_zero(&P2));
		xDBLv2(&Q2, &Q2, &A24);
		assert(ec_is_zero(&Q2));
		xDBLv2(&PmQ2, &PmQ2, &A24);
		assert(ec_is_zero(&PmQ2));
		xDBLv2(&PpQ2, &PpQ2, &A24);
		assert(ec_is_zero(&PpQ2));
	}
	printf("[%2d%%] Tested 3^g-isogenies:\t\tNo errors!\n", 100);

	// ----------------- TEST FOR ODD-DEGREE ISOGENIES ----------------- //

	for(int iter = 0; iter < TEST_LOOPS/10; iter++){

		printf("[%3d%%] Testing odd-degree isogenies", 100 * iter / ((int)TEST_LOOPS/10));
		fflush(stdout);
		printf("\r\x1b[K");

		// Initial curve with A = 0
		fp2_set(&E.A, 0);
		fp_mont_setone(E.C.re);
		fp_set(E.C.im, 0);

		// Curve coefficient A24=(A+2C:4C)
		ec_point_t A24;
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Compute a 2^f-basis
		ec_basis_t B2;
		ec_point_t R;
		ec_curve_to_basis_2(&B2, &E);

		// Compute kernel points
		ec_point_t P, Q, PQ, R_plus, R_minus;
		fp_t k;
		fp2_tomont(&P.x, &xPA);
		fp2_tomont(&Q.x, &xQA);
		fp2_tomont(&PQ.x, &xPQA);
		fp_mont_setone(P.z.re);
		fp_set(P.z.im, 0);
		fp2_copy(&Q.z, &P.z);
		fp2_copy(&PQ.z, &P.z);
		random_scalar(k);
		ladder3pt(&R_plus, k, &P, &Q, &PQ, &A24);
		fp2_tomont(&P.x, &xPB);
		fp2_tomont(&Q.x, &xQB);
		fp2_tomont(&PQ.x, &xPQB);
		random_scalar(k);
		ladder3pt(&R_minus, k, &P, &Q, &PQ, &A24);

		// Evaluate an odd-degree isogeny
		ec_isog_odd_t isog;
		copy_point(&isog.ker_plus, &R_plus);
		copy_point(&isog.ker_minus, &R_minus);
		fp2_copy(&isog.curve.A, &E.A);
		fp2_copy(&isog.curve.C, &E.C);
		for(int i = 0; i < P_LEN; i++){
			isog.degree[i] = rand() % (TORSION_ODD_POWERS[i]+1);
			int j = isog.degree[i];
			while(j < TORSION_ODD_POWERS[i]){
				xMULv2(&isog.ker_plus, &isog.ker_plus, &TORSION_ODD_PRIMES[i], p_plus_minus_bitlength[i], &A24);
				j++;
			}
		}
		for(int i = P_LEN; i < P_LEN+M_LEN; i++){
			isog.degree[i] = rand() % (TORSION_ODD_POWERS[i]+1);
			int j = isog.degree[i];
			while(j < TORSION_ODD_POWERS[i]){
				xMULv2(&isog.ker_minus, &isog.ker_minus, &TORSION_ODD_PRIMES[i], p_plus_minus_bitlength[i], &A24);
				j++;
			}
		}

		ec_eval_odd_basis(&E, &isog, &B2, 1);
		ec_eval_odd(&E, &isog, &R_plus, 1);
		ec_eval_odd(&E, &isog, &R_minus, 1);

		// Update A24
    	fp2_add(&A24.z, &E.C, &E.C);
    	fp2_add(&A24.x, &E.A, &A24.z);
    	fp2_add(&A24.z, &A24.z, &A24.z);

		// Check order of the pushed 2^f-basis
		ec_point_t P2, Q2, PmQ2, PpQ, PpQ2;
		xADD(&PpQ, &B2.P, &B2.Q, &B2.PmQ);
		copy_point(&P2, &B2.P);
		copy_point(&Q2, &B2.Q);
		copy_point(&PmQ2, &B2.PmQ);
		copy_point(&PpQ2, &PpQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			assert(!ec_is_zero(&P2));
			xDBLv2(&Q2, &Q2, &A24);
			assert(!ec_is_zero(&Q2));
			xDBLv2(&PmQ2, &PmQ2, &A24);
			assert(!ec_is_zero(&PmQ2));
			xDBLv2(&PpQ2, &PpQ2, &A24);
			assert(!ec_is_zero(&PpQ2));
		}
		assert(is_point_equal(&PmQ2, &PpQ2));
		xDBLv2(&P2, &P2, &A24);
		assert(ec_is_zero(&P2));
		xDBLv2(&Q2, &Q2, &A24);
		assert(ec_is_zero(&Q2));
		xDBLv2(&PmQ2, &PmQ2, &A24);
		assert(ec_is_zero(&PmQ2));
		xDBLv2(&PpQ2, &PpQ2, &A24);
		assert(ec_is_zero(&PpQ2));

		// Check order of the pushed R_plus point
		int last = -1;
		for(int i = 0; i < P_LEN; i++){
			int j = isog.degree[i];
			while(j < TORSION_ODD_POWERS[i]){
				xMULv2(&R_plus, &R_plus, &TORSION_ODD_PRIMES[i], p_plus_minus_bitlength[i], &A24);
				if(ec_is_zero(&R_plus)){
					last = i;
					break;
				}
				j++;
			}
		}
		assert(last >= 0);
		for(int i = last+1; i < P_LEN; i++)
			assert(isog.degree[i] == TORSION_ODD_POWERS[i]);

		// Check order of the pushed R_minus point
		last = -1;
		for(int i = P_LEN; i < P_LEN+M_LEN; i++){
			int j = isog.degree[i];
			while(j < TORSION_ODD_POWERS[i]){
				xMULv2(&R_minus, &R_minus, &TORSION_ODD_PRIMES[i], p_plus_minus_bitlength[i], &A24);
				if(ec_is_zero(&R_minus)){
					last = i;
					break;
				}
				j++;
			}
		}
		assert(last >= 0);
		for(int i = last+1; i < P_LEN+M_LEN; i++)
			assert(isog.degree[i] == TORSION_ODD_POWERS[i]);
	}
	printf("[%2d%%] Tested odd-degree isogenies:\tNo errors!\n", 100);

	// ----------------- TEST FOR ISOMORPHISMS ----------------- //

	// Initial curve with A = 0
	fp2_set(&E.A, 0);
	fp_mont_setone(E.C.re);
	fp_set(E.C.im, 0);

	for(int iter = 0; iter < TEST_LOOPS; iter++){

		printf("[%3d%%] Testing odd-degree isogenies", 100 * iter / (int)TEST_LOOPS);
		fflush(stdout);
		printf("\r\x1b[K");

		// Compute jinv
		fp2_t j;
		ec_j_inv(&j, &E);

		// Normalize the curve
		ec_curve_t Enorm;
		ec_isom_t isom;
		ec_curve_normalize(&Enorm, &isom, &E);

		// Compare j invariants
		fp2_t jnorm;
		ec_j_inv(&jnorm, &Enorm);
		assert(fp2_is_equal(&j, &jnorm));

		// Check that normalization is consistent
		ec_curve_t Enorm2;
		ec_isom_t isom2;
		ec_curve_normalize(&Enorm2, &isom2, &Enorm);
		assert(curve_equal(&Enorm, &Enorm2));

		// Compute a basis
		ec_basis_t B;
		ec_curve_to_basis_2(&B, &E);

		// Push it through isomorphism
		ec_basis_t Bnorm;
		copy_point(&Bnorm.P, &B.P);
		ec_iso_eval(&Bnorm.P, &isom);
		copy_point(&Bnorm.Q, &B.Q);
		ec_iso_eval(&Bnorm.Q, &isom);
		copy_point(&Bnorm.PmQ, &B.PmQ);
		ec_iso_eval(&Bnorm.PmQ, &isom);

		// Check that the new basis is rational
		assert(ec_is_on_curve(&Enorm, &Bnorm.P));
		assert(ec_is_on_curve(&Enorm, &Bnorm.Q));
		assert(ec_is_on_curve(&Enorm, &Bnorm.PmQ));

		// Compute P+Q
		ec_point_t PpQ;
		xADD(&PpQ, &Bnorm.P, &Bnorm.Q, &Bnorm.PmQ);

		// // Check the order
		ec_point_t P2, Q2, PmQ2, PpQ2, R, A24;
		fp2_add(&A24.z, &Enorm.C, &Enorm.C);
		fp2_add(&A24.x, &Enorm.A, &A24.z);
		fp2_add(&A24.z, &A24.z, &A24.z);
		copy_point(&P2, &Bnorm.P);
		copy_point(&Q2, &Bnorm.Q);
		copy_point(&PmQ2, &Bnorm.PmQ);
		copy_point(&PpQ2, &PpQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			assert(!ec_is_zero(&P2));
			xDBLv2(&Q2, &Q2, &A24);
			assert(!ec_is_zero(&Q2));
			xDBLv2(&PmQ2, &PmQ2, &A24);
			assert(!ec_is_zero(&PmQ2));
			xDBLv2(&PpQ2, &PpQ2, &A24);
			assert(!ec_is_zero(&PpQ2));
		}
		assert(is_point_equal(&PmQ2, &PpQ2));
		xDBLv2(&P2, &P2, &A24);
		assert(ec_is_zero(&P2));
		xDBLv2(&Q2, &Q2, &A24);
		assert(ec_is_zero(&Q2));
		xDBLv2(&PmQ2, &PmQ2, &A24);
		assert(ec_is_zero(&PmQ2));
		xDBLv2(&PpQ2, &PpQ2, &A24);
		assert(ec_is_zero(&PpQ2));

		// Compute a 2^e-basis with 2^(e-1)*Q=(0,0)
		fp2_add(&A24.z, &E.C, &E.C);
		fp2_add(&A24.x, &E.A, &A24.z);
		fp2_add(&A24.z, &A24.z, &A24.z);
		copy_point(&P2, &B.P);
		copy_point(&Q2, &B.Q);
		copy_point(&PmQ2, &B.PmQ);
		for(int i = 0; i < POWER_OF_2 - 1; i++){
			xDBLv2(&P2, &P2, &A24);
			xDBLv2(&Q2, &Q2, &A24);
			xDBLv2(&PmQ2, &PmQ2, &A24);
		}
		if(fp2_is_zero(&P2.x)){
			copy_point(&Q2, &B.Q);
			copy_point(&B.Q, &B.P);
			copy_point(&B.P, &Q2);
		}
		else if(fp2_is_zero(&PmQ2.x)){
			copy_point(&Q2, &B.Q);
			copy_point(&B.Q, &B.PmQ);
			copy_point(&B.PmQ, &Q2);
		}

		// Compute a 2^e-isogeny
		ec_isog_even_t isog;
		fp_t k;
		random_scalar(k);
		ladder3pt(&R, k, &B.P, &B.Q, &B.PmQ, &A24);
		fp2_copy(&isog.curve.A, &E.A);
		fp2_copy(&isog.curve.C, &E.C);
		copy_point(&isog.kernel, &R);
		isog.length = POWER_OF_2;
		ec_eval_even_nonzero(&E, &isog, &P2, 1);
	}
	printf("[%2d%%] Tested curve normalization+isom:\tNo errors!\n", 100);

	printf("-- All tests passed!\n");
	return 0;
}
