#include "isog.h"


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

int ec_is_on_curve(const ec_curve_t* curve, const ec_point_t* P){

    fp2_t t0, t1, t2;

    // Check if xz*(C^2x^2+zACx+z^2C^2) is a square
    fp2_mul(&t0, &curve->C, &P->x); 
    fp2_mul(&t1, &t0, &P->z);       
    fp2_mul(&t1, &t1, &curve->A);   
    fp2_mul(&t2, &curve->C, &P->z); 
    fp2_sqr(&t0, &t0);              
    fp2_sqr(&t2, &t2);              
    fp2_add(&t0, &t0, &t1);
    fp2_add(&t0, &t0, &t2);
    fp2_mul(&t0, &t0, &P->x);
    fp2_mul(&t0, &t0, &P->z);
    return fp2_is_square(&t0);
}

static void difference_point(ec_point_t* PQ, const ec_point_t* P, const ec_point_t* Q, const ec_curve_t* curve){
    // Given P,Q in affine x-only, computes a deterministic choice for (P-Q)
    // The points must be normalized to z=1 and the curve to C=1

    fp2_t t0, t1, t2, t3;
    
    fp2_sub(&PQ->z, &P->x, &Q->x);  // P - Q
    fp2_mul(&t2, &P->x, &Q->x);     // P*Q
    fp_mont_setone(t1.re);
    fp_set(t1.im, 0);
    fp2_sub(&t3, &t2, &t1);         // P*Q-1
    fp2_mul(&t0, &PQ->z, &t3);      // (P-Q)*(P*Q-1)
    fp2_sqr(&PQ->z, &PQ->z);        // (P-Q)^2
    fp2_sqr(&t0, &t0);              // (P-Q)^2*(P*Q-1)^2
    fp2_add(&t1, &t2, &t1);         // P*Q+1
    fp2_add(&t3, &P->x, &Q->x);     // P+Q
    fp2_mul(&t1, &t1, &t3);         // (P+Q)*(P*Q+1)
    fp2_mul(&t2, &t2, &curve->A);   // A*P*Q
    fp2_add(&t2, &t2, &t2);         // 2*A*P*Q
    fp2_add(&t1, &t1, &t2);         // (P+Q)*(P*Q+1) + 2*A*P*Q
    fp2_sqr(&t2, &t1);              // ((P+Q)*(P*Q+1) + 2*A*P*Q)^2
    fp2_sub(&t0, &t2, &t0);         // ((P+Q)*(P*Q+1) + 2*A*P*Q)^2 - (P-Q)^2*(P*Q-1)^2
    fp2_sqrt(&t0);
    fp2_add(&PQ->x, &t0, &t1);
}

void ec_curve_to_basis_2(ec_basis_t *PQ2, const ec_curve_t *curve){
    fp2_t x, t0, t1, t2;
    ec_point_t P, Q, Q2, P2, A24;

    // Curve coefficient in the form A24 = (A+2C:4C)
    fp2_add(&A24.z, &curve->C, &curve->C);
    fp2_add(&A24.x, &curve->A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);

    fp_mont_setone(x.re);
    fp_set(x.im, 0);

    // Find P
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&P.x, &x);
            fp_mont_setone(P.z.re);
            fp_set(P.z.im, 0);
        }
        else
            continue;

        // Clear odd factors from the order
        xMULv2(&P, &P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, &A24);

        // Check if point has order 2^f
        copy_point(&P2, &P);
        for(int i = 0; i < POWER_OF_2 - 1; i++)
            xDBLv2(&P2, &P2, &A24);
        if(ec_is_zero(&P2))
            continue;
        else
            break;
    }
    
    // Find Q
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&Q.x, &x);
            fp_mont_setone(Q.z.re);
            fp_set(Q.z.im, 0);
        }
        else
            continue;

        // Clear odd factors from the order
        xMULv2(&Q, &Q, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, &A24);

        // Check if point has order 2^f
        copy_point(&Q2, &Q);
        for(int i = 0; i < POWER_OF_2 - 1; i++)
            xDBLv2(&Q2, &Q2, &A24);
        if(ec_is_zero(&Q2))
            continue;

        // Check if point is orthogonal to P
        if(is_point_equal(&P2, &Q2))
            continue;
        else
            break;
    }

    // Normalize points
    ec_curve_t E;
    fp2_mul(&t0, &P.z, &Q.z);
    fp2_mul(&t1, &t0, &curve->C);
    fp2_inv(&t1);
    fp2_mul(&P.x, &P.x, &t1);
    fp2_mul(&Q.x, &Q.x, &t1);
    fp2_mul(&E.A, &curve->A, &t1);
    fp2_mul(&P.x, &P.x, &Q.z);
    fp2_mul(&P.x, &P.x, &curve->C);
    fp2_mul(&Q.x, &Q.x, &P.z);
    fp2_mul(&Q.x, &Q.x, &curve->C);
    fp2_mul(&E.A, &E.A, &t0);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    fp2_copy(&Q.z, &P.z);
    fp2_copy(&E.C, &P.z);

    // Compute P-Q
    difference_point(&PQ2->PmQ, &P, &Q, &E);
    copy_point(&PQ2->P, &P);
    copy_point(&PQ2->Q, &Q);
}


void ec_complete_basis_2(ec_basis_t* PQ2, const ec_curve_t* curve, const ec_point_t* P){

    fp2_t x, t0, t1, t2;
    ec_point_t Q, Q2, P2, A24;

    // Curve coefficient in the form A24 = (A+2C:4C)
    fp2_add(&A24.z, &curve->C, &curve->C);
    fp2_add(&A24.x, &curve->A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);

    // Point of order 2 generated by P
    copy_point(&P2, P);
    for(int i = 0; i < POWER_OF_2 - 1; i++)
        xDBLv2(&P2, &P2, &A24);

    // Find Q
    fp_mont_setone(x.re);
    fp_set(x.im, 0);
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&Q.x, &x);
            fp_mont_setone(Q.z.re);
            fp_set(Q.z.im, 0);
        }
        else
            continue;

        // Clear odd factors from the order
        xMULv2(&Q, &Q, p_cofactor_for_2f, (int)P_COFACTOR_FOR_2F_BITLENGTH, &A24);

        // Check if point has order 2^f
        copy_point(&Q2, &Q);
        for(int i = 0; i < POWER_OF_2 - 1; i++)
            xDBLv2(&Q2, &Q2, &A24);
        if(ec_is_zero(&Q2))
            continue;

        // Check if point is orthogonal to P
        if(is_point_equal(&P2, &Q2))
            continue;
        else
            break;
    }

    // Normalize points
    ec_curve_t E;
    ec_point_t PP;
    fp2_mul(&t0, &P->z, &Q.z);
    fp2_mul(&t1, &t0, &curve->C);
    fp2_inv(&t1);
    fp2_mul(&PP.x, &P->x, &t1);
    fp2_mul(&Q.x, &Q.x, &t1);
    fp2_mul(&E.A, &curve->A, &t1);
    fp2_mul(&PP.x, &PP.x, &Q.z);
    fp2_mul(&PP.x, &PP.x, &curve->C);
    fp2_mul(&Q.x, &Q.x, &P->z);
    fp2_mul(&Q.x, &Q.x, &curve->C);
    fp2_mul(&E.A, &E.A, &t0);
    fp_mont_setone(PP.z.re);
    fp_set(PP.z.im, 0);
    fp2_copy(&Q.z, &PP.z);
    fp2_copy(&E.C, &PP.z);

    // Compute P-Q
    difference_point(&PQ2->PmQ, &PP, &Q, &E);
    copy_point(&PQ2->P, &PP);
    copy_point(&PQ2->Q, &Q);
}

void ec_curve_to_basis_3(ec_basis_t* PQ3, const ec_curve_t* curve){

    fp2_t x, t0, t1, t2;
    ec_point_t P, Q, Q3, P3, A24, A3;

    // Curve coefficient in the form A24 = (A+2C:4C)
    fp2_add(&A24.z, &curve->C, &curve->C);
    fp2_add(&A24.x, &curve->A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);

    // Curve coefficient in the form A3 = (A+2C:A-2C)
    fp2_sub(&A3.z, &A24.x, &A24.z);
    fp2_copy(&A3.x, &A24.x);

    fp_mont_setone(x.re);
    fp_set(x.im, 0);

    // Find P
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&P.x, &x);
            fp_mont_setone(P.z.re);
            fp_set(P.z.im, 0);
        }
        else
            continue;

        // Clear non-3 factors from the order
        xMULv2(&P, &P, p_cofactor_for_3g, (int)P_COFACTOR_FOR_3G_BITLENGTH, &A24);

        // Check if point has order 3^g
        copy_point(&P3, &P);
        for(int i = 0; i < POWER_OF_3 - 1; i++)
            xTPL(&P3, &P3, &A3);
        if(ec_is_zero(&P3))
            continue;
        else
            break;
    }
    
    // Find Q
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&Q.x, &x);
            fp_mont_setone(Q.z.re);
            fp_set(Q.z.im, 0);
        }
        else
            continue;

        // Clear non-3 factors from the order
        xMULv2(&Q, &Q, p_cofactor_for_3g, (int)P_COFACTOR_FOR_3G_BITLENGTH, &A24);

        // Check if point has order 3^g
        copy_point(&Q3, &Q);
        for(int i = 0; i < POWER_OF_3 - 1; i++)
            xTPL(&Q3, &Q3, &A3);
        if(ec_is_zero(&Q3))
            continue;

        // Check if point is orthogonal to P
        if(is_point_equal(&P3, &Q3))
            continue;
        xDBLv2(&P3, &P3, &A24);
        if(is_point_equal(&P3, &Q3))
            continue;
        else
            break;
    }

    // Normalize points
    ec_curve_t E;
    fp2_mul(&t0, &P.z, &Q.z);
    fp2_mul(&t1, &t0, &curve->C);
    fp2_inv(&t1);
    fp2_mul(&P.x, &P.x, &t1);
    fp2_mul(&Q.x, &Q.x, &t1);
    fp2_mul(&E.A, &curve->A, &t1);
    fp2_mul(&P.x, &P.x, &Q.z);
    fp2_mul(&P.x, &P.x, &curve->C);
    fp2_mul(&Q.x, &Q.x, &P.z);
    fp2_mul(&Q.x, &Q.x, &curve->C);
    fp2_mul(&E.A, &E.A, &t0);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    fp2_copy(&Q.z, &P.z);
    fp2_copy(&E.C, &P.z);

    // Compute P-Q
    difference_point(&PQ3->PmQ, &P, &Q, &E);
    copy_point(&PQ3->P, &P);
    copy_point(&PQ3->Q, &Q);
}

void ec_curve_to_basis_6(ec_basis_t* PQ6, const ec_curve_t* curve){

    fp2_t x, t0, t1, t2;
    ec_point_t P, Q, Q6, P6, R, T, A24, A3;

    // Curve coefficient in the form A24 = (A+2C:4C)
    fp2_add(&A24.z, &curve->C, &curve->C);
    fp2_add(&A24.x, &curve->A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);

    // Curve coefficient in the form A3 = (A+2C:A-2C)
    fp2_sub(&A3.z, &A24.x, &A24.z);
    fp2_copy(&A3.x, &A24.x);

    fp_mont_setone(x.re);
    fp_set(x.im, 0);

    // Find P
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&P.x, &x);
            fp_mont_setone(P.z.re);
            fp_set(P.z.im, 0);
        }
        else
            continue;

        // Clear non-2 factors and non-3 factors from the order
        xMULv2(&P, &P, p_cofactor_for_6fg, (int)P_COFACTOR_FOR_6FG_BITLENGTH, &A24);

        // Check if point has order 2^f*3^g
        copy_point(&P6, &P);
        for(int i = 0; i < POWER_OF_2 - 1; i++)
            xDBLv2(&P6, &P6, &A24);
        for(int i = 0; i < POWER_OF_3 - 1; i++)
            xTPL(&P6, &P6, &A3);
        if(ec_is_zero(&P6))
            continue;
        xDBLv2(&T, &P6, &A24);
        if (ec_is_zero(&T))
            continue;
        xTPL(&T, &P6, &A3);
        if (ec_is_zero(&T))
            continue;
        break;
    }

    // Find Q
    while(1){
        fp_add(x.im, x.re, x.im);

        // Check if point is rational
        fp2_sqr(&t0, &curve->C);
        fp2_mul(&t1, &t0, &x);
        fp2_mul(&t2, &curve->A, &curve->C);
        fp2_add(&t1, &t1, &t2);
        fp2_mul(&t1, &t1, &x);
        fp2_add(&t1, &t1, &t0);
        fp2_mul(&t1, &t1, &x);
        if(fp2_is_square(&t1)){
            fp2_copy(&Q.x, &x);
            fp_mont_setone(Q.z.re);
            fp_set(Q.z.im, 0);
        }
        else
            continue;

        // Clear non-6 factors from the order
        xMULv2(&Q, &Q, p_cofactor_for_6fg, (int)P_COFACTOR_FOR_6FG_BITLENGTH, &A24);

        // Check first if point has order 2^f*3^g
        copy_point(&Q6, &Q);
        for(int i = 0; i < POWER_OF_2 - 1; i++)
            xDBLv2(&Q6, &Q6, &A24);
        for(int i = 0; i < POWER_OF_3 - 1; i++)
            xTPL(&Q6, &Q6, &A3);
        if(ec_is_zero(&Q6))
            continue;
        xDBLv2(&T, &Q6, &A24);
        if (ec_is_zero(&T))
            continue;
        xTPL(&T, &Q6, &A3);
        if (ec_is_zero(&T))
            continue;

        // Check if point P is independent from point Q
        xTPL(&R, &P6, &A3);
        xTPL(&T, &Q6, &A3);
        if(is_point_equal(&R, &T))
            continue;
        xDBLv2(&R, &P6, &A24);
        xDBLv2(&T, &Q6, &A24);
        if(is_point_equal(&R, &T))
            continue;
        break;
    }

    // Normalize points
    ec_curve_t E;
    fp2_mul(&t0, &P.z, &Q.z);
    fp2_mul(&t1, &t0, &curve->C);
    fp2_inv(&t1);
    fp2_mul(&P.x, &P.x, &t1);
    fp2_mul(&Q.x, &Q.x, &t1);
    fp2_mul(&E.A, &curve->A, &t1);
    fp2_mul(&P.x, &P.x, &Q.z);
    fp2_mul(&P.x, &P.x, &curve->C);
    fp2_mul(&Q.x, &Q.x, &P.z);
    fp2_mul(&Q.x, &Q.x, &curve->C);
    fp2_mul(&E.A, &E.A, &t0);
    fp_mont_setone(P.z.re);
    fp_set(P.z.im, 0);
    fp2_copy(&Q.z, &P.z);
    fp2_copy(&E.C, &P.z);

    // Compute P-Q
    difference_point(&PQ6->PmQ, &P, &Q, &E);
    copy_point(&PQ6->P, &P);
    copy_point(&PQ6->Q, &Q);
}
