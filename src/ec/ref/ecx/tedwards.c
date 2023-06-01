#include <tedwards.h>
#include <assert.h>

// a*x^2+y^2=1+d*x^2*y^2
// a = A.x/A.z + 2, d = A.x/A.z - 2

void ted_init(ted_point_t* P)
{ // Initialize point as identity element (X:Y:Z:T) <- (0:1:1:0)
    fp_t one = {0};

    memset((digit_t*)P, 0, NWORDS_FIELD*RADIX*8/8);
    one[0] = 1;
    fp_tomont(P->x.re, one);
}

void copy_ted_point(ted_point_t* P, ted_point_t const* Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->y), &(Q->y));
    fp2_copy(&(P->z), &(Q->z));
    fp2_copy(&(P->t), &(Q->t));
}

void ted_dbl(ted_point_t *Q, ted_point_t const *P, ec_curve_t const* E) 
{
    // A = X1^2
    // B = Y1^2
    // C = 2*Z1^2
    // D = a*A
    // K = (X1+Y1)^2-A-B
    // G = D+B
    // F = G-C
    // H = D-B
    // X3 = K*F
    // Y3 = G*H
    // T3 = K*H
    // Z3 = F*G

    // TODO: neutral element
    fp2_t A, B, C, D, K, G, F, H;

    fp2_sqr(&A, &P->x);
    fp2_sqr(&B, &P->y);
    fp2_sqr(&C, &P->z);
    fp2_add(&C, &C, &C);
    fp2_mul(&D, &A, &E->A);
    fp2_add(&K, &P->x, &P->y);
    fp2_sqr(&K, &K);
    fp2_sub(&K, &K, &A);
    fp2_sub(&K, &K, &B);
    fp2_add(&G, &D, &B);
    fp2_sub(&F, &G, &C);
    fp2_sub(&H, &D, &B);
    fp2_mul(&Q->x, &K, &F);
    fp2_mul(&Q->y, &G, &H);
    fp2_mul(&Q->t, &K, &H);
    fp2_mul(&Q->z, &F, &G);
}

void ted_add(ted_point_t* S, ted_point_t const* P, ted_point_t const* Q, ec_curve_t const* E)
{
    // A = X1*X2
    // B = Y1*Y2
    // C = Z1*T2
    // D = T1*Z2
    // K = D+C
    // F = (X1-Y1)*(X2+Y2)+B-A
    // G = B+a*A
    // H = D-C
    // X3 = K*F
    // Y3 = G*H
    // T3 = K*H
    // Z3 = F*G

    // TODO: neutral element

    ted_point_t res;

    if (is_ted_equal(P, Q)) {
      ted_dbl(S, P, E);
      return;
    }
    //assert(!is_ted_equal(P, Q));
    
    ted_neg(&res, P);
    if (is_ted_equal(&res, Q)) {
       ted_init(S);
       return;
    }
    // assert(!ted_equal(&res,Q));
    fp2_t A, B, C, D, K, F, G, H, tmp;

    fp2_mul(&A, &P->x, &Q->x);
    fp2_mul(&B, &P->y, &Q->y);
    fp2_mul(&C, &P->z, &Q->t);
    fp2_mul(&D, &P->t, &Q->z);
    fp2_add(&K, &D, &C);
    fp2_add(&F, &Q->x, &Q->y);
    fp2_sub(&tmp, &P->x, &P->y);
    fp2_mul(&F, &F, &tmp);
    fp2_add(&F, &F, &B);
    fp2_sub(&F, &F, &A);
    fp2_mul(&G, &A, &E->A);
    fp2_add(&G, &G, &B);
    fp2_sub(&H, &D, &C);
    fp2_mul(&res.x, &K, &F);
    fp2_mul(&res.y, &G, &H);
    fp2_mul(&res.t, &K, &H);
    fp2_mul(&res.z, &F, &G);

    if (fp2_is_zero(&res.x) && fp2_is_zero(&res.y) && fp2_is_zero(&res.z)) {
        ted_dbl(S, P, E);
    } else {
        copy_ted_point(S, &res);
    }
}

void ted_neg(ted_point_t* Q, ted_point_t const* P)
{
    fp2_neg(&Q->x, &P->x);
    fp2_copy(&Q->y, &P->y);
    fp2_copy(&Q->z, &P->z);
    fp2_neg(&Q->t, &P->t);
}

static bool xLIFT(fp2_t* y, const ec_point_t* P, const ec_curve_t* curve)
{ // Returns false if it is on the curve, true if it is on the twist
    fp2_t z2, tmp1, tmp2, y2;

    if (fp2_is_zero(&P->z)) return false;

    // (X^2 + Z^2) C
    fp2_sqr(&tmp1, &P->x);
    fp2_sqr(&z2, &P->z);
    fp2_add(&tmp1, &tmp1, &z2);
    fp2_mul(&tmp1, &tmp1, &curve->C);

    // X^2C + AXZ + Z^2C
    fp2_mul(&tmp2, &P->x, &P->z);
    fp2_mul(&tmp2, &tmp2, &curve->A);
    fp2_add(&tmp1, &tmp1, &tmp2);

    // X^3C + AX^2Z + XZ^2C = Z^3(Cx^3 + Ax^2 + Cx) = Z^3 C (B*y^2) = Z C (B*Y^2) // x = X/Z
    fp2_mul(&tmp1, &tmp1, &P->x);
    // (ZC)^(-1)
    fp2_mul(&tmp2, &curve->C, &P->z);

    assert(!fp2_is_zero(&tmp2));
    
    fp2_inv(&tmp2);    
    fp2_mul(&y2, &tmp1, &tmp2);    // (B*Y^2)
    fp2_copy(y, &y2);

    if (fp2_is_square(&y2)) {  // on the curve
        fp2_sqrt(y);
        return false;
    } else { // on the twist
        fp2_t tmp = fp2_non_residue();
        fp2_mul(y, y, &tmp);
        fp2_sqrt(y);
        return true;
    }
}

//void mont_to_ted(ec_point_t* E, ec_point_t const* A, bool twist)
void mont_to_ted(ec_curve_t* ted_curve, ec_curve_t const* curve)
{
    fp2_t tmp, two;

    // A : y^2 = x^3 + (a/c)x^2 + x
    fp2_copy(&tmp, &curve->C);         
    fp2_inv(&tmp);                    // 1/c
    fp2_mul(&tmp, &tmp, &curve->A);   // a/c
    fp2_set(&two, 2);
    fp2_tomont(&two, &two);
    fp2_add(&ted_curve->A, &tmp, &two);       // a/c + 2
    fp2_sub(&ted_curve->C, &tmp, &two);       // a/c - 2
    //if (twist) {
        // B = Fp2_inv(fp2_non_residue)
    //    tmp = fp2_non_residue();
    //    fp2_mul2(&E->x,&tmp);
    //    fp2_mul2(&E->z,&tmp);
    //}
}

void mont_to_ted_point(ted_point_t* Q, ec_point_t const* P, ec_curve_t const* curve)
{
    if (fp2_is_zero(&P->z)) {
        fp2_set(&Q->x, 0);
        fp2_set(&Q->y, 1);
        fp2_set(&Q->z, 1);
        fp2_set(&Q->t, 0);
        fp_tomont(Q->y.re, Q->y.re);
        fp_tomont(Q->z.re, Q->z.re);
    } else {
        fp2_t tmp, y;

        xLIFT(&y, P, curve);
        fp2_add(&tmp, &P->x, &P->z);
        fp2_mul(&Q->x, &P->x, &tmp);
        fp2_sub(&Q->y, &P->x, &P->z);
        fp2_mul(&Q->y, &Q->y, &y);
        fp2_mul(&Q->z, &tmp, &y);
        fp2_copy(&Q->t, &Q->z);
        fp2_inv(&Q->t);
        fp2_mul(&Q->t, &Q->t, &Q->x);
        fp2_mul(&Q->t, &Q->t, &Q->y);
    }
}

void ted_to_mont_point(ec_point_t* Q, ted_point_t const* P)
{
    fp2_add(&Q->x, &P->z, &P->y);
    fp2_sub(&Q->z, &P->z, &P->y);
}

bool is_ted_equal(ted_point_t const* P1, ted_point_t const* P2)
{
    fp2_t x1z2, y1z2;
    fp2_t y2z1, x2z1;
    fp2_t t1y2, t2y1;

    fp2_mul(&x1z2, &P1->x, &P2->z);
    fp2_mul(&y1z2, &P1->y, &P2->z);
    fp2_mul(&y2z1, &P2->y, &P1->z);
    fp2_mul(&x2z1, &P2->x, &P1->z);
    fp2_mul(&t1y2, &P1->t, &P2->y);
    fp2_mul(&t2y1, &P2->t, &P1->y);

    return fp2_is_equal(&x1z2, &x2z1) && fp2_is_equal(&y1z2, &y2z1) && fp2_is_equal(&t1y2, &t2y1);
}