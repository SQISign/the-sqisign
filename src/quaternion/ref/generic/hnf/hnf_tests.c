
#include "hnf_internal.h"
#include "quaternion_tests.h"
#include <rng.h>

// test helper for xgcd_not_0

// round a/b to closest integer
void
ibz_rounded_div(ibz_t *q, const ibz_t *a, const ibz_t *b)
{
    ibz_t r, sign_q, abs_b;
    ibz_init(&r);
    ibz_init(&sign_q);
    ibz_init(&abs_b);

    // assumed to round towards 0
    ibz_abs(&abs_b, b);
    // q is of same sign as a*b (and 0 if a is 0)
    ibz_mul(&sign_q, a, b);
    ibz_div(q, &r, a, b);
    ibz_abs(&r, &r);
    ibz_add(&r, &r, &r);
    ibz_set(&sign_q, (1 - 2 * (ibz_cmp(&sign_q, &ibz_const_zero) < 0)) * (ibz_cmp(&r, &abs_b) > 0));
    ibz_add(q, q, &sign_q);
    ibz_finalize(&r);
    ibz_finalize(&sign_q);
    ibz_finalize(&abs_b);
}

// old HNF version, used only in tests of the new HNF
// Algorithm used is the one at number 2.4.5 in Henri Cohen's "A Course in Computational Algebraic
// Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993
//  assumes ibz_xgcd outputs u,v which are small in absolute value (as described in the
//  book)
void
ibz_mat_4xn_hnf_core(ibz_mat_4x4_t *hnf, int generator_number, const ibz_vec_4_t *generators)
{
    assert(generator_number > 3);
    int n = generator_number;
    int i = 3;
    int j = n - 1;
    int k = n - 1;
    ibz_t b, u, v, d, zero, coeff_1, coeff_2, r;
    ibz_vec_4_t c;
    ibz_vec_4_t a[n];
    ibz_init(&b);
    ibz_init(&d);
    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&r);
    ibz_init(&coeff_1);
    ibz_init(&coeff_2);
    ibz_init(&zero);
    ibz_set(&zero, 0);
    ibz_vec_4_init(&c);
    for (int h = 0; h < n; h++) {
        ibz_vec_4_init(&(a[h]));
        ibz_copy(&(a[h][0]), &(generators[h][0]));
        ibz_copy(&(a[h][1]), &(generators[h][1]));
        ibz_copy(&(a[h][2]), &(generators[h][2]));
        ibz_copy(&(a[h][3]), &(generators[h][3]));
    }
    while (i != -1) {
        while (j != 0) {
            j = j - 1;
            if (!ibz_is_zero(&(a[j][i]))) {
                // assumtion that ibz_xgcd outputs u,v which are small in absolute
                // value is needed here
                ibz_xgcd(&d, &u, &v, &(a[k][i]), &(a[j][i]));
                // also, needs u non 0, but v can be 0 if needed
                if (ibz_is_zero(&u)) {
                    ibz_div(&v, &r, &(a[k][i]), &(a[j][i]));
                    ibz_set(&u, 1);
                    ibz_sub(&v, &u, &v);
                }
                ibz_vec_4_linear_combination(&c, &u, &(a[k]), &v, &(a[j]));
                ibz_div(&coeff_1, &r, &(a[k][i]), &d);
                ibz_div(&coeff_2, &r, &(a[j][i]), &d);
                ibz_neg(&coeff_2, &coeff_2);
                ibz_vec_4_linear_combination(&(a[j]), &coeff_1, &(a[j]), &coeff_2, &(a[k]));
                ibz_vec_4_copy(&(a[k]), &c);
            }
        }
        ibz_copy(&b, &(a[k][i]));
        if (ibz_cmp(&b, &zero) < 0) {
            ibz_vec_4_negate(&(a[k]), &(a[k]));
            ibz_neg(&b, &b);
        }
        if (ibz_is_zero(&b)) {
            k = k + 1;
        } else {
            for (j = k + 1; j < n; j++) {
                ibz_div(&d, &r, &(a[j][i]), &b);
                if (ibz_cmp(&r, &zero) < 0) {
                    ibz_set(&r, 1);
                    ibz_sub(&d, &d, &r);
                }
                ibz_set(&r, 1);
                ibz_neg(&d, &d);
                ibz_vec_4_linear_combination(&(a[j]), &r, &(a[j]), &d, &(a[k]));
            }
        }
        if (i != 0) {
            k = k - 1;
            j = k;
        }
        i = i - 1;
    }
    for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) {
            ibz_copy(&((*hnf)[i][j]), &(a[n - 4 + j][i]));
        }
    }

    ibz_finalize(&b);
    ibz_finalize(&d);
    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&r);
    ibz_finalize(&coeff_1);
    ibz_finalize(&coeff_2);
    ibz_finalize(&zero);
    ibz_vec_4_finalize(&c);
    for (int h = 0; h < n; h++) {
        ibz_vec_4_finalize(&(a[h]));
    }
}

// integer hnf helpers

// void ibz_mod_not_zero(ibz_t *res, const ibz_t *x, const ibz_t *mod);
int
quat_test_ibz_mod_not_zero(void)
{
    int res = 0;
    ibz_t m, a, r, d;
    int s, x;
    ibz_init(&m);
    ibz_init(&a);
    ibz_init(&r);
    ibz_init(&d);
    s = 5;
    ibz_set(&m, s);
    for (x = -20; x < 20; x++) {
        ibz_set(&a, x);
        ibz_mod_not_zero(&r, &a, &m);
        ibz_sub(&d, &r, &a);
        res = res || !ibz_divides(&d, &m);
        res = res || !(ibz_cmp(&r, &m) <= 0);
        res = res || !(ibz_cmp(&ibz_const_one, &r) <= 0);
        ibz_mod_not_zero(&a, &a, &m);
        res = res || !(ibz_cmp(&a, &r) == 0);
        res = res || !(ibz_cmp(&r, &m) <= 0);
        res = res || !(ibz_cmp(&ibz_const_one, &r) <= 0);
    }

    if (res != 0) {
        printf("Quaternion unit test ibz_mod_not_zero failed\n");
    }
    ibz_finalize(&m);
    ibz_finalize(&a);
    ibz_finalize(&r);
    ibz_finalize(&d);
    return (res);
}

// void ibz_centered_mod(ibz_t *remainder, const ibz_t *a, const ibz_t *mod);
int
quat_test_ibz_centered_mod(void)
{
    int res = 0;
    ibz_t m, a, r, d, h;
    int s, x;
    ibz_init(&m);
    ibz_init(&a);
    ibz_init(&r);
    ibz_init(&d);
    ibz_init(&h);
    s = 5;
    ibz_set(&m, s);
    for (x = -20; x < 20; x++) {
        ibz_set(&a, x);
        ibz_centered_mod(&r, &a, &m);
        ibz_sub(&d, &r, &a);
        res = res || !ibz_divides(&d, &m);
        ibz_mul(&h, &r, &ibz_const_two);
        res = res || !(ibz_cmp(&h, &m) <= 0);
        ibz_neg(&h, &h);
        res = res || !(ibz_cmp(&h, &m) < 0);
        ibz_centered_mod(&a, &a, &m);
        res = res || !(ibz_cmp(&a, &r) == 0);
    }

    if (res != 0) {
        printf("Quaternion unit test ibz_centered_mod failed\n");
    }
    ibz_finalize(&m);
    ibz_finalize(&a);
    ibz_finalize(&r);
    ibz_finalize(&d);
    ibz_finalize(&h);
    return (res);
}

// void ibz_conditional_assign(ibz_t *res, const ibz_t *x, const ibz_t *y, int c);
int
quat_test_ibz_conditional_assign(void)
{
    int res = 0;
    ibz_t x, y, r;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&r);
    ibz_set(&x, 5);
    ibz_set(&y, -9);
    ibz_conditional_assign(&r, &x, &y, 1);
    res = res || !(ibz_cmp(&x, &r) == 0);
    ibz_conditional_assign(&r, &x, &y, 0);
    res = res || !(ibz_cmp(&y, &r) == 0);
    ibz_conditional_assign(&x, &x, &y, 0);
    res = res || !(ibz_cmp(&y, &x) == 0);
    ibz_set(&x, -5);
    ibz_set(&y, -0);
    ibz_conditional_assign(&y, &x, &y, 1);
    res = res || !(ibz_cmp(&y, &x) == 0);

    if (res != 0) {
        printf("Quaternion unit test ibz_conditional_assign failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&r);
    return (res);
}

// This tests that our xgcd_with_u_non_0 respects all criteria required by our HNF
int
quat_test_helper_ibz_xgcd_with_u_not_0(const ibz_t *gcd, const ibz_t *u, const ibz_t *v, const ibz_t *x, const ibz_t *y)
{
    int res = 0;
    ibz_t sum, prod, test, cmp;
    ibz_init(&sum);
    ibz_init(&prod);
    ibz_init(&cmp);
    ibz_init(&test);
    // sign correct
    res = res | !(ibz_cmp(gcd, &ibz_const_zero) >= 0);
    if (ibz_is_zero(x) && ibz_is_zero(y)) {
        res = res | !(ibz_is_zero(v) && ibz_is_one(u) && ibz_is_one(gcd));
    } else {
        if (!ibz_is_zero(x) && !ibz_is_zero(y)) {
            // GCD divides x
            ibz_div(&sum, &prod, x, gcd);
            res = res | !ibz_is_zero(&prod);
            // Small enough
            ibz_mul(&prod, x, u);
            res = res | !(ibz_cmp(&prod, &ibz_const_zero) > 0);
            ibz_mul(&sum, &sum, y);
            ibz_abs(&sum, &sum);
            res = res | !(ibz_cmp(&prod, &sum) <= 0);

            // GCD divides y
            ibz_div(&sum, &prod, y, gcd);
            res = res | !ibz_is_zero(&prod);
            // Small enough
            ibz_mul(&prod, y, v);
            res = res | !(ibz_cmp(&prod, &ibz_const_zero) <= 0);
            ibz_mul(&sum, &sum, x);
            ibz_abs(&sum, &sum);
            res = res | !(ibz_cmp(&prod, &sum) < 0);
        } else {
            // GCD divides x
            ibz_div(&sum, &prod, x, gcd);
            res = res | !ibz_is_zero(&prod);
            // GCD divides y
            ibz_div(&sum, &prod, y, gcd);
            res = res | !ibz_is_zero(&prod);
            if (ibz_is_zero(x) && !ibz_is_zero(y)) {
                ibz_abs(&prod, v);
                res = res | !(ibz_is_one(&prod));
                res = res | !(ibz_is_one(u));
            } else {
                ibz_abs(&prod, u);
                res = res | !(ibz_is_one(&prod));
                res = res | !(ibz_is_zero(v));
            }
        }
        // Bezout coeffs
        ibz_mul(&sum, x, u);
        ibz_mul(&prod, y, v);
        ibz_add(&sum, &sum, &prod);
        res = res | !(ibz_cmp(&sum, gcd) == 0);
    }
    assert(!res);
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    ibz_finalize(&cmp);
    ibz_finalize(&test);

    return (res);
}
// void ibz_xgcd_with_u_not_0(ibz_t *d, ibz_t *u, ibz_t *v, const ibz_t *x, const ibz_t
// *y);
int
quat_test_ibz_xgcd_with_u_not_0(void)
{
    int res = 0;
    ibz_t x, y, u, v, d, s, p;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&d);
    ibz_init(&s);
    ibz_init(&p);

    ibz_set(&x, 75);
    ibz_set(&y, 50);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 25) != 0);

    ibz_set(&x, -75);
    ibz_set(&y, 50);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 25) != 0);

    ibz_set(&x, -75);
    ibz_set(&y, -50);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 25) != 0);

    ibz_set(&x, 75);
    ibz_set(&y, -50);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 25) != 0);

    ibz_set(&x, 50);
    ibz_set(&y, 50);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 50) != 0);

    ibz_set(&x, 0);
    ibz_set(&y, -50);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 50) != 0);

    ibz_set(&x, -50);
    ibz_set(&y, 0);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | (ibz_cmp_int32(&d, 50) != 0);

    ibz_set(&x, 0);
    ibz_set(&y, 0);
    ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
    res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);

    for (int i = -20; i <= 20; i++) {
        for (int j = -20; j <= 20; j++) {
            ibz_set(&x, i);
            ibz_set(&y, j);
            ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
            res = res | quat_test_helper_ibz_xgcd_with_u_not_0(&d, &u, &v, &x, &y);
        }
    }

    if (res != 0) {
        printf("Quaternion unit test ibz_xgcd_with_u_not_0 failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&d);
    ibz_finalize(&s);
    ibz_finalize(&p);
    return (res);
}

// void ibz_rounded_div(ibz_t *q, const ibz_t *a, const ibz_t *b);
int
quat_test_ibz_rounded_div(void)
{
    int res = 0;
    ibz_t q, a, b;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&q);

    // basic tests
    ibz_set(&a, 15);
    ibz_set(&b, 3);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 5) == 0);
    ibz_set(&a, 16);
    ibz_set(&b, 3);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 5) == 0);
    ibz_set(&a, 17);
    ibz_set(&b, 3);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 6) == 0);
    ibz_set(&a, 37);
    ibz_set(&b, 5);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 7) == 0);
    // test sign combination
    ibz_set(&a, 149);
    ibz_set(&b, 12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 12) == 0);
    ibz_set(&a, 149);
    ibz_set(&b, -12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -12) == 0);
    ibz_set(&a, -149);
    ibz_set(&b, -12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 12) == 0);
    ibz_set(&a, -149);
    ibz_set(&b, 12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -12) == 0);
    ibz_set(&a, 151);
    ibz_set(&b, 12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 13) == 0);
    ibz_set(&a, -151);
    ibz_set(&b, -12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 13) == 0);
    ibz_set(&a, 151);
    ibz_set(&b, -12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -13) == 0);
    ibz_set(&a, -151);
    ibz_set(&b, 12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -13) == 0);
    // divisibles with sign
    ibz_set(&a, 144);
    ibz_set(&b, 12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 12) == 0);
    ibz_set(&a, -144);
    ibz_set(&b, -12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 12) == 0);
    ibz_set(&a, 144);
    ibz_set(&b, -12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -12) == 0);
    ibz_set(&a, -144);
    ibz_set(&b, 12);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -12) == 0);
    // tests close to 0
    ibz_set(&a, -12);
    ibz_set(&b, -25);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 0) == 0);
    ibz_set(&a, 12);
    ibz_set(&b, 25);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 0) == 0);
    ibz_set(&a, -12);
    ibz_set(&b, 25);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 0) == 0);
    ibz_set(&a, 12);
    ibz_set(&b, -25);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 0) == 0);
    ibz_set(&a, -12);
    ibz_set(&b, -23);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 1) == 0);
    ibz_set(&a, 12);
    ibz_set(&b, 23);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, 1) == 0);
    ibz_set(&a, -12);
    ibz_set(&b, 23);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -1) == 0);
    ibz_set(&a, 12);
    ibz_set(&b, -23);
    ibz_rounded_div(&q, &a, &b);
    res = res || !(ibz_cmp_int32(&q, -1) == 0);
    // test output equal input
    ibz_set(&a, -151);
    ibz_set(&b, 12);
    ibz_rounded_div(&a, &a, &b);
    res = res || !(ibz_cmp_int32(&a, -13) == 0);
    ibz_set(&a, -151);
    ibz_set(&b, 12);
    ibz_rounded_div(&b, &a, &b);
    res = res || !(ibz_cmp_int32(&b, -13) == 0);
    // test for cmp not returning 1 or -1 or 0
    ibz_set(&a, 2046874253);
    ibz_set(&b, 412308266);
    ibz_rounded_div(&b, &a, &b);
    res = res || !(ibz_cmp_int32(&b, 5) == 0);

    if (res != 0) {
        printf("Quaternion unit test ibz_rounded_div failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&q);
    return (res);
}

// int ibz_mat_4x4_is_hnf(const ibz_mat_4x4_t *mat);
int
quat_test_ibz_mat_4x4_is_hnf(void)
{
    int res = 0;
    ibz_mat_4x4_t mat;
    ibz_mat_4x4_init(&mat);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(mat[i][j]), 0);
        }
    }
    res = res || (!ibz_mat_4x4_is_hnf(&mat));
    ibz_set(&(mat[0][0]), 7);
    ibz_set(&(mat[0][1]), 6);
    ibz_set(&(mat[0][2]), 5);
    ibz_set(&(mat[0][3]), 4);
    ibz_set(&(mat[1][1]), 6);
    ibz_set(&(mat[1][2]), 5);
    ibz_set(&(mat[1][3]), 4);
    ibz_set(&(mat[2][2]), 5);
    ibz_set(&(mat[2][3]), 4);
    ibz_set(&(mat[3][3]), 4);
    res = res || (!ibz_mat_4x4_is_hnf(&mat));

    ibz_set(&(mat[0][0]), 7);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 5);
    ibz_set(&(mat[0][3]), 4);
    ibz_set(&(mat[1][1]), 0);
    ibz_set(&(mat[1][2]), 0);
    ibz_set(&(mat[1][3]), 0);
    ibz_set(&(mat[2][2]), 5);
    ibz_set(&(mat[2][3]), 4);
    ibz_set(&(mat[3][3]), 4);
    res = res || (!ibz_mat_4x4_is_hnf(&mat));

    // negative tests
    ibz_set(&(mat[0][0]), 7);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 5);
    ibz_set(&(mat[0][3]), 4);
    ibz_set(&(mat[1][1]), 1);
    ibz_set(&(mat[1][2]), 5);
    ibz_set(&(mat[1][3]), 9);
    ibz_set(&(mat[2][2]), 5);
    ibz_set(&(mat[2][3]), 4);
    ibz_set(&(mat[3][3]), 4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));

    ibz_set(&(mat[0][0]), 7);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 5);
    ibz_set(&(mat[0][3]), 4);
    ibz_set(&(mat[1][1]), 1);
    ibz_set(&(mat[1][2]), -5);
    ibz_set(&(mat[1][3]), 1);
    ibz_set(&(mat[2][2]), 5);
    ibz_set(&(mat[2][3]), 4);
    ibz_set(&(mat[3][3]), 4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));

    ibz_set(&(mat[0][0]), 7);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 5);
    ibz_set(&(mat[0][3]), 4);
    ibz_set(&(mat[1][0]), 2);
    ibz_set(&(mat[1][1]), 3);
    ibz_set(&(mat[1][2]), 1);
    ibz_set(&(mat[1][3]), 1);
    ibz_set(&(mat[2][2]), 5);
    ibz_set(&(mat[2][3]), 4);
    ibz_set(&(mat[3][3]), 4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));

    ibz_set(&(mat[0][0]), 7);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 5);
    ibz_set(&(mat[0][3]), 4);
    ibz_set(&(mat[1][0]), 2);
    ibz_set(&(mat[1][1]), 3);
    ibz_set(&(mat[1][2]), -1);
    ibz_set(&(mat[1][3]), 7);
    ibz_set(&(mat[2][2]), 0);
    ibz_set(&(mat[2][3]), 0);
    ibz_set(&(mat[3][3]), 4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));

    if (res != 0) {
        printf("Quaternion unit test hnf_ibz_mat_4x4_is_hnf failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    return (res);
}

// void ibz_mat_4xn_hnf_core(ibz_mat_4x4_t *hnf, int generator_number, const
// ibz_vec_4_t
// (*generators)[generator_number]);
int
quat_test_ibz_mat_4xn_hnf_core(void)
{
    int res = 0;
    ibz_vec_4_t generators[8];
    ibz_mat_4x4_t hnf, cmp;
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&cmp);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_set(&(generators[j][i]), 0);
        }
    }
    ibz_mat_4xn_hnf_core(&hnf, 8, generators);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));
    // also should test that they generate the same lattice. Since HNF is unique, copute the HNF for
    // test vectors might be ok

    ibz_set(&(generators[2][0]), 2);
    ibz_set(&(generators[3][1]), 3);
    ibz_set(&(generators[4][0]), 4);
    ibz_set(&(generators[2][3]), 5);
    ibz_set(&(generators[7][3]), 6);
    ibz_set(&(generators[7][1]), 7);
    ibz_set(&(generators[3][1]), 8);
    ibz_set(&(generators[1][1]), 9);
    ibz_set(&(generators[6][0]), 10);
    ibz_set(&(generators[5][0]), 11);
    ibz_set(&(generators[0][0]), 12);
    ibz_mat_4xn_hnf_core(&hnf, 8, generators);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));

    ibz_set(&(generators[5][2]), 1);
    ibz_set(&(generators[0][2]), 2);
    ibz_mat_4xn_hnf_core(&hnf, 8, generators);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));

    // test equality of result to a known hnf
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_set(&(generators[j][i]), 0);
        }
    }
    ibz_set(&(generators[0][0]), 4);
    ibz_set(&(generators[2][0]), 3);
    ibz_set(&(generators[4][0]), 1);
    ibz_set(&(generators[7][0]), -1);
    ibz_set(&(generators[1][1]), 5);
    ibz_set(&(generators[5][1]), -2);
    ibz_set(&(generators[2][2]), 3);
    ibz_set(&(generators[6][2]), 1);
    ibz_set(&(generators[5][2]), 1);
    ibz_set(&(generators[3][3]), 7);
    ibz_set(&(generators[7][3]), -3);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(cmp[i][j]), 0);
        }
        ibz_set(&(cmp[i][i]), 1);
    }
    ibz_mat_4xn_hnf_core(&hnf, 8, generators);
    res = res || (!ibz_mat_4x4_equal(&cmp, &hnf));

    // test known hnf encountered in
    // https://github.com/SQISign/sqisign-nist/issues/38#issuecomment-1554585079
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_set(&(generators[j][i]), 0);
        }
    }
    ibz_set(&(generators[4][0]), 438);
    ibz_set(&(generators[4][1]), 400);
    ibz_set(&(generators[4][2]), 156);
    ibz_set(&(generators[4][3]), -2);
    ibz_set(&(generators[5][0]), -400);
    ibz_set(&(generators[5][1]), 438);
    ibz_set(&(generators[5][2]), 2);
    ibz_set(&(generators[5][3]), 156);
    ibz_set(&(generators[6][0]), -28826);
    ibz_set(&(generators[6][1]), -148);
    ibz_set(&(generators[6][2]), 220);
    ibz_set(&(generators[6][3]), -122);
    ibz_set(&(generators[7][0]), 586);
    ibz_set(&(generators[7][1]), -28426);
    ibz_set(&(generators[7][2]), 278);
    ibz_set(&(generators[7][3]), 218);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(cmp[i][j]), 0);
        }
    }
    ibz_set(&(cmp[0][0]), 2321156);
    ibz_set(&(cmp[1][1]), 2321156);
    ibz_set(&(cmp[0][2]), 620252);
    ibz_set(&(cmp[1][2]), 365058);
    ibz_set(&(cmp[2][2]), 2);
    ibz_set(&(cmp[0][3]), 1956098);
    ibz_set(&(cmp[1][3]), 620252);
    ibz_set(&(cmp[3][3]), 2);
    ibz_mat_4xn_hnf_core(&hnf, 8, generators);
    res = res || (!ibz_mat_4x4_equal(&cmp, &hnf));

    if (res != 0) {
        printf("Quaternion unit test dim4_ibz_mat_4xn_hnf_core failed\n");
    }
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&cmp);
    return (res);
}

// void ibz_mat_4xn_hnf_mod_core(ibz_mat_4x4_t *hnf, int generator_number, const ibz_vec_4_t
// (*generators)[generator_number], const ibz_t *mod);
int
quat_test_ibz_mat_4xn_hnf_mod_core(void)
{
    int res = 0;
    ibz_t det;
    ibz_vec_4_t generators[8];
    ibz_mat_4x4_t hnf, cmp;
    ibz_init(&det);
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&cmp);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));

    ibz_set(&(generators[2][0]), 2);
    ibz_set(&(generators[3][1]), 3);
    ibz_set(&(generators[4][0]), 4);
    ibz_set(&(generators[2][3]), 5);
    ibz_set(&(generators[7][3]), 6);
    ibz_set(&(generators[7][1]), 7);
    ibz_set(&(generators[3][1]), 8);
    ibz_set(&(generators[1][1]), 9);
    ibz_set(&(generators[6][0]), 10);
    ibz_set(&(generators[5][0]), 11);
    ibz_set(&(generators[0][0]), 12);
    ibz_set(&(generators[5][2]), 1);
    ibz_set(&(generators[0][2]), 2);
    ibz_set(&det, 4);
    ibz_mat_4xn_hnf_mod_core(&hnf, 8, generators, &det);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));

    // test equality of result to a known hnf
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_set(&(generators[j][i]), 0);
        }
    }
    ibz_set(&(generators[0][0]), 4);
    ibz_set(&(generators[2][0]), 3);
    ibz_set(&(generators[4][0]), 1);
    ibz_set(&(generators[7][0]), -1);
    ibz_set(&(generators[1][1]), 5);
    ibz_set(&(generators[5][1]), -2);
    ibz_set(&(generators[2][2]), 3);
    ibz_set(&(generators[6][2]), 1);
    ibz_set(&(generators[5][2]), 1);
    ibz_set(&(generators[3][3]), 7);
    ibz_set(&(generators[7][3]), -3);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(cmp[i][j]), 0);
        }
        ibz_set(&(cmp[i][i]), 1);
    }
    ibz_set(&det, 1);
    ibz_mat_4xn_hnf_mod_core(&hnf, 8, generators, &det);
    res = res || (!ibz_mat_4x4_equal(&cmp, &hnf));

    // test known hnf encountered in
    // https://github.com/SQISign/sqisign-nist/issues/38#issuecomment-1554585079
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 8; j++) {
            ibz_set(&(generators[j][i]), 0);
        }
    }
    ibz_set(&(generators[4][0]), 438);
    ibz_set(&(generators[4][1]), 400);
    ibz_set(&(generators[4][2]), 156);
    ibz_set(&(generators[4][3]), -2);
    ibz_set(&(generators[5][0]), -400);
    ibz_set(&(generators[5][1]), 438);
    ibz_set(&(generators[5][2]), 2);
    ibz_set(&(generators[5][3]), 156);
    ibz_set(&(generators[6][0]), -28826);
    ibz_set(&(generators[6][1]), -148);
    ibz_set(&(generators[6][2]), 220);
    ibz_set(&(generators[6][3]), -122);
    ibz_set(&(generators[7][0]), 586);
    ibz_set(&(generators[7][1]), -28426);
    ibz_set(&(generators[7][2]), 278);
    ibz_set(&(generators[7][3]), 218);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(cmp[i][j]), 0);
        }
    }
    ibz_set(&(cmp[0][0]), 2321156);
    ibz_set(&(cmp[1][1]), 2321156);
    ibz_set(&(cmp[0][2]), 620252);
    ibz_set(&(cmp[1][2]), 365058);
    ibz_set(&(cmp[2][2]), 2);
    ibz_set(&(cmp[0][3]), 1956098);
    ibz_set(&(cmp[1][3]), 620252);
    ibz_set(&(cmp[3][3]), 2);
    ibz_set_from_str(&det, "21551060705344", 10);
    ibz_mat_4xn_hnf_mod_core(&hnf, 8, generators, &det);
    res = res || (!ibz_mat_4x4_equal(&cmp, &hnf));

    // use non-modular hnf version to test
    ibz_set(&(generators[4][0]), 438);
    ibz_set(&(generators[4][0]), 438);
    ibz_set(&(generators[4][1]), 400);
    ibz_set(&(generators[4][2]), 156);
    ibz_set(&(generators[5][3]), -2);
    ibz_set(&(generators[5][0]), -40);
    ibz_set(&(generators[5][1]), 438);
    ibz_set(&(generators[5][2]), 20);
    ibz_set(&(generators[5][3]), 156);
    ibz_set(&(generators[6][0]), -28826);
    ibz_set(&(generators[6][1]), -148);
    ibz_set(&(generators[6][2]), 220);
    ibz_set(&(generators[6][3]), -122);
    ibz_set(&(generators[7][0]), 586);
    ibz_set(&(generators[7][1]), -28426);
    ibz_set(&(generators[7][2]), 278);
    ibz_set(&(generators[7][3]), 218);
    ibz_mat_4xn_hnf_core(&cmp, 8, generators);
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &cmp);
    ibz_abs(&det, &det);
    ibz_mat_4xn_hnf_mod_core(&hnf, 8, generators, &det);
    res = res || (!ibz_mat_4x4_equal(&cmp, &hnf));

    if (res != 0) {
        printf("Quaternion unit test dim4_ibz_mat_4x8_hnf_mod_core failed\n");
    }
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&cmp);
    ibz_finalize(&det);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    return (res);
}

// void ibz_mat_4xn_hnf_core(ibz_mat_4x4_t *hnf, int generator_number, const
// ibz_vec_4_t
// (*generators)[generator_number]);
int
quat_test_ibz_mat_4xn_hnf_core_randomized(void)
{
    // only partial test, since lattice equality cannot be tested without hnf, so only one inclusion
    // is checked
    int res = 0;
    ibz_vec_4_t generators[8];
    quat_alg_elem_t vec;
    quat_lattice_t lat;
    ibz_mat_4x4_t hnf, cmp;
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));
    quat_lattice_init(&lat);
    quat_alg_elem_init(&vec);
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&cmp);
    int32_t rand[8][4];
    int det_non_0;

    for (int iter = 0; iter < 10; iter++) {
        int randret = randombytes((unsigned char *)rand, sizeof(rand));
        if (randret != 0)
            return 1;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                ibz_set(&(generators[j][i]), rand[j][i]);
            }
        }
        ibz_mat_4xn_hnf_core(&hnf, 8, generators);
        res = res || (!ibz_mat_4x4_is_hnf(&hnf));
        // also should test that they generate the same lattice. However, can only check one
        // inclusion efficiently (and this only if full rank), so do so
        det_non_0 = 1;
        for (int i = 0; i < 4; i++) {
            det_non_0 = det_non_0 && ibz_is_zero(&(hnf[i][i]));
        }
        if (det_non_0) {
            ibz_mat_4x4_copy(&(lat.basis), &hnf);
            ibz_set(&(lat.denom), 1);
            ibz_set(&(vec.denom), 1);
            for (int i = 0; i < 8; i++) {
                ibz_vec_4_copy(&(vec.coord), &(generators[i]));
                res = res || !quat_lattice_contains(NULL, &lat, &vec);
            }
        }
    }

    if (res != 0) {
        printf("Quaternion unit test with randomization for ibz_mat_4xn_hnf_core "
               "failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&vec);
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&cmp);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    return (res);
}

// void ibz_mat_4xn_hnf_mod_core(ibz_mat_4x4_t *hnf, int generator_number, const ibz_vec_4_t
// (*generators)[generator_number],const ibz_t *mod);
int
quat_test_ibz_mat_4xn_hnf_mod_core_randomized(void)
{
    // only partial test, since lattice equality cannot be tested without hnf, so only one inclusion
    // is checked
    int res = 0;
    int32_t rand[8][4];
    int32_t rand_m;
    int randret;
    quat_alg_elem_t vec;
    quat_lattice_t lat;
    ibz_vec_4_t generators[8];
    ibz_t det, m;
    ibz_mat_4x4_t hnf, cmp;
    ibz_init(&det);
    ibz_init(&m);
    quat_lattice_init(&lat);
    quat_alg_elem_init(&vec);
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&cmp);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));

    for (int iter = 0; iter < 10; iter++) {
        randret = randombytes((unsigned char *)rand, sizeof(rand));
        if (randret != 0) {
            res = 1;
            goto end;
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                ibz_set(&(generators[j][i]), rand[j][i]);
            }
        }
        rand_m = 0;
        while (rand_m <= 0) {
            randret = randombytes((unsigned char *)&rand_m, sizeof(int32_t));
            if (randret != 0) {
                res = 1;
                goto end;
            }
        }
        ibz_set(&m, rand_m);

        ibz_mat_4xn_hnf_core(&cmp, 8, generators);
        res = res || (!ibz_mat_4x4_is_hnf(&cmp));
        ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &cmp);
        if (!ibz_is_zero(&det)) {
            ibz_mul(&det, &det, &m);
            ibz_mat_4xn_hnf_mod_core(&hnf, 8, generators, &det);
            res = res || (!ibz_mat_4x4_equal(&hnf, &cmp));
        }
    }

end:;
    if (res != 0) {
        printf("Quaternion unit test with randomization for ibz_mat_4xn_hnf_mod_core failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&vec);
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&cmp);
    ibz_finalize(&det);
    ibz_finalize(&m);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    return (res);
}

// run all previous tests
int
quat_test_hnf(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of matrices, vectors and quadratic forms in dimension 4\n");
    res = res | quat_test_ibz_mod_not_zero();
    res = res | quat_test_ibz_centered_mod();
    res = res | quat_test_ibz_conditional_assign();
    res = res | quat_test_ibz_xgcd_with_u_not_0();
    res = res | quat_test_ibz_rounded_div();
    res = res | quat_test_ibz_mat_4x4_is_hnf();
    res = res | quat_test_ibz_mat_4xn_hnf_core();
    res = res | quat_test_ibz_mat_4xn_hnf_mod_core();
    res = res | quat_test_ibz_mat_4xn_hnf_core_randomized();
    res = res | quat_test_ibz_mat_4xn_hnf_mod_core_randomized();
    return (res);
}
