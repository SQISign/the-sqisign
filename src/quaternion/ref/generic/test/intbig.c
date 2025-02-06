#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <rng.h>
#include "intbig_internal.h"

// void ibz_init(ibz_t *x);
// void ibz_finalize(ibz_t *x);
// int ibz_cmp(const ibz_t *a, const ibz_t *b);
// int ibz_is_zero(const ibz_t *x);
// int ibz_is_one(const ibz_t *x);
// int ibz_cmp_int32(const ibz_t *x, int32_t y);
// int ibz_is_even(const ibz_t *x);
// int ibz_is_odd(const ibz_t *x);
// void ibz_set(ibz_t *i, int32_t x);
// void ibz_copy(ibz_t *target, const ibz_t *value);
// void ibz_swap(ibz_t *a, ibz_t *b);
// int32_t ibz_get(const ibz_t *i);
// int ibz_bitsize(const ibz_t *a);
int
ibz_test_init_set_cmp()
{
    int res = 0;
    ibz_t a, b, c;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);

    res = res | !ibz_is_zero(&a);
    ibz_set(&a, 1);
    res = res | !ibz_is_one(&a);
    res = res | ibz_is_zero(&a);
    res = res | (0 != ibz_cmp(&b, &c));
    res = res | (0 == ibz_cmp(&a, &c));
    res = res | ibz_is_even(&a);
    res = res | !ibz_is_odd(&a);
    res = res | !ibz_is_even(&b);
    res = res | ibz_is_odd(&b);
    ibz_copy(&b, &a);
    res = res | !ibz_is_one(&a);
    res = res | !ibz_is_one(&b);
    res = res | (0 != ibz_cmp(&a, &b));
    res = res | (0 == ibz_cmp(&c, &b));
    ibz_swap(&b, &c);
    res = res | (0 == ibz_cmp(&a, &b));
    res = res | (0 == ibz_cmp(&c, &b));
    res = res | (0 != ibz_cmp(&a, &c));
    res = res | (ibz_bitsize(&a) != 1);
    res = res | (ibz_get(&a) != 1);
    res = res | (ibz_cmp_int32(&a, 1) != 0);
    res = res | (ibz_bitsize(&b) > 1);
    res = res | (ibz_get(&b) != 0);
    res = res | (ibz_cmp_int32(&b, 0) != 0);

    ibz_set(&a, -1);
    res = res | !(ibz_cmp(&a, &b) < 0);
    res = res | !(ibz_cmp(&b, &c) < 0);
    res = res | !(ibz_cmp(&b, &a) > 0);
    res = res | !(ibz_cmp(&c, &b) > 0);
    ibz_copy(&b, &a);
    res = res | (0 != ibz_cmp(&a, &b));
    res = res | (0 == ibz_cmp(&c, &b));
    ibz_swap(&b, &c);
    res = res | (0 == ibz_cmp(&a, &b));
    res = res | (0 != ibz_cmp(&c, &a));

    ibz_set_from_str(&a, "-10000000000000000011111100000001", 10);
    res = res | !(ibz_cmp(&a, &b) < 0);
    res = res | !(ibz_cmp(&b, &a) > 0);
    ibz_copy(&b, &a);
    res = res | (0 != ibz_cmp(&a, &b));
    res = res | (0 == ibz_cmp(&c, &b));
    ibz_swap(&b, &c);
    res = res | (0 == ibz_cmp(&a, &b));
    res = res | (0 != ibz_cmp(&c, &a));

    ibz_set_from_str(&a, "1aaaa00000000000000000123", 16);
    res = res | !(ibz_cmp(&a, &c) > 0);
    res = res | !(ibz_cmp(&c, &a) < 0);
    res = res | (ibz_bitsize(&a) != 4 * 24 + 1);
    res = res | (ibz_get(&a) != 16 * 18 + 3);
    if (res) {
        printf("Quaternion module integer group test ibz_test_init_set_cmp failed\n");
    }

    ibz_set_from_str(&a, "deadbeef12345678", 16);
    res = res | (ibz_get(&a) != 0x12345678);

    // Test INT32_MAX and INT32_MIN
    ibz_set_from_str(&a, "-2147483648", 10);
    ibz_set(&b, INT32_MIN);
    res = res | (ibz_get(&a) != -2147483648);
    res = res | (ibz_cmp_int32(&a, -2147483648) != 0);
    res = res | (0 != ibz_cmp(&a, &b));
    ibz_set_from_str(&a, "2147483647", 10);
    ibz_set(&b, INT32_MAX);
    res = res | (ibz_get(&a) != 2147483647);
    res = res | (ibz_cmp_int32(&a, 2147483647) != 0);
    res = res | (0 != ibz_cmp(&a, &b));

    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    return (res);
}

// void ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b);
// void ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b);
// void ibz_neg(ibz_t *neg, const ibz_t *a);
// void ibz_abs(ibz_t *abs, const ibz_t *a);
int
ibz_test_add_sub_neg_abs()
{
    int res = 0;
    ibz_t a, b, c, d, r, q, m;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&r);
    ibz_init(&q);
    ibz_init(&m);

    ibz_set_from_str(&a, "10000111100002222", 16);
    ibz_copy(&c, &a);
    ibz_add(&d, &a, &b);
    res = res | (0 != ibz_cmp(&b, &q));
    res = res | (0 != ibz_cmp(&a, &c));
    res = res | (0 != ibz_cmp(&d, &c));
    ibz_add(&d, &b, &a);
    res = res | (0 != ibz_cmp(&b, &q));
    res = res | (0 != ibz_cmp(&a, &c));
    res = res | (0 != ibz_cmp(&d, &c));
    // add test for adding non-0, sub, mul, ...
    ibz_set_from_str(&b, "20000111100002223", 16);
    ibz_set_from_str(&c, "30000222200004445", 16);
    ibz_add(&d, &a, &b);
    res = res | (0 != ibz_cmp(&d, &c));
    ibz_neg(&m, &a);
    ibz_add(&q, &m, &d);
    res = res | (0 != ibz_cmp(&q, &b));
    ibz_neg(&m, &m);
    res = res | (0 != ibz_cmp(&m, &a));
    ibz_add(&q, &m, &b);
    res = res | (0 != ibz_cmp(&q, &c));
    ibz_neg(&m, &m);
    ibz_sub(&q, &b, &m);
    res = res | (0 != ibz_cmp(&q, &c));
    ibz_neg(&m, &m);
    res = res | (0 != ibz_cmp(&a, &m));
    ibz_sub(&q, &d, &a);
    res = res | (0 != ibz_cmp(&q, &b));
    ibz_neg(&d, &d);
    ibz_neg(&m, &a);
    ibz_neg(&c, &b);
    ibz_sub(&q, &d, &m);
    res = res | (0 != ibz_cmp(&q, &c));
    ibz_sub(&q, &m, &b);
    res = res | (0 != ibz_cmp(&q, &d));
    ibz_abs(&r, &r);
    res = res | !ibz_is_zero(&r);
    ibz_abs(&r, &a);
    res = res | (0 != ibz_cmp(&a, &r));
    ibz_abs(&r, &m);
    res = res | (0 != ibz_cmp(&a, &r));
    ibz_neg(&m, &m);
    res = res | (0 != ibz_cmp(&a, &m));

    if (res) {
        printf("Quaternion module integer group test ibz_test_add_sub_neg_abs failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&r);
    ibz_finalize(&q);
    ibz_finalize(&m);
    return (res);
}

// void ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b);
// void ibz_sqrt_floor(ibz_t *sqrt, const ibz_t *a);
int
ibz_test_mul_sqrt()
{
    int res = 0;
    ibz_t a, b, c, d, r, q, m;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&r);
    ibz_init(&q);
    ibz_init(&m);

    // zero
    ibz_set_from_str(&a, "2113309833171849999003363", 10);
    ibz_copy(&c, &a);
    ibz_set(&b, 0);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &b));
    res = res | (0 != ibz_cmp(&a, &c));
    ibz_mul(&m, &b, &a);
    res = res | (0 != ibz_cmp(&m, &b));
    res = res | (0 != ibz_cmp(&a, &c));
    // one
    ibz_set(&b, 1);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &a));
    res = res | (0 != ibz_cmp(&a, &c));
    ibz_mul(&m, &b, &a);
    res = res | (0 != ibz_cmp(&m, &a));
    res = res | (0 != ibz_cmp(&a, &c));
    // -1
    ibz_neg(&b, &b);
    ibz_neg(&d, &a);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &d));
    res = res | (0 != ibz_cmp(&a, &c));
    ibz_mul(&m, &b, &a);
    res = res | (0 != ibz_cmp(&m, &d));
    res = res | (0 != ibz_cmp(&a, &c));
    // larger
    ibz_set_from_str(&b, "34575345632322576567896", 10);
    ibz_set_from_str(&c, "73068417910102676801285574959599851857101834248", 10);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &c));
    ibz_neg(&b, &b);
    ibz_neg(&a, &a);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &c));
    ibz_neg(&a, &a);
    ibz_neg(&c, &c);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &c));
    ibz_neg(&b, &b);
    ibz_neg(&a, &a);
    ibz_mul(&m, &a, &b);
    res = res | (0 != ibz_cmp(&m, &c));

    // sqrt tests
    ibz_abs(&b, &b);
    ibz_abs(&a, &a);
    ibz_mul(&m, &a, &a);
    ibz_sqrt_floor(&d, &m);
    res = res | (0 != ibz_cmp(&d, &a));
    ibz_mul(&m, &b, &b);
    ibz_sqrt_floor(&d, &m);
    res = res | (0 != ibz_cmp(&d, &b));
    ibz_set(&d, 1);
    ibz_mul(&m, &b, &b);
    ibz_sub(&m, &m, &d);
    ibz_sub(&c, &b, &d);
    ibz_sqrt_floor(&d, &m);
    res = res | (0 != ibz_cmp(&d, &c));

    if (res) {
        printf("Quaternion module integer group test ibz_test_mul_sqrt failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&r);
    ibz_finalize(&q);
    ibz_finalize(&m);
    return (res);
}

// void ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b);
// int ibz_divides(const ibz_t *a, const ibz_t *b);
int
ibz_test_div()
{
    int res = 0;
    ibz_t a, b, c, d, r, q, m;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&r);
    ibz_init(&q);
    ibz_init(&m);

    // one
    ibz_set_from_str(&a, "2113309833171849999003363", 10);
    ibz_copy(&c, &a);
    ibz_set(&b, 1);
    ibz_copy(&d, &b);
    res = res | !ibz_divides(&a, &b);
    ibz_div(&q, &r, &a, &b);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &a));
    res = res | (0 != ibz_cmp(&c, &a));
    res = res | (0 != ibz_cmp(&b, &d));
    // not one, zero remainder
    ibz_set_from_str(&b, "15678200126527887351125", 10);
    ibz_mul(&d, &a, &b);
    ibz_copy(&c, &d);
    res = res | !ibz_divides(&d, &b);
    res = res | !ibz_divides(&d, &a);
    res = res | !ibz_divides(&d, &d);
    ibz_div(&q, &r, &d, &b);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &a));
    ibz_div(&q, &r, &d, &a);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &b));
    // flipping signs
    ibz_neg(&a, &a);
    ibz_neg(&b, &b);
    ibz_div(&q, &r, &d, &b);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &a));
    ibz_div(&q, &r, &d, &a);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &b));
    ibz_neg(&a, &a);
    ibz_neg(&d, &d);
    ibz_div(&q, &r, &d, &b);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &a));
    ibz_div(&q, &r, &d, &a);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &b));
    ibz_neg(&a, &a);
    ibz_neg(&b, &b);
    ibz_div(&q, &r, &d, &b);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &a));
    ibz_div(&q, &r, &d, &a);
    res = res | !ibz_is_zero(&r);
    res = res | (0 != ibz_cmp(&q, &b));
    // non-zero remainder
    ibz_neg(&a, &a);
    ibz_neg(&d, &d);
    ibz_set_from_str(&c, "8678205677345432110000", 10);
    ibz_add(&d, &d, &c);
    ibz_div(&q, &r, &d, &b);
    ibz_mul(&m, &q, &b);
    ibz_add(&m, &m, &r);
    res = res | (0 != ibz_cmp(&d, &m));
    ibz_abs(&m, &r);
    ibz_abs(&q, &b);
    res = res | (0 <= ibz_cmp(&m, &q));
    ibz_set(&q, 0);
    res =
        res | !(((ibz_cmp(&q, &r) <= 0) && (ibz_cmp(&q, &d) < 0)) || ((ibz_cmp(&q, &r) >= 0) && (ibz_cmp(&r, &d) > 0)));
    // flip signs
    ibz_neg(&d, &d);
    ibz_div(&q, &r, &d, &b);
    ibz_mul(&m, &q, &b);
    ibz_add(&m, &m, &r);
    res = res | (0 != ibz_cmp(&d, &m));
    ibz_abs(&m, &r);
    ibz_abs(&q, &b);
    res = res | (0 <= ibz_cmp(&m, &q));
    ibz_set(&q, 0);
    res =
        res | !(((ibz_cmp(&q, &r) <= 0) && (ibz_cmp(&q, &d) < 0)) || ((ibz_cmp(&q, &r) >= 0) && (ibz_cmp(&r, &d) > 0)));
    ibz_neg(&b, &b);
    ibz_div(&q, &r, &d, &b);
    ibz_mul(&m, &q, &b);
    ibz_add(&m, &m, &r);
    res = res | (0 != ibz_cmp(&d, &m));
    ibz_abs(&m, &r);
    ibz_abs(&q, &b);
    res = res | (0 <= ibz_cmp(&m, &q));
    ibz_set(&q, 0);
    res =
        res | !(((ibz_cmp(&q, &r) <= 0) && (ibz_cmp(&q, &d) < 0)) || ((ibz_cmp(&q, &r) >= 0) && (ibz_cmp(&r, &d) > 0)));
    ibz_neg(&d, &d);
    ibz_div(&q, &r, &d, &b);
    ibz_mul(&m, &q, &b);
    ibz_add(&m, &m, &r);
    res = res | (0 != ibz_cmp(&d, &m));
    ibz_abs(&m, &r);
    ibz_abs(&q, &b);
    res = res | (0 <= ibz_cmp(&m, &q));
    ibz_set(&q, 0);
    res =
        res | !(((ibz_cmp(&q, &r) <= 0) && (ibz_cmp(&q, &d) < 0)) || ((ibz_cmp(&q, &r) >= 0) && (ibz_cmp(&r, &d) > 0)));

    if (res) {
        printf("Quaternion module integer group test ibz_test_div failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&r);
    ibz_finalize(&q);
    ibz_finalize(&m);
    return (res);
}

// void ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b);
// unsigned long int ibz_mod_ui(const mpz_t *n, unsigned long int d);
int
ibz_test_mod()
{
    int res = 0;
    ibz_t a, b, c, r, m;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&r);
    ibz_init(&m);

    ibz_set_from_str(&a, "2113309833171849999003363", 10);
    res = res | (ibz_mod_ui(&a, 3) != 0);
    res = res | (ibz_mod_ui(&a, 2) != 1);
    ibz_set_from_str(&m, "2113309833171840000000000", 10);
    ibz_mod(&r, &a, &m);
    ibz_add(&c, &r, &m);
    res = res | (0 != ibz_cmp(&c, &a));
    ibz_neg(&b, &a);
    res = res | (ibz_mod_ui(&b, 3) != 0);
    res = res | (ibz_mod_ui(&b, 2) != 1);
    ibz_mod(&r, &b, &m);
    ibz_sub(&c, &r, &m);
    ibz_sub(&c, &c, &m);
    res = res | (0 != ibz_cmp(&c, &b));

    if (res) {
        printf("Quaternion module integer group test ibz_test_mod failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&r);
    ibz_finalize(&m);
    return (res);
}

// void ibz_pow(ibz_t *pow, const ibz_t *x, uint32_t e);
// void ibz_pow_mod(ibz_t *pow, const ibz_t *x, const ibz_t *e, const ibz_t *m);
// void ibz_div_2exp(ibz_t *quotient, const ibz_t *a, uint32_t exp);
int
ibz_test_pow()
{
    int res = 0;
    ibz_t a, b, c, d, e, m;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&e);
    ibz_init(&m);

    ibz_set_from_str(&a, "aaaaaaaabbbbbbbb2222221111", 16);
    ibz_copy(&c, &a);
    int exp = 10;
    ibz_set(&b, 1);
    for (int i = 1; i < exp + 1; i++)
        ibz_mul(&b, &b, &a);
    ibz_pow(&d, &a, exp);
    res = res | (0 != ibz_cmp(&c, &a));
    res = res | (0 != ibz_cmp(&d, &b));
    ibz_neg(&a, &a);
    ibz_copy(&c, &a);
    ibz_pow(&d, &a, exp);
    res = res | (0 != ibz_cmp(&c, &a));
    res = res | (0 != ibz_cmp(&d, &b));
    exp = 9;
    ibz_set(&b, 1);
    for (int i = 1; i < exp + 1; i++)
        ibz_mul(&b, &b, &a);
    ibz_pow(&d, &a, exp);
    res = res | (0 != ibz_cmp(&c, &a));
    res = res | (0 != ibz_cmp(&d, &b));

    ibz_set_from_str(&a, "aaaaaaaabbbbbbbb2222221111", 16);
    ibz_set_from_str(&m, "cdabde24864122912341", 16);
    exp = 10;
    ibz_set(&m, exp);
    ibz_copy(&c, &a);
    ibz_set(&b, 1);
    for (int i = 1; i < exp + 1; i++) {
        ibz_mul(&b, &b, &a);
        ibz_mod(&b, &b, &m);
    }
    ibz_pow_mod(&d, &a, &e, &m);
    res = res | (0 != ibz_cmp(&c, &a));
    res = res | (0 != ibz_cmp(&d, &b));

    exp = 23;
    ibz_set(&b, 2);
    ibz_pow(&b, &b, exp);
    ibz_div(&m, &b, &a, &b);
    ibz_div_2exp(&d, &a, exp);
    res = res | (0 != ibz_cmp(&c, &a));
    res = res | (0 != ibz_cmp(&m, &d));

    if (res) {
        printf("Quaternion module integer group test ibz_test_pow failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&e);
    ibz_finalize(&m);
    return (res);
}

// void ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b);
// int ibz_invmod(ibz_t *inv, const ibz_t *a, const ibz_t *mod);ibz_test_mod()
int
ibz_test_gcd()
{
    int res = 0;
    ibz_t a, b, c, d, ac, bc, m;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&ac);
    ibz_init(&bc);
    ibz_init(&m);

    // fixed large test
    ibz_set_from_str(&c, "25791357069084", 10);
    ibz_set_from_str(&a, "6173271838293993987767", 10);
    ibz_set_from_str(&b, "89882267321617266071838286", 10);
    ibz_mul(&ac, &a, &c);
    ibz_mul(&bc, &b, &c);
    ibz_gcd(&d, &ac, &bc);
    res = res | (0 != ibz_cmp(&c, &d));
    ibz_mul(&m, &a, &c);
    res = res | (0 != ibz_cmp(&ac, &m));
    ibz_mul(&m, &b, &c);
    res = res | (0 != ibz_cmp(&bc, &m));
    // gcd is positif whatever the inputs are
    ibz_neg(&ac, &ac);
    ibz_gcd(&d, &ac, &bc);
    res = res | (0 != ibz_cmp(&c, &d));
    ibz_neg(&bc, &bc);
    ibz_gcd(&d, &ac, &bc);
    res = res | (0 != ibz_cmp(&c, &d));
    ibz_neg(&ac, &ac);
    ibz_gcd(&d, &ac, &bc);
    res = res | (0 != ibz_cmp(&c, &d));
    // different sizes do not matter
    ibz_set(&d, 2);
    ibz_gcd(&m, &ac, &d);
    res = res | (0 != ibz_cmp(&m, &d));

    // invmod
    ibz_invmod(&d, &a, &b);
    res = res | (0 >= ibz_cmp(&b, &d));
    ibz_set(&m, 0);
    res = res | (0 <= ibz_cmp(&m, &d));
    ibz_mul(&m, &d, &a);
    ibz_mod(&m, &m, &b);
    res = res | !ibz_is_one(&m);
    // negative changes nothing
    ibz_neg(&a, &a);
    ibz_invmod(&d, &a, &b);
    res = res | (0 >= ibz_cmp(&b, &d));
    ibz_set(&m, 0);
    res = res | (0 <= ibz_cmp(&m, &d));
    ibz_mul(&m, &d, &a);
    ibz_mod(&m, &m, &b);
    res = res | !ibz_is_one(&m);

    if (res) {
        printf("Quaternion module integer group test ibz_test_gcd failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&ac);
    ibz_finalize(&bc);
    ibz_finalize(&m);
    return (res);
}

// Tests ibz_sqrt_mod_p
// Allows to provide the number of repetitions and the bit-size of the primes.
int
ibz_test_sqrt_mod_p(int reps, int prime_n)
{
    int res = 0;

    // Initialize GMP variables
    ibz_t prime, a, prime_minus_a, asq, sqrt, tmp;
    ibz_t prime_p4m3, prime_p5m8, prime_p1m8;
    ibz_t prime_p4m3_x2, prime_p5m8_x2, prime_p1m8_x2;
    ibz_t two_to_the_n_minus_1;
    ibz_init(&prime);
    ibz_init(&prime_p4m3);
    ibz_init(&prime_p5m8);
    ibz_init(&prime_p1m8);
    ibz_init(&prime_p4m3_x2);
    ibz_init(&prime_p5m8_x2);
    ibz_init(&prime_p1m8_x2);
    ibz_init(&a);
    ibz_init(&prime_minus_a);
    ibz_init(&asq);
    ibz_init(&sqrt);
    ibz_init(&tmp);
    ibz_init(&two_to_the_n_minus_1);

    // Generate random prime number
    int n = prime_n; // Number of bits
    ibz_set(&two_to_the_n_minus_1, 1);
    mpz_mul_2exp(two_to_the_n_minus_1, two_to_the_n_minus_1, n);
    ibz_sub(&two_to_the_n_minus_1, &two_to_the_n_minus_1, &ibz_const_one);

    for (int r = 0; r < reps; ++r) {
        ibz_rand_interval(&prime, &ibz_const_zero, &two_to_the_n_minus_1);
        if (ibz_is_even(&prime))
            ibz_add(&prime, &prime, &ibz_const_one);

        int p4m3 = 0, p5m8 = 0, p1m8 = 0;

        while (p4m3 == 0 || p5m8 == 0 || p1m8 == 0) {
            do {
                ibz_add(&prime, &prime, &ibz_const_two);
            } while (!mpz_probab_prime_p(prime, 25));

            if (mpz_mod_ui(tmp, prime, 4) == 3) {
                ibz_copy(&prime_p4m3, &prime);
                p4m3 = 1;
            } else if (mpz_mod_ui(tmp, prime, 8) == 5) {
                ibz_copy(&prime_p5m8, &prime);
                p5m8 = 1;
            } else if (mpz_mod_ui(tmp, prime, 8) == 1) {
                ibz_copy(&prime_p1m8, &prime);
                p1m8 = 1;
            } else {
                res = 1;
                goto err;
            }
        }

        ibz_t *primes[] = { &prime_p4m3, &prime_p5m8, &prime_p1m8 };
        ibz_t *primes_x2[] = { &prime_p4m3_x2, &prime_p5m8_x2, &prime_p1m8_x2 };
        for (int i = 0; i < 3; ++i) // 2p
            mpz_mul_2exp(*primes_x2[i], *primes[i], 1);

        // Test sqrt mod p
        for (int i = 0; i < 3; ++i) {
            ibz_sub(&tmp, primes[i], &ibz_const_one);
            ibz_rand_interval(&a, &ibz_const_zero, &tmp);
            ibz_sub(&prime_minus_a, (primes[i]), &a);
            mpz_powm_ui(asq, a, 2, *primes[i]);

            int no_sqrt = !ibz_sqrt_mod_p(&sqrt, &asq, primes[i]);
            mpz_powm_ui(tmp, sqrt, 2, *primes[i]);

            if (no_sqrt || (ibz_cmp(&sqrt, &a) && ibz_cmp(&sqrt, &prime_minus_a))) {
                res = 1;
                goto err;
            }
        }
    }

err:
    if (res) {
        printf("Quaternion module integer test ibz_test_sqrt_mod_p failed\n");
    }
    ibz_finalize(&prime);
    ibz_finalize(&prime_p4m3);
    ibz_finalize(&prime_p5m8);
    ibz_finalize(&prime_p1m8);
    ibz_finalize(&prime_p4m3_x2);
    ibz_finalize(&prime_p5m8_x2);
    ibz_finalize(&prime_p1m8_x2);
    ibz_finalize(&a);
    ibz_finalize(&prime_minus_a);
    ibz_finalize(&asq);
    ibz_finalize(&sqrt);
    ibz_finalize(&tmp);
    ibz_finalize(&two_to_the_n_minus_1);
    return res;
}

int
ibz_test_rand_interval(int reps)
{
    int res = 0;
    ibz_t low, high, rand;
    ibz_init(&low);
    ibz_init(&high);
    ibz_init(&rand);
    ibz_set_from_str(&low, "ffa", 16);
    ibz_set_from_str(&high, "eeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeef", 16);

    for (int i = 0; i < reps; ++i) {
        res = ibz_rand_interval(&rand, &low, &high);
        if (res != 1) {
            res = 1;
            goto err;
        } else {
            res = 0;
        }

        if (ibz_cmp(&rand, &low) < 0 || ibz_cmp(&rand, &high) > 0) {
            res = 1;
            goto err;
        }
    }
err:
    if (res) {
        printf("Quaternion module integer test ibz_test_rand_interval failed\n");
    }
    ibz_finalize(&low);
    ibz_finalize(&high);
    ibz_finalize(&rand);
    return res;
}

int
ibz_test_rand_interval_i(int reps)
{
    int res = 0;
    int32_t low, high;
    ibz_t rand, cmp;
    ibz_init(&rand);
    ibz_init(&cmp);

    for (int i = 0; i < reps; ++i) {
        randombytes((unsigned char *)&low, sizeof(int32_t));
        randombytes((unsigned char *)&high, sizeof(int32_t));
        if (low < 0)
            low = -low;
        if (high < 0)
            high = -high;
        // The function requires a < b, thus a != b
        if (low == high) {
            continue;
        }
        if (low > high) {
            int32_t tmp = low;
            low = high;
            high = tmp;
        }

        res = ibz_rand_interval_i(&rand, low, high);
        if (res != 1) {
            res = 1;
            goto err;
        } else {
            res = 0;
        }
        ibz_set(&cmp, low);
        if (ibz_cmp(&rand, &cmp) < 0) {
            res = 1;
            goto err;
        }
        ibz_set(&cmp, high);
        if (ibz_cmp(&rand, &cmp) > 0) {
            res = 1;
            goto err;
        }
    }

err:
    if (res) {
        printf("Quaternion module integer test ibz_test_rand_interval_i failed\n");
    }
    ibz_finalize(&rand);
    ibz_finalize(&cmp);
    return res;
}

int
ibz_test_rand_interval_minm_m(int reps)
{
    int res = 0;
    int32_t m;
    ibz_t rand;
    ibz_init(&rand);

    for (int i = 0; i < reps; ++i) {
        randombytes((unsigned char *)&m, sizeof(int32_t));
        if (m < 0)
            m = -m;
        m >>= 1; // less than 32 bit

        if (m < 0) {
            res = 1;
            goto err;
        }

        res = ibz_rand_interval_minm_m(&rand, m);
        if (res != 1) {
            res = 1;
            goto err;
        } else {
            res = 0;
        }
        if (ibz_cmp_int32(&rand, -m) < 0) {
            res = 1;
            goto err;
        }
        if (ibz_cmp_int32(&rand, m) > 0) {
            res = 1;
            goto err;
        }
    }

err:
    if (res) {
        printf("Quaternion module integer test ibz_test_rand_interval_minm_m failed\n");
    }
    ibz_finalize(&rand);
    return res;
}

int
ibz_test_copy_digits(void)
{
    int res = 0;
    const digit_t d1[] = { 0x12345678 };
    const digit_t d2[] = { 2, 1 };

    const char d1str[] = "12345678";
#if RADIX == 32
    const char d2str[] = "100000002";
#elif RADIX == 64
    const char d2str[] = "10000000000000002";
#endif

    char d1_intbig_str[80] = { 0 };
    char d2_intbig_str[80] = { 0 };

    ibz_t d1_intbig, d2_intbig;
    ibz_init(&d1_intbig);
    ibz_init(&d2_intbig);

    ibz_copy_digits(&d1_intbig, d1, 1);
    ibz_copy_digits(&d2_intbig, d2, 2);

    ibz_convert_to_str(&d1_intbig, d1_intbig_str, 16);
    ibz_convert_to_str(&d2_intbig, d2_intbig_str, 16);

    if (memcmp(d1str, d1_intbig_str, sizeof(d1str))) {
        res = 1;
        goto err;
    }
    if (memcmp(d2str, d2_intbig_str, sizeof(d2str))) {
        res = 1;
        goto err;
    }
err:
    if (res) {
        printf("Quaternion module integer test ibz_test_copy_digits failed\n");
    }
    ibz_finalize(&d1_intbig);
    ibz_finalize(&d2_intbig);
    return res;
}

int
ibz_test_to_digits(void)
{
    int res = 0;
    ibz_t d1_intbig, d2_intbig, zero_intbig;
    ibz_t d1_intbig_rec, d2_intbig_rec, zero_intbig_rec, cof, cof2;
    const char d1str[] = "12345678";
    const char d2str[] = "10000000000000002";
    ibz_init(&d1_intbig);
    ibz_init(&d2_intbig);

    ibz_set_from_str(&d1_intbig, d1str, 16);
    ibz_set_from_str(&d2_intbig, d2str, 16);
    ibz_init(&zero_intbig);

    ibz_init(&d1_intbig_rec);
    ibz_init(&d2_intbig_rec);
    ibz_init(&zero_intbig_rec);
    ibz_init(&cof);
    ibz_init(&cof2);

    size_t d1_digits = (mpz_sizeinbase(d1_intbig, 2) + sizeof(digit_t) * 8 - 1) / (sizeof(digit_t) * 8);
    size_t d2_digits = (mpz_sizeinbase(d2_intbig, 2) + sizeof(digit_t) * 8 - 1) / (sizeof(digit_t) * 8);

    digit_t d1[d1_digits];
    digit_t d2[d2_digits];
    digit_t zero[1];

    ibz_to_digits(d1, &d1_intbig);
    ibz_to_digits(d2, &d2_intbig);
    ibz_to_digits(zero, &zero_intbig);

    // A lazy test, but we know that this conversion should be correct from the previous test

    ibz_copy_digits(&d1_intbig_rec, d1, d1_digits);
    ibz_copy_digits(&d2_intbig_rec, d2, d2_digits);
    ibz_copy_digits(&zero_intbig_rec, zero, 1);

    if (ibz_cmp(&d1_intbig, &d1_intbig_rec) || ibz_cmp(&zero_intbig, &zero_intbig_rec)) {
        res = 1;
        goto err;
    }

#if RADIX == 64
    const digit_t p_cofactor_for_3g[5] = {
        0x0000000000000000, 0x74f9dace0d9ec800, 0x63a25b437f655001, 0x0000000000000019, 0
    };
    digit_t p_cofactor_for_3g_rec[10] = { 0 };
#elif RADIX == 32
    const digit_t p_cofactor_for_3g[10] = { 0x00000000, 0x00000000, 0x0d9ec800, 0x74f9dace, 0x7f655001,
                                            0x63a25b43, 0x00000019, 0x00000000, 0,          0 };
    digit_t p_cofactor_for_3g_rec[5] = { 0 };
#endif
    ibz_copy_digits(&cof, p_cofactor_for_3g, 5);
    ibz_to_digits(p_cofactor_for_3g_rec, &cof);
    ibz_copy_digits(&cof2, p_cofactor_for_3g_rec, 5);

    if (ibz_cmp(&cof, &cof2)) {
        res = 1;
        goto err;
    }

#if RADIX == 32
    digit_t da[3] = { 0, 0, 0 };
#elif RADIX == 64
    digit_t da[2] = { 0, 0 };
#endif

    ibz_t strval, strval_check;
    ibz_init(&strval_check);
    ibz_init(&strval);
    ibz_set_from_str(&strval, "1617406613339667622221321", 10);
    ibz_to_digits(da, &strval);
    ibz_copy_digits(&strval_check, da, sizeof(da) / sizeof(digit_t));
    // ibz_printf("strval:       %Zd\nstrval_check: %Zd\n", strval, strval_check);
    if (ibz_cmp(&strval, &strval_check)) {
        res = 1;
        goto err;
    }

err:
    if (res) {
        printf("Quaternion module integer test ibz_test_to_digits failed\n");
    }
    ibz_finalize(&d1_intbig);
    ibz_finalize(&d2_intbig);
    ibz_finalize(&zero_intbig);
    ibz_finalize(&d1_intbig_rec);
    ibz_finalize(&d2_intbig_rec);
    ibz_finalize(&zero_intbig_rec);
    ibz_finalize(&cof);
    ibz_finalize(&cof2);
    ibz_finalize(&strval);
    ibz_finalize(&strval_check);
    return res;
}

int
ibz_test_intbig()
{
    int reps = 2;
    int prime_n = 103;
    int res = 0;
    printf("\nRunning quaternion tests of gmp-based integer functions\n");
    res = res | ibz_test_init_set_cmp();
    res = res | ibz_test_add_sub_neg_abs();
    res = res | ibz_test_mul_sqrt();
    res = res | ibz_test_div();
    res = res | ibz_test_mod();
    res = res | ibz_test_pow();
    res = res | ibz_test_gcd();
    res = res | ibz_test_sqrt_mod_p(reps, prime_n);
    res = res | ibz_test_rand_interval(reps);
    res = res | ibz_test_rand_interval_i(reps);
    res = res | ibz_test_rand_interval_minm_m(reps);
    res = res | ibz_test_copy_digits();
    res = res | ibz_test_to_digits();
    return res;
}
