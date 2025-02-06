#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include "mini-gmp-extra.h"
#include "intbig_internal.h"
#include "rng.h"

#define RANDOM_TEST_ITERS 1000

int
mini_gmp_test_mpz_legendre(void)
{
    int res = 0;
    const int levels = 3;
    const int cofactor[] = { 5, 65, 27 };
    const int e[] = { 248, 376, 500 };

    int as[] = { -2, -1, 0, 1, 2 };
    // clang-format off
    int legendre[3][5] = {
        { -1, -1, 0, 1, 1, },
        { -1, -1, 0, 1, 1, },
        { -1, -1, 0, 1, 1, },
    };
    // clang-format on

    mpz_t a, p;
    mpz_init(a);
    mpz_init(p);

    for (int i = 0; i < levels; i++) {
        // build cofactor*2^e - 1 for the i-th level
        mpz_set_ui(p, 1);
        mpz_mul_2exp(p, p, e[i]);
        mpz_mul_ui(p, p, cofactor[i]);
        mpz_sub_ui(p, p, 1);

        for (unsigned long j = 0; j < sizeof(as) / sizeof(as[0]); j++) {
            mpz_set_si(a, as[j]);
            res = res | (mini_mpz_legendre(a, p) != legendre[i][j]);
        }

#if defined(__GMP_H__)
        for (int j = 0; j < RANDOM_TEST_ITERS; j++) {
            ibz_rand_interval(&a, &ibz_const_zero, &p);
            // Compare against the full GMP implementation
            res = res | (mini_mpz_legendre(a, p) != mpz_legendre(a, p));
        }
#endif
    }

    mpz_clear(a);
    mpz_clear(p);

    if (res) {
        printf("mini-gmp test mpz_legendre failed\n");
    }

    return res;
}

typedef union
{
    double d;
    int64_t s64;
    unsigned char uc[sizeof(int64_t)];
} double_int64_uchar_t;

int
mini_gmp_test_mpz_get_d_2exp(void)
{
    int res = 0;
    signed long int e, e2 UNUSED;
    double d, d2 UNUSED;
    mpz_t op;

    mpz_init(op);

    // Test 0
    mpz_set_si(op, 0);
    d = mini_mpz_get_d_2exp(&e, op);
    res = res | (e != 0);
    res = res | (d != 0.0); // exact floating point comparison

    // Test 1
    mpz_set_si(op, 1);
    d = mini_mpz_get_d_2exp(&e, op);
    res = res | (e != 1);
    res = res | (d != 0.5); // exact floating point comparison

    // Test -1
    mpz_set_si(op, -1);
    d = mini_mpz_get_d_2exp(&e, op);
    res = res | (e != 1);
    res = res | (d != -0.5); // exact floating point comparison

    // Test a few powers of 2: 2^1, 2^2, 2^4, 2^8, ..., 2^65536, and their negatives
    for (int i = 0; i <= 16; i++) {
        mpz_set_ui(op, 1);
        mpz_mul_2exp(op, op, 1 << i);

        d = mini_mpz_get_d_2exp(&e, op);
        res = res | (e != (1 << i) + 1);
        res = res | (d != 0.5); // exact floating point comparison

        mpz_neg(op, op);
        d = mini_mpz_get_d_2exp(&e, op);
        res = res | (e != (1 << i) + 1);
        res = res | (d != -0.5); // exact floating point comparison
    }

    // Test random doubles with random exponent
    double_int64_uchar_t dd;
    volatile uint16_t ee;

    for (uint32_t i = 0; i < RANDOM_TEST_ITERS; i++) {
        randombytes(dd.uc, sizeof(dd.uc));
        randombytes((unsigned char *)&ee, sizeof(ee));

#ifdef TARGET_BIG_ENDIAN
        // Ensure reproducibility in big-endian systems
        dd.s64 = BSWAP64(dd.s64);
        ee = BSWAP16(ee);
#endif

        dd.s64 &= 0x800fffffffffffff; // clear exponent
        dd.s64 |= 0x3fe0000000000000; // set exponent to -1

        if (ee >= DBL_MANT_DIG) {
            mpz_set_d(op, dd.d * (INT64_C(1) << DBL_MANT_DIG));
            mpz_mul_2exp(op, op, ee - DBL_MANT_DIG);
        } else {
            // Since it fits in a double, round it first to ensure it's an integer
            dd.d = round(dd.d * (INT64_C(1) << ee));
            mpz_set_d(op, dd.d);
            dd.d /= INT64_C(1) << ee;
            // These cases (-1, 0, 1) were already tested, and +/- 1 would require special-casing below
            if (fabs(dd.d) <= 1.0) {
                continue;
            }
        }

        d = mini_mpz_get_d_2exp(&e, op);

        res = res | (e != ee);
        res = res | (d != dd.d);

#if defined(__GMP_H__)
        // Compare against the full GMP implementation
        d2 = mpz_get_d_2exp(&e2, op);

        res = res | (e != e2);
        res = res | (d != d2);
#endif
    }

#if defined(__GMP_H__)
    for (int exp2 = 1023; exp2 <= 1025; exp2++) {
        for (int sign_outer = -1; sign_outer <= 1; sign_outer += 2) {
            for (int sign_inner = -1; sign_inner <= 1; sign_inner += 2) {
                for (int i = 1; i < 15; i++) {
                    mpz_set_si(op, i);
                    mpz_mul_2exp(op, op, 55);
                    if (sign_inner > 0)
                        mpz_add_ui(op, op, 1);
                    else
                        mpz_sub_ui(op, op, 1);
                    mpz_mul_2exp(op, op, exp2 - 55);
                    mpz_mul_si(op, op, sign_outer);

                    d = mini_mpz_get_d_2exp(&e, op);
                    d2 = mpz_get_d_2exp(&e2, op);

                    res = res | (e != e2);
                    res = res | (d != d2);
                }
            }
        }
    }

    // Test random integers
    for (uint32_t i = 0; i < RANDOM_TEST_ITERS; i++) {
        ibz_rand_interval_bits(&op, 8 * (i + 1));

        d = mini_mpz_get_d_2exp(&e, op);
        d2 = mpz_get_d_2exp(&e2, op);

        res = res | (e != e2);
        res = res | (d != d2);
    }
#endif

    if (res) {
        printf("mini-gmp test mpz_get_d_2exp failed\n");
    }

    mpz_clear(op);

    return res;
}

int
mini_gmp_test(void)
{
    int res = 0;
    printf("\nRunning tests for implementations of GMP functions missing from the mini-GMP API\n");
    res = res | mini_gmp_test_mpz_legendre();
    res = res | mini_gmp_test_mpz_get_d_2exp();
    return res;
}
