#include <assert.h>
#include <stdio.h>
#include <inttypes.h>

#include "test_extras.h"
#include <ec.h>
#include <isog.h>
#include <rng.h>
#include <bench_test_arguments.h>

/******************************
Test functions
******************************/

int
test_xDBL_xADD(const ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_point_t P, Q, PQ, R1, R2;

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);
        projective_difference_point(&PQ, &P, &Q, curve);

        // 2(P + Q) = 2P + 2Q
        xADD(&R1, &P, &Q, &PQ);
        ec_dbl(&R1, &R1, curve);
        ec_dbl(&P, &P, curve);
        ec_dbl(&Q, &Q, curve);
        ec_dbl(&PQ, &PQ, curve);
        xADD(&R2, &P, &Q, &PQ);
        if (!ec_is_equal(&R1, &R2)) {
            printf("Failed 2(P + Q) = 2P + 2Q\n");
            return 1;
        }

        // (P+Q) + (P-Q) = 2P
        xADD(&R1, &P, &Q, &PQ);
        ec_dbl(&Q, &Q, curve);
        xADD(&R1, &R1, &PQ, &Q);
        ec_dbl(&P, &P, curve);
        ec_dbl(&PQ, &PQ, curve);
        if (!ec_is_equal(&R1, &P)) {
            printf("Failed (P+Q) + (P-Q) = 2P\n");
            return 1;
        }
    }

    return 0;
}

int
test_xDBLADD(const ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_point_t P, Q, PQ, R1, R2;

    ec_point_t A24;
    AC_to_A24(&A24, curve);

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);
        projective_difference_point(&PQ, &P, &Q, curve);

        xDBLADD(&R1, &R2, &P, &Q, &PQ, &A24, false);
        xADD(&PQ, &P, &Q, &PQ);
        if (!ec_is_equal(&R2, &PQ)) {
            printf("Failed addition in xDBLADD\n");
            return 1;
        }
        ec_dbl(&P, &P, curve);
        if (!ec_is_equal(&R1, &P)) {
            printf("Failed doubling in xDBLADD\n");
            return 1;
        }
    }
    return 0;
}

int
test_xDBL_variants(ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;
    ec_curve_t E;
    ec_point_t P, R1, R2, R3, R4;
    ec_point_t A24, A24norm;
    fp2_t z;

    AC_to_A24(&A24, curve);
    copy_point(&A24norm, &A24);
    ec_normalize_point(&A24norm);

    // Randomize projective representation
    copy_curve(&E, curve);
    fp2_random_test(&z);
    fp2_mul(&(E.A24.x), &(A24.x), &z);
    fp2_mul(&(E.A24.z), &(A24.z), &z);
    E.is_A24_computed_and_normalized = false;

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        xDBL(&R1, &P, (const ec_point_t *)curve);
        xDBL_A24(&R2, &P, &(E.A24), false);
        xDBL_A24(&R3, &P, &A24norm, true);
        xDBL_E0(&R4, &P);
        if (!ec_is_equal(&R1, &R2)) {
            printf("xDBL and xDBL_A24 dont match\n");
            return 1;
        }
        if (!ec_is_equal(&R1, &R3)) {
            printf("xDBL and xDBL_A24 normalized dont match\n");
            return 1;
        }
        if (!ec_is_equal(&R1, &R4)) {
            printf("xDBL and xDBL_E0 dont match\n");
            return 1;
        }
    }
    return 0;
}

int
test_zero_identities(ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_point_t P, Q, R, ec_zero;

    fp2_set_one(&(P.x));
    fp2_set_zero(&(P.z));

    fp2_set_one(&(ec_zero.x));
    fp2_set_zero(&(ec_zero.z));

    assert(ec_is_zero(&P));

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);

        xADD(&R, &ec_zero, &ec_zero, &ec_zero);
        if (!ec_is_zero(&R)) {
            printf("Failed 0 + 0 = 0\n");
            return 1;
        }

        ec_dbl(&R, &P, curve);
        xADD(&R, &P, &P, &R);
        if (!ec_is_zero(&R)) {
            printf("Failed P - P = 0\n");
            return 1;
        }

        ec_dbl(&R, &ec_zero, curve);
        if (!ec_is_zero(&R)) {
            printf("Failed 2*0 = 0\n");
            return 1;
        }

        xADD(&R, &P, &ec_zero, &P);
        if (!ec_is_equal(&R, &P)) {
            printf("Failed P + 0 = P\n");
            return 1;
        }
        xADD(&R, &ec_zero, &P, &P);
        if (!ec_is_equal(&R, &P)) {
            printf("Failed P + 0 = P\n");
            return 1;
        }

        xDBLADD(&R, &Q, &P, &ec_zero, &P, &curve->A24, false);
        if (!ec_is_equal(&Q, &P)) {
            printf("Failed P + 0 = P in xDBLADD\n");
            return 1;
        }
        xDBLADD(&R, &Q, &ec_zero, &P, &P, &curve->A24, false);
        if (!ec_is_equal(&Q, &P)) {
            printf("Failed P + 0 = P in xDBLADD\n");
            return 1;
        }
        if (!ec_is_zero(&R)) {
            printf("Failed 2*0 = 0 in xDBLADD\n");
            return 1;
        }
    }
    return 0;
}

int
test_jacobian(const ec_curve_t *curve, unsigned int Ntest)
{
    unsigned int i;

    ec_point_t P, Q;
    jac_point_t R, S, T, U, jac_zero;
    fp2_t t0, t1;

    jac_init(&jac_zero);

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_normalize_point(&P);
        ec_random_test(&Q, curve);
        ec_normalize_point(&Q);

        /* Convert to Jacobian coordinates. */
        fp2_copy(&(S.x), &(P.x));
        ec_recover_y(&(S.y), &(S.x), curve);
        fp2_set_one(&(S.z));
        fp2_copy(&(T.x), &(Q.x));
        ec_recover_y(&(T.y), &(T.x), curve);
        fp2_set_one(&(T.z));

        ADD(&R, &jac_zero, &jac_zero, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed 0 + 0 = 0 in jac\n");
            return 1;
        }

        DBL(&R, &jac_zero, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed 2*0 = 0 in jac\n");
            return 1;
        }

        jac_neg(&R, &S);
        ADD(&R, &S, &R, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed P - P = 0 in jac\n");
            return 1;
        }

        ADD(&R, &S, &jac_zero, curve);
        if (!jac_is_equal(&R, &S)) {
            printf("Failed P + 0 = P in jac\n");
            return 1;
        }
        ADD(&R, &jac_zero, &S, curve);
        if (!jac_is_equal(&R, &S)) {
            printf("Failed P + 0 = P in jac\n");
            return 1;
        }
        ADD(&R, &S, &jac_zero, curve);
        if (!jac_is_equal(&R, &S)) {
            printf("Failed 0 + P = P in jac\n");
            return 1;
        }

        DBL(&R, &S, curve);
        ADD(&U, &S, &S, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed P + P = 2*P in jac\n");
            return 1;
        }

        ADD(&R, &T, &S, curve);
        ADD(&T, &S, &T, curve);
        if (!jac_is_equal(&R, &T)) {
            printf("Failed P + Q = Q + P in jac\n");
            return 1;
        }
        ADD(&R, &T, &S, curve);
        ADD(&U, &S, &T, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed P + Q = Q + P in jac\n");
            return 1;
        }

        // Double R to make it different than (T + S).
        DBL(&R, &R, curve);
        ADD(&U, &S, &T, curve);
        ADD(&U, &U, &R, curve);
        ADD(&R, &R, &T, curve);
        ADD(&R, &R, &S, curve);
        if (!jac_is_equal(&R, &U)) {
            printf("Failed (P + Q) + R = P + (Q + R) in jac\n");
            return 1;
        }

        jac_to_ws(&R, &t0, &t1, &jac_zero, curve);
        jac_from_ws(&R, &R, &t1, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed converting to Weierstrass\n");
            return 1;
        }

        jac_to_ws(&R, &t0, &t1, &S, curve);
        jac_from_ws(&R, &R, &t1, curve);
        if (!jac_is_equal(&S, &R)) {
            printf("Failed converting to Weierstrass\n");
            return 1;
        }
        DBL(&S, &S, curve);
        jac_to_ws(&R, &t0, &t1, &S, curve);
        jac_from_ws(&R, &R, &t1, curve);
        if (!jac_is_equal(&S, &R)) {
            printf("Failed converting to Weierstrass\n");
            return 1;
        }

        jac_to_ws(&R, &t0, &t1, &jac_zero, curve);
        DBLW(&R, &t0, &R, &t0);
        jac_from_ws(&R, &R, &t1, curve);
        if (!jac_is_equal(&R, &jac_zero)) {
            printf("Failed 2*0 = 0 in Weierstrass\n");
            return 1;
        }
        jac_to_ws(&R, &t0, &t1, &S, curve);
        DBLW(&R, &t0, &R, &t0);
        jac_from_ws(&R, &R, &t1, curve);
        DBL(&S, &S, curve);
        if (!jac_is_equal(&S, &R)) {
            printf("Failed doubling in Weierstrass\n");
            return 1;
        }
    }
    return 0;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 100 * SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;
    int res = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for testing; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    // Curve A=6
    ec_curve_t curve;
    ec_curve_init(&curve);
    fp2_set_small(&(curve.A), 0);
    fp2_set_small(&(curve.C), 1);
    // fp2_random_test(&(curve.C));
    // fp2_mul(&(curve.A), &(curve.A), &(curve.C));
    ec_curve_normalize_A24(&curve);

    res |= test_xDBL_xADD(&curve, iterations);
    res |= test_xDBLADD(&curve, iterations);
    res |= test_xDBL_variants(&curve, iterations);
    res |= test_zero_identities(&curve, iterations);
    res |= test_jacobian(&curve, iterations);

    fp2_random_test(&(curve.C));
    fp2_mul(&(curve.A), &(curve.A), &(curve.C));
    ec_curve_normalize_A24(&curve);

    res |= test_xDBL_xADD(&curve, iterations);
    res |= test_xDBLADD(&curve, iterations);
    res |= test_xDBL_variants(&curve, iterations);
    res |= test_zero_identities(&curve, iterations);
    res |= test_jacobian(&curve, iterations);

    if (res) {
        printf("Tests failed!\n");
    } else {
        printf("All ec arithmetic tests passed.\n");
    }

    return res;
}
