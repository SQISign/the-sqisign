#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <tools.h>
#include <mp.h>
#include "biextension.h"
#include <rng.h>
#include <bench_test_arguments.h>

void
fp2_exp_2e(fp2_t *r, uint32_t e, const fp2_t *x)
{
    fp2_copy(r, x);
    for (uint32_t i = 0; i < e; i++) {
        fp2_sqr(r, r);
    }
}

void
biextension_test()
{
    clock_t t;
    ec_curve_t curve;
    ec_basis_t even_torsion;
    uint32_t e = TORSION_EVEN_POWER;
    fp2_t one, r1, rr1, rrr1, r2, r3, tp;
    ec_point_t P, Q, PmQ, A24;
    ec_point_t tmp, tmp2, PQ, PP, QQ, PPQ, PQQ, PPP, QQQ, PPPQ, PQQQ;

    // Get constants form curve E6 : y^2 = x^3 + 6*x^2 + x
    ec_curve_init(&curve);
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);
    copy_point(&A24, &curve.A24);

    // Compute 2^e torsion on curve
    (void)ec_curve_to_basis_2f_to_hint(&even_torsion, &curve, e);
    copy_point(&P, &even_torsion.P);
    copy_point(&Q, &even_torsion.Q);
    copy_point(&PmQ, &even_torsion.PmQ);

    printf("Testing order of points\n");
    t = tic();
    ec_dbl_iter(&tmp, e, &P, &curve);
    TOC_clock(t, "Doublings");
    assert(ec_is_zero(&tmp));
    ec_dbl_iter(&tmp, e, &Q, &curve);
    assert(ec_is_zero(&tmp));
    ec_dbl_iter(&tmp, e, &PmQ, &curve);
    assert(ec_is_zero(&tmp));

    printf("Computing Weil pairing\n");
    xADD(&PQ, &P, &Q, &PmQ);
    t = tic();

    weil(&r1, e, &P, &Q, &PQ, &curve);
    TOC_clock(t, "Weil pairing");

    printf("Computing Tate pairing\n");
    t = tic();

    reduced_tate(&tp, e, &P, &Q, &PQ, &curve);
    TOC_clock(t, "Tate pairing");

    printf("Testing order of Weil pairing\n");
    fp2_set_one(&one);
    fp2_exp_2e(&r2, e - 1, &r1);
    assert(!fp2_is_equal(&r2, &one));
    fp2_exp_2e(&r2, e, &r1);
    assert(fp2_is_equal(&r2, &one));

    printf("Testing order of Tate pairing\n");
    fp2_set_one(&one);
    fp2_exp_2e(&r2, e - 1, &tp);
    assert(!fp2_is_equal(&r2, &one));
    fp2_exp_2e(&r2, e, &tp);
    assert(fp2_is_equal(&r2, &one));

    printf("Bilinearity tests\n");
    weil(&r2, e, &P, &Q, &PmQ, &curve);
    fp2_inv(&r2);
    assert(fp2_is_equal(&r1, &r2));

    xDBL_A24(&PP, &P, &A24, false);
    xDBL_A24(&QQ, &Q, &A24, false);
    xADD(&PPQ, &PQ, &P, &Q);
    xADD(&PQQ, &PQ, &Q, &P);

    weil(&r2, e, &PP, &Q, &PPQ, &curve);
    weil(&r3, e, &P, &QQ, &PQQ, &curve);
    assert(fp2_is_equal(&r2, &r3));
    fp2_sqr(&rr1, &r1);
    assert(fp2_is_equal(&rr1, &r2));

    xADD(&PPP, &PP, &P, &P);
    xADD(&QQQ, &QQ, &Q, &Q);
    xADD(&PPPQ, &PPQ, &P, &PQ);
    xADD(&PQQQ, &PQQ, &Q, &PQ);
    weil(&r2, e, &PPP, &Q, &PPPQ, &curve);
    weil(&r3, e, &P, &QQQ, &PQQQ, &curve);
    assert(fp2_is_equal(&r2, &r3));
    fp2_mul(&rrr1, &rr1, &r1);
    assert(fp2_is_equal(&rrr1, &r2));

    printf("dlog tests\n");
    ec_basis_t BPQ, BRS;
    digit_t scal_r1[NWORDS_ORDER] = { 0 };
    digit_t scal_r2[NWORDS_ORDER] = { 0 };
    digit_t scal_s1[NWORDS_ORDER] = { 0 };
    digit_t scal_s2[NWORDS_ORDER] = { 0 };
    digit_t scal_d1[NWORDS_ORDER] = { 0 };
    digit_t scal_d2[NWORDS_ORDER] = { 0 };

    // original even torsion
    BPQ = even_torsion;
    BRS = even_torsion;

    // alternative torsion, just mix the points up a little...
    // not filling top word so the addition below can overflow into it
    // so the scalars are "random enough" but we still keep the difference
    // scal_d1 and scal_d2 required to compute the right multiple of RmS
    randombytes((unsigned char *)scal_d1, (NWORDS_ORDER - 1) * sizeof(digit_t));
    randombytes((unsigned char *)scal_d2, (NWORDS_ORDER - 1) * sizeof(digit_t));
    randombytes((unsigned char *)scal_s1, (NWORDS_ORDER - 1) * sizeof(digit_t));
    randombytes((unsigned char *)scal_s2, (NWORDS_ORDER - 1) * sizeof(digit_t));

    // Ensure that r1*s2 - r2*s1 is odd such that the matrix
    // [[r1, r2], [s1, s2]] is invertible
    scal_s1[0] = (scal_s1[0] & ((digit_t)(-1) - 1)) + 1; // s1 needs to be odd
    scal_d1[0] = (scal_d1[0] & ((digit_t)(-1) - 1));     // d1 needs to be even to make r1 odd
    scal_s2[0] = (scal_s2[0] & ((digit_t)(-1) - 1)) + 1; // s2 needs to be odd
    scal_d2[0] = (scal_d2[0] & ((digit_t)(-1) - 1)) + 1; // d2 needs to be odd to make r2 even

    // Compute r1 and r2 from the difference di = ri - si
    mp_add(scal_r1, scal_d1, scal_s1, NWORDS_ORDER);
    mp_add(scal_r2, scal_d2, scal_s2, NWORDS_ORDER);

    ec_biscalar_mul(&BRS.P, scal_r1, scal_r2, e, &BPQ, &curve);
    ec_biscalar_mul(&BRS.Q, scal_s1, scal_s2, e, &BPQ, &curve);
    ec_biscalar_mul(&BRS.PmQ, scal_d1, scal_d2, e, &BPQ, &curve);

    printf("mixed\n");

    // Now solve the discrete log
    ec_dlog_2_weil(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e);

    // assert everything matches
    // R = [r1]P + [r2]Q
    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.P));

    // S = [s1]P + [s2]Q
    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.Q));

    printf("weil solved\n");

    // now repeat using the tate pairing
    ec_dlog_2_tate(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e);

    // assert everything matches
    // R = [r1]P + [r2]Q
    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.P));

    // S = [s1]P + [s2]Q
    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BPQ, &curve);
    assert(ec_is_equal(&tmp, &BRS.Q));

    printf("tate solved\n");

    // now we try with bases for partial torsion E[2^e] with e < e_full
    int e_full = TORSION_EVEN_POWER;
    int e_partial = 126;

    ec_dbl_iter(&BRS.P, e_full - e_partial, &BRS.P, &curve);
    ec_dbl_iter(&BRS.Q, e_full - e_partial, &BRS.Q, &curve);
    ec_dbl_iter(&BRS.PmQ, e_full - e_partial, &BRS.PmQ, &curve);

    ec_dlog_2_tate(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e_partial);

    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BPQ, &curve);
    ec_dbl_iter(&tmp, e_full - e_partial, &tmp, &curve);
    assert(ec_is_equal(&tmp, &BRS.P));

    // S = [s1]P + [s2]Q
    // then S = [2^e_diff] S
    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BPQ, &curve);
    ec_dbl_iter(&tmp, e_full - e_partial, &tmp, &curve);
    assert(ec_is_equal(&tmp, &BRS.Q));

    printf("tate from full basis solved\n");

    ec_dlog_2_tate(scal_r1, scal_r2, scal_s1, scal_s2, &BPQ, &BRS, &curve, e_partial);
    mp_invert_matrix(scal_r1, scal_r2, scal_s1, scal_s2, e_partial, NWORDS_ORDER);

    // assert everything matches
    ec_biscalar_mul(&tmp, scal_r1, scal_r2, e, &BRS, &curve);
    ec_dbl_iter(&tmp2, e_full - e_partial, &BPQ.P, &curve);
    assert(ec_is_equal(&tmp, &tmp2));

    ec_biscalar_mul(&tmp, scal_s1, scal_s2, e, &BRS, &curve);
    ec_dbl_iter(&tmp2, e_full - e_partial, &BPQ.Q, &curve);
    assert(ec_is_equal(&tmp, &tmp2));

    printf("tate to full basis solved\n");
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }
    }

    if (help) {
        printf("Usage: %s [--seed=<seed>]\n", argv[0]);
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

    printf("Running biextension unit tests\n");

    biextension_test();

    // Failures will be caught by asserts in biextension_test
    printf("\nAll tests passed!\n");

    return 0;
}
