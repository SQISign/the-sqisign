#include <inttypes.h>
#include "id2iso_tests.h"
#include <hd.h>
#include <tools.h>

static int
ec_is_on_curve(ec_point_t *P, ec_curve_t *E)
{
    fp2_t y_sq, X, tmp;

    fp2_copy(&X, &P->z);
    fp2_inv(&X);
    fp2_mul(&X, &X, &P->x);
    fp2_mul(&y_sq, &X, &X);
    fp2_mul(&y_sq, &y_sq, &X); // X^3
    fp2_copy(&tmp, &E->C);
    fp2_inv(&tmp);
    fp2_mul(&tmp, &tmp, &E->A);
    fp2_mul(&tmp, &tmp, &X);
    fp2_mul(&tmp, &tmp, &X);
    fp2_add(&y_sq, &y_sq, &tmp); // X^3 + A X^2
    fp2_add(&y_sq, &y_sq, &X);   // X^3 + A X^2 + X

    if (fp2_is_square(&y_sq)) {
        return 1;
    } else {
        return 0;
    }
}

int
dim2id2iso_test_fixed_degree_isogeny(void)
{
    int ret;
    ibz_t u, two_pow;
    ibz_t tmp;
    quat_left_ideal_t lideal;
    ibz_init(&u);
    ibz_init(&tmp);
    ibz_init(&two_pow);
    quat_left_ideal_init(&lideal);

    for (int i = 0; i <= NUM_ALTERNATE_EXTREMAL_ORDERS; i++) {
        if (EXTREMAL_ORDERS[i].q % 4 == 1) {
            // u is a random prime
            ibz_generate_random_prime(&u, 1, 110, QUAT_represent_integer_params.primality_test_iterations);

            // now we check that we get a correct codomain in the end
            // by evaluating some point and checking the pairing
            theta_couple_curve_t E00;
            ec_curve_t E0;
            copy_curve(&E0, &CURVES_WITH_ENDOMORPHISMS[i].curve);
            ec_curve_normalize_A24(&E0);
            copy_curve(&E00.E1, &E0);
            copy_curve(&E00.E2, &E0);
            ec_basis_t bas = CURVES_WITH_ENDOMORPHISMS[i].basis_even;

            assert(test_point_order_twof(&bas.P, &E00.E1, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&bas.Q, &E00.E2, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&bas.PmQ, &E00.E1, TORSION_EVEN_POWER));
            assert(ec_is_on_curve(&bas.P, &E00.E1));

            // creating points for translation
            theta_couple_point_t T1, T2, T3;
            copy_point(&T1.P1, &bas.P);
            copy_point(&T2.P1, &bas.Q);
            copy_point(&T3.P1, &bas.PmQ);
            ec_point_init(&T1.P2);
            ec_point_init(&T2.P2);
            ec_point_init(&T3.P2);

            assert(ec_is_on_curve(&T2.P1, &E00.E1));
            assert(ec_is_on_curve(&T1.P1, &E00.E1));
            assert(ec_is_on_curve(&T3.P1, &E00.E1));
            assert(ec_is_zero(&T1.P2));
            assert(ec_is_zero(&T2.P2));
            assert(ec_is_zero(&T3.P2));

            theta_couple_curve_t codomain;
            theta_couple_point_t pushed_points[3];
            theta_couple_point_t *const Tev1 = pushed_points + 0, *const Tev2 = pushed_points + 1,
                                        *const Tev3 = pushed_points + 2;
            pushed_points[0] = T1;
            pushed_points[1] = T2;
            pushed_points[2] = T3;

            // Compute the isogeny and evaluation
            ret = fixed_degree_isogeny_and_eval(
                &lideal, &u, 0, &codomain, pushed_points, sizeof(pushed_points) / sizeof(*pushed_points), i);
            assert(ret);
            unsigned int length = (unsigned int)ret;

            assert(ec_is_on_curve(&Tev1->P2, &codomain.E2));
            assert(ec_is_on_curve(&Tev2->P1, &codomain.E1));
            assert(ec_is_on_curve(&Tev1->P1, &codomain.E1));
            assert(ec_is_on_curve(&Tev3->P1, &codomain.E1));
            assert(ec_is_on_curve(&Tev2->P2, &codomain.E2));
            assert(ec_is_on_curve(&Tev3->P2, &codomain.E2));

            assert(test_point_order_twof(&Tev1->P1, &codomain.E1, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&Tev1->P2, &codomain.E2, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&Tev2->P1, &codomain.E1, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&Tev2->P2, &codomain.E2, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&Tev3->P1, &codomain.E1, TORSION_EVEN_POWER));
            assert(test_point_order_twof(&Tev3->P2, &codomain.E2, TORSION_EVEN_POWER));
            fp2_t w0, w1, w2;

            weil(&w0, TORSION_EVEN_POWER, &bas.P, &bas.Q, &bas.PmQ, &E0);
            weil(&w1, TORSION_EVEN_POWER, &Tev1->P1, &Tev2->P1, &Tev3->P1, &codomain.E1);
            weil(&w2, TORSION_EVEN_POWER, &Tev1->P2, &Tev2->P2, &Tev3->P2, &codomain.E2);
            ibz_pow(&two_pow, &ibz_const_two, length);
            ibz_sub(&two_pow, &two_pow, &u);
            // now we are checking that one of the two is equal to the correct value
            digit_t digit_u[NWORDS_ORDER] = { 0 };
            ibz_to_digit_array(digit_u, &u);
            fp2_t test_pow;
            fp2_pow_vartime(&test_pow, &w0, digit_u, NWORDS_ORDER);

            // The way we set up the 2d isogeny we should always get the first curve
            assert(fp2_is_equal(&test_pow, &w1));
            ibz_to_digit_array(digit_u, &two_pow);
            fp2_pow_vartime(&test_pow, &w0, digit_u, NWORDS_ORDER);
            assert(fp2_is_equal(&test_pow, &w2));
        }
    }

    ibz_finalize(&tmp);
    ibz_finalize(&u);
    ibz_finalize(&two_pow);
    quat_left_ideal_finalize(&lideal);

    return 1;
}

int
dim2id2iso_test_find_uv(void)
{
    // var dec
    int found = 1;
    ibz_t temp, remainder, n1, n2;
    ibz_t norm_d;
    quat_alg_elem_t gen;

    quat_left_ideal_t lideal_small;
    quat_lattice_t right_order;
    ibz_mat_4x4_t reduced, gram;
    ibz_vec_4_t coeffs;
    quat_alg_elem_t beta1, beta2;
    ibz_t u, v, d1, d2, target;
    // var init
    ibz_init(&target);
    ibz_init(&norm_d);
    ibz_init(&temp);
    ibz_init(&remainder);
    ibz_init(&n1);
    ibz_init(&n2);
    quat_alg_elem_init(&gen);
    quat_left_ideal_init(&lideal_small);
    quat_lattice_init(&right_order);
    ibz_mat_4x4_init(&reduced);
    ibz_mat_4x4_init(&gram);
    ibz_vec_4_init(&coeffs);

    quat_alg_elem_init(&beta1);
    quat_alg_elem_init(&beta2);

    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&d1);
    ibz_init(&d2);

    // computation of lideal_small
    ibz_generate_random_prime(
        &n1, 1, ibz_bitsize(&QUATALG_PINFTY.p), QUAT_represent_integer_params.primality_test_iterations);
    ibz_generate_random_prime(
        &n2, 1, ibz_bitsize(&QUATALG_PINFTY.p), QUAT_represent_integer_params.primality_test_iterations);
    ibz_mul(&temp, &n1, &n2);
    found = found && quat_represent_integer(&gen, &temp, 0, &QUAT_represent_integer_params);
    assert(found);
    quat_lideal_create(&lideal_small, &gen, &n1, &STANDARD_EXTREMAL_ORDER.order, &QUATALG_PINFTY);

    int exp = TORSION_EVEN_POWER;
    ibz_pow(&target, &ibz_const_two, exp);

    found = 0;
    int num_rerun = 0;
    int max_num_rerun = 2;
    int index_alternate_order_1;
    int index_alternate_order_2;
    while (!found && num_rerun < max_num_rerun) {
        found = find_uv(&u,
                        &v,
                        &beta1,
                        &beta2,
                        &d1,
                        &d2,
                        &index_alternate_order_1,
                        &index_alternate_order_2,
                        &target,
                        &lideal_small,
                        &QUATALG_PINFTY,
                        NUM_ALTERNATE_EXTREMAL_ORDERS);
        if (num_rerun > 0) {
            printf("alternate rerun ! \n");
        }
        num_rerun++;
    }
    // TOC_clock(t,"alternate find_uv");
    // printf("\n\n");
    // assert(found);

    if (!found) {
        printf("alternate failed \n");
    }

    ibz_finalize(&norm_d);
    ibz_finalize(&temp);
    ibz_finalize(&remainder);
    ibz_finalize(&n1);
    ibz_finalize(&n2);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal_small);
    quat_lattice_finalize(&right_order);
    ibz_mat_4x4_finalize(&reduced);
    ibz_mat_4x4_finalize(&gram);
    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&d1);
    ibz_finalize(&d2);
    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&target);

    quat_alg_elem_finalize(&beta1);
    quat_alg_elem_finalize(&beta2);

    return found;
}

int
dim2id2iso_test_dimid2iso(void)
{

    // var dec
    int found = 1;
    ibz_t temp, remainder, n1, n2;
    quat_alg_elem_t gen;

    quat_left_ideal_t lideal_small;
    quat_lattice_t right_order;
    ibz_mat_4x4_t reduced, gram;
    quat_alg_elem_t beta1, beta2;
    ibz_t u, v, d1, d2;

    // var init
    ibz_init(&temp);
    ibz_init(&remainder);
    ibz_init(&n1);
    ibz_init(&n2);
    quat_alg_elem_init(&gen);
    quat_left_ideal_init(&lideal_small);
    quat_lattice_init(&right_order);
    ibz_mat_4x4_init(&reduced);
    ibz_mat_4x4_init(&gram);

    quat_alg_elem_init(&beta1);
    quat_alg_elem_init(&beta2);

    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&d1);
    ibz_init(&d2);

    // we start by checking that the precomputation is correct, while checking at the same time our
    // algorithm
    for (int i = 0; i < NUM_ALTERNATE_EXTREMAL_ORDERS; i++) {
        ec_basis_t bas_check;
        ec_curve_t codom_check;
        quat_left_ideal_t id;
        quat_left_ideal_init(&id);
        quat_lideal_copy(&id, &ALTERNATE_CONNECTING_IDEALS[i]);
        found = found && dim2id2iso_arbitrary_isogeny_evaluation(&bas_check, &codom_check, &id);
        ec_isom_t isom;
        ec_isomorphism(&isom, &codom_check, &ALTERNATE_STARTING_CURVES[i].curve);
        ec_iso_eval(&bas_check.P, &isom);
        ec_iso_eval(&bas_check.Q, &isom);
        ec_iso_eval(&bas_check.PmQ, &isom);
        assert(test_point_order_twof(&bas_check.P, &ALTERNATE_STARTING_CURVES[i].curve, TORSION_EVEN_POWER));
        assert(test_point_order_twof(&bas_check.Q, &ALTERNATE_STARTING_CURVES[i].curve, TORSION_EVEN_POWER));
        assert(test_point_order_twof(&bas_check.PmQ, &ALTERNATE_STARTING_CURVES[i].curve, TORSION_EVEN_POWER));
        assert(ec_is_equal(&bas_check.P, &ALTERNATE_STARTING_CURVES[i].basis_even.P));
        assert(ec_is_equal(&bas_check.Q, &ALTERNATE_STARTING_CURVES[i].basis_even.Q));
        assert(ec_is_equal(&bas_check.PmQ, &ALTERNATE_STARTING_CURVES[i].basis_even.PmQ));
        quat_left_ideal_finalize(&id);
    }

    // computation of lideal_small
    ibz_generate_random_prime(
        &n1, 1, ibz_bitsize(&QUATALG_PINFTY.p), QUAT_represent_integer_params.primality_test_iterations);
    ibz_generate_random_prime(
        &n2, 1, ibz_bitsize(&QUATALG_PINFTY.p), QUAT_represent_integer_params.primality_test_iterations);
    ibz_mul(&temp, &n1, &n2);
    found = found && quat_represent_integer(&gen, &temp, 0, &QUAT_represent_integer_params);
    assert(found);
    quat_lideal_create(&lideal_small, &gen, &n1, &STANDARD_EXTREMAL_ORDER.order, &QUATALG_PINFTY);
    ec_basis_t bas_end;
    ec_curve_t codom;
    ec_curve_init(&codom);

    found = dim2id2iso_ideal_to_isogeny_clapotis(
        &beta1, &beta2, &u, &v, &d1, &d2, &codom, &bas_end, &lideal_small, &QUATALG_PINFTY);

    for (int i = 0; i < 10; i++) {
        ibz_generate_random_prime(
            &n1, 1, ibz_bitsize(&QUATALG_PINFTY.p), QUAT_represent_integer_params.primality_test_iterations);
        ibz_generate_random_prime(
            &n2, 1, ibz_bitsize(&QUATALG_PINFTY.p), QUAT_represent_integer_params.primality_test_iterations);
        ibz_mul(&temp, &n1, &n2);
        found = found && quat_represent_integer(&gen, &temp, 0, &QUAT_represent_integer_params);
        assert(found);
        quat_lideal_create(&lideal_small, &gen, &n1, &STANDARD_EXTREMAL_ORDER.order, &QUATALG_PINFTY);
        found = dim2id2iso_ideal_to_isogeny_clapotis(
            &beta1, &beta2, &u, &v, &d1, &d2, &codom, &bas_end, &lideal_small, &QUATALG_PINFTY);
    }

    ibz_finalize(&temp);
    ibz_finalize(&remainder);
    ibz_finalize(&n1);
    ibz_finalize(&n2);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal_small);
    quat_lattice_finalize(&right_order);
    ibz_mat_4x4_finalize(&reduced);
    ibz_mat_4x4_finalize(&gram);
    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&d1);
    ibz_finalize(&d2);

    quat_alg_elem_finalize(&beta1);
    quat_alg_elem_finalize(&beta2);

    return found;
}
