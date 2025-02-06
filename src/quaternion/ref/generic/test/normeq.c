#include "quaternion_tests.h"
#include <stdlib.h>
#include <assert.h>

// helpers for setting parameters

// test if parameters are such that represent_inter is likely to find a solution
int
quat_test_input_for_repres_bound(const ibz_t *p, const ibz_t *M, int q)
{
    ibz_t c, r;
    ibz_init(&c);
    ibz_init(&r);
    ibz_set(&r, q);
    ibz_mul(&r, &r, p);
    ibz_div(&c, &r, M, &r);
    ibz_set(&r, 1 << 20);
    int res = (ibz_cmp(&r, &c) <= 0);
    ibz_finalize(&c);
    ibz_finalize(&r);
    return (res);
}

// 1 if ok, 0 if error
int
quat_test_special_extremal_setup(quat_represent_integer_params_t *params, const quat_alg_t *alg)
{ // check the order is maximal and i,j generate a suborder
    int res = 1;
    ibz_vec_4_t ij;
    quat_lattice_t lat;
    quat_lattice_t test;
    ibz_vec_4_init(&ij);
    quat_lattice_init(&lat);
    quat_lattice_init(&test);
    quat_alg_coord_mul(&ij, &(params->order->z.coord), &(params->order->t.coord), alg);
    ibz_copy(&(lat.denom), &(params->order->z.denom));
    ibz_copy(&(lat.basis[0][0]), &(lat.denom));
    for (int i = 0; i < 4; i++) {
        ibz_copy(&(lat.basis[i][3]), &(ij[i]));
        ibz_mul(&(lat.basis[i][2]), &(params->order->t.coord[i]), &(lat.denom));
        ibz_copy(&(lat.basis[i][1]), &(params->order->z.coord[i]));
    }
    quat_lattice_hnf(&lat);
    quat_lattice_mul(&test, &lat, &lat, alg);
    res = res && quat_lattice_inclusion(&test, &lat);
    quat_lattice_mul(&test, &(params->order->order), &(params->order->order), alg);
    res = res && quat_lattice_inclusion(&test, &(params->order->order));
    res = res && quat_order_is_maximal(&(params->order->order), alg);
    res = res && quat_lattice_inclusion(&lat, &(params->order->order));
    ibz_vec_4_finalize(&ij);
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&test);
    return (res);
}

void
quat_test_set_params_standard(quat_represent_integer_params_t *params,
                              quat_p_extremal_maximal_order_t *order,
                              const quat_alg_t *alg)
{
    quat_lattice_O0_set_extremal(order);

    params->algebra = alg;
    params->order = order;
    assert(quat_test_special_extremal_setup(params, alg));
}

void
quat_test_set_params_non_standard(quat_represent_integer_params_t *params,
                                  quat_p_extremal_maximal_order_t *order,
                                  const quat_alg_t *alg)
{
    ibz_set_from_str(&order->z.coord[1], "214764116738303679745780048598183569015", 10);
    ibz_set(&order->z.coord[2], 0);
    ibz_set_from_str(&order->z.coord[3], "1", 10);
    ibz_set_from_str(&order->z.denom, "73403150722045993989123427738005336972", 10);
    ibz_set(&order->t.coord[2], 1);
    ibz_set(&order->t.denom, 1);
    order->q = 13;
    ibz_set_from_str(&(order->order.basis[0][0]), "73403150722045993989123427738005336972", 10);
    ibz_set_from_str(&(order->order.basis[1][1]),
                     "2694011267961700664357934052637599020390646823337886018360381743577635064392",
                     10);
    ibz_set_from_str(&(order->order.basis[2][2]), "36701575361022996994561713869002668486", 10);
    ibz_set(&(order->order.basis[3][3]), 1);
    ibz_set_from_str(&(order->order.basis[0][2]), "36701575361022996994561713869002668486", 10);
    ibz_set_from_str(&(order->order.basis[0][3]), "0", 10);
    ibz_set_from_str(&(order->order.basis[1][2]), "0", 10);
    ibz_set_from_str(&(order->order.basis[1][3]), "214764116738303679745780048598183569015", 10);
    ibz_set_from_str(&(order->order.denom), "73403150722045993989123427738005336972", 10);

    params->algebra = alg;
    params->order = order;

    assert(quat_test_special_extremal_setup(params, alg));
}

// void quat_order_elem_create(quat_alg_elem_t *elem, const quat_p_extremal_maximal_order_t *order,
// const ibz_vec_4_t *coeffs, const quat_alg_t *Bpoo);
int
quat_test_order_elem_create()
{
    int res = 0;
    quat_p_extremal_maximal_order_t O0;
    ibz_vec_4_t vec, cmp;
    quat_alg_t alg;
    quat_alg_elem_t elem;
    ibz_vec_4_init(&vec);
    ibz_vec_4_init(&cmp);
    quat_alg_elem_init(&elem);
    quat_alg_elem_init(&(O0.z));
    quat_alg_elem_init(&(O0.t));
    quat_lattice_init(&(O0.order));
    quat_alg_init_set_ui(&alg, 103);

    quat_lattice_O0_set_extremal(&O0);
    ibz_vec_4_set(&vec, 1, 7, 2, -2);
    ibz_vec_4_copy(&cmp, &vec);
    ibz_neg(&(cmp[3]), &(cmp[3]));
    quat_order_elem_create(&elem, &O0, &vec, &alg);
    res = res | ibz_cmp(&(elem.denom), &ibz_const_one);
    for (int i = 0; i < 4; i++)
        res = res | ibz_cmp(&(elem.coord[i]), &cmp[i]);

    if (res) {
        printf("Quaternion unit test order_elem_create failed\n");
    }
    ibz_vec_4_finalize(&vec);
    ibz_vec_4_finalize(&cmp);
    quat_alg_elem_finalize(&elem);
    quat_alg_elem_finalize(&(O0.z));
    quat_alg_elem_finalize(&(O0.t));
    quat_lattice_finalize(&(O0.order));
    quat_alg_finalize(&alg);
    return (res);
}

// int quat_sampling_random_ideal_O0_given_norm(quat_left_ideal_t *lideal,const ibz_t *norm,int
//  is_prime,const quat_represent_integer_params_t *params,int prime_sampling_attempts);
int
quat_test_sampling_random_ideal_O0_given_norm()
{
    int res = 0;
    quat_left_ideal_t lideal;
    quat_lattice_t test;
    quat_represent_integer_params_t params;
    quat_p_extremal_maximal_order_t O0;
    quat_alg_t alg;
    ibz_t p, norm, coprime;
    quat_lattice_init(&test);
    ibz_init(&norm);
    ibz_init(&p);
    ibz_init(&coprime);
    quat_left_ideal_init(&lideal);
    quat_alg_elem_init(&(O0.z));
    quat_alg_elem_init(&(O0.t));
    quat_lattice_init(&(O0.order));
    params.primality_test_iterations = 30;
    ibz_set(&p, 1533069323);
    quat_alg_init_set(&alg, &p);
    quat_test_set_params_standard(&params, &O0, &alg);
    ibz_set(&p, 1338708463);
    res = res || !quat_sampling_random_ideal_O0_given_norm(&lideal, &p, 1, &params, NULL);
    res = res || (ibz_cmp(&(lideal.norm), &p) != 0);
    quat_lideal_norm(&lideal);
    res = res || (ibz_cmp(&(lideal.norm), &p) != 0);
    res = res || !ibz_mat_4x4_equal(&(O0.order.basis), &(lideal.parent_order->basis));
    res = res || (0 != ibz_cmp(&(O0.order.denom), &(lideal.parent_order->denom)));
    quat_lattice_mul(&test, &(lideal.lattice), &(lideal.lattice), &alg);
    res = res || !quat_lattice_inclusion(&test, &(lideal.lattice));

    ibz_set(&p, 3093 * 59471);
    ibz_set(&coprime, 1533069337);
    res = res || !quat_sampling_random_ideal_O0_given_norm(&lideal, &p, 0, &params, &coprime);
    res = res || (ibz_cmp(&(lideal.norm), &p) != 0);
    quat_lideal_norm(&lideal);
    res = res || (ibz_cmp(&(lideal.norm), &p) != 0);
    res = res || !ibz_mat_4x4_equal(&(O0.order.basis), &(lideal.parent_order->basis));
    res = res || (0 != ibz_cmp(&(O0.order.denom), &(lideal.parent_order->denom)));
    quat_lattice_mul(&test, &(lideal.lattice), &(lideal.lattice), &alg);
    res = res || !quat_lattice_inclusion(&test, &(lideal.lattice));

    if (res) {
        printf("Quaternion unit test sampling_random_ideal_O0_given_norm failed\n");
    }
    quat_left_ideal_finalize(&lideal);
    quat_alg_elem_finalize(&(O0.z));
    quat_alg_elem_finalize(&(O0.t));
    quat_lattice_finalize(&(O0.order));
    ibz_finalize(&norm);
    ibz_finalize(&p);
    ibz_finalize(&coprime);
    quat_lattice_finalize(&test);
    quat_alg_finalize(&alg);
    return (res);
}

// void quat_change_to_O0_basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el);
int
quat_test_change_to_O0_basis()
{
    int res = 0;
    quat_alg_elem_t cmp, elem;
    ibz_vec_4_t out;
    quat_lattice_t O0;
    quat_alg_elem_init(&elem);
    quat_alg_elem_init(&cmp);
    ibz_vec_4_init(&out);
    quat_lattice_init(&O0);
    quat_lattice_O0_set(&O0);

    quat_alg_elem_set(&cmp, 2, 2, 7, 1, -4);
    quat_alg_elem_copy(&elem, &cmp);
    quat_change_to_O0_basis(&out, &elem);
    res = res || !quat_alg_elem_equal(&elem, &cmp);
    quat_alg_elem_set(&elem, 1, 0, 0, 0, 0);
    ibz_mat_4x4_eval(&(elem.coord), &(O0.basis), &out);
    ibz_copy(&(elem.denom), &(O0.denom));
    res = res || !quat_alg_elem_equal(&elem, &cmp);

    quat_alg_elem_set(&cmp, 2, 1, 0, 6, -3);
    quat_alg_elem_copy(&elem, &cmp);
    quat_change_to_O0_basis(&out, &elem);
    res = res || !quat_alg_elem_equal(&elem, &cmp);
    quat_alg_elem_set(&elem, 1, 0, 0, 0, 0);
    ibz_mat_4x4_eval(&(elem.coord), &(O0.basis), &out);
    ibz_copy(&(elem.denom), &(O0.denom));
    res = res || !quat_alg_elem_equal(&elem, &cmp);

    quat_alg_elem_set(&cmp, 1, -8, 2, 1, 3);
    quat_alg_elem_copy(&elem, &cmp);
    quat_change_to_O0_basis(&out, &elem);
    res = res || !quat_alg_elem_equal(&elem, &cmp);
    quat_alg_elem_set(&elem, 1, 0, 0, 0, 0);
    ibz_mat_4x4_eval(&(elem.coord), &(O0.basis), &out);
    ibz_copy(&(elem.denom), &(O0.denom));
    res = res || !quat_alg_elem_equal(&elem, &cmp);

    if (res) {
        printf("Quaternion unit test change_to_O0_basis failed");
    }
    quat_alg_elem_finalize(&elem);
    quat_alg_elem_finalize(&cmp);
    ibz_vec_4_finalize(&out);
    quat_lattice_finalize(&O0);
    return (res);
}

int
quat_test_represent_integer_internal(int gamma_iter,
                                     int prim_iter,
                                     int nbits,
                                     int randomized,
                                     int standard,
                                     int non_diag)
{
    int tested = 0;
    int rand_ret = 1;
    int res = 0;
    ibz_t p, norm_n, M, norm_d;
    ibz_vec_4_t coord;
    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    // must initialize special extremal order, product and primes
    quat_represent_integer_params_t params;
    quat_p_extremal_maximal_order_t order;
    quat_lattice_init(&order.order);
    quat_alg_elem_init(&order.z);
    quat_alg_elem_init(&order.t);
    quat_alg_elem_init(&gamma);
    ibz_vec_4_init(&coord);
    ibz_init(&M);
    ibz_init(&p);
    ibz_init(&norm_d);
    ibz_init(&norm_n);
    quat_alg_init_set_ui(&Bpoo, 101);

    // setup
    params.primality_test_iterations = gamma_iter;

    while (!tested) {
        if (standard) {
            ibz_set(&p, 0);
            rand_ret = !ibz_generate_random_prime(&p, 1, 2 * nbits, 32);
            rand_ret = rand_ret || !ibz_rand_interval_bits(&M, 2 * nbits + 24 + 100 * non_diag);
            ibz_copy(&(Bpoo.p), &p);
            quat_test_set_params_standard(&params, &order, &Bpoo);
        } else {
            ibz_set_from_str(&p, "23920667128620486487914848107166358953830561597426178123910317653495243603967", 10);
            int real_nbits = ibz_bitsize(&p);
            rand_ret = !ibz_rand_interval_bits(&M, real_nbits + 24 + 100 * non_diag);
            ibz_copy(&(Bpoo.p), &p);
            quat_test_set_params_non_standard(&params, &order, &Bpoo);
        }
        // make test more realistic by keeping some properties always given in SQIsign
        tested = (ibz_get(&p) % 8 == 7);
        tested = tested && (ibz_get(&M) % 2 == 1);
        tested = tested && quat_test_input_for_repres_bound(&(Bpoo.p), &M, params.order->q);
        if (rand_ret)
            break;
    }

    if (!rand_ret) {
        res = !quat_represent_integer(&gamma, &M, non_diag, &params);
        if (!res) {
            quat_alg_norm(&norm_n, &norm_d, &gamma, &Bpoo);
            res = res || !ibz_is_one(&norm_d);
            res = res || !(ibz_cmp(&norm_n, &M) == 0);
            res = res || !quat_lattice_contains(NULL, &(params.order->order), &gamma);
            if (non_diag) {
                if (standard) {
                    ibz_mul(&norm_n, &ibz_const_two, &ibz_const_two);
                    // add not sub since basis in quat_order_elem_create is 1,i,j,-ij for O0
                    ibz_add(&norm_d, &(gamma.coord[0]), &(gamma.coord[3]));
                    ibz_mod(&norm_d, &norm_d, &norm_n);
                    res = res || (0 != ibz_cmp(&ibz_const_two, &norm_d));
                    ibz_sub(&norm_d, &(gamma.coord[1]), &(gamma.coord[2]));
                    ibz_mod(&norm_d, &norm_d, &norm_n);
                    res = res || (0 != ibz_cmp(&ibz_const_two, &norm_d));
                } else {
                    quat_lattice_contains(&coord, &(params.order->order), &gamma);
                    ibz_gcd(&norm_d, &(coord[1]), &(coord[2]));
                    ibz_gcd(&norm_d, &norm_d, &(coord[3]));
                    res = res || ibz_is_even(&norm_d);
                }
            }
        }
    } else {
        printf("Randomness failure in quat_represent_integer test\n");
    }

    quat_alg_finalize(&Bpoo);
    ibz_finalize(&p);
    ibz_finalize(&M);
    ibz_finalize(&norm_n);
    ibz_finalize(&norm_d);
    ibz_vec_4_finalize(&coord);
    quat_lattice_finalize(&order.order);
    quat_alg_elem_finalize(&order.z);
    quat_alg_elem_finalize(&order.t);
    quat_alg_elem_finalize(&gamma);
    return res | (rand_ret);
}

// int quat_represent_integer(quat_alg_elem_t *gamma, const ibz_t *n_gamma, const
// quat_represent_integer_params_t *params);
int
quat_test_represent_integer(void)
{
    int res = 0;
    int prim_iter = 32;
    int gamma_iter = 16384;

    res = res | quat_test_represent_integer_internal(prim_iter, gamma_iter, 64, 1, 1, 0);
    res = res | quat_test_represent_integer_internal(prim_iter, gamma_iter, 64, 1, 1, 1);

    gamma_iter = 32768;
    res = res | quat_test_represent_integer_internal(prim_iter, gamma_iter, 128, 1, 1, 0);
    res = res | quat_test_represent_integer_internal(prim_iter, gamma_iter, 128, 1, 1, 1);

    gamma_iter = 100000;
    res = res | quat_test_represent_integer_internal(prim_iter, gamma_iter, 128, 1, 0, 0);
    res = res | quat_test_represent_integer_internal(prim_iter, gamma_iter, 128, 1, 0, 1);

    if (res) {
        printf("Quaternion unit test represent_integer failed\n");
    }
    return (res);
}

int
quat_test_normeq(void)
{

    int res = 0;
    printf("\nRunning quaternion tests of functions for special extremal orders\n");

    res = res | quat_test_change_to_O0_basis();
    res = res | quat_test_order_elem_create();
    res = res | quat_test_sampling_random_ideal_O0_given_norm();
    res = res | quat_test_represent_integer();

    return res;
}
