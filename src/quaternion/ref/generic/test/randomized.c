#include "quaternion_tests.h"
#include <rng.h>

// int ibz_mat_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);
int
quat_test_randomized_ibz_mat_2x2_inv_mod(int bitsize_matrix, int bitsize_modulus, int iterations)
{
    int randret = 0;
    int res = 0;
    ibz_t m, det, gcd;
    ibz_mat_2x2_t a, inv, id, prod;
    ibz_init(&m);
    ibz_init(&det);
    ibz_init(&gcd);
    ibz_mat_2x2_init(&a);
    ibz_mat_2x2_init(&inv);
    ibz_mat_2x2_init(&prod);
    ibz_mat_2x2_init(&id);
    ibz_mat_2x2_set(&id, 1, 0, 0, 1);

    for (int iter = 0; iter < iterations; iter++) {
        // generate random matrix and modulo, with modulo larger than 2
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                randret = randret | !ibz_rand_interval_bits(&(a[i][j]), bitsize_matrix);
            }
        }
        randret = randret | !ibz_rand_interval_bits(&m, bitsize_modulus);
        ibz_abs(&m, &m);
        ibz_add(&m, &m, &ibz_const_two);
        if (randret != 0)
            goto fin;

        // compute det
        ibz_mat_2x2_det_from_ibz(&det, &(a[0][0]), &(a[0][1]), &(a[1][0]), &(a[1][1]));
        // is it prime to mod
        ibz_gcd(&gcd, &det, &m);
        if (ibz_is_one(&gcd)) {
            // matrix should be invertible mod m
            if (ibz_mat_2x2_inv_mod(&inv, &a, &m)) {
                // ibz_2x2_mul_mod(&prod,&a,&inv, &m);
                ibz_2x2_mul_mod(&prod, &inv, &a, &m);
                res = res || ibz_cmp(&(prod[0][0]), &(id[0][0]));
                res = res || ibz_cmp(&(prod[0][1]), &(id[0][1]));
                res = res || ibz_cmp(&(prod[1][0]), &(id[1][0]));
                res = res || ibz_cmp(&(prod[1][1]), &(id[1][1]));
            } else {
                res = 1;
            }
        } else {
            res = res || ibz_mat_2x2_inv_mod(&inv, &a, &m);
        }
    }

fin:;
    if (randret != 0) {
        printf("Randomness failed in quaternion unit test with randomization for "
               "ibz_mat_2x2_inv_mod\n");
    }
    if (res != 0) {
        printf("Quaternion unit test with randomization for ibz_mat_2x2_inv_mod failed\n");
    }
    ibz_mat_2x2_finalize(&a);
    ibz_mat_2x2_finalize(&inv);
    ibz_mat_2x2_finalize(&prod);
    ibz_mat_2x2_finalize(&id);
    ibz_finalize(&m);
    ibz_finalize(&det);
    ibz_finalize(&gcd);
    return (res);
}

// int ibz_4x4_inv_with_det_as_denom(ibz_mat_4x4_t *inv, ibz_t *det, const ibz_mat_4x4_t mat);
int
quat_test_randomized_ibz_mat_4x4_inv_with_det_as_denom(int matrix_bitsize, int iterations)
{
    int res = 0;
    int randret = 0;
    ibz_t det;
    ibz_mat_4x4_t mat, inv;
    ibz_mat_4x4_t det_id, prod;
    ibz_mat_4x4_init(&det_id);
    ibz_mat_4x4_init(&prod);
    ibz_init(&det);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&inv);

    for (int r = 0; r < iterations; r++) {
        do {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++) {
                    randret = randret | !ibz_rand_interval_bits(&mat[i][j], matrix_bitsize);
                    if (randret != 0)
                        goto fin;
                }
        } while (!ibz_mat_4x4_inv_with_det_as_denom(&inv, &det, &mat));
        ibz_mat_4x4_identity(&det_id);
        ibz_mat_4x4_scalar_mul(&det_id, &det, &det_id);
        ibz_mat_4x4_mul(&prod, &inv, &mat);
        res = res || !ibz_mat_4x4_equal(&det_id, &prod);
        ibz_mat_4x4_mul(&prod, &mat, &inv);
        res = res || !ibz_mat_4x4_equal(&det_id, &prod);
    }

fin:;
    if (randret != 0) {
        printf("Randomness failed in quaternion unit test with randomization for "
               "ibz_mat_2x2_inv_mod\n");
    }
    if (res != 0) {
        printf("Quaternion unit test with randomization for ibz_mat_4x4_inv_with_det_as_denom "
               "failed\n");
    }
    ibz_mat_4x4_finalize(&det_id);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&inv);
    ibz_finalize(&det);
    return res;
}

// int quat_lattice_contains(ibz_vec_4_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t
// *x);
int
quat_test_randomized_lattice_contains(int lattice_bitsize, int coord_bitsize, int iterations)
{
    // only tests the case where the element is in the lattice
    int res = 0;
    int randret = 0;
    ibz_t det;
    quat_alg_elem_t x, cmp;
    ibz_vec_4_t coord, set_coord;
    quat_lattice_t lat;
    ibz_init(&det);
    ibz_vec_4_init(&coord);
    ibz_vec_4_init(&set_coord);
    quat_alg_elem_init(&cmp);
    quat_alg_elem_init(&x);
    quat_lattice_init(&lat);

    for (int iter = 0; iter < iterations; iter++) {
        randret = quat_test_input_random_lattice_generation(&lat, lattice_bitsize, 1, 1);
        for (int i = 0; i < 4; i++) {
            randret = randret | !ibz_rand_interval_bits(&(set_coord[i]), coord_bitsize);
        }
        if (randret != 0)
            goto fin;

        ibz_mat_4x4_eval(&(x.coord), &(lat.basis), &set_coord);
        ibz_copy(&(x.denom), &(lat.denom));
        ibz_vec_4_set(&coord, 1, 0, 1, 0);
        if (quat_lattice_contains(&coord, &lat, &x)) {
            ibz_mat_4x4_eval(&(cmp.coord), &(lat.basis), &coord);
            ibz_copy(&(cmp.denom), &(lat.denom));
            quat_alg_sub(&cmp, &x, &cmp);
            res = res || !quat_alg_elem_is_zero(&cmp);
            ibz_vec_4_sub(&coord, &coord, &set_coord);
            res = res || !ibz_vec_4_is_zero(&coord);
        } else {
            res = 1;
        }
    }
fin:;
    if (randret != 0) {
        printf("Randomness failed in quaternion unit test with randomization for lattice_contains\n");
    }
    if (res != 0) {
        printf("Quaternion unit test with randomization for lattice_contains failed\n");
    }
    ibz_finalize(&det);
    ibz_vec_4_finalize(&coord);
    ibz_vec_4_finalize(&set_coord);
    quat_alg_elem_finalize(&x);
    quat_alg_elem_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return (res);
}

// void quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
// *lat2)
int
quat_test_randomized_lattice_add(int lattice_bitsize, int iterations)
{
    int res = 0;
    int randret = 0;
    ibz_t det;
    quat_lattice_t lat1, lat2, sum;
    ibz_init(&det);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&sum);

    for (int iter = 0; iter < iterations; iter++) {
        randret = quat_test_input_random_lattice_generation(&lat1, lattice_bitsize, 1, 1);
        randret = randret | quat_test_input_random_lattice_generation(&lat2, lattice_bitsize, 1, 1);
        if (!randret)
            goto fin;
        quat_lattice_add(&sum, &lat1, &lat2);
        res = res | !quat_lattice_inclusion(&lat1, &sum);
        res = res | !quat_lattice_inclusion(&lat2, &sum);
    }

fin:;
    if (randret != 0) {
        printf("Randomness failed in quaternion unit test with randomization for lattice_add\n");
    }
    if (res != 0) {
        printf("Quaternion unit test with randomization for lattice_add failed\n");
    }
    ibz_finalize(&det);
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&sum);
    return (res);
}

// int ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p);
int
quat_test_randomized_ibz_cornacchia_prime(int bitsize, int n_bound, int iterations)
{
    int res = 0;
    ibz_t x, y, n, prod, c_res, p;
    int32_t rand_fact;
    int randret = 0;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&n);
    ibz_init(&p);
    ibz_init(&prod);
    ibz_init(&c_res);
    for (int iter = 0; iter < iterations; iter++) {
        // Sample small n for cornacchia
        rand_fact = 0;
        while (rand_fact < 1) {
            randret = randret | randombytes((unsigned char *)&rand_fact, sizeof(int32_t));
            if (randret != 0)
                goto fin;
            if (rand_fact < 0)
                rand_fact = -rand_fact;
            rand_fact = rand_fact % n_bound;
            ibz_set(&n, rand_fact);
        }
        randret = randret | !ibz_generate_random_prime(&p, 0, bitsize, 32);
        if (randret != 0)
            goto fin;
        // If the legendre symbol is ok, Cornacchia should sometimes be able to solve
        ibz_neg(&prod, &n);
        ibz_mod(&prod, &prod, &p);
        if (ibz_legendre(&prod, &p) > -1) {
            //  If there is output, check the output is correct
            if (ibz_cornacchia_prime(&x, &y, &n, &p)) {
                ibz_mul(&c_res, &x, &x);
                ibz_mul(&prod, &y, &y);
                ibz_mul(&prod, &prod, &n);
                ibz_add(&c_res, &c_res, &prod);
                res = res || (0 != ibz_cmp(&p, &c_res));
            }
        } else {
            // Otherwise Cornacchia should fail
            res = res || (ibz_cornacchia_prime(&x, &y, &n, &p));
        }
    }
fin:;

    if (randret != 0) {
        printf("Randomness failed in quaternion unit test with randomization for "
               "ibz_cornacchia_prime\n");
    }
    if (res != 0) {
        printf("Quaternion unit test with randomization for ibz_cornacchia_prime failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&n);
    ibz_finalize(&p);
    ibz_finalize(&prod);
    ibz_finalize(&c_res);
    return res;
}

// run all previous tests
int
quat_test_with_randomization(void)
{
    int res = 0;
    printf("\nRunning randomized tests from quaternion module\n");
    res = res | quat_test_randomized_ibz_mat_2x2_inv_mod(370, 270, 100);
    res = res | quat_test_randomized_ibz_mat_4x4_inv_with_det_as_denom(1500, 10);
    res = res | quat_test_randomized_lattice_contains(250, 250, 10);
    res = res | quat_test_randomized_lattice_add(700, 100);
    res = res | quat_test_randomized_ibz_cornacchia_prime(128, 6, 10);
    return (res);
}
