#include "lll_internals.h"
#include "quaternion_tests.h"
#include <rng.h>

int
quat_test_lll_ibq_consts(void)
{
    int ret = 0;
    ibq_t t;
    ibz_t tmp1, tmp2, tmp3;
    ibz_init(&tmp1);
    ibz_init(&tmp2);
    ibz_init(&tmp3);
    ibq_init(&t);

    ibz_set(&tmp1, 123);
    ibz_set(&tmp2, -123);
    if (!ibq_set(&t, &tmp1, &tmp2)) {
        ret = -1;
        goto err;
    }

    if (ibq_is_one(&t)) {
        ret = -1;
        goto err;
    }

    if (!ibq_is_ibz(&t)) {
        ret = -1;
        goto err;
    }

    if (!ibq_to_ibz(&tmp3, &t)) {
        ret = -1;
        goto err;
    }

    if (ibz_is_one(&tmp3)) {
        ret = -1;
        goto err;
    }

    ibz_set(&tmp2, 123);
    if (!ibq_set(&t, &tmp1, &tmp2)) {
        ret = -1;
        goto err;
    }

    if (!ibq_is_one(&t)) {
        ret = -1;
        goto err;
    }

    if (!ibq_is_ibz(&t)) {
        ret = -1;
        goto err;
    }

    if (!ibq_to_ibz(&tmp3, &t)) {
        ret = -1;
        goto err;
    }

    if (!ibz_is_one(&tmp3)) {
        ret = -1;
        goto err;
    }

    ibz_set(&tmp1, 0);
    ibq_set(&t, &tmp1, &tmp2);

    if (!ibq_is_zero(&t)) {
        ret = -1;
        goto err;
    }

    if (!ibq_is_ibz(&t)) {
        ret = -1;
        goto err;
    }

    if (!ibq_to_ibz(&tmp3, &t)) {
        ret = -1;
        goto err;
    }

    if (!ibz_is_zero(&tmp3)) {
        ret = -1;
        goto err;
    }

err:
    ibq_finalize(&t);
    ibz_finalize(&tmp1);
    ibz_finalize(&tmp2);
    ibz_finalize(&tmp3);
    return ret;
}

// test for lll verification
// void ibq_vec_4_copy_ibz(ibq_vec_4_t *vec, const ibz_t *coeff0, const ibz_t *coeff1,const ibz_t
// *coeff2,const ibz_t *coeff3);
int
quat_test_lll_ibq_vec_4_copy_ibz(void)
{
    int res = 0;
    ibq_vec_4_t vec;
    ibz_vec_4_t vec_z;
    ibz_vec_4_init(&vec_z);
    ibq_vec_4_init(&vec);
    ibz_vec_4_set(&vec_z, 2, 3, 4, 5);
    ibq_vec_4_copy_ibz(&vec, &(vec_z[0]), &(vec_z[1]), &(vec_z[2]), &(vec_z[3]));
    for (int i = 0; i < 4; i++) {
        ibq_to_ibz(&(vec_z[i]), &(vec[i]));
        res = res || (ibz_cmp_int32(&(vec_z[i]), i + 2) != 0);
    }

    if (res != 0) {
        printf("Quaternion unit test lll_ibq_vec_4_copy_ibz failed\n");
    }
    ibz_vec_4_finalize(&vec_z);
    ibq_vec_4_finalize(&vec);
    return (res);
}

// void quat_lll_bilinear(ibq_t *b, const ibq_vec_4_t *vec0, const ibq_vec_4_t *vec1, const
// ibz_t *q);
int
quat_test_lll_bilinear(void)
{
    int res = 0;
    ibz_vec_4_t init_helper;
    ibq_vec_4_t vec0, vec1;
    ibz_t q;
    ibq_t cmp, b;
    ibz_vec_4_init(&init_helper);
    ibq_init(&cmp);
    ibq_init(&b);
    ibz_init(&q);
    ibq_vec_4_init(&vec0);
    ibq_vec_4_init(&vec1);
    ibz_vec_4_set(&init_helper, 1, 2, 3, 4);
    ibq_vec_4_copy_ibz(&vec0, &(init_helper[0]), &(init_helper[1]), &(init_helper[2]), &(init_helper[3]));
    ibz_vec_4_set(&init_helper, 9, -8, 7, -6);
    ibq_vec_4_copy_ibz(&vec1, &(init_helper[0]), &(init_helper[1]), &(init_helper[2]), &(init_helper[3]));
    for (int i = 0; i < 4; i++) {
        ibq_inv(&(vec0[i]), &(vec0[i]));
    }
    ibz_set(&q, 3);
    ibz_vec_4_set(&init_helper, 15, 2, 0, 0);
    ibq_set(&cmp, &(init_helper[0]), &(init_helper[1]));
    quat_lll_bilinear(&b, &vec0, &vec1, &q);
    res = res || (ibq_cmp(&b, &cmp));

    if (res != 0) {
        printf("Quaternion unit test quat_lll_bilinear failed\n");
    }
    ibq_finalize(&cmp);
    ibq_finalize(&b);
    ibz_finalize(&q);
    ibz_vec_4_finalize(&init_helper);
    ibq_vec_4_finalize(&vec0);
    ibq_vec_4_finalize(&vec1);
    return (res);
}

// void quat_lll_gram_schmidt_transposed_with_ibq(ibq_mat_4x4_t *orthogonalised_transposed, const
// ibz_mat_4x4_t *mat, const ibz_t *q);
int
quat_test_lll_gram_schmidt_transposed_with_ibq(void)
{
    int res = 0;
    int zero;
    ibq_mat_4x4_t ot, cmp;
    ibz_mat_4x4_t mat;
    ibz_t q, num, denom;
    ibq_t b;
    ibz_init(&q);
    ibz_init(&num);
    ibz_init(&denom);
    ibq_init(&b);
    ibz_mat_4x4_init(&mat);
    ibq_mat_4x4_init(&ot);
    ibq_mat_4x4_init(&cmp);

    ibz_mat_4x4_zero(&mat);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(mat[i][j]), i * i + (j + 5) * j - 2 + (i == j));
        }
    }
    ibz_set(&q, 3);
    quat_lll_gram_schmidt_transposed_with_ibq(&ot, &mat, &q);
    // test orthogonality
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            quat_lll_bilinear(&b, &(ot[i]), &(ot[j]), &q);
            res = res || !ibq_is_zero(&b);
        }
    }
    // test first vector is identical to mat
    for (int i = 0; i < 4; i++) {
        ibq_to_ibz(&q, &(ot[0][i]));
        res = res || ibz_cmp(&q, &(mat[i][0]));
    }
    // test no zero vector
    for (int i = 0; i < 4; i++) {
        zero = 1;
        for (int j = 0; j < 4; j++) {
            zero = zero && ibq_is_zero(&(ot[i][j]));
        }
        res = res || zero;
    }

    ibz_set(&(mat[0][0]), 1);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 1);
    ibz_set(&(mat[0][3]), 0);
    ibz_set(&(mat[1][0]), 0);
    ibz_set(&(mat[1][1]), 1);
    ibz_set(&(mat[1][2]), 0);
    ibz_set(&(mat[1][3]), 1);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), 0);
    ibz_set(&(mat[2][2]), 2);
    ibz_set(&(mat[2][3]), 0);
    ibz_set(&(mat[3][0]), 0);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 2);
    ibz_set(&denom, 1);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibq_set(&(cmp[i][j]), &(mat[j][i]), &denom);
        }
    }
    ibz_set(&denom, 3);
    ibz_set(&num, -2);
    ibq_set(&(cmp[2][0]), &num, &denom);
    ibq_set(&(cmp[3][1]), &num, &denom);
    ibz_set(&num, 1);
    ibq_set(&(cmp[2][2]), &num, &denom);
    ibq_set(&(cmp[3][3]), &num, &denom);
    ibz_set(&q, 2);
    quat_lll_gram_schmidt_transposed_with_ibq(&ot, &mat, &q);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res || ibq_cmp(&(cmp[i][j]), &(ot[i][j]));
        }
    }

    ibz_set(&(mat[0][0]), 1);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 1);
    ibz_set(&(mat[0][3]), 0);
    ibz_set(&(mat[1][0]), 0);
    ibz_set(&(mat[1][1]), 1);
    ibz_set(&(mat[1][2]), 0);
    ibz_set(&(mat[1][3]), 1);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), 0);
    ibz_set(&(mat[2][2]), 2);
    ibz_set(&(mat[2][3]), 1);
    ibz_set(&(mat[3][0]), 0);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 2);
    ibz_set(&denom, 1);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibq_set(&(cmp[i][j]), &(mat[j][i]), &denom);
        }
    }
    ibz_set(&denom, 3);
    ibz_set(&num, -2);
    ibq_set(&(cmp[2][0]), &num, &denom);
    ibq_set(&(cmp[3][1]), &num, &denom);
    ibz_set(&num, 1);
    ibq_set(&(cmp[2][2]), &num, &denom);
    ibq_set(&(cmp[3][3]), &num, &denom);
    ibz_set(&num, 0);
    ibq_set(&(cmp[3][0]), &num, &denom);
    ibq_set(&(cmp[3][2]), &num, &denom);
    ibz_set(&q, 2);
    quat_lll_gram_schmidt_transposed_with_ibq(&ot, &mat, &q);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res || ibq_cmp(&(cmp[i][j]), &(ot[i][j]));
        }
    }

    if (res != 0) {
        printf("Quaternion unit test dim4_gram_schmidt_transposed_with_ibq failed\n");
    }
    ibz_finalize(&q);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibq_finalize(&b);
    ibz_mat_4x4_finalize(&mat);
    ibq_mat_4x4_finalize(&ot);
    ibq_mat_4x4_finalize(&cmp);
    return (res);
}

// int quat_lll_verify(const ibz_mat_4x4_t *mat, const ibq_t *delta, const ibq_t *eta, const
// quat_alg_t *alg);
int
quat_test_lll_verify(void)
{
    int res = 0;
    ibz_mat_4x4_t mat;
    ibz_t q, coeff_num, coeff_denom;
    ibq_t eta, delta;
    quat_alg_t alg;
    ibz_mat_4x4_init(&mat);
    ibz_init(&q);
    ibq_init(&delta);
    ibq_init(&eta);
    ibz_init(&coeff_num);
    ibz_init(&coeff_denom);

    // reduced: non-1 norm
    ibz_set(&q, 3);
    quat_alg_init_set(&alg, &q);
    ibq_set(&eta, &ibz_const_one, &ibz_const_two);
    ibq_set(&delta, &ibz_const_three, &ibz_const_two);
    ibq_mul(&delta, &delta, &eta);
    ibz_set(&(mat[0][0]), 0);
    ibz_set(&(mat[0][1]), 2);
    ibz_set(&(mat[0][2]), 3);
    ibz_set(&(mat[0][3]), -14);
    ibz_set(&(mat[1][0]), 2);
    ibz_set(&(mat[1][1]), -1);
    ibz_set(&(mat[1][2]), -4);
    ibz_set(&(mat[1][3]), -8);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), -2);
    ibz_set(&(mat[2][2]), 1);
    ibz_set(&(mat[2][3]), 0);
    ibz_set(&(mat[3][0]), 1);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 7);
    res = res || !quat_lll_verify(&mat, &delta, &eta, &alg);
    quat_alg_finalize(&alg);

    // reduced: non-1 norm
    ibz_set(&q, 103);
    quat_alg_init_set(&alg, &q);
    ibz_set(&coeff_num, 99);
    ibz_set(&coeff_denom, 100);
    ibq_set(&delta, &coeff_num, &coeff_denom);
    ibz_set(&(mat[0][0]), 3);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 90);
    ibz_set(&(mat[0][3]), -86);
    ibz_set(&(mat[1][0]), 11);
    ibz_set(&(mat[1][1]), 15);
    ibz_set(&(mat[1][2]), 12);
    ibz_set(&(mat[1][3]), 50);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), -2);
    ibz_set(&(mat[2][2]), 0);
    ibz_set(&(mat[2][3]), 3);
    ibz_set(&(mat[3][0]), -1);
    ibz_set(&(mat[3][1]), 0);
    ibz_set(&(mat[3][2]), 5);
    ibz_set(&(mat[3][3]), 5);
    res = res || !quat_lll_verify(&mat, &delta, &eta, &alg);
    quat_alg_finalize(&alg);

    if (res != 0) {
        printf("Quaternion unit test quat_lll_verify failed\n");
    }
    ibz_finalize(&q);
    ibq_finalize(&delta);
    ibq_finalize(&eta);
    ibz_finalize(&coeff_num);
    ibz_finalize(&coeff_denom);
    ibz_mat_4x4_finalize(&mat);
    return (res);
}

// int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int
// precision);
int
quat_test_lll_lattice_lll(void)
{
    int res = 0;
    quat_lattice_t lat, test;
    ibz_mat_4x4_t red;
    ibz_t num, denom, q;
    ibq_t eta, delta;
    quat_alg_t alg;
    ibz_init(&num);
    ibz_init(&denom);
    ibz_init(&q);
    ibq_init(&delta);
    ibq_init(&eta);
    ibz_mat_4x4_init(&red);
    quat_lattice_init(&lat);
    quat_lattice_init(&test);
    ibz_set(&q, 103);
    quat_alg_init_set(&alg, &q);

    // set lattice
    ibz_set(&lat.denom, 60);
    ibz_mat_4x4_zero(&(lat.basis));
    ibz_set(&lat.basis[0][0], 3);
    ibz_set(&lat.basis[1][0], 7);
    ibz_set(&lat.basis[0][1], 1);
    ibz_set(&lat.basis[3][1], -6);
    ibz_set(&lat.basis[1][2], 12);
    ibz_set(&lat.basis[2][2], 5);
    ibz_set(&lat.basis[0][3], -19);
    ibz_set(&lat.basis[3][3], 3);

    quat_lattice_hnf(&lat);

    res = res || quat_lattice_lll(&red, &lat, &alg);
    // test lll reduced
    quat_lll_set_ibq_parameters(&delta, &eta);
    res = res || !quat_lll_verify(&red, &delta, &eta, &alg);
    // test lattice equality
    ibz_copy(&(test.denom), &(lat.denom));
    ibz_mat_4x4_copy(&(test.basis), &(red));
    quat_lattice_hnf(&test);
    res = res || !quat_lattice_equal(&test, &lat);

    if (res != 0) {
        printf("Quaternion unit test lll_lattice_lll failed\n");
    }
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibz_finalize(&q);
    ibq_finalize(&eta);
    ibq_finalize(&delta);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&test);
    quat_alg_finalize(&alg);
    return (res);
}

// int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const quat_alg_t *alg);
int
quat_test_lll_randomized_lattice_lll(void)
{
    int res = 0;
    quat_lattice_t lat, test;
    ibz_mat_4x4_t red;
    ibz_t q, det;
    ibq_t delta, eta;
    int32_t rand[4][4];
    int32_t rand_denom;
    uint32_t rand_q;
    ibz_init(&q);
    ibz_init(&det);
    ibq_init(&eta);
    ibq_init(&delta);
    ibz_mat_4x4_init(&red);
    quat_lattice_init(&lat);
    quat_lattice_init(&test);
    quat_lll_set_ibq_parameters(&delta, &eta);

    for (int iter = 0; iter < 20; iter++) {
        quat_alg_t alg;
        rand_denom = 0;
        while (rand_denom <= 0) {
            int randret = randombytes((unsigned char *)&rand_denom, sizeof(int32_t));
            if (randret != 0)
                return 1;
        }
        int randret = randombytes((unsigned char *)&rand_q, sizeof(uint32_t));
        if (randret != 0)
            return 1;
        // generate random invertible matrix
        ibz_set(&det, 0);
        while (ibz_is_zero(&det)) {
            randret = randombytes((unsigned char *)rand, sizeof(rand));
            if (randret != 0)
                return 1;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    ibz_set(&(lat.basis[i][j]), rand[j][i]);
                }
            }
            ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &(lat.basis));
        }

        // set lattice
        ibz_set(&lat.denom, rand_denom);
        quat_lattice_hnf(&lat);
        // set algebra
        ibz_set(&q, 1 + (rand_q % 1023));
        quat_alg_init_set(&alg, &q);
        // reduce
        res = res || quat_lattice_lll(&red, &lat, &alg);
        // test lll reduced
        res = res || !quat_lll_verify(&red, &delta, &eta, &alg);
        // test lattice equality
        ibz_copy(&(test.denom), &(lat.denom));
        ibz_mat_4x4_copy(&(test.basis), &(red));
        quat_lattice_hnf(&test);
        res = res || !quat_lattice_equal(&test, &lat);
        quat_alg_finalize(&alg);
    }

    if (res != 0) {
        printf("Quaternion unit test of lll with randomization for lattice_lll failed\n");
    }
    ibz_finalize(&q);
    ibz_finalize(&det);
    ibq_finalize(&delta);
    ibq_finalize(&eta);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&test);
    return (res);
}

// int quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced, ibz_mat_4x4_t *gram, const
// quat_left_ideal_t *lideal, const quat_alg_t *alg);
int
quat_test_lideal_reduce_basis()
{
    int res = 0;
    ibz_mat_4x4_t red, gram, prod, gram_norm;
    ibz_vec_4_t vec;
    quat_left_ideal_t lideal;
    quat_lattice_t test;
    quat_alg_t alg;
    quat_alg_elem_t init_helper;
    quat_lattice_t order;
    ibq_t delta, eta;
    ibz_t num, denom, norm, test_norm;
    ibz_init(&num);
    ibz_init(&denom);
    ibz_init(&norm);
    ibz_init(&test_norm);
    ibq_init(&delta);
    ibq_init(&eta);
    ibz_vec_4_init(&vec);
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&gram_norm);
    quat_lattice_init(&test);
    quat_alg_elem_init(&init_helper);
    quat_lattice_init(&order);
    quat_alg_init_set_ui(&alg, 19);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);
    quat_left_ideal_init(&lideal);
    quat_lattice_O0_set(&order);
    quat_alg_elem_set(&init_helper, 1, 1, 2, 8, 8);
    quat_lideal_create_principal(&lideal, &init_helper, &order, &alg);
    quat_lattice_reduce_denom(&(lideal.lattice), &(lideal.lattice));
    quat_lideal_reduce_basis(&red, &gram, &lideal, &alg);
    quat_lll_set_ibq_parameters(&delta, &eta);
    res = res || !quat_lll_verify(&red, &delta, &eta, &alg);
    // test reduced and lideal generate same lattice
    ibz_mat_4x4_copy(&(test.basis), &red);
    ibz_copy(&(test.denom), &(lideal.lattice.denom));
    quat_lattice_hnf(&test);
    res = res || !quat_lattice_equal(&(lideal.lattice), &test);
    // test gram matrix is gram matrix
    ibz_mat_4x4_identity(&gram_norm);
    ibz_copy(&(gram_norm[2][2]), &(alg.p));
    ibz_copy(&(gram_norm[3][3]), &(alg.p));
    ibz_mat_4x4_transpose(&prod, &red);
    ibz_mat_4x4_mul(&prod, &prod, &gram_norm);
    ibz_mat_4x4_mul(&prod, &prod, &red);
    for (int i = 0; i < 4; i++) {
        ibz_vec_4_set(&vec, (i == 0), (i == 1), (i == 2), (i == 3));
        quat_qf_eval(&norm, &gram, &vec);
        quat_qf_eval(&test_norm, &prod, &vec);
        ibz_mul(&norm, &(lideal.norm), &norm);
        res = res || !(ibz_cmp(&norm, &test_norm) == 0);
    }

    if (res != 0) {
        printf("Quaternion unit test lideal_reduce_basis failed\n");
    }
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibz_finalize(&norm);
    ibz_finalize(&test_norm);
    ibq_finalize(&delta);
    ibq_finalize(&eta);
    quat_lattice_finalize(&test);
    quat_alg_elem_finalize(&init_helper);
    ibz_mat_4x4_finalize(&prod);
    quat_lattice_finalize(&order);
    quat_alg_finalize(&alg);
    ibz_vec_4_finalize(&vec);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);
    ibz_mat_4x4_finalize(&gram_norm);
    quat_left_ideal_finalize(&lideal);
    return (res);
}

// int quat_lideal_lideal_mul_reduced(quat_left_ideal_t *prod, ibz_mat_4x4_t *gram, const
// quat_left_ideal_t *lideal1,const quat_left_ideal_t *lideal2, const quat_alg_t *alg);
int
quat_test_lll_lideal_lideal_mul_reduced()
{
    int res = 0;
    ibz_t n, norm, test_norm;
    ibq_t delta, eta;
    ibz_vec_4_t vec;
    quat_alg_t alg;
    quat_alg_elem_t gen;
    quat_lattice_t order, ro;
    quat_left_ideal_t lideal1, lideal2, prod, i1, i2;
    ibz_mat_4x4_t gram;
    ibz_mat_4x4_t gram_test;

    ibz_init(&n);
    ibz_init(&norm);
    ibz_init(&test_norm);
    ibq_init(&delta);
    ibq_init(&eta);
    ibz_vec_4_init(&vec);
    quat_alg_elem_init(&gen);
    quat_left_ideal_init(&lideal1);
    quat_left_ideal_init(&lideal2);
    quat_left_ideal_init(&i1);
    quat_left_ideal_init(&i2);
    quat_left_ideal_init(&prod);
    quat_lattice_init(&order);
    quat_lattice_init(&ro);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&gram_test);

    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_O0_set(&order);
    quat_lll_set_ibq_parameters(&delta, &eta);

    ibz_set(&n, 113);
    quat_alg_elem_set(&gen, 1, 10, 0, 1, 3);
    quat_lideal_create(&lideal1, &gen, &n, &order, &alg);

    ibz_set(&n, 89);
    quat_alg_elem_set(&gen, 2, 2, 5, 1, 4);
    quat_lideal_create(&lideal2, &gen, &n, &order, &alg);

    quat_lideal_copy(&i1, &lideal1);
    quat_lideal_copy(&i2, &lideal2);

    quat_lideal_lideal_mul_reduced(&prod, &gram, &lideal1, &lideal2, &alg);
    res = res || !quat_lll_verify(&(prod.lattice.basis), &delta, &eta, &alg);
    ibz_mat_4x4_identity(&(gram_test));
    ibz_copy(&(gram_test[2][2]), &(alg.p));
    ibz_copy(&(gram_test[3][3]), &(alg.p));
    ibz_mat_4x4_mul(&(gram_test), &(gram_test), &(prod.lattice.basis));
    ibz_mat_4x4_transpose(&(gram_test), &(gram_test));
    ibz_mat_4x4_mul(&(gram_test), &(gram_test), &(prod.lattice.basis));
    for (int i = 0; i < 4; i++) {
        ibz_vec_4_set(&vec, (i == 0), (i == 1), (i == 2), (i == 3));
        quat_qf_eval(&norm, &gram, &vec);
        quat_qf_eval(&test_norm, &gram_test, &vec);
        ibz_mul(&norm, &(prod.norm), &norm);
        res = res || !(ibz_cmp(&norm, &test_norm) == 0);
    }
    quat_lattice_hnf(&(prod.lattice));

    res = res || !quat_lideal_equals(&i1, &lideal1, &alg);
    res = res || !quat_lideal_equals(&i2, &lideal2, &alg);
    quat_lattice_mul(&i1.lattice, &i1.lattice, &i2.lattice, &alg);
    res = res || !quat_lattice_equal(&i1.lattice, &prod.lattice);
    res = res || !(prod.parent_order == lideal1.parent_order);
    i1.parent_order = lideal1.parent_order;
    quat_lideal_norm(&i1);
    res = res || !quat_lideal_equals(&i1, &prod, &alg);

    if (res != 0) {
        printf("Quaternion unit test lideal_lideal_mul_reduced failed\n");
    }
    ibz_finalize(&n);
    ibz_finalize(&norm);
    ibz_finalize(&test_norm);
    ibq_finalize(&delta);
    ibq_finalize(&eta);
    ibz_vec_4_finalize(&vec);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal1);
    quat_left_ideal_finalize(&lideal2);
    quat_left_ideal_finalize(&i1);
    quat_left_ideal_finalize(&i2);
    quat_left_ideal_finalize(&prod);
    quat_lattice_finalize(&order);
    quat_lattice_finalize(&ro);
    quat_alg_finalize(&alg);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&gram_test);
    return (res);
}

// int quat_lideal_prime_norm_reduced_equivalent(quat_left_ideal_t *lideal, const quat_alg_t *alg,
// const int primality_num_iter, const int equiv_bound_coeff, const int equiv_num_iter);
int
quat_test_lll_lideal_prime_norm_reduced_equivalent()
{
    int res = 0;
    ibz_t n, d;
    quat_alg_t alg;
    quat_alg_elem_t gen;
    ibz_mat_4x4_t red, gram;
    quat_lattice_t order, ro, ro2;
    quat_left_ideal_t lideal1, lideal2, i1;

    ibz_init(&n);
    ibz_init(&d);
    quat_alg_elem_init(&gen);
    quat_left_ideal_init(&lideal1);
    quat_left_ideal_init(&lideal2);
    quat_left_ideal_init(&i1);
    quat_lattice_init(&order);
    quat_lattice_init(&ro);
    quat_lattice_init(&ro2);
    ibz_mat_4x4_init(&red);
    ibz_mat_4x4_init(&gram);

    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_O0_set(&order);

    ibz_set(&n, 113);
    quat_alg_elem_set(&gen, 1, 10, 0, 1, 3);
    quat_lideal_create(&lideal1, &gen, &n, &order, &alg);

    quat_lideal_copy(&i1, &lideal1);
    quat_lideal_right_order(&ro, &lideal1, &alg);

    quat_lideal_prime_norm_reduced_equivalent(&lideal1, &alg, 20, 20);

    // test norm correctness
    quat_lattice_hnf(&(lideal1.lattice));
    ibz_copy(&n, &(lideal1.norm));
    quat_lideal_norm(&lideal1);
    res = res || (0 != ibz_cmp(&n, &(lideal1.norm)));

    // test norm primality
    res = res || !ibz_probab_prime(&n, 20);

    // test equivalence
    quat_lideal_right_order(&ro2, &lideal1, &alg);
    quat_lattice_mul(&(lideal2.lattice), &ro, &ro2, &alg);
    ibz_set(&(lideal2.lattice.denom), 1);
    lideal2.parent_order = &ro;
    quat_lattice_hnf(&(lideal2.lattice));
    quat_lideal_norm(&lideal2);
    // now lideal2 is a connecting idea of ro and ro2
    quat_lideal_reduce_basis(&red, &gram, &lideal2, &alg);
    quat_alg_elem_copy_ibz(&gen, &(lideal2.lattice.denom), &(red[0][0]), &(red[1][0]), &(red[2][0]), &(red[3][0]));
    quat_alg_norm(&n, &d, &gen, &alg);
    assert(ibz_is_one(&d));
    res = res || (0 != ibz_cmp(&n, &(lideal2.norm)));

    if (res != 0) {
        printf("Quaternion unit test lideal_prime_norm_reduced_equivalent failed\n");
    }
    ibz_finalize(&n);
    ibz_finalize(&d);
    ibz_mat_4x4_finalize(&red);
    ibz_mat_4x4_finalize(&gram);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&lideal1);
    quat_left_ideal_finalize(&lideal2);
    quat_left_ideal_finalize(&i1);
    quat_lattice_finalize(&order);
    quat_lattice_finalize(&ro);
    quat_lattice_finalize(&ro2);
    quat_alg_finalize(&alg);
    return (res);
}

// run all previous tests
int
quat_test_lll(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of lll and its subfunctions\n");
    res = res | quat_test_lll_ibq_consts();
    res = res | quat_test_lll_ibq_vec_4_copy_ibz();
    res = res | quat_test_lll_bilinear();
    res = res | quat_test_lll_gram_schmidt_transposed_with_ibq();
    res = res | quat_test_lll_verify();
    res = res | quat_test_lll_lattice_lll();
    res = res | quat_test_lll_randomized_lattice_lll();
    res = res | quat_test_lideal_reduce_basis();
    res = res | quat_test_lll_lideal_lideal_mul_reduced();
    res = res | quat_test_lll_lideal_prime_norm_reduced_equivalent();
    return (res);
}
