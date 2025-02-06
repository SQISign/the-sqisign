#include "quaternion_tests.h"

// helper functions

// int quat_lattice_equal(const quat_lattice_t *lat1, const quat_lattice_t *lat2);
int
quat_test_lattice_equal(void)
{
    int res = 0;
    quat_lattice_t lat, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&cmp);

    ibz_mat_4x4_identity(&(lat.basis));
    ibz_mat_4x4_identity(&(cmp.basis));
    res = res || !quat_lattice_equal(&lat, &cmp);
    ibz_set(&(lat.denom), 5);
    ibz_set(&(cmp.denom), 4);
    res = res || quat_lattice_equal(&lat, &cmp);
    ibz_set(&(lat.denom), 1);
    ibz_set(&(cmp.denom), -1);
    res = res || !quat_lattice_equal(&lat, &cmp);
    ibz_set(&(lat.denom), 3);
    ibz_set(&(cmp.denom), 3);
    res = res || !quat_lattice_equal(&lat, &cmp);
    ibz_set(&(lat.basis[0][0]), 1);
    ibz_set(&(lat.basis[0][3]), -1);
    ibz_set(&(lat.basis[1][1]), -2);
    ibz_set(&(lat.basis[2][2]), 1);
    ibz_set(&(lat.basis[2][1]), 1);
    ibz_set(&(lat.basis[3][3]), -3);
    ibz_set(&(lat.denom), 6);
    quat_lattice_hnf(&lat);
    ibz_mat_4x4_copy(&(cmp.basis), &(lat.basis));
    ibz_set(&(cmp.denom), 6);
    res = res || !quat_lattice_equal(&lat, &cmp);
    ibz_set(&(cmp.denom), -7);
    res = res || quat_lattice_equal(&lat, &cmp);
    ibz_set(&(cmp.denom), 6);
    ibz_set(&(cmp.basis[3][3]), 165);
    res = res || quat_lattice_equal(&lat, &cmp);

    if (res != 0) {
        printf("Quaternion unit test lattice_equal failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&cmp);
    return (res);
}

// int quat_lattice_inclusion(const quat_lattice_t *sublat, const quat_lattice_t *overlat)
int
quat_test_lattice_inclusion(void)
{
    int res = 0;
    quat_lattice_t lat, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&cmp);

    ibz_mat_4x4_identity(&(lat.basis));
    ibz_mat_4x4_identity(&(cmp.basis));
    res = res || !quat_lattice_inclusion(&lat, &cmp);
    ibz_set(&(lat.denom), 5);
    ibz_set(&(cmp.denom), 4);
    res = res || quat_lattice_inclusion(&lat, &cmp);
    ibz_set(&(lat.denom), 1);
    ibz_set(&(cmp.denom), 3);
    res = res || !quat_lattice_inclusion(&lat, &cmp);
    ibz_set(&(lat.denom), 3);
    ibz_set(&(cmp.denom), 3);
    res = res || !quat_lattice_inclusion(&lat, &cmp);
    ibz_set(&(lat.basis[0][0]), 1);
    ibz_set(&(lat.basis[0][3]), -1);
    ibz_set(&(lat.basis[1][1]), -2);
    ibz_set(&(lat.basis[2][2]), 1);
    ibz_set(&(lat.basis[2][1]), 1);
    ibz_set(&(lat.basis[3][3]), -3);
    ibz_set(&(lat.denom), 6);
    quat_lattice_hnf(&lat);
    ibz_mat_4x4_copy(&(cmp.basis), &(lat.basis));
    ibz_set(&(cmp.denom), 6);
    res = res || !quat_lattice_inclusion(&lat, &cmp);
    ibz_set(&(cmp.denom), 12);
    res = res || !quat_lattice_inclusion(&lat, &cmp);
    ibz_set(&(cmp.denom), 6);
    ibz_set(&(cmp.basis[3][3]), 165);
    res = res || quat_lattice_inclusion(&lat, &cmp);

    if (res != 0) {
        printf("Quaternion unit test lattice_inclusion failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&cmp);
    return (res);
}

// void quat_lattice_reduce_denom(quat_lattice_t *reduced, const quat_lattice_t *lat);
int
quat_test_lattice_reduce_denom(void)
{
    int res = 0;
    int s;
    quat_lattice_t red, lat, cmp;
    quat_lattice_init(&red);
    quat_lattice_init(&cmp);
    quat_lattice_init(&lat);

    s = 15;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat.basis[i][j]), (i + j) * s);
            ibz_set(&(cmp.basis[i][j]), (i + j));
        }
    }
    ibz_set(&(lat.denom), 4 * s);
    ibz_set(&(cmp.denom), 4);

    quat_lattice_reduce_denom(&red, &lat);
    res = res || (!ibz_mat_4x4_equal(&(red.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(red.denom), &(cmp.denom));

    quat_lattice_reduce_denom(&lat, &lat);
    res = res || (!ibz_mat_4x4_equal(&(lat.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(lat.denom), &(cmp.denom));

    if (res != 0) {
        printf("Quaternion unit test lattice_reduce_denom failed\n");
    }
    quat_lattice_finalize(&red);
    quat_lattice_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return (res);
}

// void quat_lattice_conjugate_without_hnf(quat_lattice_t *conj, const quat_lattice_t *lat);
int
quat_test_lattice_conjugate_without_hnf(void)
{
    int res = 0;
    quat_lattice_t lat, conj, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&conj);
    quat_lattice_init(&cmp);
    // set lattice
    ibz_mat_4x4_zero(&(lat.basis));
    ibz_set(&(lat.basis[0][0]), 4);
    ibz_set(&(lat.basis[0][3]), 1);
    ibz_set(&(lat.basis[1][1]), -2);
    ibz_set(&(lat.basis[2][2]), -1);
    ibz_set(&(lat.basis[2][1]), -1);
    ibz_set(&(lat.basis[3][3]), -3);
    ibz_set(&(lat.denom), 6);
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_set(&(cmp.basis[0][0]), 4);
    ibz_set(&(cmp.basis[0][3]), 1);
    ibz_set(&(cmp.basis[1][1]), 2);
    ibz_set(&(cmp.basis[2][2]), 1);
    ibz_set(&(cmp.basis[2][1]), 1);
    ibz_set(&(cmp.basis[3][3]), 3);
    ibz_set(&(cmp.denom), 6);
    quat_lattice_hnf(&lat);
    quat_lattice_conjugate_without_hnf(&conj, &lat);
    quat_lattice_hnf(&conj);
    quat_lattice_hnf(&cmp);
    res = res || !quat_lattice_equal(&conj, &cmp);
    // test whether coj of conj is original lattice
    quat_lattice_conjugate_without_hnf(&conj, &conj);
    quat_lattice_hnf(&conj);
    res = res || !quat_lattice_equal(&conj, &lat);

    if (res != 0) {
        printf("Quaternion unit test lattice_conjugate_without_hnf failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&conj);
    quat_lattice_finalize(&cmp);
    return (res);
}

// void quat_lattice_dual_without_hnf(quat_lattice_t *dual, const quat_lattice_t *lat);
int
quat_test_lattice_dual_without_hnf(void)
{
    int res = 0;
    quat_lattice_t lat, dual, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&dual);
    quat_lattice_init(&cmp);
    // set lattice
    ibz_mat_4x4_zero(&(lat.basis));
    ibz_set(&(lat.basis[0][0]), 1);
    ibz_set(&(lat.basis[0][3]), -1);
    ibz_set(&(lat.basis[1][1]), -2);
    ibz_set(&(lat.basis[2][2]), 1);
    ibz_set(&(lat.basis[2][1]), 1);
    ibz_set(&(lat.basis[3][3]), -3);
    ibz_set(&(lat.denom), 6);
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_set(&(cmp.basis[0][0]), 6);
    ibz_set(&(cmp.basis[1][1]), 3);
    ibz_set(&(cmp.basis[2][2]), 6);
    ibz_set(&(cmp.basis[3][3]), 2);
    ibz_set(&(cmp.denom), 1);
    quat_lattice_hnf(&lat);
    // test whether dual of dual is original lattice, but dual is not.
    quat_lattice_dual_without_hnf(&dual, &lat);
    quat_lattice_hnf(&dual);
    quat_lattice_hnf(&cmp);
    res = res || !quat_lattice_equal(&dual, &cmp);
    res = res || quat_lattice_equal(&dual, &lat);
    quat_lattice_dual_without_hnf(&dual, &dual);
    quat_lattice_hnf(&dual);
    res = res || !quat_lattice_equal(&dual, &lat);

    if (res != 0) {
        printf("Quaternion unit test lattice_dual_without_hnf failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&dual);
    quat_lattice_finalize(&cmp);
    return (res);
}

// void quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
// *lat2);
int
quat_test_lattice_add(void)
{
    int res = 0;
    quat_lattice_t lat1, lat2, cmp, sum;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&sum);
    quat_lattice_init(&cmp);
    ibz_mat_4x4_zero(&(lat1.basis));
    ibz_mat_4x4_zero(&(lat2.basis));
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_set(&(lat1.basis[0][0]), 44);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat1.basis[0][3]), 32);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(cmp.basis[0][0]), 2);
    ibz_set(&(cmp.basis[0][2]), 1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), 2);
    ibz_set(&(cmp.basis[1][1]), 1);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(cmp.basis[2][2]), 1);
    ibz_set(&(lat1.basis[3][3]), 1);
    ibz_set(&(lat2.basis[3][3]), 3);
    ibz_set(&(cmp.basis[3][3]), 3);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);
    ibz_set(&(cmp.denom), 12);

    quat_lattice_add(&sum, &lat1, &lat2);
    res = res || (!ibz_mat_4x4_equal(&(sum.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(sum.denom), &(cmp.denom));

    // same lattices but not under hnf
    ibz_mat_4x4_zero(&(lat1.basis));
    ibz_mat_4x4_zero(&(lat2.basis));
    ibz_set(&(lat1.basis[0][0]), 4);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(lat2.basis[0][3]), -1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), -2);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(lat2.basis[2][1]), 1);
    ibz_set(&(lat1.basis[3][3]), 7);
    ibz_set(&(lat2.basis[3][3]), -3);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);

    quat_lattice_add(&sum, &lat1, &lat2);
    res = res || (!ibz_mat_4x4_equal(&(sum.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(sum.denom), &(cmp.denom));

    // double in place gives hnf
    ibz_mat_4x4_copy(&(cmp.basis), &lat2.basis);
    ibz_copy(&(cmp.denom), &(lat2.denom));
    quat_lattice_hnf(&cmp);
    quat_lattice_add(&lat2, &lat2, &lat2);
    res = res || (!ibz_mat_4x4_equal(&(lat2.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(lat2.denom), &(cmp.denom));

    if (res != 0) {
        printf("Quaternion unit test lattice_add failed\n");
    }
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&sum);
    quat_lattice_finalize(&cmp);
    return (res);
}

// void quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
// *lat2);
int
quat_test_lattice_intersect(void)
{
    int res = 0;
    quat_lattice_t lat1, lat2, inter, cmp;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&inter);
    quat_lattice_init(&cmp);
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_mat_4x4_zero(&(lat1.basis));
    ibz_mat_4x4_zero(&(lat2.basis));
    ibz_set(&(lat1.basis[0][0]), 4);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(lat2.basis[0][3]), -1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), -2);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(lat2.basis[2][1]), 1);
    ibz_set(&(lat1.basis[3][3]), 7);
    ibz_set(&(lat2.basis[3][3]), -3);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);
    quat_lattice_hnf(&lat1);
    quat_lattice_hnf(&lat2);

    ibz_set(&(cmp.basis[0][0]), 2);
    ibz_set(&(cmp.basis[0][2]), 1);
    ibz_set(&(cmp.basis[1][1]), 10);
    ibz_set(&(cmp.basis[2][2]), 3);
    ibz_set(&(cmp.basis[3][3]), 7);
    ibz_set(&(cmp.denom), 2);
    quat_lattice_intersect(&inter, &lat1, &lat2);

    res = res || !quat_lattice_equal(&inter, &cmp);
    quat_lattice_intersect(&lat2, &lat1, &lat2);
    res = res || !quat_lattice_equal(&lat2, &cmp);
    ibz_mat_4x4_copy(&(cmp.basis), &(lat1.basis));
    ibz_copy(&(cmp.denom), &(lat1.denom));
    quat_lattice_intersect(&lat1, &lat1, &lat1);
    res = res || !quat_lattice_equal(&lat1, &cmp);

    if (res != 0) {
        printf("Quaternion unit test lattice_intersect failed\n");
    }
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&inter);
    quat_lattice_finalize(&cmp);
    return (res);
}

// void quat_lattice_mat_alg_coord_mul_without_hnf(ibz_mat_4x4_t *prod, const ibz_mat_4x4_t *lat,
// const ibz_vec_4_t *coord, const quat_alg_t *alg);
int
quat_test_lattice_mat_alg_coord_mul_without_hnf(void)
{
    int res = 0;
    ibz_mat_4x4_t prod, cmp, lat;
    ibz_vec_4_t elem;
    quat_alg_t alg;
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_init(&lat);
    ibz_vec_4_init(&elem);
    quat_alg_init_set_ui(&alg, 23);
    ibz_vec_4_set(&elem, 3, 4, -1, 0);
    ibz_set(&(lat[0][0]), 11);
    ibz_set(&(lat[1][1]), 13);
    ibz_set(&(lat[2][2]), 15);
    ibz_set(&(lat[3][3]), 4);
    ibz_set(&(lat[0][1]), 9);
    ibz_set(&(lat[1][3]), 1);
    quat_lattice_mat_alg_coord_mul_without_hnf(&prod, &lat, &elem, &alg);
    ibz_set(&(cmp[0][0]), 33);
    ibz_set(&(cmp[1][0]), 44);
    ibz_set(&(cmp[2][0]), -11);
    ibz_set(&(cmp[3][0]), 0);
    ibz_set(&(cmp[0][1]), 27 - 4 * 13);
    ibz_set(&(cmp[1][1]), 36 + 3 * 13);
    ibz_set(&(cmp[2][1]), -9 + 0);
    ibz_set(&(cmp[3][1]), 0 - 13);
    ibz_set(&(cmp[0][2]), 15 * 23);
    ibz_set(&(cmp[1][2]), 0);
    ibz_set(&(cmp[2][2]), 45);
    ibz_set(&(cmp[3][2]), -60);
    ibz_set(&(cmp[0][3]), -4);
    ibz_set(&(cmp[1][3]), 3 + 23 * 4);
    ibz_set(&(cmp[2][3]), 0 + 4 * 4);
    ibz_set(&(cmp[3][3]), -1 + 4 * 3);
    res = res || !ibz_mat_4x4_equal(&cmp, &prod);
    quat_lattice_mat_alg_coord_mul_without_hnf(&lat, &lat, &elem, &alg);
    res = res || !ibz_mat_4x4_equal(&cmp, &lat);

    if (res != 0) {
        printf("Quaternion unit test lattice_mat_alg_coord_mul_without_hnf failed\n");
    }
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&cmp);
    ibz_mat_4x4_finalize(&lat);
    ibz_vec_4_finalize(&elem);
    quat_alg_finalize(&alg);
    return (res);
}

// void quat_lattice_alg_elem_mul(quat_lattice_t *prod, const quat_lattice_t *lat, const
// quat_alg_elem_t *elem, const quat_alg_t *alg);
int
quat_test_lattice_alg_elem_mul(void)
{
    int res = 0;
    quat_lattice_t prod, cmp, lat;
    quat_alg_elem_t elem;
    quat_alg_t alg;
    quat_lattice_init(&prod);
    quat_lattice_init(&cmp);
    quat_lattice_init(&lat);
    quat_alg_elem_init(&elem);
    quat_alg_init_set_ui(&alg, 23);
    quat_alg_elem_set(&elem, 2, 3, 4, -1, 0);
    ibz_set(&(lat.basis[0][0]), 11);
    ibz_set(&(lat.basis[1][1]), -13);
    ibz_set(&(lat.basis[2][2]), 15);
    ibz_set(&(lat.basis[3][3]), -4);
    ibz_set(&(lat.basis[0][1]), 2);
    ibz_set(&(lat.basis[1][3]), -1);
    ibz_set(&(lat.denom), 5);
    quat_lattice_hnf(&lat);
    quat_lattice_alg_elem_mul(&prod, &lat, &elem, &alg);
    ibz_set(&(cmp.basis[0][0]), 33);
    ibz_set(&(cmp.basis[1][0]), 44);
    ibz_set(&(cmp.basis[2][0]), -11);
    ibz_set(&(cmp.basis[3][0]), 0);
    ibz_set(&(cmp.basis[0][1]), 27 - 4 * 13);
    ibz_set(&(cmp.basis[1][1]), 36 + 3 * 13);
    ibz_set(&(cmp.basis[2][1]), -9 + 0);
    ibz_set(&(cmp.basis[3][1]), 0 - 13);
    ibz_set(&(cmp.basis[0][2]), 15 * 23);
    ibz_set(&(cmp.basis[1][2]), 0);
    ibz_set(&(cmp.basis[2][2]), 45);
    ibz_set(&(cmp.basis[3][2]), -60);
    ibz_set(&(cmp.basis[0][3]), -4);
    ibz_set(&(cmp.basis[1][3]), 3 + 23 * 4);
    ibz_set(&(cmp.basis[2][3]), 0 + 4 * 4);
    ibz_set(&(cmp.basis[3][3]), -1 + 4 * 3);
    ibz_set(&(cmp.denom), 10);
    quat_lattice_hnf(&cmp);
    res = res || !quat_lattice_equal(&cmp, &prod);
    quat_lattice_alg_elem_mul(&lat, &lat, &elem, &alg);
    res = res || !quat_lattice_equal(&cmp, &lat);

    if (res != 0) {
        printf("Quaternion unit test lattice_alg_elem_mul failed\n");
    }
    quat_lattice_finalize(&prod);
    quat_lattice_finalize(&cmp);
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&elem);
    quat_alg_finalize(&alg);
    return (res);
}

// void quat_lattice_mul(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
// *lat2, const quat_alg_t *alg);
int
quat_test_lattice_mul(void)
{
    int res = 0;
    quat_lattice_t lat1, lat2, cmp, prod;
    quat_alg_t alg;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&prod);
    quat_lattice_init(&cmp);
    quat_alg_init_set_ui(&alg, 19);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat1.basis[i][j]), 0);
            ibz_set(&(lat2.basis[i][j]), 0);
            ibz_set(&(cmp.basis[i][j]), 0);
        }
    }

    ibz_set(&(lat1.basis[0][0]), 44);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat1.basis[0][3]), 32);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(cmp.basis[0][0]), 1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), 2);
    ibz_set(&(cmp.basis[1][1]), 1);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(cmp.basis[2][2]), 1);
    ibz_set(&(lat1.basis[3][3]), 1);
    ibz_set(&(lat2.basis[3][3]), 3);
    ibz_set(&(cmp.basis[3][3]), 1);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);
    ibz_set(&(cmp.denom), 24);

    quat_lattice_mul(&prod, &lat1, &lat2, &alg);
    res = res || (!ibz_mat_4x4_equal(&(prod.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(prod.denom), &(cmp.denom));

    // same lattices but not under hnf
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat1.basis[i][j]), 0);
            ibz_set(&(lat2.basis[i][j]), 0);
        }
    }
    ibz_set(&(lat1.basis[0][0]), 4);
    ibz_set(&(lat1.basis[0][2]), 3);
    ibz_set(&(lat2.basis[0][0]), 1);
    ibz_set(&(lat2.basis[0][3]), -1);
    ibz_set(&(lat1.basis[1][1]), 5);
    ibz_set(&(lat2.basis[1][1]), -2);
    ibz_set(&(lat1.basis[2][2]), 3);
    ibz_set(&(lat2.basis[2][2]), 1);
    ibz_set(&(lat2.basis[2][1]), 1);
    ibz_set(&(lat1.basis[3][3]), 7);
    ibz_set(&(lat2.basis[3][3]), -3);
    ibz_set(&(lat1.denom), 4);
    ibz_set(&(lat2.denom), 6);

    quat_lattice_mul(&prod, &lat1, &lat2, &alg);
    res = res || (!ibz_mat_4x4_equal(&(prod.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(prod.denom), &(cmp.denom));

    // double in place gives hnf
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(cmp.basis[i][j]), 0);
        }
    }
    ibz_set(&(cmp.basis[0][0]), 1);
    ibz_set(&(cmp.basis[1][1]), 1);
    ibz_set(&(cmp.basis[2][2]), 1);
    ibz_set(&(cmp.basis[3][3]), 1);
    ibz_set(&(cmp.denom), 36);
    quat_lattice_mul(&lat2, &lat2, &lat2, &alg);
    res = res || (!ibz_mat_4x4_equal(&(lat2.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(lat2.denom), &(cmp.denom));

    if (res != 0) {
        printf("Quaternion unit test lattice_mul failed\n");
    }
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&prod);
    quat_lattice_finalize(&cmp);
    quat_alg_finalize(&alg);
    return (res);
}

// int quat_lattice_contains(quat_alg_coord_t *coord, const quat_lattice_t *lat, const
// quat_alg_elem_t *x, const quat_alg_t *alg);
int
quat_test_lattice_contains(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_alg_elem_t x;
    ibz_vec_4_t coord, cmp;
    quat_lattice_t lat;
    quat_alg_init_set_ui(&alg, 103);
    quat_alg_elem_init(&x);
    ibz_vec_4_init(&coord);
    ibz_vec_4_init(&cmp);
    quat_lattice_init(&lat);

    // lattice 1
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat.basis[i][j]), 0);
        }
    }
    ibz_set(&(lat.basis[0][0]), 4);
    ibz_set(&(lat.basis[0][2]), 3);
    ibz_set(&(lat.basis[1][1]), 5);
    ibz_set(&(lat.basis[2][2]), 3);
    ibz_set(&(lat.basis[3][3]), 7);
    ibz_set(&(lat.denom), 4);

    // x 1, should fail
    ibz_set(&(x.denom), 3);
    ibz_set(&(x.coord[0]), 1);
    ibz_set(&(x.coord[1]), -2);
    ibz_set(&(x.coord[2]), 26);
    ibz_set(&(x.coord[3]), 9);

    res = res || quat_lattice_contains(&coord, &lat, &x);
    // again, but with NULL
    res = res || quat_lattice_contains(NULL, &lat, &x);

    // lattice 2
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat.basis[i][j]), 0);
        }
    }
    ibz_set(&(lat.basis[0][0]), 1);
    ibz_set(&(lat.basis[0][3]), -1);
    ibz_set(&(lat.basis[1][1]), -2);
    ibz_set(&(lat.basis[2][2]), 1);
    ibz_set(&(lat.basis[2][1]), 1);
    ibz_set(&(lat.basis[3][3]), -3);
    ibz_set(&(lat.denom), 6);
    quat_lattice_hnf(&lat);
    // x 1, should succeed
    ibz_set(&(x.denom), 3);
    ibz_set(&(x.coord[0]), 1);
    ibz_set(&(x.coord[1]), -2);
    ibz_set(&(x.coord[2]), 26);
    ibz_set(&(x.coord[3]), 9);
    ibz_set(&(cmp[0]), 2);
    ibz_set(&(cmp[1]), -2);
    ibz_set(&(cmp[2]), 52);
    ibz_set(&(cmp[3]), 6);

    res = res || (0 == quat_lattice_contains(&coord, &lat, &x));

    res = res || ibz_cmp(&(coord[0]), &(cmp[0]));
    res = res || ibz_cmp(&(coord[1]), &(cmp[1]));
    res = res || ibz_cmp(&(coord[2]), &(cmp[2]));
    res = res || ibz_cmp(&(coord[3]), &(cmp[3]));
    // again, but with NULL
    res = res || (0 == quat_lattice_contains(NULL, &lat, &x));

    if (res != 0) {
        printf("Quaternion unit test lattice_contains failed\n");
    }
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&x);
    ibz_vec_4_finalize(&coord);
    ibz_vec_4_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return (res);
}

// void quat_lattice_index(ibz_t *index, const quat_lattice_t *sublat, const quat_lattice_t
// *overlat);
int
quat_test_lattice_index(void)
{
    int res = 0;
    quat_lattice_t sublat, overlat;
    ibz_t index;
    ibz_init(&index);
    quat_lattice_init(&sublat);
    quat_lattice_init(&overlat);

    ibz_mat_4x4_zero(&(sublat.basis));
    ibz_mat_4x4_identity(&(overlat.basis));
    ibz_set(&(overlat.denom), 2);
    ibz_set(&(sublat.basis[0][0]), 2);
    ibz_set(&(sublat.basis[0][1]), 0);
    ibz_set(&(sublat.basis[0][2]), 1);
    ibz_set(&(sublat.basis[0][3]), 0);
    ibz_set(&(sublat.basis[1][0]), 0);
    ibz_set(&(sublat.basis[1][1]), 4);
    ibz_set(&(sublat.basis[1][2]), 2);
    ibz_set(&(sublat.basis[1][3]), 3);
    ibz_set(&(sublat.basis[2][0]), 0);
    ibz_set(&(sublat.basis[2][1]), 0);
    ibz_set(&(sublat.basis[2][2]), 1);
    ibz_set(&(sublat.basis[2][3]), 0);
    ibz_set(&(sublat.basis[3][0]), 0);
    ibz_set(&(sublat.basis[3][1]), 0);
    ibz_set(&(sublat.basis[3][2]), 0);
    ibz_set(&(sublat.basis[3][3]), 1);
    ibz_set(&(sublat.denom), 2);
    quat_lattice_index(&index, &sublat, &overlat);

    res = res || !(ibz_cmp_int32(&index, 8) == 0);

    if (res != 0) {
        printf("Quaternion unit test lattice_index failed\n");
    }
    quat_lattice_finalize(&sublat);
    quat_lattice_finalize(&overlat);
    ibz_finalize(&index);
    return (res);
}

// void quat_lattice_hnf(quat_lattice_t *lat);
int
quat_test_lattice_hnf(void)
{
    int res = 0;
    quat_lattice_t lat, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&cmp);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat.basis[i][j]), 0);
            ibz_set(&(cmp.basis[i][j]), 0);
        }
    }
    ibz_set(&(lat.basis[0][0]), 1);
    ibz_set(&(lat.basis[0][3]), -1);
    ibz_set(&(lat.basis[1][1]), -2);
    ibz_set(&(lat.basis[2][2]), 1);
    ibz_set(&(lat.basis[2][1]), 1);
    ibz_set(&(lat.basis[3][3]), -3);
    ibz_set(&(cmp.basis[0][0]), 1);
    ibz_set(&(cmp.basis[1][1]), 2);
    ibz_set(&(cmp.basis[2][2]), 1);
    ibz_set(&(cmp.basis[3][3]), 3);
    ibz_set(&(cmp.denom), 6);
    ibz_set(&(lat.denom), 6);

    quat_lattice_hnf(&lat);
    res = res || (!ibz_mat_4x4_equal(&(lat.basis), &(cmp.basis)));
    res = res || ibz_cmp(&(lat.denom), &(cmp.denom));

    if (res != 0) {
        printf("Quaternion unit test lattice_hnf failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&cmp);
    return (res);
}

int
quat_test_lattice_gram()
{
    int res = 0;

    quat_lattice_t lattice;
    ibz_mat_4x4_t gram;
    ibz_t test, norm1, norm2;
    ibz_vec_4_t vec1, vec2;
    quat_alg_elem_t elem1, elem2;
    quat_alg_t alg;
    quat_lattice_init(&lattice);
    ibz_mat_4x4_init(&gram);
    ibz_init(&test);
    ibz_init(&norm1);
    ibz_init(&norm2);
    ibz_vec_4_init(&vec1);
    ibz_vec_4_init(&vec2);
    quat_alg_elem_init(&elem1);
    quat_alg_elem_init(&elem2);
    quat_alg_init_set_ui(&alg, 103);

    quat_lattice_O0_set(&lattice);
    quat_lattice_gram(&gram, &lattice, &alg);
    quat_alg_elem_set(&elem1, 1, 2, 3, 4, 1);
    quat_alg_elem_set(&elem2, 1, 2, 4, 4, 1);
    quat_lattice_contains(&vec1, &lattice, &elem1);
    quat_lattice_contains(&vec2, &lattice, &elem2);
    quat_alg_conj(&elem2, &elem2);
    quat_alg_mul(&elem1, &elem1, &elem2, &alg);
    ibz_mul(&norm1, &(elem1.coord[0]), &ibz_const_two);
    ibz_div(&norm1, &test, &norm1, &(elem1.denom));

    ibz_mat_4x4_eval(&vec1, &gram, &vec1);
    ibz_set(&norm2, 0);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&test, &(vec1[i]), &(vec2[i]));
        ibz_add(&norm2, &norm2, &test);
    }
    ibz_div(&norm2, &test, &norm2, &(lattice.denom));
    ibz_div(&norm2, &test, &norm2, &(lattice.denom));
    res = res | !(ibz_cmp(&norm1, &norm2) == 0);

    ibz_mat_4x4_zero(&(lattice.basis));
    ibz_set(&(lattice.basis[0][0]), 202);
    ibz_set(&(lattice.basis[1][1]), 202);
    ibz_set(&(lattice.basis[2][2]), 1);
    ibz_set(&(lattice.basis[3][3]), 1);
    ibz_set(&(lattice.basis[0][2]), 158);
    ibz_set(&(lattice.basis[0][3]), 53);
    ibz_set(&(lattice.basis[1][2]), 149);
    ibz_set(&(lattice.basis[1][3]), 158);
    ibz_set(&(lattice.denom), 2);
    quat_lattice_gram(&gram, &lattice, &alg);

    quat_alg_elem_set(&elem1, 2, 360, 149, 1, 0);
    quat_alg_elem_set(&elem2, 2, 53, 360, 0, 1);
    int ok = quat_lattice_contains(&vec1, &lattice, &elem1);
    ok = ok && quat_lattice_contains(&vec2, &lattice, &elem2);
    assert(ok);
    quat_alg_conj(&elem2, &elem2);
    quat_alg_mul(&elem1, &elem1, &elem2, &alg);
    ibz_mul(&norm1, &(elem1.coord[0]), &ibz_const_two);
    ibz_div(&norm1, &test, &norm1, &(elem1.denom));

    ibz_mat_4x4_eval(&vec1, &gram, &vec1);
    ibz_set(&norm2, 0);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&test, &(vec1[i]), &(vec2[i]));
        ibz_add(&norm2, &norm2, &test);
    }
    ibz_div(&norm2, &test, &norm2, &(lattice.denom));
    ibz_div(&norm2, &test, &norm2, &(lattice.denom));
    res = res | !(ibz_cmp(&norm1, &norm2) == 0);

    if (res != 0) {
        printf("Quaternion unit test lattice_gram failed\n");
    }
    quat_lattice_finalize(&lattice);
    ibz_mat_4x4_finalize(&gram);
    ibz_finalize(&test);
    ibz_finalize(&norm1);
    ibz_finalize(&norm2);
    ibz_vec_4_finalize(&vec1);
    ibz_vec_4_finalize(&vec2);
    quat_alg_elem_finalize(&elem1);
    quat_alg_elem_finalize(&elem2);
    quat_alg_finalize(&alg);
    return (res);
}

// run all previous tests
int
quat_test_lattice(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of lattice functions\n");
    res = res | quat_test_lattice_equal();
    res = res | quat_test_lattice_inclusion();
    res = res | quat_test_lattice_reduce_denom();
    res = res | quat_test_lattice_conjugate_without_hnf();
    res = res | quat_test_lattice_dual_without_hnf();
    res = res | quat_test_lattice_add();
    res = res | quat_test_lattice_intersect();
    res = res | quat_test_lattice_mat_alg_coord_mul_without_hnf();
    res = res | quat_test_lattice_alg_elem_mul();
    res = res | quat_test_lattice_mul();
    res = res | quat_test_lattice_contains();
    res = res | quat_test_lattice_index();
    res = res | quat_test_lattice_hnf();
    res = res | quat_test_lattice_gram();
    return (res);
}
