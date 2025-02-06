#include "quaternion_tests.h"

// void ibz_vec_2_set(ibz_vec_2_t *vec, int a0, int a1);
int
quat_test_dim2_ibz_vec_2_set(void)
{
    int res = 0;
    ibz_vec_2_t vec;
    ibz_vec_2_init(&vec);
    ibz_vec_2_set(&vec, 2, 5);
    res = res || (ibz_cmp_int32(&(vec[0]), 2) != 0);
    res = res || (ibz_cmp_int32(&(vec[1]), 5) != 0);
    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_vec_2_set failed\n");
    }
    ibz_vec_2_finalize(&vec);
    return (res);
}

// void ibz_mat_2x2_set(ibz_mat_2x2_t *mat, int a00, int a01, int a10, int a11);
int
quat_test_dim2_ibz_mat_2x2_set(void)
{
    int res = 0;
    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);
    ibz_mat_2x2_set(&mat, 2, 7, -1, 5);
    res = res || (ibz_cmp_int32(&(mat[0][0]), 2) != 0);
    res = res || (ibz_cmp_int32(&(mat[0][1]), 7) != 0);
    res = res || (ibz_cmp_int32(&(mat[1][0]), -1) != 0);
    res = res || (ibz_cmp_int32(&(mat[1][1]), 5) != 0);
    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_mat_2x2_set failed\n");
    }
    ibz_mat_2x2_finalize(&mat);
    return (res);
}

// void ibz_mat_2x2_copy(ibz_vec_2_t *copy, const ibz_mat_2x2_t *copied);
int
quat_test_dim2_ibz_mat_2x2_copy(void)
{
    int res = 0;
    ibz_mat_2x2_t mat, copy;
    ibz_mat_2x2_init(&mat);
    ibz_mat_2x2_init(&copy);

    ibz_mat_2x2_set(&mat, 1, -1, 2, 4);
    res = res || (0 == ibz_cmp_int32(&(copy[0][0]), 1));
    res = res || (0 == ibz_cmp_int32(&(copy[0][1]), -1));
    res = res || (0 == ibz_cmp_int32(&(copy[1][0]), 2));
    res = res || (0 == ibz_cmp_int32(&(copy[1][1]), 4));
    res = res || (0 != ibz_cmp_int32(&(mat[0][0]), 1));
    res = res || (0 != ibz_cmp_int32(&(mat[0][1]), -1));
    res = res || (0 != ibz_cmp_int32(&(mat[1][0]), 2));
    res = res || (0 != ibz_cmp_int32(&(mat[1][1]), 4));
    ibz_mat_2x2_copy(&copy, &mat);
    res = res || (0 != ibz_cmp_int32(&(copy[0][0]), 1));
    res = res || (0 != ibz_cmp_int32(&(copy[0][1]), -1));
    res = res || (0 != ibz_cmp_int32(&(copy[1][0]), 2));
    res = res || (0 != ibz_cmp_int32(&(copy[1][1]), 4));
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res = res || (0 != ibz_cmp(&(mat[i][j]), &(copy[i][j])));
        }
    }

    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_mat_2x2_copy failed\n");
    }
    ibz_mat_2x2_finalize(&mat);
    ibz_mat_2x2_finalize(&copy);
    return (res);
}

// void ibz_mat_2x2_det_from_ibz(ibz_t *det, const ibz_t *a11, const ibz_t *a12, const ibz_t *a21,
// const ibz_t *a22);
int
quat_test_dim2_ibz_mat_2x2_det_from_ibz(void)
{
    int res = 0;
    ibz_t det, cmp, a11, a12, a21, a22;
    ibz_init(&a11);
    ibz_init(&a12);
    ibz_init(&a21);
    ibz_init(&a22);
    ibz_init(&det);
    ibz_init(&cmp);

    ibz_set(&a11, 1);
    ibz_set(&a12, 0);
    ibz_set(&a21, 0);
    ibz_set(&a22, 1);
    ibz_set(&cmp, 1);
    ibz_mat_2x2_det_from_ibz(&det, &a11, &a12, &a21, &a22);
    res = res || ibz_cmp(&cmp, &det);

    ibz_set(&a11, 2);
    ibz_set(&a12, 3);
    ibz_set(&a21, 1);
    ibz_set(&a22, -2);
    ibz_set(&cmp, -7);
    ibz_mat_2x2_det_from_ibz(&det, &a11, &a12, &a21, &a22);
    res = res || ibz_cmp(&cmp, &det);

    ibz_set(&a11, 0);
    ibz_set(&a12, 3);
    ibz_set(&a21, -1);
    ibz_set(&a22, 0);
    ibz_set(&cmp, 3);
    ibz_mat_2x2_det_from_ibz(&det, &a11, &a12, &a21, &a22);
    res = res || ibz_cmp(&cmp, &det);

    ibz_set(&a11, 2);
    ibz_set(&cmp, 0);
    ibz_mat_2x2_det_from_ibz(&a11, &a11, &a11, &a11, &a11);
    res = res || ibz_cmp(&cmp, &a11);

    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_mat_2x2_det_from_ibz failed\n");
    }
    ibz_finalize(&a11);
    ibz_finalize(&a12);
    ibz_finalize(&a21);
    ibz_finalize(&a22);
    ibz_finalize(&cmp);
    ibz_finalize(&det);
    return res;
}

// void ibz_mat_2x2_eval(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec);
int
quat_test_dim2_ibz_mat_2x2_eval(void)
{
    int res = 0;
    ibz_vec_2_t vec, ret, cmp;
    ibz_mat_2x2_t mat;
    ibz_vec_2_init(&vec);
    ibz_mat_2x2_init(&mat);
    ibz_vec_2_init(&ret);
    ibz_vec_2_init(&cmp);

    ibz_mat_2x2_set(&mat, 1, -1, 2, 4);
    ibz_vec_2_set(&vec, 1, -1);
    ibz_vec_2_set(&cmp, 2, -2);
    ibz_mat_2x2_eval(&ret, &mat, &vec);
    res = res || ibz_cmp(&(ret[0]), &(cmp[0]));
    res = res || ibz_cmp(&(ret[1]), &(cmp[1]));

    ibz_mat_2x2_set(&mat, 2, -2, 1, 3);
    ibz_vec_2_set(&vec, 2, 4);
    ibz_vec_2_set(&cmp, -4, 14);
    ibz_mat_2x2_eval(&vec, &mat, &vec);
    res = res || ibz_cmp(&(vec[0]), &(cmp[0]));
    res = res || ibz_cmp(&(vec[1]), &(cmp[1]));

    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_mat_2x2_eval failed\n");
    }
    ibz_vec_2_finalize(&vec);
    ibz_mat_2x2_finalize(&mat);
    ibz_vec_2_finalize(&ret);
    ibz_vec_2_finalize(&cmp);
    return (res);
}

// modular 2x2 operations

// void ibz_2x2_mul_mod(ibz_mat_2x2_t *prod, const ibz_mat_2x2_t *mat_a, const ibz_mat_2x2_t *mat_b,
// const ibz_t *m);
int
quat_test_dim2_ibz_mat_2x2_mul_mod(void)
{
    int res = 0;
    ibz_t m;
    ibz_mat_2x2_t a, b, cmp, prod;
    ibz_init(&m);
    ibz_mat_2x2_init(&a);
    ibz_mat_2x2_init(&b);
    ibz_mat_2x2_init(&prod);
    ibz_mat_2x2_init(&cmp);

    ibz_set(&m, 7);
    ibz_mat_2x2_set(&a, 2, -2, 1, 3);
    ibz_mat_2x2_set(&b, 5, 3, 4, 1);
    ibz_mat_2x2_set(&cmp, 2, 4, 3, 6);
    ibz_2x2_mul_mod(&prod, &a, &b, &m);
    res = res || ibz_cmp(&(prod[0][0]), &(cmp[0][0]));
    res = res || ibz_cmp(&(prod[0][1]), &(cmp[0][1]));
    res = res || ibz_cmp(&(prod[1][0]), &(cmp[1][0]));
    res = res || ibz_cmp(&(prod[1][1]), &(cmp[1][1]));
    ibz_mat_2x2_set(&cmp, 6, 6, 2, 2);
    ibz_2x2_mul_mod(&prod, &b, &a, &m);
    res = res || ibz_cmp(&(prod[0][0]), &(cmp[0][0]));
    res = res || ibz_cmp(&(prod[0][1]), &(cmp[0][1]));
    res = res || ibz_cmp(&(prod[1][0]), &(cmp[1][0]));
    res = res || ibz_cmp(&(prod[1][1]), &(cmp[1][1]));

    ibz_set(&m, 12);
    ibz_mat_2x2_set(&a, 2, 7, 1, -2);
    ibz_mat_2x2_set(&cmp, 11, 0, 0, 11);
    ibz_2x2_mul_mod(&a, &a, &a, &m);
    res = res || ibz_cmp(&(a[0][0]), &(cmp[0][0]));
    res = res || ibz_cmp(&(a[0][1]), &(cmp[0][1]));
    res = res || ibz_cmp(&(a[1][0]), &(cmp[1][0]));
    res = res || ibz_cmp(&(a[1][1]), &(cmp[1][1]));

    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_mat_2x2_mul_mod failed\n");
    }
    ibz_mat_2x2_finalize(&a);
    ibz_mat_2x2_finalize(&b);
    ibz_mat_2x2_finalize(&prod);
    ibz_mat_2x2_finalize(&cmp);
    ibz_finalize(&m);
    return (res);
}

// int ibz_mat_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);
int
quat_test_dim2_ibz_mat_2x2_inv_mod(void)
{
    int res = 0;
    ibz_t m;
    ibz_mat_2x2_t a, inv, id, prod;
    ibz_init(&m);
    ibz_mat_2x2_init(&a);
    ibz_mat_2x2_init(&inv);
    ibz_mat_2x2_init(&prod);
    ibz_mat_2x2_init(&id);
    ibz_mat_2x2_set(&id, 1, 0, 0, 1);

    // inverse exists
    ibz_set(&m, 7);
    ibz_mat_2x2_set(&a, 2, -3, 1, 3);
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
    ibz_set(&m, 12);
    ibz_mat_2x2_set(&a, 2, 7, 1, -2);
    ibz_mat_2x2_set(&inv, 2, 7, 1, -2);
    if (ibz_mat_2x2_inv_mod(&inv, &inv, &m)) {
        ibz_2x2_mul_mod(&prod, &a, &inv, &m);
        res = res || ibz_cmp(&(prod[0][0]), &(id[0][0]));
        res = res || ibz_cmp(&(prod[0][1]), &(id[0][1]));
        res = res || ibz_cmp(&(prod[1][0]), &(id[1][0]));
        res = res || ibz_cmp(&(prod[1][1]), &(id[1][1]));
    } else {
        res = 1;
    }

    // no inverse
    ibz_set(&m, 25);
    ibz_mat_2x2_set(&a, 2, -2, -1, 1);
    res = res || ibz_mat_2x2_inv_mod(&inv, &a, &m);
    ibz_set(&m, 7);
    ibz_mat_2x2_set(&a, 2, 3, 1, -2);
    res = res || ibz_mat_2x2_inv_mod(&inv, &a, &m);
    ibz_set(&m, 25);
    ibz_mat_2x2_set(&a, 2, 1, 1, -2);
    res = res || ibz_mat_2x2_inv_mod(&inv, &a, &m);

    if (res != 0) {
        printf("Quaternion unit test dim2_ibz_mat_2x2_inv_mod failed\n");
    }
    ibz_mat_2x2_finalize(&a);
    ibz_mat_2x2_finalize(&inv);
    ibz_mat_2x2_finalize(&prod);
    ibz_mat_2x2_finalize(&id);
    ibz_finalize(&m);
    return (res);
}

// run all previous tests
int
quat_test_dim2(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of functions for matrices, vectors and lattices in "
           "dimension 2\n");
    res = res | quat_test_dim2_ibz_vec_2_set();
    res = res | quat_test_dim2_ibz_mat_2x2_set();
    res = res | quat_test_dim2_ibz_mat_2x2_copy();
    res = res | quat_test_dim2_ibz_mat_2x2_det_from_ibz();
    res = res | quat_test_dim2_ibz_mat_2x2_eval();
    res = res | quat_test_dim2_ibz_mat_2x2_mul_mod();
    res = res | quat_test_dim2_ibz_mat_2x2_inv_mod();
    return (res);
}
