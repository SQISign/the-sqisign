#include "quaternion_tests.h"

// Schema of tests: initialize structure, assign values, finalize

// void quat_alg_init(quat_alg_t *alg);
// void quat_alg_finalize(quat_alg_t *alg);
int
quat_test_finit_alg(void)
{
    int res = 0;
    quat_alg_t alg;
    ibz_t p;
    ibz_init(&p);
    ibz_set(&p, 7);
    quat_alg_init_set(&alg, &p);
    res = res || ibz_cmp(&(alg.p), &p);
    if (res != 0) {
        printf("Quaternion unit test finit_alg failed\n");
    }
    ibz_finalize(&p);
    quat_alg_finalize(&alg);
    return res;
}

// void quat_alg_elem_init(quat_alg_elem_t *elem);
// void quat_alg_elem_finalize(quat_alg_elem_t *elem);
int
quat_test_finit_alg_elem(void)
{
    quat_alg_elem_t elem;
    int res;
    quat_alg_elem_init(&elem);
    ibz_set(&(elem.coord[0]), 0);
    ibz_set(&(elem.coord[1]), 1);
    ibz_set(&(elem.coord[2]), 2);
    ibz_set(&(elem.coord[3]), 3);
    ibz_set(&(elem.denom), 1);
    res = 1 - (1 == ibz_is_one(&(elem.denom)));
    for (int i = 0; i < 4; i++) {
        res = res || (ibz_cmp_int32(&(elem.coord[i]), i) != 0);
    }
    if (res != 0) {
        printf("Quaternion unit test finit_alg_elem failed\n");
    }
    quat_alg_elem_finalize(&elem);
    return res;
}

// void ibz_vec_2_init(ibz_vec_2_t *vec);
// void ibz_vec_2_finalize(ibz_vec_2_t *vec);
int
quat_test_finit_ibz_vec_2(void)
{
    ibz_vec_2_t vec;
    int res = 0;
    ibz_vec_2_init(&vec);
    for (int i = 0; i < 2; i++) {
        ibz_set(&(vec[i]), i);
    }
    for (int i = 0; i < 2; i++) {
        res = res || (ibz_cmp_int32(&(vec[i]), i) != 0);
    }
    if (res != 0) {
        printf("Quaternion unit test finit_ibz_vec_2 failed\n");
    }
    ibz_vec_2_finalize(&vec);
    return res;
}

// void ibz_vec_4_init(ibz_vec_4_t *vec);
// void ibz_vec_4_finalize(ibz_vec_4_t *vec);
int
quat_test_finit_ibz_vec_4(void)
{
    ibz_vec_4_t vec;
    int res = 0;
    ibz_vec_4_init(&vec);
    for (int i = 0; i < 4; i++) {
        ibz_set(&(vec[i]), i);
    }
    for (int i = 0; i < 4; i++) {
        res = res || (ibz_cmp_int32(&(vec[i]), i) != 0);
    }
    if (res != 0) {
        printf("Quaternion unit test finit_ibz_vec_4 failed\n");
    }
    ibz_vec_4_finalize(&vec);
    return res;
}

// void ibz_mat_2x2_init(ibz_mat_2x2_t *mat);
// void ibz_mat_2x2_finalize(ibz_mat_2x2_t *mat);
int
quat_test_finit_ibz_mat_2x2(void)
{
    ibz_mat_2x2_t mat;
    int res = 0;
    ibz_mat_2x2_init(&mat);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_set(&(mat[i][j]), i + j);
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            res = res || (ibz_cmp_int32(&(mat[i][j]), i + j) != 0);
        }
    }
    if (res != 0) {
        printf("Quaternion unit test finit_ibz_mat_2x2 failed\n");
    }
    ibz_mat_2x2_finalize(&mat);
    return res;
}

// void ibz_mat_4x4_init(ibz_mat_4x4_t *mat);
// void ibz_mat_4x4_finalize(ibz_mat_4x4_t *mat);
int
quat_test_finit_ibz_mat_4x4(void)
{
    ibz_mat_4x4_t mat;
    int res = 0;
    ibz_mat_4x4_init(&mat);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(mat[i][j]), i + j);
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res || (ibz_cmp_int32(&(mat[i][j]), i + j) != 0);
        }
    }
    if (res != 0) {
        printf("Quaternion unit test finit_ibz_mat_4x4 failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    return res;
}

// void quat_lattice_init(quat_lattice_t *lat);
// void quat_lattice_finalize(quat_lattice_t *lat);
int
quat_test_finit_lattice(void)
{
    quat_lattice_t lat;
    int res = 0;
    quat_lattice_init(&lat);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lat.basis[i][j]), i + j);
        }
    }
    ibz_set(&(lat.denom), 1);
    res = 1 - (1 == ibz_is_one(&(lat.denom)));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res || (ibz_cmp_int32(&(lat.basis[i][j]), i + j) != 0);
        }
    }
    if (res != 0) {
        printf("Quaternion unit test finit_alg_lattice failed\n");
    }
    quat_lattice_finalize(&lat);
    return res;
}

// void quat_left_ideal_init(quat_left_ideal_t *lideal);
// void quat_left_ideal_finalize(quat_left_ideal_t *lideal);
int
quat_test_finit_lideal(void)
{
    quat_left_ideal_t lideal;
    int res = 0;
    quat_left_ideal_init(&lideal);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(lideal.lattice.basis[i][j]), i + j);
        }
    }
    ibz_set(&(lideal.lattice.denom), 1);
    ibz_set(&(lideal.norm), 5);
    lideal.parent_order = NULL;
    res = 1 - (1 == ibz_is_one(&(lideal.lattice.denom)));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res || (ibz_cmp_int32(&(lideal.lattice.basis[i][j]), i + j) != 0);
        }
    }
    res = res || (ibz_cmp_int32(&(lideal.norm), 5) != 0);
    if (res != 0) {
        printf("Quaternion unit test finit_alg_lideal failed\n");
    }
    quat_left_ideal_finalize(&lideal);
    return res;
}

// run all previous tests
int
quat_test_finit(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of initializers and finalizers\n");
    res = res | quat_test_finit_alg();
    res = res | quat_test_finit_alg_elem();
    res = res | quat_test_finit_ibz_vec_2();
    res = res | quat_test_finit_ibz_vec_4();
    res = res | quat_test_finit_ibz_mat_2x2();
    res = res | quat_test_finit_ibz_mat_4x4();
    res = res | quat_test_finit_lattice();
    res = res | quat_test_finit_lideal();
    return (res);
}
