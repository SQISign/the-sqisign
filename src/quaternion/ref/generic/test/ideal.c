#include "quaternion_tests.h"

// void quat_lideal_norm(quat_left_ideal_t *lideal);
int
quat_test_lideal_norm(void)
{
    int res = 0;
    quat_left_ideal_t lideal;
    quat_lattice_t order;
    quat_alg_t alg;
    quat_alg_elem_t gen;
    ibz_t norm;
    ibz_init(&norm);
    quat_left_ideal_init(&lideal);
    quat_lattice_init(&order);
    quat_alg_elem_init(&gen);
    quat_alg_init_set_ui(&alg, 23);
    quat_lattice_O0_set(&order);
    quat_alg_elem_set(&gen, 1, 3, 2, 1, 4);
    assert(quat_lattice_contains(NULL, &order, &gen));
    ibz_set(&norm, 17);
    quat_lideal_create(&lideal, &gen, &norm, &order, &alg);
    ibz_copy(&norm, &(lideal.norm));
    ibz_set(&(lideal.norm), 0);
    quat_lideal_norm(&lideal);
    res = res || ibz_cmp(&(lideal.norm), &norm);

    if (res != 0) {
        printf("Quaternion unit test lideal_norm failed\n");
    }
    ibz_finalize(&norm);
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&gen);
    quat_alg_finalize(&alg);
    return (res);
}

// void quat_lideal_copy(quat_left_ideal_t *copy, const quat_left_ideal_t *copied);
int
quat_test_lideal_copy()
{
    int res = 0;
    quat_left_ideal_t copy, copied;
    quat_alg_elem_t elem;
    quat_alg_t alg;
    quat_lattice_t order;
    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_init(&order);
    quat_alg_elem_init(&elem);
    quat_left_ideal_init(&copy);
    quat_left_ideal_init(&copied);
    quat_lattice_O0_set(&order);
    quat_alg_elem_set(&elem, 1, 4, 2, 9, -1);
    quat_lideal_create_principal(&copied, &elem, &order, &alg);
    quat_lideal_copy(&copy, &copied);
    res = res | !(quat_lideal_equals(&copy, &copied, &alg));
    ibz_set(&(elem.coord[0]), 23);
    quat_lideal_create_principal(&copied, &elem, &order, &alg);
    res = res | (quat_lideal_equals(&copy, &copied, &alg));
    quat_lideal_copy(&copy, &copied);
    res = res | !(quat_lideal_equals(&copy, &copied, &alg));
    if (res) {
        printf("Quaternion unit test lideal_copy failed\n");
    }
    quat_alg_elem_finalize(&elem);
    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order);
    quat_left_ideal_finalize(&copy);
    quat_left_ideal_finalize(&copied);
    return (res);
}

// void quat_lideal_create_principal(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const
// quat_lattice_t *order, const quat_alg_t *alg);
int
quat_test_lideal_create_principal(void)
{
    int res = 0;

    quat_alg_t alg;
    quat_lattice_t lat;
    quat_alg_elem_t gamma;
    quat_left_ideal_t I;
    quat_alg_init_set_ui(&alg, 367);
    quat_lattice_init(&lat);
    quat_alg_elem_init(&gamma);
    quat_left_ideal_init(&I);
    quat_lattice_O0_set(&lat);
    ibz_set(&gamma.coord[0], 219);
    ibz_set(&gamma.coord[1], 200);
    ibz_set(&gamma.coord[2], 78);
    ibz_set(&gamma.coord[3], -1);

    quat_lideal_create_principal(&I, &gamma, &lat, &alg);

    res |= I.parent_order != &lat;
    res |= ibz_cmp_int32(&I.norm, 2321156);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_one);

    res |= ibz_cmp_int32(&I.lattice.basis[0][0], 1160578);
    res |= ibz_cmp_int32(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[1][1], 1160578);
    res |= ibz_cmp_int32(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][2], 310126);
    res |= ibz_cmp_int32(&I.lattice.basis[1][2], 182529);
    res |= ibz_cmp_int32(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_int32(&I.lattice.basis[3][2], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][3], 978049);
    res |= ibz_cmp_int32(&I.lattice.basis[1][3], 310126);
    res |= ibz_cmp_int32(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][3], 1);

    // same test, just with gamma not reduced
    ibz_set(&gamma.coord[0], 438);
    ibz_set(&gamma.coord[1], 400);
    ibz_set(&gamma.coord[2], 156);
    ibz_set(&gamma.coord[3], -2);
    ibz_set(&gamma.denom, 2);

    quat_lideal_create_principal(&I, &gamma, &lat, &alg);

    res |= I.parent_order != &lat;
    res |= ibz_cmp_int32(&I.norm, 2321156);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_one);

    res |= ibz_cmp_int32(&I.lattice.basis[0][0], 1160578);
    res |= ibz_cmp_int32(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[1][1], 1160578);
    res |= ibz_cmp_int32(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][2], 310126);
    res |= ibz_cmp_int32(&I.lattice.basis[1][2], 182529);
    res |= ibz_cmp_int32(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_int32(&I.lattice.basis[3][2], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][3], 978049);
    res |= ibz_cmp_int32(&I.lattice.basis[1][3], 310126);
    res |= ibz_cmp_int32(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][3], 1);

    // same test, just with gamma and basis not reduced
    ibz_set(&lat.denom, 6);
    ibz_set(&lat.basis[0][0], 6);
    ibz_set(&lat.basis[1][1], 6);
    ibz_set(&lat.basis[1][2], 3);
    ibz_set(&lat.basis[2][2], 3);
    ibz_set(&lat.basis[3][3], 3);
    ibz_set(&lat.basis[0][3], 3);
    ibz_set(&gamma.coord[0], 438);
    ibz_set(&gamma.coord[1], 400);
    ibz_set(&gamma.coord[2], 156);
    ibz_set(&gamma.coord[3], -2);
    ibz_set(&gamma.denom, 2);

    quat_lideal_create_principal(&I, &gamma, &lat, &alg);

    res |= I.parent_order != &lat;
    res |= ibz_cmp_int32(&I.norm, 2321156);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_one);

    res |= ibz_cmp_int32(&I.lattice.basis[0][0], 1160578);
    res |= ibz_cmp_int32(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[1][1], 1160578);
    res |= ibz_cmp_int32(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][2], 310126);
    res |= ibz_cmp_int32(&I.lattice.basis[1][2], 182529);
    res |= ibz_cmp_int32(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_int32(&I.lattice.basis[3][2], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][3], 978049);
    res |= ibz_cmp_int32(&I.lattice.basis[1][3], 310126);
    res |= ibz_cmp_int32(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][3], 1);

    quat_alg_finalize(&alg);
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&gamma);
    quat_left_ideal_finalize(&I);

    if (res != 0) {
        printf("Quaternion unit test lideal_create_principal failed\n");
    }
    return (res);
}

// void quat_lideal_create(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const
// ibz_t *N, const quat_lattice_t *order, const quat_alg_t *alg);
int
quat_test_lideal_create_from_primitive(void)
{
    int res = 0;

    quat_alg_t alg;
    quat_lattice_t lat;
    quat_alg_elem_t gamma;
    ibz_t N;
    quat_left_ideal_t I;
    quat_alg_init_set_ui(&alg, 367);
    quat_lattice_init(&lat);
    quat_alg_elem_init(&gamma);
    ibz_init(&N);
    quat_left_ideal_init(&I);
    quat_lattice_O0_set(&lat);
    ibz_set(&gamma.coord[0], 219);
    ibz_set(&gamma.coord[1], 200);
    ibz_set(&gamma.coord[2], 78);
    ibz_set(&gamma.coord[3], -1);
    ibz_set(&N, 31);

    quat_lideal_create(&I, &gamma, &N, &lat, &alg);

    res |= I.parent_order != &lat;
    res |= ibz_cmp(&I.norm, &N);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_two);

    res |= ibz_cmp_int32(&I.lattice.basis[0][0], 62);
    res |= ibz_cmp_int32(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[1][1], 62);
    res |= ibz_cmp_int32(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][2], 2);
    res |= ibz_cmp_int32(&I.lattice.basis[1][2], 1);
    res |= ibz_cmp_int32(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_int32(&I.lattice.basis[3][2], 0);

    res |= ibz_cmp_int32(&I.lattice.basis[0][3], 61);
    res |= ibz_cmp_int32(&I.lattice.basis[1][3], 2);
    res |= ibz_cmp_int32(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_int32(&I.lattice.basis[3][3], 1);

    quat_alg_finalize(&alg);
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&gamma);
    ibz_finalize(&N);
    quat_left_ideal_finalize(&I);

    if (res != 0) {
        printf("Quaternion unit test lideal_create_from_primitive failed\n");
    }
    return (res);
}

// int quat_lideal_generator(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t
// *alg)
int
quat_test_lideal_generator(void)
{
    int res = 0;

    quat_alg_t alg;
    quat_lattice_t order;
    quat_alg_elem_t gen;
    ibz_t N;
    quat_left_ideal_t lideal, lideal2;
    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_init(&order);
    quat_alg_elem_init(&gen);
    ibz_init(&N);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal2);
    quat_lattice_O0_set(&order);

    ibz_set(&gen.coord[0], 3);
    ibz_set(&gen.coord[1], 5);
    ibz_set(&gen.coord[2], 7);
    ibz_set(&gen.coord[3], 11);
    ibz_set(&N, 17);

    quat_lideal_create(&lideal, &gen, &N, &order, &alg);

    // Try to regenerate the same ideal
    for (int i = 0; i <= 10; i++) {
        res |= !quat_lideal_generator(&gen, &lideal, &alg);
        quat_lideal_create(&lideal2, &gen, &N, &order, &alg);
        res |= !quat_lideal_equals(&lideal, &lideal2, &alg);
    }

    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&N);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal2);

    if (res != 0) {
        printf("Quaternion unit test lideal_generator failed\n");
    }
    return (res);
}

// void quat_lideal_mul(quat_left_ideal_t *product, const quat_left_ideal_t *lideal, const
// quat_alg_elem_t *alpha, const quat_alg_t *alg);
int
quat_test_lideal_mul(void)
{
    int res = 0;

    quat_alg_t alg;
    quat_lattice_t order;
    quat_alg_elem_t gen1, gen2, gen_prod;
    quat_left_ideal_t lideal, lideal2;
    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_init(&order);
    quat_alg_elem_init(&gen1);
    quat_alg_elem_init(&gen2);
    quat_alg_elem_init(&gen_prod);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal2);
    quat_lattice_O0_set(&order);

    ibz_set(&gen1.coord[0], 3);
    ibz_set(&gen1.coord[1], 5);
    ibz_set(&gen1.coord[2], 7);
    ibz_set(&gen1.coord[3], 11);
    ibz_set(&gen2.coord[0], -2);
    ibz_set(&gen2.coord[1], 13);
    ibz_set(&gen2.coord[2], -17);
    ibz_set(&gen2.coord[3], 19);

    // Check that (O·gen1)·gen2 == O·(gen1·gen2)
    quat_lideal_create_principal(&lideal, &gen1, &order, &alg);
    quat_lideal_mul(&lideal, &lideal, &gen2, &alg);
    quat_alg_mul(&gen_prod, &gen1, &gen2, &alg);
    quat_lideal_create_principal(&lideal2, &gen_prod, &order, &alg);
    res |= !quat_lideal_equals(&lideal, &lideal2, &alg);

    // Same test, but with non-integral gen2
    ibz_set(&(gen2.denom), 2);
    quat_lideal_create_principal(&lideal, &gen1, &order, &alg);
    quat_lideal_mul(&lideal, &lideal, &gen2, &alg);
    quat_alg_mul(&gen_prod, &gen1, &gen2, &alg);
    quat_lideal_create_principal(&lideal2, &gen_prod, &order, &alg);
    res |= !quat_lideal_equals(&lideal, &lideal2, &alg);

    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&gen1);
    quat_alg_elem_finalize(&gen2);
    quat_alg_elem_finalize(&gen_prod);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal2);

    if (res != 0) {
        printf("Quaternion unit test lideal_mul failed\n");
    }
    return (res);
}

// void quat_lideal_conjugate_without_hnf(quat_left_ideal_t *conj, quat_lattice_t *new_parent_order,
// const quat_left_ideal_t *lideal, const quat_alg_t *alg);
int
quat_test_lideal_conj_without_hnf(void)
{
    int res = 0;
    ibz_t n;
    quat_lattice_t o, ro, rro, test, cmp;
    quat_left_ideal_t lideal, conj;
    quat_alg_elem_t gen;
    quat_alg_t alg;
    ibz_init(&n);
    quat_alg_elem_init(&gen);
    quat_lattice_init(&o);
    quat_lattice_init(&ro);
    quat_lattice_init(&rro);
    quat_lattice_init(&test);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&conj);
    quat_lattice_init(&cmp);
    quat_alg_init_set_ui(&alg, 103);
    ibz_set(&n, 25);
    quat_alg_elem_set(&gen, 1, 2, 3, 0, 2);
    quat_lattice_O0_set(&o);
    quat_lideal_create(&lideal, &gen, &n, &o, &alg);

    quat_lideal_conjugate_without_hnf(&conj, &ro, &lideal, &alg);
    quat_lattice_hnf(&(conj.lattice));
    quat_lattice_mul(&test, &(lideal.lattice), &(conj.lattice), &alg);
    ibz_copy(&(cmp.denom), &(o.denom));
    ibz_mat_4x4_scalar_mul(&(cmp.basis), &n, &(o.basis));
    res = res || !quat_lattice_equal(&test, &cmp);
    quat_lattice_mul(&test, &(conj.lattice), &(lideal.lattice), &alg);
    ibz_copy(&(cmp.denom), &(ro.denom));
    ibz_mat_4x4_scalar_mul(&(cmp.basis), &n, &(ro.basis));
    res = res || !quat_lattice_equal(&test, &cmp);
    quat_lideal_conjugate_without_hnf(&conj, &rro, &conj, &alg);
    res = res || !quat_lattice_equal(&rro, &o);
    conj.parent_order = &o;
    quat_lattice_hnf(&(conj.lattice));
    res = res || !quat_lideal_equals(&conj, &lideal, &alg);

    if (res != 0) {
        printf("Quaternion unit test lideal_conj_without_hnf failed\n");
    }
    ibz_finalize(&n);
    quat_alg_elem_finalize(&gen);
    quat_lattice_finalize(&o);
    quat_lattice_finalize(&ro);
    quat_lattice_finalize(&rro);
    quat_lattice_finalize(&test);
    quat_lattice_finalize(&cmp);
    quat_alg_finalize(&alg);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&conj);
    return (res);
}

// void quat_lideal_add(quat_left_ideal_t *sum, const quat_left_ideal_t *I1, const quat_left_ideal_t
// *I2, const quat_alg_t *alg); void quat_lideal_inter(quat_left_ideal_t *intersection, const
// quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg); int
// quat_lideal_equals(const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t
// *alg);
int
quat_test_lideal_add_intersect_equals(void)
{
    int res = 0;

    quat_alg_t alg;
    quat_lattice_t order;
    quat_alg_elem_t gen1, gen2, gen3;
    ibz_t N1, N2, N3;
    quat_left_ideal_t lideal1, lideal2, lideal3, lideal4, lideal5;
    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_init(&order);
    quat_alg_elem_init(&gen1);
    quat_alg_elem_init(&gen2);
    quat_alg_elem_init(&gen3);
    ibz_init(&N1);
    ibz_init(&N2);
    ibz_init(&N3);
    quat_left_ideal_init(&lideal1);
    quat_left_ideal_init(&lideal2);
    quat_left_ideal_init(&lideal3);
    quat_left_ideal_init(&lideal4);
    quat_left_ideal_init(&lideal5);

    quat_lattice_O0_set(&order);

    ibz_set(&gen1.coord[0], 3);
    ibz_set(&gen1.coord[1], 5);
    ibz_set(&gen1.coord[2], 7);
    ibz_set(&gen1.coord[3], 11);
    ibz_set(&N1, 17);
    quat_lideal_create(&lideal1, &gen1, &N1, &order, &alg);

    ibz_set(&gen2.coord[0], -2);
    ibz_set(&gen2.coord[1], 13);
    ibz_set(&gen2.coord[2], -17);
    ibz_set(&gen2.coord[3], 19);
    ibz_set(&N2, 43);
    quat_lideal_create(&lideal2, &gen2, &N2, &order, &alg);

    quat_alg_mul(&gen3, &gen2, &gen1, &alg);
    quat_lideal_create_principal(&lideal3, &gen3, &order, &alg);

    // Union should be the whole ring
    quat_lideal_add(&lideal4, &lideal1, &lideal2, &alg);
    res |= !ibz_is_one(&lideal4.norm);
    res |= !quat_lattice_equal(&lideal4.lattice, &order);

    // Self-intersection should be stable
    quat_lideal_inter(&lideal4, &lideal1, &lideal1, &alg);
    res |= !quat_lideal_equals(&lideal4, &lideal1, &alg);

    // Self-union should be stable
    quat_lideal_add(&lideal4, &lideal1, &lideal1, &alg);
    res |= !quat_lideal_equals(&lideal4, &lideal1, &alg);

    // lideal3 ⊂ lideal1
    quat_lideal_add(&lideal4, &lideal1, &lideal3, &alg);
    res |= !quat_lideal_equals(&lideal4, &lideal1, &alg);
    quat_lideal_inter(&lideal4, &lideal1, &lideal3, &alg);
    res |= !quat_lideal_equals(&lideal4, &lideal3, &alg);

    // Intersection then union should be stable
    quat_lideal_inter(&lideal4, &lideal1, &lideal2, &alg);
    quat_lideal_add(&lideal4, &lideal4, &lideal2, &alg);
    res |= !quat_lideal_equals(&lideal4, &lideal2, &alg);

    // (A ∩ B) ∪ (A ∩ C) == A ∩ (B ∪ C)
    quat_lideal_inter(&lideal4, &lideal1, &lideal2, &alg);
    quat_lideal_inter(&lideal5, &lideal1, &lideal3, &alg);
    quat_lideal_add(&lideal4, &lideal4, &lideal5, &alg);
    quat_lideal_add(&lideal5, &lideal2, &lideal3, &alg);
    quat_lideal_inter(&lideal5, &lideal1, &lideal5, &alg);
    res |= !quat_lideal_equals(&lideal4, &lideal5, &alg);
    res |= ibz_cmp_int32(&lideal4.norm, 17);

    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&gen1);
    quat_alg_elem_finalize(&gen2);
    quat_alg_elem_finalize(&gen3);
    ibz_finalize(&N1);
    ibz_finalize(&N2);
    ibz_finalize(&N3);
    quat_left_ideal_finalize(&lideal1);
    quat_left_ideal_finalize(&lideal2);
    quat_left_ideal_finalize(&lideal3);
    quat_left_ideal_finalize(&lideal4);
    quat_left_ideal_finalize(&lideal5);

    if (res != 0) {
        printf("Quaternion unit test lideal_add_intersect_equals failed\n");
    }
    return (res);
}

// void quat_lideal_inverse_lattice_without_hnf(quat_lattice_t *inv, const quat_left_ideal_t
// *lideal, const quat_alg_t *alg);
int
quat_test_lideal_inverse_lattice_without_hnf(void)
{
    int res = 0;
    quat_left_ideal_t lideal;
    quat_lattice_t inv, prod;
    ibz_t norm;
    quat_alg_elem_t init_helper;
    quat_lattice_t order;
    quat_alg_t alg;
    ibz_init(&norm);
    quat_left_ideal_init(&lideal);
    quat_lattice_init(&inv);
    quat_lattice_init(&prod);
    quat_alg_init_set_ui(&alg, 19);
    quat_lattice_init(&order);
    quat_alg_elem_init(&init_helper);

    ibz_set(&(order.basis[0][0]), 4);
    ibz_set(&(order.basis[0][1]), 0);
    ibz_set(&(order.basis[0][2]), 2);
    ibz_set(&(order.basis[0][3]), 2);
    ibz_set(&(order.basis[1][0]), 0);
    ibz_set(&(order.basis[1][1]), 8);
    ibz_set(&(order.basis[1][2]), 4);
    ibz_set(&(order.basis[1][3]), 3);
    ibz_set(&(order.basis[2][0]), 0);
    ibz_set(&(order.basis[2][1]), 0);
    ibz_set(&(order.basis[2][2]), 2);
    ibz_set(&(order.basis[2][3]), 0);
    ibz_set(&(order.basis[3][0]), 0);
    ibz_set(&(order.basis[3][1]), 0);
    ibz_set(&(order.basis[3][2]), 0);
    ibz_set(&(order.basis[3][3]), 1);
    ibz_set(&(order.denom), 4);
    quat_alg_elem_set(&init_helper, 1, 2, 3, 0, 1);
    ibz_set(&norm, 15);
    quat_lideal_create(&lideal, &init_helper, &norm, &order, &alg);
    quat_lideal_inverse_lattice_without_hnf(&inv, &lideal, &alg);
    quat_lattice_mul(&prod, &(lideal.lattice), &inv, &alg);
    res = res || !quat_lattice_equal(&prod, &order);

    if (res != 0) {
        printf("Quaternion unit test lideal_inverse_lattice_without_hnf failed\n");
    }
    ibz_finalize(&norm);
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&inv);
    quat_lattice_finalize(&prod);
    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&init_helper);
    return (res);
}

// void quat_lideal_right_transporter(quat_lattice_t *trans, const quat_left_ideal_t *lideal1, const
// quat_left_ideal_t *lideal2, const quat_alg_t *alg);
int
quat_test_lideal_right_transporter(void)
{
    int res = 0;
    quat_left_ideal_t lideal1, lideal2;
    quat_lattice_t trans, prod;
    ibz_t norm;
    quat_alg_elem_t init_helper;
    quat_lattice_t order;
    quat_alg_t alg;
    ibz_init(&norm);
    quat_left_ideal_init(&lideal1);
    quat_left_ideal_init(&lideal2);
    quat_lattice_init(&trans);
    quat_lattice_init(&prod);
    quat_alg_init_set_ui(&alg, 19);
    quat_lattice_init(&order);
    quat_alg_elem_init(&init_helper);

    ibz_set(&(order.basis[0][0]), 4);
    ibz_set(&(order.basis[0][1]), 0);
    ibz_set(&(order.basis[0][2]), 2);
    ibz_set(&(order.basis[0][3]), 2);
    ibz_set(&(order.basis[1][0]), 0);
    ibz_set(&(order.basis[1][1]), 8);
    ibz_set(&(order.basis[1][2]), 4);
    ibz_set(&(order.basis[1][3]), 3);
    ibz_set(&(order.basis[2][0]), 0);
    ibz_set(&(order.basis[2][1]), 0);
    ibz_set(&(order.basis[2][2]), 2);
    ibz_set(&(order.basis[2][3]), 0);
    ibz_set(&(order.basis[3][0]), 0);
    ibz_set(&(order.basis[3][1]), 0);
    ibz_set(&(order.basis[3][2]), 0);
    ibz_set(&(order.basis[3][3]), 1);
    ibz_set(&(order.denom), 4);
    quat_alg_elem_set(&init_helper, 1, 2, 3, 0, 1);
    ibz_set(&norm, 15);
    quat_lideal_create(&lideal1, &init_helper, &norm, &order, &alg);
    quat_alg_elem_set(&init_helper, 1, 4, 3, 0, 1);
    ibz_set(&norm, 11);
    quat_lideal_create(&lideal2, &init_helper, &norm, &order, &alg);
    quat_lideal_right_transporter(&trans, &lideal1, &lideal2, &alg);
    quat_lattice_mul(&prod, &(lideal1.lattice), &trans, &alg);
    res = res || !quat_lattice_inclusion(&prod, &(lideal2.lattice));
    quat_lideal_right_transporter(&trans, &lideal2, &lideal1, &alg);
    quat_lattice_mul(&prod, &(lideal2.lattice), &trans, &alg);
    res = res || !quat_lattice_inclusion(&prod, &(lideal1.lattice));
    quat_lideal_right_transporter(&(lideal1.lattice), &lideal2, &lideal1, &alg);
    res = res || !quat_lattice_equal(&trans, &(lideal1.lattice));
    quat_alg_elem_set(&init_helper, 1, 2, 3, 0, 1);
    ibz_set(&norm, 15);
    quat_lideal_create(&lideal1, &init_helper, &norm, &order, &alg);
    quat_lideal_right_transporter(&(lideal2.lattice), &lideal2, &lideal1, &alg);
    res = res || !quat_lattice_equal(&trans, &(lideal2.lattice));
    if (res != 0) {
        printf("Quaternion unit test lideal_right_transporter failed\n");
    }
    ibz_finalize(&norm);
    quat_left_ideal_finalize(&lideal1);
    quat_left_ideal_finalize(&lideal2);
    quat_lattice_finalize(&trans);
    quat_lattice_finalize(&prod);
    quat_alg_finalize(&alg);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&init_helper);
    return (res);
}

// void quat_lideal_right_order(quat_lattice_t *order, const quat_left_ideal_t *lideal, const
// quat_alg_t *alg);
int
quat_test_lideal_right_order(void)
{
    int res = 0;
    ibz_t norm;
    quat_alg_t alg;
    quat_alg_elem_t gen;
    quat_lattice_t order, rorder;
    quat_lattice_t prod;
    quat_alg_elem_t test;
    quat_left_ideal_t lideal;
    quat_lattice_init(&order);
    quat_lattice_init(&rorder);
    quat_lattice_init(&prod);
    quat_alg_init_set_ui(&alg, 19);
    quat_left_ideal_init(&lideal);
    quat_alg_elem_init(&test);
    quat_alg_elem_init(&gen);
    ibz_init(&norm);

    ibz_set(&(order.basis[0][0]), 4);
    ibz_set(&(order.basis[0][1]), 0);
    ibz_set(&(order.basis[0][2]), 2);
    ibz_set(&(order.basis[0][3]), 2);
    ibz_set(&(order.basis[1][0]), 0);
    ibz_set(&(order.basis[1][1]), 8);
    ibz_set(&(order.basis[1][2]), 4);
    ibz_set(&(order.basis[1][3]), 3);
    ibz_set(&(order.basis[2][0]), 0);
    ibz_set(&(order.basis[2][1]), 0);
    ibz_set(&(order.basis[2][2]), 2);
    ibz_set(&(order.basis[2][3]), 0);
    ibz_set(&(order.basis[3][0]), 0);
    ibz_set(&(order.basis[3][1]), 0);
    ibz_set(&(order.basis[3][2]), 0);
    ibz_set(&(order.basis[3][3]), 1);
    ibz_set(&(order.denom), 4);
    quat_alg_elem_set(&gen, 1, 3, 3, 0, 1);
    ibz_set(&norm, 15);
    quat_lideal_create(&lideal, &gen, &norm, &order, &alg);

    quat_lideal_right_order(&rorder, &lideal, &alg);
    // test order is in HNF
    res = res || !ibz_mat_4x4_is_hnf(&(rorder.basis));
    // test order is of dimension 4 (assuming HNF)
    for (int i = 0; i < 4; i++) {
        res = res || ibz_is_zero(&(rorder.basis[i][i]));
    }
    // test order contains 1
    quat_alg_elem_set(&test, 1, 1, 0, 0, 0);
    res = res || !quat_lattice_contains(NULL, &rorder, &test);
    // test it is stable by product
    quat_lattice_mul(&prod, &rorder, &rorder, &alg);
    res = res || !quat_lattice_inclusion(&prod, &rorder);
    // test it is right order of ideal
    quat_lattice_mul(&prod, &(lideal.lattice), &rorder, &alg);
    res = res || !quat_lattice_inclusion(&prod, &(lideal.lattice));

    quat_lattice_O0_set(&order);

    quat_alg_elem_set(&gen, 1, 1, 1, 1, 0);
    ibz_set(&norm, 35);
    quat_lideal_create(&lideal, &gen, &norm, &order, &alg);

    quat_lideal_right_order(&rorder, &lideal, &alg);
    quat_lattice_mul(&prod, &(lideal.lattice), &rorder, &alg);
    // test order is in HNF
    res = res || !ibz_mat_4x4_is_hnf(&(rorder.basis));
    // test order is of dimension 4 (assuming HNF)
    for (int i = 0; i < 4; i++) {
        res = res || ibz_is_zero(&(rorder.basis[i][i]));
    }
    // test order contains 1
    quat_alg_elem_set(&test, 1, 1, 0, 0, 0);
    res = res || !quat_lattice_contains(NULL, &rorder, &test);
    // test it is stable by product
    quat_lattice_mul(&prod, &rorder, &rorder, &alg);
    res = res || !quat_lattice_inclusion(&prod, &rorder);
    // test it is right order of ideal
    quat_lattice_mul(&prod, &(lideal.lattice), &rorder, &alg);
    res = res || !quat_lattice_inclusion(&prod, &(lideal.lattice));

    if (res != 0) {
        printf("Quaternion unit test lideal_right_order failed\n");
    }
    quat_alg_elem_finalize(&test);
    quat_lattice_finalize(&order);
    quat_lattice_finalize(&rorder);
    quat_lattice_finalize(&prod);
    quat_left_ideal_finalize(&lideal);
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&norm);
    return (res);
}

/********************* Functions from quaternion_tools.c **********************/

// int quat_order_discriminant(ibz_t *disc, const quat_lattice_t *order, const quat_alg_t *alg);
int
quat_test_lideal_order_discriminant(void)
{
    int res = 0;
    ibz_t disc;
    quat_alg_t alg;
    quat_lattice_t order;
    ibz_init(&disc);
    quat_lattice_init(&order);
    quat_alg_init_set_ui(&alg, 43);
    quat_lattice_O0_set(&order);
    quat_order_discriminant(&disc, &order, &alg);
    res = res | ibz_cmp(&(alg.p), &disc);
    if (res != 0) {
        printf("Quaternion unit test lideal_order_discriminant failed\n");
    }
    ibz_finalize(&disc);
    quat_lattice_finalize(&order);
    quat_alg_finalize(&alg);
    return (res);
}

// int quat_order_is_maximal(const quat_lattice_t *order, const quat_alg_t *alg);
int
quat_test_lideal_order_is_maximal(void)
{
    int res = 0;
    quat_alg_t alg;
    quat_lattice_t order;
    quat_lattice_init(&order);
    quat_alg_init_set_ui(&alg, 43);
    quat_lattice_O0_set(&order);
    res = res | !quat_order_is_maximal(&order, &alg);
    ibz_mat_4x4_identity(&(order.basis));
    ibz_set(&(order.denom), 1);
    res = res | quat_order_is_maximal(&order, &alg);
    if (res != 0) {
        printf("Quaternion unit test lideal__order_is_maximal failed\n");
    }
    quat_lattice_finalize(&order);
    quat_alg_finalize(&alg);
    return (res);
}

// void quat_lideal_class_gram(ibz_mat_4x4_t *G, const quat_left_ideal_t *lideal, const quat_alg_t
// *alg);
int
quat_test_lideal_class_gram()
{
    int res = 0;
    ibz_t norm, norm1, norm2, cmp;
    ibz_mat_4x4_t gram, gram_l;
    quat_left_ideal_t lideal;
    quat_lattice_t order;
    quat_alg_elem_t elem, elem1, elem2;
    ibz_vec_4_t vec, vec1, vec2;
    quat_alg_t alg;
    ibz_init(&norm);
    ibz_init(&norm1);
    ibz_init(&norm2);
    ibz_init(&cmp);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&gram_l);
    ibz_vec_4_init(&vec);
    ibz_vec_4_init(&vec1);
    ibz_vec_4_init(&vec2);
    quat_left_ideal_init(&lideal);
    quat_lattice_init(&order);
    quat_alg_elem_init(&elem);
    quat_alg_elem_init(&elem1);
    quat_alg_elem_init(&elem2);
    quat_alg_init_set_ui(&alg, 103);

    // setup
    quat_lattice_O0_set(&order);
    ibz_set(&cmp, 101);
    quat_alg_elem_set(&elem, 1, 9, 4, 1, 1);
    quat_lideal_create(&lideal, &elem, &cmp, &order, &alg);

    // test relation to lattice_gram output
    quat_lideal_class_gram(&gram, &lideal, &alg);
    ibz_mul(&cmp, &(lideal.lattice.denom), &(lideal.lattice.denom));
    ibz_mul(&cmp, &cmp, &(lideal.norm));
    ibz_mat_4x4_scalar_mul(&gram, &cmp, &gram);
    quat_lattice_gram(&gram_l, &(lideal.lattice), &alg);
    res = res | !ibz_mat_4x4_equal(&gram, &gram_l);

    // repeat lattice_gram test
    quat_alg_elem_set(&elem1, 2, 360, 149, 1, 0);
    quat_alg_elem_set(&elem2, 2, 53, 360, 0, 1);
    int ok = quat_lattice_contains(&vec1, &(lideal.lattice), &elem1);
    ok = ok && quat_lattice_contains(&vec2, &(lideal.lattice), &elem2);
    assert(ok);
    quat_alg_conj(&elem2, &elem2);
    quat_alg_mul(&elem1, &elem1, &elem2, &alg);
    ibz_mul(&norm1, &(elem1.coord[0]), &ibz_const_two);
    ibz_div(&norm1, &cmp, &norm1, &(elem1.denom));

    ibz_mat_4x4_eval(&vec1, &gram, &vec1);
    ibz_set(&norm2, 0);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&cmp, &(vec1[i]), &(vec2[i]));
        ibz_add(&norm2, &norm2, &cmp);
    }
    ibz_div(&norm2, &cmp, &norm2, &(lideal.lattice.denom));
    ibz_div(&norm2, &cmp, &norm2, &(lideal.lattice.denom));

    res = res | !(0 == ibz_cmp(&norm1, &norm2));

    // test for norm
    quat_lideal_class_gram(&gram, &lideal, &alg);
    quat_lattice_contains(&vec, &(lideal.lattice), &elem);
    quat_alg_norm(&norm, &cmp, &elem, &alg);
    assert(ibz_is_one(&cmp));
    ibz_mul(&norm, &norm, &ibz_const_two);
    ibz_div(&norm, &cmp, &norm, &(lideal.norm));

    assert(ibz_is_zero(&cmp));
    quat_qf_eval(&cmp, &gram, &vec);
    res = res | !(0 == ibz_cmp(&cmp, &norm));

    if (res != 0) {
        printf("Quaternion unit test lideal_class_gram failed\n");
    }
    ibz_finalize(&norm);
    ibz_finalize(&norm1);
    ibz_finalize(&norm2);
    ibz_finalize(&cmp);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&gram_l);
    ibz_vec_4_finalize(&vec);
    ibz_vec_4_finalize(&vec1);
    ibz_vec_4_finalize(&vec2);
    quat_left_ideal_finalize(&lideal);
    quat_lattice_finalize(&order);
    quat_alg_elem_finalize(&elem);
    quat_alg_elem_finalize(&elem1);
    quat_alg_elem_finalize(&elem2);
    quat_alg_finalize(&alg);
    return (res);
}

// run all previous tests
int
quat_test_lideal(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of ideal and order functions\n");
    res = res | quat_test_lideal_norm();
    res = res | quat_test_lideal_copy();
    res = res | quat_test_lideal_create_principal();
    res = res | quat_test_lideal_create_from_primitive();
    res = res | quat_test_lideal_generator();
    res = res | quat_test_lideal_mul();
    res = res | quat_test_lideal_conj_without_hnf();
    res = res | quat_test_lideal_add_intersect_equals();
    res = res | quat_test_lideal_inverse_lattice_without_hnf();
    res = res | quat_test_lideal_right_transporter();
    res = res | quat_test_lideal_right_order();
    res = res | quat_test_lideal_order_discriminant();
    res = res | quat_test_lideal_order_is_maximal();
    res = res | quat_test_lideal_class_gram();
    return (res);
}
