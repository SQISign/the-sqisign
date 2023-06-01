#include "quaternion_tests.h"

//void quat_lideal_create_principal(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const quat_order_t *order, const quat_alg_t *alg);
int quat_test_lideal_create_principal(){
    int res = 0;
    
    // Example from https://github.com/SQISign/sqisign-nist/issues/38#issuecomment-1554585079
    quat_alg_t alg;
    quat_order_t lat;
    quat_alg_elem_t gamma;
    quat_left_ideal_t I;
    quat_alg_init_set_ui(&alg, 367);
    quat_order_init(&lat);
    quat_alg_elem_init(&gamma);
    quat_left_ideal_init(&I);
    ibz_set(&lat.denom, 2);
    ibz_set(&lat.basis[0][0], 2);
    ibz_set(&lat.basis[1][1], 2);
    ibz_set(&lat.basis[1][2], 1);
    ibz_set(&lat.basis[2][2], 1);
    ibz_set(&lat.basis[3][3], 1);
    ibz_set(&lat.basis[0][3], 1);
    ibz_set(&gamma.coord[0], 219);
    ibz_set(&gamma.coord[1], 200);
    ibz_set(&gamma.coord[2], 78);
    ibz_set(&gamma.coord[3], -1);

    quat_lideal_create_principal(&I, &gamma, &lat, &alg);
    
    res |= I.parent_order != &lat;
    res |= ibz_cmp_si(&I.norm, 2321156);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_one);
    
    res |= ibz_cmp_si(&I.lattice.basis[0][0], 1160578);
    res |= ibz_cmp_si(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[1][1], 1160578);
    res |= ibz_cmp_si(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][2], 310126);
    res |= ibz_cmp_si(&I.lattice.basis[1][2], 182529);
    res |= ibz_cmp_si(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_si(&I.lattice.basis[3][2], 0);
    
    res |= ibz_cmp_si(&I.lattice.basis[0][3], 978049);
    res |= ibz_cmp_si(&I.lattice.basis[1][3], 310126);
    res |= ibz_cmp_si(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][3], 1);

    // same test, just with gamma not reduced
    ibz_set(&gamma.coord[0], 438);
    ibz_set(&gamma.coord[1], 400);
    ibz_set(&gamma.coord[2], 156);
    ibz_set(&gamma.coord[3], -2);
    ibz_set(&gamma.denom, 2);

    quat_lideal_create_principal(&I, &gamma, &lat, &alg);

    res |= I.parent_order != &lat;
    res |= ibz_cmp_si(&I.norm, 2321156);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_one);

    res |= ibz_cmp_si(&I.lattice.basis[0][0], 1160578);
    res |= ibz_cmp_si(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[1][1], 1160578);
    res |= ibz_cmp_si(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][2], 310126);
    res |= ibz_cmp_si(&I.lattice.basis[1][2], 182529);
    res |= ibz_cmp_si(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_si(&I.lattice.basis[3][2], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][3], 978049);
    res |= ibz_cmp_si(&I.lattice.basis[1][3], 310126);
    res |= ibz_cmp_si(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][3], 1);

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
    res |= ibz_cmp_si(&I.norm, 2321156);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_one);
    
    res |= ibz_cmp_si(&I.lattice.basis[0][0], 1160578);
    res |= ibz_cmp_si(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[1][1], 1160578);
    res |= ibz_cmp_si(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][1], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][2], 310126);
    res |= ibz_cmp_si(&I.lattice.basis[1][2], 182529);
    res |= ibz_cmp_si(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_si(&I.lattice.basis[3][2], 0);
    
    res |= ibz_cmp_si(&I.lattice.basis[0][3], 978049);
    res |= ibz_cmp_si(&I.lattice.basis[1][3], 310126);
    res |= ibz_cmp_si(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][3], 1);


    quat_alg_finalize(&alg);
    quat_order_finalize(&lat);
    quat_alg_elem_finalize(&gamma);
    quat_left_ideal_finalize(&I);

    if (res != 0){
        printf("Quaternion unit test lideal_create_principal failed\n");
    }
    return(res);
}

//void quat_lideal_create_from_primitive(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const ibz_t *N, const quat_order_t *order, const quat_alg_t *alg); 
int quat_test_lideal_create_from_primitive(){
    int res = 0;

    // Example from https://github.com/SQISign/sqisign-nist/issues/38#issuecomment-1554585079
    quat_alg_t alg;
    quat_order_t lat;
    quat_alg_elem_t gamma;
    ibz_t N;
    quat_left_ideal_t I;
    quat_alg_init_set_ui(&alg, 367);
    quat_order_init(&lat);
    quat_alg_elem_init(&gamma);
    ibz_init(&N);
    quat_left_ideal_init(&I);
    ibz_set(&lat.denom, 2);
    ibz_set(&lat.basis[0][0], 2);
    ibz_set(&lat.basis[1][1], 2);
    ibz_set(&lat.basis[1][2], 1);
    ibz_set(&lat.basis[2][2], 1);
    ibz_set(&lat.basis[3][3], 1);
    ibz_set(&lat.basis[0][3], 1);
    ibz_set(&gamma.coord[0], 219);
    ibz_set(&gamma.coord[1], 200);
    ibz_set(&gamma.coord[2], 78);
    ibz_set(&gamma.coord[3], -1);
    ibz_set(&N, 31);

    quat_lideal_create_from_primitive(&I, &gamma, &N, &lat, &alg);
    
    res |= I.parent_order != &lat;
    res |= ibz_cmp(&I.norm, &N);
    res |= ibz_cmp(&I.lattice.denom, &ibz_const_two);

    res |= ibz_cmp_si(&I.lattice.basis[0][0], 62);
    res |= ibz_cmp_si(&I.lattice.basis[1][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[2][0], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][0], 0);

    res |= ibz_cmp_si(&I.lattice.basis[0][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[1][1], 62);
    res |= ibz_cmp_si(&I.lattice.basis[2][1], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][1], 0);
    
    res |= ibz_cmp_si(&I.lattice.basis[0][2], 2);
    res |= ibz_cmp_si(&I.lattice.basis[1][2], 1);
    res |= ibz_cmp_si(&I.lattice.basis[2][2], 1);
    res |= ibz_cmp_si(&I.lattice.basis[3][2], 0);
    
    res |= ibz_cmp_si(&I.lattice.basis[0][3], 61);
    res |= ibz_cmp_si(&I.lattice.basis[1][3], 2);
    res |= ibz_cmp_si(&I.lattice.basis[2][3], 0);
    res |= ibz_cmp_si(&I.lattice.basis[3][3], 1);

    quat_alg_finalize(&alg);
    quat_order_finalize(&lat);
    quat_alg_elem_finalize(&gamma);
    ibz_finalize(&N);
    quat_left_ideal_finalize(&I);
                   
    if (res != 0){
        printf("Quaternion unit test lideal_create_from_primitive failed\n");
    }
    return(res);
}

//void quat_lideal_make_primitive_then_create(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const ibz_t *N, const quat_order_t *order, const quat_alg_t *alg); 
int quat_test_lideal_make_primitive_then_create(){
    int res = 0;
    ibz_t n, cnt;
    quat_alg_t alg;
    quat_alg_elem_t gen, x;
    quat_alg_coord_t prim;
    quat_order_t order;
    quat_left_ideal_t lideal,cmp;
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&cmp);
    quat_alg_init_set_ui(&alg, 103);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&x);
    quat_alg_coord_init(&prim);
    quat_order_init(&order);
    ibz_init(&n);
    ibz_init(&cnt);

    
    ibz_set(&order.denom, 2);
    ibz_set(&order.basis[0][0], 2);
    ibz_set(&order.basis[1][1], 2);
    ibz_set(&order.basis[1][2], 1);
    ibz_set(&order.basis[2][2], 1);
    ibz_set(&order.basis[3][3], 1);
    ibz_set(&order.basis[0][3], 1);

    ibz_set(&x.coord[0], 6);
    ibz_set(&x.coord[1], 10);
    ibz_set(&x.coord[2], 14);
    ibz_set(&x.coord[3], 22);
    ibz_set(&n, 17);

    quat_lideal_make_primitive_then_create(&lideal,&x,&n,&order,&alg);
    quat_alg_make_primitive(&prim,&cnt,&x,&order,&alg);
    ibz_mat_4x4_eval(&(gen.coord),&(order.basis),&prim);
    ibz_copy(&gen.denom,&order.denom);
    quat_lideal_create_from_primitive(&cmp,&gen,&n,&order,&alg);
    res = res ||  !quat_lideal_equals(&lideal,&cmp,&alg);

    if (res != 0){
        printf("Quaternion unit test lideal_make_primitive_then_create failed\n");
    }
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&x);
    quat_alg_coord_finalize(&prim);
    quat_order_finalize(&order);
    ibz_finalize(&n);
    ibz_finalize(&cnt);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&cmp);
    return(res);
}

//void quat_lideal_random_2e(quat_left_ideal_t *lideal, const quat_order_t *order, const quat_alg_t *alg, int64_t e, unsigned char n);
int quat_test_lideal_random_2e(){
    int res = 0;

    quat_alg_t alg;
    quat_order_t order;
    quat_left_ideal_t lideal;
    quat_alg_init_set_ui(&alg, 103);
    quat_order_init(&order);
    quat_left_ideal_init(&lideal);

    ibz_set(&order.denom, 2);
    ibz_set(&order.basis[0][0], 2);
    ibz_set(&order.basis[1][1], 2);
    ibz_set(&order.basis[1][2], 1);
    ibz_set(&order.basis[2][2], 1);
    ibz_set(&order.basis[3][3], 1);
    ibz_set(&order.basis[0][3], 1);

    // Pretty much self-testing (thanks to embedded asserts)
    res |= quat_lideal_random_2e(&lideal, &order, &alg, 100, 10) == 0;
    
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
    quat_left_ideal_finalize(&lideal);
    
    if (res != 0){
        printf("Quaternion unit test lideal_random_2e failed\n");
    }
    return(res);
}

//int quat_lideal_generator(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t *alg, int bound) 
int quat_test_lideal_generator(){
    int res = 0;

    quat_alg_t alg;
    quat_order_t order;
    quat_alg_elem_t gen;
    ibz_t N;
    quat_left_ideal_t lideal, lideal2;
    quat_alg_init_set_ui(&alg, 103);
    quat_order_init(&order);
    quat_alg_elem_init(&gen);
    ibz_init(&N);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal2);

    ibz_set(&order.denom, 2);
    ibz_set(&order.basis[0][0], 2);
    ibz_set(&order.basis[1][1], 2);
    ibz_set(&order.basis[1][2], 1);
    ibz_set(&order.basis[2][2], 1);
    ibz_set(&order.basis[3][3], 1);
    ibz_set(&order.basis[0][3], 1);

    ibz_set(&gen.coord[0], 3);
    ibz_set(&gen.coord[1], 5);
    ibz_set(&gen.coord[2], 7);
    ibz_set(&gen.coord[3], 11);
    ibz_set(&N, 17);

    quat_lideal_create_from_primitive(&lideal, &gen, &N, &order, &alg);
    
    // Try to regenerate the same ideal
    for (int i = 0; i <= 10; i++) {
        res |= !quat_lideal_generator(&gen, &lideal, &alg,0);
        quat_lideal_create_from_primitive(&lideal2, &gen, &N, &order, &alg);
        res |= !quat_lideal_equals(&lideal, &lideal2, &alg);
    }
    
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&N);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal2);
    
    if (res != 0){
        printf("Quaternion unit test lideal_generator failed\n");
    }
    return(res);
}

//int quat_lideal_generator_coprime(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const ibz_t *n, const quat_alg_t *alg, int bound);
int quat_test_lideal_generator_coprime(){
    int res = 0;

    quat_alg_t alg;
    quat_order_t order, order2;
    quat_alg_elem_t gen;
    ibz_t N, M, gcd, norm_int, prod;
    ibq_t norm;
    quat_left_ideal_t lideal, lideal2;
    quat_alg_init_set_ui(&alg, 103);
    quat_order_init(&order);
    quat_order_init(&order2);
    quat_alg_elem_init(&gen);
    ibz_init(&N);
    ibz_init(&M);
    ibz_init(&norm_int);
    ibz_init(&prod);
    ibz_init(&gcd);
    ibq_init(&norm);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal2);

    ibz_set(&order.denom, 2);
    ibz_set(&order.basis[0][0], 2);
    ibz_set(&order.basis[1][1], 2);
    ibz_set(&order.basis[1][2], 1);
    ibz_set(&order.basis[2][2], 1);
    ibz_set(&order.basis[3][3], 1);
    ibz_set(&order.basis[0][3], 1);

    ibz_set(&gen.coord[0], 3);
    ibz_set(&gen.coord[1], 5);
    ibz_set(&gen.coord[2], 7);
    ibz_set(&gen.coord[3], 11);
    ibz_set(&N, 5);
    // product of primes < 50
    ibz_set(&M, 614889782588491410l);

    quat_lideal_create_from_primitive(&lideal, &gen, &N, &order, &alg);
    
    res |= !quat_lideal_generator_coprime(&gen, &lideal, &M, &alg,0);
    quat_alg_norm(&norm, &gen, &alg);
    res |= !ibq_to_ibz(&norm_int, &norm);
    // case with coprimality
    ibz_gcd(&gcd, &norm_int, &M);
    res |= !ibz_is_one(&gcd);
    // test other formula
    ibz_mul(&prod,&M,&M);
    ibz_gcd(&prod,&prod,&norm_int);
    ibz_gcd(&gcd, &(lideal.norm), &M);
    res |= ibz_cmp(&gcd,&prod);

    quat_lideal_create_from_primitive(&lideal2, &gen, &N, &order, &alg);
    res |= !quat_lideal_equals(&lideal, &lideal2, &alg);
    res |= !quat_alg_is_primitive(&gen,&(lideal.lattice), &alg);
    
    

    ibz_set(&(order.denom),1);
    quat_alg_elem_set(&gen,1,2,4,2,6);

    quat_lideal_create_from_primitive(&lideal, &gen, &N, &order, &alg);
    
    ibz_set(&M,15);
    res |= !quat_lideal_generator_coprime(&gen, &lideal, &M, &alg,0);
    quat_alg_norm(&norm, &gen, &alg);
    res |= !ibq_to_ibz(&norm_int, &norm);
    // test other formula
    ibz_mul(&prod,&M,&M);
    ibz_gcd(&prod,&prod,&norm_int);
    ibz_gcd(&gcd, &(lideal.norm), &M);
    res |= ibz_cmp(&gcd,&prod);

    quat_lideal_create_from_primitive(&lideal2, &gen, &N, &order, &alg);
    res |= !quat_lideal_equals(&lideal, &lideal2, &alg);
    res |= !quat_alg_is_primitive(&gen,&(lideal.lattice), &alg);
    
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
    quat_order_finalize(&order2);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&N);
    ibz_finalize(&M);
    ibz_finalize(&norm_int);
    ibz_finalize(&prod);
    ibz_finalize(&gcd);
    ibq_finalize(&norm);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal2);
    
    if (res != 0){
        printf("Quaternion unit test lideal_generator_coprime failed\n");
    }
    return(res);
}

//int quat_lideal_mul(quat_left_ideal_t *product, const quat_left_ideal_t *lideal, const quat_alg_elem_t *alpha, const quat_alg_t *alg, int bound); 
int quat_test_lideal_mul(){
    int res = 0;

    quat_alg_t alg;
    quat_order_t order;
    quat_alg_elem_t gen1, gen2, gen_prod;
    quat_left_ideal_t lideal, lideal2;
    quat_alg_init_set_ui(&alg, 103);
    quat_order_init(&order);
    quat_alg_elem_init(&gen1);
    quat_alg_elem_init(&gen2);
    quat_alg_elem_init(&gen_prod);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal2);

    ibz_set(&order.denom, 2);
    ibz_set(&order.basis[0][0], 2);
    ibz_set(&order.basis[1][1], 2);
    ibz_set(&order.basis[1][2], 1);
    ibz_set(&order.basis[2][2], 1);
    ibz_set(&order.basis[3][3], 1);
    ibz_set(&order.basis[0][3], 1);

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
    res |= !quat_lideal_mul(&lideal, &lideal, &gen2, &alg,0);
    quat_alg_mul(&gen_prod, &gen1, &gen2, &alg);
    quat_lideal_create_principal(&lideal2, &gen_prod, &order, &alg);
    res |= !quat_lideal_equals(&lideal, &lideal2, &alg);

    // Same test, but with non-integral gen2
    ibz_set(&(gen2.denom),2);
    quat_lideal_create_principal(&lideal, &gen1, &order, &alg);
    res |= !quat_lideal_mul(&lideal, &lideal, &gen2, &alg,0);
    quat_alg_mul(&gen_prod, &gen1, &gen2, &alg);
    quat_lideal_create_principal(&lideal2, &gen_prod, &order, &alg);
    res |= !quat_lideal_equals(&lideal, &lideal2, &alg);
    
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
    quat_alg_elem_finalize(&gen1);
    quat_alg_elem_finalize(&gen2);
    quat_alg_elem_finalize(&gen_prod);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal2);
    
    if (res != 0){
        printf("Quaternion unit test lideal_mul failed\n");
    }
    return(res);
}

//void quat_lideal_add(quat_left_ideal_t *sum, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
//void quat_lideal_inter(quat_left_ideal_t *intersection, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
//int quat_lideal_equals(const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
int quat_test_lideal_add_intersect_equals(){
    int res = 0;

    quat_alg_t alg;
    quat_order_t order;
    quat_alg_elem_t gen1, gen2, gen3;
    ibz_t N1, N2, N3;
    quat_left_ideal_t lideal1, lideal2, lideal3, lideal4, lideal5;
    quat_alg_init_set_ui(&alg, 103);
    quat_order_init(&order);
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

    ibz_set(&order.denom, 2);
    ibz_set(&order.basis[0][0], 2);
    ibz_set(&order.basis[1][1], 2);
    ibz_set(&order.basis[1][2], 1);
    ibz_set(&order.basis[2][2], 1);
    ibz_set(&order.basis[3][3], 1);
    ibz_set(&order.basis[0][3], 1);

    ibz_set(&gen1.coord[0], 3);
    ibz_set(&gen1.coord[1], 5);
    ibz_set(&gen1.coord[2], 7);
    ibz_set(&gen1.coord[3], 11);
    ibz_set(&N1, 17);
    quat_lideal_create_from_primitive(&lideal1, &gen1, &N1, &order, &alg);

    ibz_set(&gen2.coord[0], -2);
    ibz_set(&gen2.coord[1], 13);
    ibz_set(&gen2.coord[2], -17);
    ibz_set(&gen2.coord[3], 19);
    ibz_set(&N2, 43);
    quat_lideal_create_from_primitive(&lideal2, &gen2, &N2, &order, &alg);

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
    res |= ibz_cmp_si(&lideal4.norm, 17);
    
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
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
    
    if (res != 0){
        printf("Quaternion unit test lideal_add_intersect_equals failed\n");
    }
    return(res);
}

//int quat_lideal_isom(quat_alg_elem_t *iso, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
int quat_test_lideal_isom(){
    int res = 0;
    quat_left_ideal_t lideal1,lideal2,prod;
    ibz_t norm;
    quat_alg_elem_t init_helper, isom;
    quat_order_t order, order2;
    quat_alg_t alg;
    ibz_init(&norm);
    quat_left_ideal_init(&lideal1);
    quat_left_ideal_init(&lideal2);
    quat_left_ideal_init(&prod);
    quat_alg_init_set_ui(&alg,19);
    quat_order_init(&order);
    quat_order_init(&order2);
    quat_alg_elem_init(&isom);
    quat_alg_elem_init(&init_helper);

    ibz_set(&(order.basis[0][0]),4);
    ibz_set(&(order.basis[0][1]),0);
    ibz_set(&(order.basis[0][2]),2);
    ibz_set(&(order.basis[0][3]),2);
    ibz_set(&(order.basis[1][0]),0);
    ibz_set(&(order.basis[1][1]),8);
    ibz_set(&(order.basis[1][2]),4);
    ibz_set(&(order.basis[1][3]),3);
    ibz_set(&(order.basis[2][0]),0);
    ibz_set(&(order.basis[2][1]),0);
    ibz_set(&(order.basis[2][2]),2);
    ibz_set(&(order.basis[2][3]),0);
    ibz_set(&(order.basis[3][0]),0);
    ibz_set(&(order.basis[3][1]),0);
    ibz_set(&(order.basis[3][2]),0);
    ibz_set(&(order.basis[3][3]),1);
    ibz_set(&(order.denom),4);
    quat_alg_elem_set(&init_helper,1,2,6,2,1);
    quat_lideal_create_principal(&lideal1,&init_helper,&order,&alg);
    quat_alg_elem_set(&init_helper,3,1,6,5,1);
    res |= !quat_lideal_mul(&lideal2,&lideal1,&init_helper,&alg,0);
    res = res || !quat_lideal_isom(&isom,&lideal1,&lideal2,&alg);
    if(!res){
        res |= !quat_lideal_mul(&prod,&lideal1,&isom,&alg,0);
        res = res || !quat_lideal_equals(&prod,&lideal2,&alg);
    }

    // not same order
    ibz_set(&(order2.basis[0][0]),2);
    ibz_set(&(order2.basis[0][1]),0);
    ibz_set(&(order2.basis[0][2]),1);
    ibz_set(&(order2.basis[0][3]),0);
    ibz_set(&(order2.basis[1][0]),0);
    ibz_set(&(order2.basis[1][1]),2);
    ibz_set(&(order2.basis[1][2]),0);
    ibz_set(&(order2.basis[1][3]),1);
    ibz_set(&(order2.basis[2][0]),0);
    ibz_set(&(order2.basis[2][1]),0);
    ibz_set(&(order2.basis[2][2]),1);
    ibz_set(&(order2.basis[2][3]),0);
    ibz_set(&(order2.basis[3][0]),0);
    ibz_set(&(order2.basis[3][1]),0);
    ibz_set(&(order2.basis[3][2]),0);
    ibz_set(&(order2.basis[3][3]),1);
    ibz_set(&(order2.denom),2);
    quat_alg_elem_set(&init_helper,1,2,6,2,1);
    quat_lideal_create_principal(&lideal2,&init_helper,&order2,&alg);
    res = res || quat_lideal_isom(&isom,&lideal1,&lideal2,&alg);

    // should add test for not isomorphic ideals
    
    if (res != 0){
        printf("Quaternion unit test lideal_isom failed\n");
    }
    quat_left_ideal_finalize(&lideal1);
    quat_left_ideal_finalize(&lideal2);
    quat_left_ideal_finalize(&prod);
    ibz_finalize(&norm);
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
    quat_order_finalize(&order2);
    quat_alg_elem_finalize(&isom);
    quat_alg_elem_finalize(&init_helper);
    return(res);
}

//void quat_lideal_right_order(quat_order_t *order, const quat_left_ideal_t *lideal, const quat_alg_t *alg);
int quat_test_lideal_right_order(){
    int res = 0;
    ibz_t norm;
    quat_alg_t alg;
    quat_alg_elem_t gen;
    quat_order_t order,rorder;
    quat_lattice_t prod;
    quat_alg_elem_t test;
    quat_left_ideal_t lideal, cmp;
    quat_order_init(&order);
    quat_order_init(&rorder);
    quat_lattice_init(&prod);
    quat_alg_init_set_ui(&alg,19);
    quat_left_ideal_init(&lideal);
    quat_alg_elem_init(&test);
    quat_alg_elem_init(&gen);

    ibz_set(&(order.basis[0][0]),4);
    ibz_set(&(order.basis[0][1]),0);
    ibz_set(&(order.basis[0][2]),2);
    ibz_set(&(order.basis[0][3]),2);
    ibz_set(&(order.basis[1][0]),0);
    ibz_set(&(order.basis[1][1]),8);
    ibz_set(&(order.basis[1][2]),4);
    ibz_set(&(order.basis[1][3]),3);
    ibz_set(&(order.basis[2][0]),0);
    ibz_set(&(order.basis[2][1]),0);
    ibz_set(&(order.basis[2][2]),2);
    ibz_set(&(order.basis[2][3]),0);
    ibz_set(&(order.basis[3][0]),0);
    ibz_set(&(order.basis[3][1]),0);
    ibz_set(&(order.basis[3][2]),0);
    ibz_set(&(order.basis[3][3]),1);
    ibz_set(&(order.denom),4);
    quat_alg_elem_set(&gen,1,2,1,8,-8);
    quat_lideal_create_principal(&lideal,&gen,&order,&alg);

    quat_lideal_right_order(&rorder,&lideal,&alg);
    // test order is in HNF
    res = res || !ibz_mat_4x4_is_hnf(&(rorder.basis));
    // test order is of dimension 4 (assuming HNF)
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(rorder.denom),&(rorder.basis[0][i]),&(rorder.basis[1][i]),&(rorder.basis[2][i]),&(rorder.basis[3][i]));
        res = res || quat_alg_elem_is_zero(&test);
    }
    // test order contains 1
    quat_alg_elem_set(&test,1,1,0,0,0);
    res = res || !quat_lattice_contains(NULL,&rorder,&test,&alg);
    // test it is right order of ideal
    quat_lattice_mul(&prod,&(lideal.lattice),&rorder,&alg);
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(prod.denom),&(prod.basis[0][i]),&(prod.basis[1][i]),&(prod.basis[2][i]),&(prod.basis[3][i]));
        res = res || !quat_lattice_contains(NULL,&(lideal.lattice),&test,&alg);
    }

    ibz_set(&(order.basis[0][0]),2);
    ibz_set(&(order.basis[0][1]),0);
    ibz_set(&(order.basis[0][2]),1);
    ibz_set(&(order.basis[0][3]),0);
    ibz_set(&(order.basis[1][0]),0);
    ibz_set(&(order.basis[1][1]),2);
    ibz_set(&(order.basis[1][2]),0);
    ibz_set(&(order.basis[1][3]),1);
    ibz_set(&(order.basis[2][0]),0);
    ibz_set(&(order.basis[2][1]),0);
    ibz_set(&(order.basis[2][2]),1);
    ibz_set(&(order.basis[2][3]),0);
    ibz_set(&(order.basis[3][0]),0);
    ibz_set(&(order.basis[3][1]),0);
    ibz_set(&(order.basis[3][2]),0);
    ibz_set(&(order.basis[3][3]),1);
    ibz_set(&(order.denom),2);
    quat_alg_elem_set(&gen,1,1,2,8,8);
    quat_lideal_create_principal(&lideal,&gen,&order,&alg);

    quat_lideal_right_order(&rorder,&lideal,&alg);
    quat_lattice_mul(&prod,&(lideal.lattice),&rorder,&alg);
    // test order is in HNF
    res = res || !ibz_mat_4x4_is_hnf(&(rorder.basis));
    // test order is of dimension 4 (assuming HNF)
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(rorder.denom),&(rorder.basis[0][i]),&(rorder.basis[1][i]),&(rorder.basis[2][i]),&(rorder.basis[3][i]));
        res = res || quat_alg_elem_is_zero(&test);
    }
    // test order contains 1
    quat_alg_elem_set(&test,1,1,0,0,0);
    res = res || !quat_lattice_contains(NULL,&rorder,&test,&alg);
    // test it is right order of ideal
    quat_lattice_mul(&prod,&(lideal.lattice),&rorder,&alg);
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(prod.denom),&(prod.basis[0][i]),&(prod.basis[1][i]),&(prod.basis[2][i]),&(prod.basis[3][i]));
        res = res || !quat_lattice_contains(NULL,&(lideal.lattice),&test,&alg);
    }

    if (res != 0){
        printf("Quaternion unit test lideal_right_order failed\n");
    }
    quat_alg_elem_finalize(&test);
    quat_order_finalize(&order);
    quat_order_finalize(&rorder);
    quat_lattice_finalize(&prod);
    quat_left_ideal_finalize(&lideal);
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&gen);
    return(res);
}

//void quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced, ibz_mat_4x4_t *gram, const quat_left_ideal_t *lideal, const quat_alg_t *alg); //replaces lideal_lll
int quat_test_lideal_reduce_basis(){
    int res = 0;
    ibz_mat_4x4_t red, gram, prod;
    quat_left_ideal_t lideal;
    quat_lattice_t test;
    quat_alg_t alg;
    quat_alg_elem_t init_helper;
    quat_order_t order;
    ibq_t coeff;
    ibz_t num,denom;
    ibz_init(&num);
    ibz_init(&denom);
    ibq_init(&coeff);
    ibz_mat_4x4_init(&prod);
    quat_lattice_init(&test);
    quat_alg_elem_init(&init_helper);
    quat_order_init(&order);
    quat_alg_init_set_ui(&alg,19);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);
    quat_left_ideal_init(&lideal);
    ibz_mat_4x4_identity(&(order.basis));
    ibz_set(&(order.denom),2);
    quat_alg_elem_set(&init_helper,1,1,2,8,8);
    quat_lideal_create_principal(&lideal,&init_helper,&order,&alg);
    quat_lattice_reduce_denom(&(lideal.lattice),&(lideal.lattice));
    quat_lideal_reduce_basis(&red,&gram,&lideal,&alg);
    // test if red is lll-reduced, assuming he used coeff was at least 3/4
    ibz_set(&num,99);
    ibz_set(&denom,100);
    ibq_set(&coeff,&num,&denom);
    res = res || !quat_dim4_lll_verify(&red,&coeff,&(alg.p));
    // test reduced and lideal generate same lattice
    ibz_mat_4x4_copy(&(test.basis),&red);
    ibz_copy(&(test.denom),&(lideal.lattice.denom));
    quat_lattice_hnf(&test);
    res = res || !quat_lattice_equal(&(lideal.lattice),&test);
    // test gram matrix is gram matrix
    ibz_mat_4x4_transpose(&prod,&red);
    ibz_mat_4x4_mul(&prod,&prod,&(alg.gram));
    ibz_mat_4x4_mul(&prod,&prod,&red);
    res = res || !ibz_mat_4x4_equal(&prod,&gram);
    if (res != 0){
        printf("Quaternion unit test lideal_reduce_basis failed\n");
    }
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibq_finalize(&coeff);
    quat_lattice_finalize(&test);
    quat_alg_elem_finalize(&init_helper);
    ibz_mat_4x4_finalize(&prod);
    quat_order_finalize(&order);
    quat_alg_finalize(&alg);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);
    quat_left_ideal_finalize(&lideal);
    return(res);
}

/***************************** Functions from quaternion_tools.c ***************************************/

//void quat_connecting_ideal(quat_left_ideal_t *connecting_ideal, const quat_order_t *O1, const quat_order_t *O2, const quat_alg_t *alg);
int quat_test_lideal_connecting_ideal(){
    int res = 0;
    ibz_t norm;
    quat_alg_t alg;
    quat_alg_elem_t gen;
    quat_order_t o1, o2;
    quat_lattice_t prod;
    quat_alg_elem_t test;
    quat_left_ideal_t lideal, cmp;
    quat_order_init(&o1);
    quat_order_init(&o2);
    quat_lattice_init(&prod);
    quat_alg_init_set_ui(&alg,19);
    quat_left_ideal_init(&lideal);
    quat_alg_elem_init(&test);

    ibz_set(&(o1.basis[0][0]),4);
    ibz_set(&(o1.basis[0][1]),0);
    ibz_set(&(o1.basis[0][2]),2);
    ibz_set(&(o1.basis[0][3]),2);
    ibz_set(&(o1.basis[1][0]),0);
    ibz_set(&(o1.basis[1][1]),8);
    ibz_set(&(o1.basis[1][2]),4);
    ibz_set(&(o1.basis[1][3]),3);
    ibz_set(&(o1.basis[2][0]),0);
    ibz_set(&(o1.basis[2][1]),0);
    ibz_set(&(o1.basis[2][2]),2);
    ibz_set(&(o1.basis[2][3]),0);
    ibz_set(&(o1.basis[3][0]),0);
    ibz_set(&(o1.basis[3][1]),0);
    ibz_set(&(o1.basis[3][2]),0);
    ibz_set(&(o1.basis[3][3]),1);
    ibz_set(&(o1.denom),4);

    ibz_set(&(o2.basis[0][0]),2);
    ibz_set(&(o2.basis[0][1]),0);
    ibz_set(&(o2.basis[0][2]),1);
    ibz_set(&(o2.basis[0][3]),0);
    ibz_set(&(o2.basis[1][0]),0);
    ibz_set(&(o2.basis[1][1]),2);
    ibz_set(&(o2.basis[1][2]),0);
    ibz_set(&(o2.basis[1][3]),1);
    ibz_set(&(o2.basis[2][0]),0);
    ibz_set(&(o2.basis[2][1]),0);
    ibz_set(&(o2.basis[2][2]),1);
    ibz_set(&(o2.basis[2][3]),0);
    ibz_set(&(o2.basis[3][0]),0);
    ibz_set(&(o2.basis[3][1]),0);
    ibz_set(&(o2.basis[3][2]),0);
    ibz_set(&(o2.basis[3][3]),1);
    ibz_set(&(o2.denom),2);

    quat_connecting_ideal(&lideal,&o1,&o2,&alg);
    quat_lattice_mul(&prod,&o1,&(lideal.lattice),&alg);
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(prod.denom),&(prod.basis[0][i]),&(prod.basis[1][i]),&(prod.basis[2][i]),&(prod.basis[3][i]));
        res = res || !quat_lattice_contains(NULL,&(lideal.lattice),&test,&alg);
    }
    quat_lattice_mul(&prod,&(lideal.lattice),&o2,&alg);
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(prod.denom),&(prod.basis[0][i]),&(prod.basis[1][i]),&(prod.basis[2][i]),&(prod.basis[3][i]));
        res = res || !quat_lattice_contains(NULL,&(lideal.lattice),&test,&alg);
    }

    quat_connecting_ideal(&lideal,&o2,&o2,&alg);
    quat_lattice_mul(&prod,&o2,&(lideal.lattice),&alg);
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(prod.denom),&(prod.basis[0][i]),&(prod.basis[1][i]),&(prod.basis[2][i]),&(prod.basis[3][i]));
        res = res || !quat_lattice_contains(NULL,&(lideal.lattice),&test,&alg);
    }
    quat_lattice_mul(&prod,&(lideal.lattice),&o2,&alg);
    for(int i = 0; i < 4; i++){
        quat_alg_elem_copy_ibz(&test,&(prod.denom),&(prod.basis[0][i]),&(prod.basis[1][i]),&(prod.basis[2][i]),&(prod.basis[3][i]));
        res = res || !quat_lattice_contains(NULL,&(lideal.lattice),&test,&alg);
    }

    if (res != 0){
        printf("Quaternion unit test lideal_connecting_ideal failed\n");
    }
    quat_alg_elem_finalize(&test);
    quat_order_finalize(&o1);
    quat_order_finalize(&o2);
    quat_lattice_finalize(&prod);
    quat_left_ideal_finalize(&lideal);
    quat_alg_finalize(&alg);
    return(res);
}

// run all previous tests
int quat_test_lideal(){
    int res = 0;
    printf("\nRunning quaternion tests of ideal and order functions\n");
    res = res | quat_test_lideal_create_principal();
    res = res | quat_test_lideal_create_from_primitive();
    res = res | quat_test_lideal_make_primitive_then_create();
    res = res | quat_test_lideal_random_2e();
    res = res | quat_test_lideal_generator();
    res = res | quat_test_lideal_generator_coprime();
    res = res | quat_test_lideal_mul();
    res = res | quat_test_lideal_add_intersect_equals();
    res = res | quat_test_lideal_isom();
    res = res | quat_test_lideal_right_order();
    res = res | quat_test_lideal_reduce_basis();
    res = res | quat_test_lideal_connecting_ideal();
    return(res);
}
