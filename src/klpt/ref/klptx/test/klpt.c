#include "klpt_tests.h"

int klpt_test_keygen_klpt() {

    int status;
    int res;
    quat_alg_elem_t gen,gamma_tau;
    ibz_t n,p,Ntau,temp;
    
    quat_alg_t Bpoo;
    quat_left_ideal_t lideal,lideal_start;
    ibq_t ibq_norm;
    ibz_t norm;
    ibq_init(&ibq_norm);
    ibz_init(&norm);
    quat_alg_elem_init(&gen);
    ibz_init(&n);ibz_init(&p);ibz_init(&Ntau);ibz_init(&temp);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal_start);
    quat_alg_elem_init(&gamma_tau);

    int exp = 128;
    res = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);

    quat_alg_coord_t coeffs;
    quat_alg_coord_init(&coeffs);

    res = res && generate_random_prime(&Ntau,1,exp/2);

    int rep_found;

    res = res && generate_random_prime(&temp,1,3*exp);
    ibz_mul(&temp,&temp,&Ntau);

    rep_found = represent_integer(&gamma_tau,&temp,&Bpoo);
    while (!rep_found) {
        rep_found = represent_integer(&gamma_tau,&temp,&Bpoo);
    }

    quat_lideal_make_primitive_then_create(&lideal_start,&gamma_tau,&Ntau,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

    res =res && klpt_keygen_klpt(&gen,&lideal_start,&Bpoo);
    if (res) {
        quat_alg_make_primitive(&coeffs,&temp,&gen,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    ibz_mul(&gen.denom,&gen.denom,&temp);
    quat_alg_normalize(&gen);



    assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen,&Bpoo));
    assert(quat_lattice_contains(&coeffs,&lideal_start.lattice,&gen,&Bpoo));

    quat_lideal_create_from_primitive(&lideal,&gen,&Ntau,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    assert(quat_lideal_equals(&lideal,&lideal_start,&Bpoo));
    
    quat_alg_conj(&gen,&gen);

    
    quat_alg_norm(&ibq_norm,&gen,&Bpoo);
    status = ibq_to_ibz(&norm,&ibq_norm);
    assert(status);

    ibz_pow(&temp,&ibz_const_two,KLPT_keygen_length);
    ibz_mul(&temp,&temp,&Ntau);
    ibz_mod(&temp,&temp,&norm);
    assert(ibz_cmp(&temp,&ibz_const_zero)==0);

    ibz_div(&n,&temp,&norm,&Ntau);

    quat_lideal_create_from_primitive(&lideal,&gen,&n,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    status = quat_lideal_isom(&gamma_tau,&lideal,&lideal_start,&Bpoo);
    assert(status);
    }


    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&n);ibz_finalize(&p);ibz_finalize(&Ntau);ibz_finalize(&temp);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal_start);
    quat_alg_elem_finalize(&gamma_tau);
    quat_alg_finalize(&Bpoo);
    quat_alg_coord_finalize(&coeffs);



    return res;


}

int klpt_test_signing_klpt() {
    
    int status;
    int res;
    quat_alg_elem_t gen,gamma_tau,delta;
    ibz_t n,p,Ntau,temp;
    ibz_mat_4x4_t reduced,gram;
    quat_order_t right_order;
    
    quat_alg_t Bpoo;
    quat_left_ideal_t lideal,lideal_start,lideal_sign,lideal_sign_equiv,lideal_equiv,lideal_long;
    ibq_t ibq_norm;
    ibz_t norm;
    ibq_init(&ibq_norm);
    ibz_init(&norm);

    quat_alg_elem_init(&gen);
    ibz_init(&n);ibz_init(&p);ibz_init(&Ntau);ibz_init(&temp);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&reduced);
    quat_order_init(&right_order);
    quat_left_ideal_init(&lideal);
    quat_left_ideal_init(&lideal_sign_equiv);
    quat_left_ideal_init(&lideal_sign);
    quat_left_ideal_init(&lideal_start);
    quat_left_ideal_init(&lideal_equiv);
    quat_left_ideal_init(&lideal_long);
    quat_alg_elem_init(&delta);quat_alg_elem_init(&gamma_tau);

    int exp = 128;
    res = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);

    quat_alg_coord_t coeffs;
    quat_alg_coord_init(&coeffs);

    generate_random_prime(&Ntau,1,exp/2);

    int rep_found;

    int cnt = 0;
    res = 0;
    while (!res & (cnt < 10)) {
        
        cnt ++;
    
        res = generate_random_prime(&temp,1,3*exp);
        ibz_mul(&temp,&temp,&Ntau);


        rep_found = represent_integer(&gamma_tau,&temp,&Bpoo);
        while (!rep_found) {
            rep_found = represent_integer(&gamma_tau,&temp,&Bpoo);
        }

        quat_lideal_make_primitive_then_create(&lideal_start,&gamma_tau,&Ntau,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

    
        quat_lideal_right_order(&right_order,&lideal_start,&Bpoo);

        quat_lideal_random_2e(&lideal_sign_equiv,&right_order,&Bpoo,2*exp,'8');

        ibz_pow(&temp,&ibz_const_two,2*exp);
        assert(ibz_cmp(&temp,&lideal_sign_equiv.norm)==0);

        int lideal_generator_ok = quat_lideal_generator(&gen,&lideal_sign_equiv,&Bpoo,0);
        assert(lideal_generator_ok);

        quat_alg_mul(&gen,&gamma_tau,&gen,&Bpoo);

        assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen,&Bpoo));

        quat_lideal_make_primitive_then_create(&lideal,&gen,&lideal_sign_equiv.norm,&STANDARD_EXTREMAL_ORDER.order,&Bpoo); 
        ibz_mul(&temp,&lideal.norm,&Ntau);
        quat_lideal_make_primitive_then_create(&lideal_long,&gen,&temp,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

        lideal_generator_ok = quat_lideal_generator(&gen,&lideal_start,&Bpoo,0);
        assert(lideal_generator_ok);
        assert(quat_lattice_contains(&coeffs,&lideal_start.lattice,&gen,&Bpoo));

        quat_lideal_reduce_basis(&reduced,&gram,&lideal,&Bpoo);

        res = res && klpt_lideal_equiv(&gen,&n,&reduced,&gram,&lideal.norm,&lideal.lattice.denom,&Bpoo);
        if (!res) {
            continue;
        }

        quat_lideal_create_from_primitive(&lideal_equiv,&gen,&n,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

        assert(ibz_cmp(&lideal_equiv.norm,&n)==0);

        quat_alg_conj(&delta,&gen);    
        assert(res);
        res = res && klpt_signing_klpt(&gen,&lideal_equiv,&lideal_start,&delta,&Bpoo);
    }

    if (res) {

        assert(quat_alg_is_primitive(&gen,&STANDARD_EXTREMAL_ORDER.order,&Bpoo));
        quat_alg_norm(&ibq_norm,&gen,&Bpoo);
        assert(ibq_to_ibz(&norm,&ibq_norm));

        ibz_pow(&temp,&ibz_const_two,KLPT_signing_klpt_length);
        ibz_mul(&temp,&temp,&n);
    
        assert(ibz_cmp(&temp,&norm)==0);

        assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen,&Bpoo));
        assert(quat_lattice_contains(&coeffs,&lideal_equiv.lattice,&gen,&Bpoo));

        quat_alg_mul(&gen,&gen,&delta,&Bpoo);
        ibz_mul(&gen.denom,&gen.denom,&lideal_equiv.norm);
        quat_alg_normalize(&gen);
        assert(quat_lattice_contains(&coeffs,&right_order,&gen,&Bpoo));
        status = quat_alg_is_primitive(&gen,&right_order,&Bpoo);
        assert(status);
        assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen,&Bpoo));
        assert(quat_lattice_contains(&coeffs,&lideal.lattice,&gen,&Bpoo));
        assert(quat_lattice_contains(&coeffs,&right_order,&gen,&Bpoo));

        // quat_alg_conj(&delta,&gen);

        assert(quat_lattice_contains(&coeffs,&lideal_sign_equiv.lattice,&gen,&Bpoo)); 

        quat_alg_conj(&gen,&gen);
        ibz_pow(&temp,&ibz_const_two,KLPT_signing_klpt_length);
        quat_lideal_create_from_primitive(&lideal_sign,&gen,&temp,&right_order,&Bpoo);
        assert(0==ibz_cmp(&temp,&lideal_sign.norm));

        status = quat_lideal_isom(&delta,&lideal_sign,&lideal_sign_equiv,&Bpoo);
        assert(status);

        assert(0==ibz_cmp(&lideal_sign.norm,&temp));
    }
    
    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&n);ibz_finalize(&p);ibz_finalize(&Ntau);ibz_finalize(&temp);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&reduced);
    quat_order_finalize(&right_order);
    quat_left_ideal_finalize(&lideal);
    quat_left_ideal_finalize(&lideal_sign_equiv);
    quat_left_ideal_finalize(&lideal_sign);
    quat_left_ideal_finalize(&lideal_start);
    quat_left_ideal_finalize(&lideal_equiv);
    quat_left_ideal_finalize(&lideal_long);
    quat_alg_elem_finalize(&delta);quat_alg_elem_finalize(&gamma_tau);
    quat_alg_coord_finalize(&coeffs);
    quat_alg_finalize(&Bpoo);

    return res;


}

int klpt_test_klpt() {
    int res = 1;
    printf("\n \nRunning klpt tests for keygen/signing klpt \n \n");

     for (int i =0; i<3;i++) {
        res = res & klpt_test_keygen_klpt();
    }
     
    if (!res) {
        printf("KLPT unit test keygen_klpt failed\n");
    }

    for (int i =0; i<10;i++) {
        res = res & klpt_test_signing_klpt();
    }
     
    if (!res) {
        printf("\nKLPT unit test signing_klpt failed\n");
    }

    return res;
}
