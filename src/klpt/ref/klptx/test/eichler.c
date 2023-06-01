#include <inttypes.h>

#include "klpt_tests.h"
#include "../eichler.c"

int test_special_cornacchia(int index) {
    int res;
    ibz_t q;
    ibz_t x,y,p;
    ibz_init(&q);
    ibz_init(&x);ibz_init(&y);ibz_init(&p);
    int qq = ALTERNATE_EXTREMAL_ORDERS[index].q;
    ibz_set(&q,qq);

    int exp = 128;
    res = 0;
    int adjust = 0;
    if (qq%8 ==7) {
        adjust =3;
    }
    else if (qq%8==3) {
        adjust =2;
    }
    while (!res) {
        generate_random_prime(&p,0,exp);
        // ibz_printf("%Zd \n",p);
        if (qq%4 ==3) {
            res = ibz_cornacchia_special_prime(&x,&y,&q,&p,adjust);
        }
        else {
            res= ibz_cornacchia_prime(&x,&y,&q,&p);
        }
        
    }


    ibz_finalize(&q);
    ibz_finalize(&x);ibz_finalize(&y);ibz_finalize(&p);
    return res;
}

int klpt_test_eichler_special_norm() {

    int status;
    int res;
    ibz_t p,n,n_beta;
    quat_alg_t Bpoo;
    quat_left_ideal_t ideal;
    quat_order_t right_order;
    quat_alg_elem_t gen,beta,gen_constraint;
    quat_alg_coord_t coeffs;
    ibq_t ibq_norm;
    ibz_t norm;

    ibz_init(&norm);ibq_init(&ibq_norm);
    ibz_init(&p);ibz_init(&n);ibz_init(&n_beta);
    ibz_copy(&p,&CHARACTERISTIC);
    quat_alg_init_set(&Bpoo,&p);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&beta);
    quat_alg_elem_init(&gen_constraint);
    quat_left_ideal_init(&ideal);
    quat_alg_coord_init(&coeffs);
    quat_order_init(&right_order);

    int exp = ibz_bitsize(&CHARACTERISTIC)/4;
    quat_lideal_random_2e(&ideal,&STANDARD_EXTREMAL_ORDER.order,&Bpoo,exp,'6');
    quat_lideal_right_order(&right_order,&ideal,&Bpoo);
    int found = 0;
    while (!found) {
        ibz_rand_interval_minm_m(&coeffs[0],20);
        ibz_rand_interval_minm_m(&coeffs[1],20);
        ibz_rand_interval_minm_m(&coeffs[2],20);
        ibz_rand_interval_minm_m(&coeffs[3],20);
        ibz_mat_4x4_eval(&gen_constraint.coord,&right_order.basis,&coeffs);
        ibz_copy(&gen_constraint.denom,&right_order.denom);

        quat_alg_norm(&ibq_norm,&gen_constraint,&Bpoo);
        status = ibq_to_ibz(&norm,&ibq_norm);
        assert(status);
        found = ibz_get(&norm)%4==2 && quat_alg_is_primitive(&gen_constraint,&right_order,&Bpoo);
    }

    assert(quat_lattice_contains(&coeffs,&right_order,&gen_constraint,&Bpoo));

    res = klpt_eichler_special_norm(&beta,&n_beta,&gen,&n,&ideal,&gen_constraint,&Bpoo);
    assert(res);

    ibz_finalize(&norm);ibq_finalize(&ibq_norm);
    ibz_finalize(&p);ibz_finalize(&n);ibz_finalize(&n_beta);
    quat_alg_finalize(&Bpoo);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&beta);
    quat_alg_elem_finalize(&gen_constraint);
    quat_left_ideal_finalize(&ideal);
    quat_alg_coord_finalize(&coeffs);
    quat_order_finalize(&right_order);
    
    return res;
}

int klpt_test_eichler_special_norm_fixed_standard() {

    int res;
    int rep_found;
    int status;
    ibz_t p,n,pow2,n_beta,temp;
    quat_alg_elem_t gen,gen_constraint,beta;
    quat_alg_t Bpoo;
    quat_alg_coord_t coeffs;
    quat_left_ideal_t ideal;
    quat_order_t right_order;
    ibq_t ibq_norm;ibz_t norm;

    ibq_init(&ibq_norm);ibz_init(&norm);
    quat_alg_coord_init(&coeffs);
    ibz_init(&p);ibz_init(&n);ibz_init(&pow2);ibz_init(&n_beta);ibz_init(&temp);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_constraint);
    quat_alg_elem_init(&beta);
    quat_left_ideal_init(&ideal);
    quat_order_init(&right_order);


    int exp = 128;
    int margin = 20;

    res = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);
    res = generate_random_prime(&n,1,exp);

    ibz_pow(&pow2,&ibz_const_two,2*exp+margin);

    ibz_mul(&temp,&n,&pow2);
    rep_found = represent_integer(&gen,&temp,&Bpoo);
    while (!rep_found) {
        rep_found = represent_integer(&gen,&temp,&Bpoo);
    }
    
    quat_lideal_make_primitive_then_create(&ideal,&gen,&n,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    quat_lideal_right_order(&right_order,&ideal,&Bpoo);

    int found = 0;
    while (!found) {
        ibz_rand_interval_minm_m(&coeffs[0],20);
        ibz_rand_interval_minm_m(&coeffs[1],20);
        ibz_rand_interval_minm_m(&coeffs[2],20);
        ibz_rand_interval_minm_m(&coeffs[3],20);
        ibz_mat_4x4_eval(&gen_constraint.coord,&right_order.basis,&coeffs);
        ibz_copy(&gen_constraint.denom,&right_order.denom);

        quat_alg_norm(&ibq_norm,&gen_constraint,&Bpoo);
        status = ibq_to_ibz(&norm,&ibq_norm);
        assert(status);

        found = ibz_get(&norm)%4==2; 
    }
    assert(quat_lattice_contains(&coeffs,&right_order,&gen_constraint,&Bpoo));
    assert(quat_lattice_contains(&coeffs,&right_order,&gen,&Bpoo));


    res = res && eichler_special_norm_fixed(&beta,&n_beta,&STANDARD_EXTREMAL_ORDER,&gen,&n,&gen_constraint,&Bpoo);
    
    quat_alg_norm(&ibq_norm,&beta,&Bpoo);
    assert(ibq_to_ibz(&norm,&ibq_norm));
    assert(ibz_cmp(&n_beta,&norm)==0);

    res = res && quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&beta,&Bpoo);
    res = res && quat_lattice_contains(&coeffs,&right_order,&beta,&Bpoo);

    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);
    quat_alg_coord_finalize(&coeffs);
    ibz_finalize(&p);ibz_finalize(&n);ibz_finalize(&pow2);ibz_finalize(&n_beta);ibz_finalize(&temp);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_constraint);
    quat_alg_elem_finalize(&beta);
    quat_left_ideal_finalize(&ideal);
    quat_order_finalize(&right_order);
    quat_alg_finalize(&Bpoo);

    return res;
}

int klpt_test_eichler_special_norm_fixed_alternate(int index) {

    int status;
    int res;
    int rep_found;
    ibz_t p,n,pow2,n_beta,temp;
    quat_alg_elem_t gen,gen_constraint,beta;
    ibz_mat_4x4_t reduced, gram;
    quat_alg_t Bpoo;
    quat_alg_coord_t coeffs;
    quat_left_ideal_t start,equiv;

    ibq_t ibq_norm;
    ibz_t norm;
    ibq_init(&ibq_norm);
    ibz_init(&norm);
    quat_alg_coord_init(&coeffs);
    ibz_init(&p);ibz_init(&n);ibz_init(&pow2);ibz_init(&n_beta);ibz_init(&temp);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_constraint);
    quat_alg_elem_init(&beta);
    quat_left_ideal_init(&start);quat_left_ideal_init(&equiv);
    ibz_mat_4x4_init(&reduced);
    ibz_mat_4x4_init(&gram);
    int exp = 128;
    int margin = 20;

    ibz_copy(&p,&CHARACTERISTIC);
    quat_alg_init_set(&Bpoo,&p);

    const quat_p_extremal_maximal_order_t *order = &ALTERNATE_EXTREMAL_ORDERS[index];

    quat_alg_norm(&ibq_norm,&order->i,&Bpoo);

    assert(ibq_to_ibz(&norm,&ibq_norm));
    assert(ibz_get(&norm)==ALTERNATE_EXTREMAL_ORDERS[index].q);
    assert(quat_lattice_contains(&coeffs,&order->order,&order->i,&Bpoo));
    assert(quat_lattice_contains(&coeffs,&order->order,&order->j,&Bpoo));


    quat_lideal_random_2e(&start,&order->order,&Bpoo,2*exp,'8');
    int lideal_generator_ok = quat_lideal_generator(&gen,&start,&Bpoo,0);
    assert(lideal_generator_ok);

    ibz_pow(&pow2,&ibz_const_two,2*exp);
    quat_lideal_create_from_primitive(&start,&gen,&pow2,&order->order,&Bpoo);

    quat_lideal_reduce_basis(&reduced,&gram,&start,&Bpoo);


    res = klpt_lideal_equiv(&gen,&n,&reduced,&gram,&start.norm,&start.lattice.denom,&Bpoo);
    assert(res);


    quat_lideal_create_from_primitive(&equiv,&gen,&n,&order->order,&Bpoo);
    assert(quat_lideal_isom(&gen_constraint,&equiv,&start,&Bpoo));

    if (order->q %4 == 3 ) {
        ibz_copy(&coeffs[0],&ibz_const_one);
        ibz_copy(&coeffs[1],&ibz_const_one);
        ibz_copy(&coeffs[2],&ibz_const_zero);
        ibz_copy(&coeffs[3],&ibz_const_zero);
        order_elem_create(&beta,order,&coeffs,&Bpoo);
        ibz_mul(&beta.denom,&beta.denom,&ibz_const_two);
        assert(quat_lattice_contains(&coeffs,&order->order,&beta,&Bpoo));
    }
    

    int found = 0; 
    quat_order_t right_order;
    quat_order_init(&right_order);

    quat_lideal_right_order(&right_order,&equiv,&Bpoo);
    while (!found) {
        ibz_rand_interval_minm_m(&coeffs[0],20);
        ibz_rand_interval_minm_m(&coeffs[1],20);
        ibz_rand_interval_minm_m(&coeffs[2],20);
        ibz_rand_interval_minm_m(&coeffs[3],20);
        ibz_mat_4x4_eval(&gen_constraint.coord,&right_order.basis,&coeffs);
        ibz_copy(&gen_constraint.denom,&right_order.denom);

        quat_alg_norm(&ibq_norm,&gen_constraint,&Bpoo);
        status = ibq_to_ibz(&norm,&ibq_norm);
        assert(status);
        found = ibz_get(&norm)%4==2; 
    }
    assert(quat_lattice_contains(&coeffs,&right_order,&gen_constraint,&Bpoo));
    assert(quat_lattice_contains(&coeffs,&right_order,&gen,&Bpoo));

    res = res && eichler_special_norm_fixed(&beta,&n_beta,order,&gen,&n,&gen_constraint,&Bpoo);
    if (!res ) {
        printf("FAIL eichler norm for q = %" PRId64 " (this is expected) \n",order->q);

    } 
    else {
        printf("SUCCESS eichler norm for q = %" PRId64 " \n",order->q);
        quat_alg_norm(&ibq_norm,&beta,&Bpoo);
        assert(ibq_to_ibz(&norm,&ibq_norm));
        assert(ibz_cmp(&n_beta,&norm)==0);

        res = res && quat_lattice_contains(&coeffs,&order->order,&beta,&Bpoo);
        res = res && quat_lattice_contains(&coeffs,&right_order,&beta,&Bpoo);
    }

    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);
    quat_alg_coord_finalize(&coeffs);
    ibz_finalize(&p);ibz_finalize(&n);ibz_finalize(&pow2);ibz_finalize(&n_beta);ibz_finalize(&temp);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_constraint);
    quat_alg_elem_finalize(&beta);
    quat_left_ideal_finalize(&start); quat_left_ideal_finalize(&equiv);
    quat_order_finalize(&right_order);
    quat_alg_finalize(&Bpoo);
    ibz_mat_4x4_finalize(&reduced);
    ibz_mat_4x4_finalize(&gram);
    return res;
}

int klpt_test_eichler() {
    int res = 1;
    printf("Running klpt tests for eichler norm \n \n");

    for (int i =0; i<3;i++) {
        res = res & klpt_test_eichler_special_norm_fixed_standard();
    }
     
    if (!res) {
        printf("KLPT unit test special eichler norm fixed for standard order failed\n \n");
    }
    for (int i =0; i<NUM_ALTERNATE_EXTREMAL_ORDERS;i++) {
        res = res & test_special_cornacchia(i);
        // if (ALTERNATE_EXTREMAL_ORDERS[i].q%4==1) {
        if (1) {
            res = res & klpt_test_eichler_special_norm_fixed_alternate(i);
        }   
    }
    if (!res) {
        printf("KLPT unit test special eichler norm fixed for alternate order failed (this is expected)\n \n");
    }
    res = 1;
    for (int i =0; i<3;i++) {
        res = res & klpt_test_eichler_special_norm();
    }
     
    if (!res) {
        printf("KLPT unit test special eichler norm full failed");
    }

    return res;
}   
