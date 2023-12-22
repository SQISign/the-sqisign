#include "klpt_tests.h"


int kltp_test_keygen_random_ideal() {

    int res;

    ibz_t n;
    quat_left_ideal_t ideal;
    ibz_init(&n);
    quat_left_ideal_init(&ideal);

    res = klpt_keygen_random_ideal(&ideal,&STANDARD_EXTREMAL_ORDER,&QUATALG_PINFTY);

    ibz_finalize(&n);
    quat_left_ideal_finalize(&ideal);

    return res;
}

int klpt_test_lideal_generator_coprime() {
    int found = 1;
    ibz_t p,temp;
    ibz_t M;
    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    quat_left_ideal_t id,id2;
    
    quat_alg_elem_init(&gamma);
    ibz_init(&M);
    ibz_init(&p);
    ibz_init(&temp); 
    quat_left_ideal_init(&id);  
    quat_left_ideal_init(&id2);   
 


    int exp = 128;
    found = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);
    
    found = found && (generate_random_prime(&M,1,exp));
    ibz_pow(&temp,&ibz_const_two,2*exp);
    ibz_mul(&temp,&M,&temp);

    found = found && represent_integer(&gamma,&temp,&Bpoo);
    assert(found);

    quat_lideal_make_primitive_then_create(&id,&gamma,&M,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

    found = found && quat_lideal_generator_coprime(&gamma,&id,&ibz_const_two,&Bpoo,0);
    assert(found);

    assert(quat_alg_is_primitive(&gamma,&STANDARD_EXTREMAL_ORDER.order,&Bpoo));

    quat_lideal_make_primitive_then_create(&id2,&gamma,&M,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

    assert(quat_lideal_equals(&id,&id2,&Bpoo));

    quat_alg_elem_finalize(&gamma);
    ibz_finalize(&M);
    ibz_finalize(&p);
    ibz_finalize(&temp); 
    quat_left_ideal_finalize(&id);  
    quat_left_ideal_finalize(&id2);
    quat_alg_finalize(&Bpoo);
    return found;

} 

int klpt_test_lideal_isom (int index) {
    
    int res =1 ;
    ibz_t p,n;
    quat_alg_elem_t gen;
    quat_alg_t Bpoo;
    quat_left_ideal_t ideal,equiv;
    quat_order_t right_order;
    ibz_mat_4x4_t gram,reduced;
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&reduced);
    ibz_init(&p);ibz_init(&n);
    ibz_copy(&p,&CHARACTERISTIC);
    quat_alg_init_set(&Bpoo,&p);
    quat_left_ideal_init(&ideal);quat_left_ideal_init(&equiv);
    const quat_p_extremal_maximal_order_t *alternate_order = &ALTERNATE_EXTREMAL_ORDERS[index];
    quat_alg_elem_init(&gen);


    for (int i=0; i<2;i++) {
        quat_connecting_ideal(&ideal,&STANDARD_EXTREMAL_ORDER.order,&alternate_order->order,&Bpoo);
        quat_lideal_reduce_basis(&reduced,&gram,&ideal,&Bpoo);
        klpt_lideal_equiv(&gen,&n,&reduced,&gram,&ideal.norm,&ideal.lattice.denom,&Bpoo);

        quat_lideal_make_primitive_then_create(&equiv,&gen,&n,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

        assert(quat_lideal_isom(&gen,&ideal,&equiv,&Bpoo));

        quat_lideal_random_2e(&ideal,&alternate_order->order,&Bpoo,128,'8');

        quat_lideal_reduce_basis(&reduced,&gram,&ideal,&Bpoo);
        klpt_lideal_equiv(&gen,&n,&reduced,&gram,&ideal.norm,&ideal.lattice.denom,&Bpoo);

        quat_lideal_make_primitive_then_create(&equiv,&gen,&n,&alternate_order->order,&Bpoo);

        assert(quat_lideal_isom(&gen,&ideal,&equiv,&Bpoo));
    }

    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&reduced);
    ibz_finalize(&p);ibz_finalize(&n);
    quat_alg_finalize(&Bpoo);
    quat_left_ideal_finalize(&ideal);quat_left_ideal_finalize(&equiv);
    quat_alg_elem_finalize(&gen);

    return res;

}

int klpt_test_connecting_ideal(int index) {
    int res;

    ibz_t p;
    quat_alg_t Bpoo;
    quat_left_ideal_t ideal,connect;
    quat_order_t right_order;
    ibz_init(&p);
    ibz_copy(&p,&CHARACTERISTIC);
    quat_alg_init_set(&Bpoo,&p);
    quat_left_ideal_init(&ideal);quat_left_ideal_init(&connect);
    quat_order_init(&right_order);
    const quat_p_extremal_maximal_order_t *alternate_order = &ALTERNATE_EXTREMAL_ORDERS[index];

    quat_lideal_random_2e(&ideal,&STANDARD_EXTREMAL_ORDER.order,&Bpoo,256,'8');

    quat_lideal_right_order(&right_order,&ideal,&Bpoo);

    quat_connecting_ideal(&connect,&STANDARD_EXTREMAL_ORDER.order,&alternate_order->order,&Bpoo);
    quat_connecting_ideal(&connect,&right_order,&alternate_order->order,&Bpoo);

    ibz_finalize(&p);
    quat_alg_finalize(&Bpoo);
    quat_left_ideal_finalize(&ideal);quat_left_ideal_finalize(&connect);
    quat_order_finalize(&right_order);

    res = 1;
    return res;
}

int klpt_test_find_linear_comb() {
    int status;
    int found = 0;
    quat_alg_t Bpoo;
    ibz_t p,M,temp,temp2;
    quat_order_t right_order;
    quat_alg_elem_t beta,gen_start,gen_end,quat_temp;
    quat_left_ideal_t ideal_two_start, ideal_two;
    quat_left_ideal_t random_ideal;
    ibz_vec_2_t C;
    ibz_t ibz_two_exp;
    ibq_t ibq_norm;
    ibz_t norm;
    quat_alg_coord_t coeffs;

    quat_alg_coord_init(&coeffs);
    quat_order_init(&right_order);
    quat_left_ideal_init(&random_ideal);
    ibq_init(&ibq_norm);
    ibz_init(&norm);
    ibz_init(&p);
    ibz_init(&M);
    ibz_init(&temp);ibz_init(&temp2);
    ibz_init(&ibz_two_exp);
    ibz_vec_2_init(&C);   
    quat_alg_elem_init(&beta);quat_alg_elem_init(&gen_start);
    quat_alg_elem_init(&gen_end);quat_alg_elem_init(&quat_temp);
    quat_left_ideal_init(&ideal_two_start);
    quat_left_ideal_init(&ideal_two);
    
    int found_rep;
    int exp = 128;
    int margin = 30;

    found = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);

    found = found && (generate_random_prime(&M,1,exp+margin));
    found = found && (generate_random_prime(&temp,1,exp+margin));
    assert(found);
    ibz_mul(&temp,&M,&temp); 
    ibz_copy(&temp2,&temp);


    found_rep = represent_integer(&beta,&temp2,&Bpoo); 
    while (!found_rep) {
        found_rep = represent_integer(&beta,&temp2,&Bpoo); 

    }
    ibz_pow(&ibz_two_exp,&ibz_const_two,exp);
    ibz_mul(&temp,&M,&ibz_two_exp);
    ibz_copy(&temp2,&temp);
    found_rep = represent_integer(&gen_start,&temp,&Bpoo); 
    while (!found_rep && quat_alg_is_primitive(&gen_start,&STANDARD_EXTREMAL_ORDER.order,&Bpoo)) {
        found_rep = represent_integer(&gen_start,&temp2,&Bpoo);
    }
    quat_lideal_make_primitive_then_create(&ideal_two_start,&gen_start,&ibz_const_two,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

    quat_alg_conj(&gen_end,&beta);
    quat_alg_mul(&gen_end,&gen_start,&gen_end,&Bpoo);


    quat_lideal_make_primitive_then_create(&ideal_two,&gen_end,&ibz_const_two,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    quat_alg_norm(&ibq_norm,&beta,&Bpoo);

    status = ibq_to_ibz(&norm,&ibq_norm);
    assert(status);

    while (!quat_lideal_equals(&ideal_two,&ideal_two_start,&Bpoo)) {

        found_rep = represent_integer(&beta,&norm,&Bpoo); 
        while (!found_rep) {
            found_rep = represent_integer(&beta,&norm,&Bpoo); 
        } 

        ibz_copy(&temp2,&temp);
        found_rep = represent_integer(&gen_start,&temp2,&Bpoo);
        while (!found_rep && quat_alg_is_primitive(&gen_start,&STANDARD_EXTREMAL_ORDER.order,&Bpoo)) {
            found_rep = represent_integer(&gen_start,&temp2,&Bpoo);
        } 
        // quat_alg_elem_print(&gen_start);
        quat_lideal_make_primitive_then_create(&ideal_two_start,&gen_start,&ibz_const_two,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

        quat_alg_conj(&gen_end,&beta);
        quat_alg_mul(&gen_end,&gen_start,&gen_end,&Bpoo);


        quat_lideal_make_primitive_then_create(&ideal_two,&gen_end,&ibz_const_two,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    }

    quat_lideal_create_from_primitive(&random_ideal,&beta,&M,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    quat_lideal_right_order(&right_order,&random_ideal,&Bpoo);

    assert(quat_lattice_contains(&coeffs,&right_order,&beta,&Bpoo));

    quat_lideal_random_2e(&ideal_two_start,&right_order,&Bpoo,exp,'8');
    
    quat_lideal_random_2e(&ideal_two,&right_order,&Bpoo,exp,'8');
    assert(!quat_lideal_equals(&ideal_two_start,&ideal_two,&Bpoo));

    found = found && quat_lideal_generator(&gen_start,&ideal_two_start,&Bpoo,0);
    assert(found);
    found = found && quat_lideal_generator(&gen_end,&ideal_two,&Bpoo,0);
    assert(found);
    assert(!quat_lattice_contains(&coeffs,&ideal_two.lattice,&gen_start,&Bpoo));

    quat_alg_conj(&quat_temp,&beta);
    quat_alg_mul(&quat_temp,&gen_start,&quat_temp,&Bpoo);
    assert(!quat_lattice_contains(&coeffs,&ideal_two_start.lattice,&quat_temp,&Bpoo));

    found = found && klpt_find_linear_comb(&C,&beta,&right_order,&ibz_two_exp,(unsigned short)exp,&gen_start,&gen_end,&Bpoo);


    assert(found);
    assert(ibz_get(&C[0])%2 !=0 || ibz_get(&C[1])%2 != 0);

    quat_left_ideal_t ideal_test;
    quat_left_ideal_init(&ideal_test);
    quat_alg_scalar(&quat_temp,&C[1],&ibz_const_one);
    quat_alg_conj(&gen_end,&beta);
    quat_alg_mul(&gen_end,&quat_temp,&gen_end,&Bpoo);
    quat_alg_scalar(&quat_temp,&C[0],&ibz_const_one);
    quat_alg_add(&gen_end,&gen_end,&quat_temp);
    quat_alg_mul(&gen_end,&gen_start,&gen_end,&Bpoo);
    assert(quat_lattice_contains(&coeffs,&ideal_two.lattice,&gen_end,&Bpoo));

    quat_alg_finalize(&Bpoo);
    quat_left_ideal_finalize(&ideal_test);
    quat_alg_coord_finalize(&coeffs);
    quat_order_finalize(&right_order);
    quat_left_ideal_finalize(&random_ideal);
    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);
    ibz_finalize(&p);
    ibz_finalize(&M);
    ibz_finalize(&temp);ibz_finalize(&temp2);
    ibz_finalize(&ibz_two_exp);
    ibz_vec_2_finalize(&C);   
    quat_alg_elem_finalize(&beta);quat_alg_elem_finalize(&gen_start);
    quat_alg_elem_finalize(&gen_end);quat_alg_elem_finalize(&quat_temp);
    quat_left_ideal_finalize(&ideal_two_start);
    quat_left_ideal_finalize(&ideal_two);
    
    return found;
}

int klpt_test_solve_combi_eichler() {
    int found = 0;
    quat_alg_t Bpoo;
    ibz_t p,M,temp;
    quat_alg_elem_t gamma1,gamma2,delta,quat_res;
    ibz_vec_2_t C;
    quat_left_ideal_t lideal;
    quat_alg_coord_t t1,t2;
    quat_order_t right_order;
    quat_alg_elem_t shift_gamma,quat_one;
    ibq_t ibq_norm;ibz_t norm;

    ibq_init(&ibq_norm);
    ibz_init(&norm);
    ibz_init(&p);
    ibz_init(&M);
    ibz_init(&temp);   
    quat_alg_elem_init(&gamma1); quat_alg_elem_init(&gamma2);quat_alg_elem_init(&delta);quat_alg_elem_init(&quat_res);
    ibz_vec_2_init(&C);
    quat_left_ideal_init(&lideal);
    quat_alg_coord_init(&t1);quat_alg_coord_init(&t2);
    quat_order_init(&right_order);
    quat_alg_elem_init(&shift_gamma);quat_alg_elem_init(&quat_one);

    quat_alg_scalar(&quat_one,&ibz_const_one,&ibz_const_one);

    int exp = 128;

    int found_rep = 1;

    found = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);

    found_rep = (generate_random_prime(&M,1,exp));
    found_rep =  (generate_random_prime(&temp,1,exp+100));
    ibz_mul(&temp,&M,&temp); 
    found = found && represent_integer(&gamma1,&temp,&Bpoo); 
    ibz_pow(&temp,&ibz_const_two,exp+100);
    ibz_mul(&temp,&M,&temp); 
    found_rep = represent_integer(&gamma2,&temp,&Bpoo); 
    while (!found_rep) {
        found_rep = represent_integer(&gamma2,&temp,&Bpoo); 
    }

    ibz_pow(&temp,&ibz_const_two,3*exp);

    found_rep = represent_integer(&shift_gamma,&temp,&Bpoo);
    while (!found_rep) {
        represent_integer(&shift_gamma,&temp,&Bpoo);
    }
    found_rep = represent_integer(&delta,&temp,&Bpoo);
    while (!found_rep) {
        found_rep = represent_integer(&delta,&temp,&Bpoo);
    }

    quat_lideal_make_primitive_then_create(&lideal,&gamma1,&M,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);

    found = found && (ibz_cmp(&M,&lideal.norm)==0);

    quat_lideal_right_order(&right_order,&lideal,&Bpoo);

    found = found && solve_combi_eichler(&C,&STANDARD_EXTREMAL_ORDER,&shift_gamma,&delta,&lideal,&Bpoo,0);
    if (!found) {
        printf("eichler solve failed \n");
    }
    ibz_copy(&t1[0],&ibz_const_zero);ibz_copy(&t1[1],&ibz_const_zero);
    ibz_copy(&t1[2],&C[0]);ibz_copy(&t1[3],&C[1]);

    order_elem_create(&quat_res,&STANDARD_EXTREMAL_ORDER,&t1,&Bpoo);
    quat_alg_mul(&quat_res,&shift_gamma,&quat_res,&Bpoo);
    quat_alg_mul(&quat_res,&quat_res,&delta,&Bpoo);

    found = found && quat_lattice_contains(&t2,&right_order,&quat_res,&Bpoo);
    if (!found) {
        printf("the solution is not in the eichler order \n");
    }
    found = found && !quat_lattice_contains(&t2,&lideal.lattice,&quat_res,&Bpoo);
    if (!found) {
        printf("the solution is in the lattice \n");
    }


    found = found && solve_combi_eichler(&C,&STANDARD_EXTREMAL_ORDER,&gamma2,&quat_one,&lideal,&Bpoo,1);
    if (!found ){
        printf("the second solve combi failed \n");
    }

    ibz_copy(&t1[2],&C[0]);ibz_copy(&t1[3],&C[1]);

    order_elem_create(&quat_res,&STANDARD_EXTREMAL_ORDER,&t1,&Bpoo);
    quat_alg_mul(&quat_res,&gamma2,&quat_res,&Bpoo);
    quat_alg_mul(&quat_res,&quat_res,&quat_one,&Bpoo);


    found = found && quat_lattice_contains(&t2,&lideal.lattice,&quat_res,&Bpoo);

    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);
    ibz_finalize(&p);
    ibz_finalize(&M);
    ibz_finalize(&temp);   
    quat_alg_elem_finalize(&gamma1); quat_alg_elem_finalize(&gamma2);quat_alg_elem_finalize(&delta);quat_alg_elem_finalize(&quat_res);
    ibz_vec_2_finalize(&C);
    quat_left_ideal_finalize(&lideal);
    quat_alg_coord_finalize(&t1);quat_alg_coord_finalize(&t2);
    quat_order_finalize(&right_order);
    quat_alg_elem_finalize(&shift_gamma);quat_alg_elem_finalize(&quat_one);
    quat_alg_finalize(&Bpoo);

     
    return found;
}

int klpt_test_represent_integer() {
    int found = 0;
    ibz_t p,temp;
    ibz_t M,M_begin;
    quat_alg_t Bpoo;
    quat_alg_elem_t gamma;
    
    quat_alg_elem_init(&gamma);
    ibz_init(&M);ibz_init(&M_begin);
    ibz_init(&p);
    ibz_init(&temp);    


    int exp = 128;
    found = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);
    
    found = found && (generate_random_prime(&M,1,exp));
    ibz_pow(&temp,&ibz_const_two,exp+15);
    ibz_mul(&M,&M,&temp);
    ibz_copy(&M_begin,&M);

    

    found = found && represent_integer(&gamma,&M,&Bpoo);
    if (!found) {
        printf("repesent integer did not find anything \n");
    }    

    ibq_t norm;
    ibq_init(&norm);
    ibz_t ibz_norm;
    ibz_init(&ibz_norm);
    quat_alg_norm(&norm,&gamma,&Bpoo);


    found = found && ibq_to_ibz(&ibz_norm,&norm);

    found = found && (ibz_cmp(&ibz_norm,&M)==0); 
    if (! found ){
        printf("unequality of norm \n");
        ibz_printf(" %Zd \n %Zd \n",ibz_norm,M);
    }

    ibz_t remainder;
    ibz_init(&remainder);
    ibz_div(&temp,&remainder,&M_begin,&M);

    found = found && (ibz_cmp(&remainder,&ibz_const_zero)==0);

    quat_alg_elem_finalize(&gamma);
    ibz_finalize(&M);
    quat_alg_finalize(&Bpoo);
    ibz_finalize(&p);
    ibz_finalize(&temp);
    ibz_finalize(&M_begin);
    ibz_finalize(&remainder); 
    ibq_finalize(&norm);
    ibz_finalize(&ibz_norm); 
    return found;
}

int klpt_test_tools() {

    int res = 1; 
    res =kltp_test_keygen_random_ideal();
    assert(res);
    for (int i =0; i<2;i++) {
        res= res && klpt_test_lideal_generator_coprime();
    }
    assert(res);
    for (int i =0; i<2;i++) {
        res= res && klpt_test_lideal_isom(i);
    }
    assert(res);
    for (int i =0; i<2;i++) {
        res= res && klpt_test_connecting_ideal(i);
    }
    assert(res);

    printf("Running klpt tests for represent integer \n \n");

    for (int i =0; i<3;i++) {
        res= res && klpt_test_represent_integer();
    }
    if (!res) {
        printf("KLPT unit test represent_integer failed\n");
    }
    printf("Running klpt tests for solve combi eichler \n \n");

    for (int i =0; i<3;i++) {
        res= res && klpt_test_solve_combi_eichler();
    }
    if (!res) {
        printf("KLPT unit test solve_combi_eichler failed\n");
    }

    printf("Running klpt tests for finding linear combination \n \n");

    for (int i =0; i<3;i++) {
        res= res && klpt_test_find_linear_comb();
    }
    if (!res) {
        printf("KLPT unit test find_linear_comb failed\n");
    }
    return res;
}
