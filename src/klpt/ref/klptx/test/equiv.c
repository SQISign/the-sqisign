#include "klpt_tests.h"

int klpt_test_equiv_eichler() {
    
    int res =1;

    ibz_t M,temp;
    quat_order_t order;
    quat_left_ideal_t id1,id2,id3;
    quat_alg_coord_t coeffs;
    quat_alg_elem_t gen,quat_temp;

    ibz_init(&M);
    quat_order_init(&order);ibz_init(&temp);
    quat_left_ideal_init(&id1);quat_left_ideal_init(&id2);
    quat_left_ideal_init(&id3);
    quat_alg_coord_init(&coeffs);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&quat_temp);

    int exp = 256;
    generate_random_prime(&M,0,exp);
    

    quat_lideal_random_2e(&id1,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY,exp,'5');
    quat_lideal_right_order(&order,&id1,&QUATALG_PINFTY);

    quat_lideal_random_2e(&id2,&order,&QUATALG_PINFTY,exp,'5');

    

    klpt_lideal_equiv_random_eichler(&gen,&M,&id2,&QUATALG_PINFTY);

    quat_alg_make_primitive(&coeffs,&temp,&gen,&order,&QUATALG_PINFTY);

    quat_alg_elem_copy(&quat_temp,&gen);
    ibz_mul(&quat_temp.denom,&quat_temp.denom,&id2.norm);

    int lideal_mul_ok = quat_lideal_mul(&id3,&id2,&quat_temp,&QUATALG_PINFTY,0);
    assert(lideal_mul_ok);

    assert(quat_lideal_isom(&quat_temp,&id3,&id2,&QUATALG_PINFTY));


    quat_alg_conj(&gen,&gen);

    assert(quat_lattice_contains(&coeffs,&id2.lattice,&gen,&QUATALG_PINFTY));

    ibz_finalize(&M);ibz_finalize(&temp);
    quat_order_finalize(&order);
    quat_left_ideal_finalize(&id1);quat_left_ideal_finalize(&id2);
    quat_alg_coord_finalize(&coeffs);
    quat_alg_elem_finalize(&gen);
    quat_left_ideal_finalize(&id3);
    quat_alg_elem_finalize(&quat_temp);

    return res;

}

int klpt_test_equiv_intern() {

    int res;
    quat_alg_elem_t gen;
    ibz_t n,p;
    ibz_mat_4x4_t reduced,gram;
    
    quat_alg_t Bpoo;
    quat_left_ideal_t lideal;

    quat_alg_elem_init(&gen);
    ibz_init(&n);ibz_init(&p);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&reduced);
    
    quat_left_ideal_init(&lideal);

    int exp = 128;
    res = generate_random_prime(&p,1,2*exp);
    quat_alg_init_set(&Bpoo,&p);
    quat_lideal_random_2e(&lideal,&STANDARD_EXTREMAL_ORDER.order,&Bpoo,2*exp,'8');

    quat_lideal_reduce_basis(&reduced,&gram,&lideal,&Bpoo);


    res = res && klpt_lideal_equiv(&gen,&n,&reduced,&gram,&lideal.norm,&lideal.lattice.denom,&Bpoo);
    assert(res);
    

    ibq_t ibq_norm;
    ibz_t norm,temp;
    ibq_init(&ibq_norm);
    ibz_init(&norm);ibz_init(&temp);


    quat_alg_norm(&ibq_norm,&gen,&Bpoo);
    res = res && ibq_to_ibz(&norm,&ibq_norm);
    if (!res) {
        printf("the norm is not an ibz \n");
    }
    ibz_div(&norm,&temp,&norm,&lideal.norm);

    res = res && ibz_cmp(&norm,&n)==0; 
    if (!res) {
        printf("the norms are not equal!\n");
        ibz_printf("Zd %Zd \n",norm,n);
    }

    quat_alg_elem_t conj;
    quat_alg_coord_t coord;
    quat_left_ideal_t lideal_test;
    quat_alg_elem_init(&conj);
    quat_alg_coord_init(&coord);
    quat_left_ideal_init(&lideal_test);

    quat_alg_conj(&conj,&gen);
    res = res && quat_lattice_contains(&coord,&lideal.lattice,&conj,&Bpoo);

    quat_lideal_create_from_primitive(&lideal_test,&conj,&lideal.norm,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    res = res && quat_lideal_equals(&lideal,&lideal_test,&Bpoo);

    quat_lideal_create_from_primitive(&lideal_test,&gen,&n,&STANDARD_EXTREMAL_ORDER.order,&Bpoo);
    res = res && quat_lideal_isom(&conj,&lideal,&lideal_test,&Bpoo);

    quat_alg_elem_finalize(&conj);
    quat_alg_coord_finalize(&coord);
    quat_left_ideal_finalize(&lideal_test);
    ibq_finalize(&ibq_norm);    
    ibz_finalize(&norm);ibz_finalize(&temp);
    quat_alg_elem_finalize(&gen);
    ibz_finalize(&n);ibz_finalize(&p);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&reduced);
    quat_alg_finalize(&Bpoo);
    quat_left_ideal_finalize(&lideal);

    return res;

}

int klpt_test_equiv() {
    int res = 1;

    printf("Running klpt tests for equivalent ideal \n \n");


    for (int i =0; i<2;i++) {
        res = res & klpt_test_equiv_intern();
        
    }
    for (int i =0; i< 2 ;i++) {
        res = res & klpt_test_equiv_eichler();
    }
     
    if (!res) {
        printf("KLPT unit test klpt_lideal_equiv failed\n");
    }

    return res;
}
