#include <inttypes.h>

#include "id2iso_tests.h"

//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.im[i]);
    printf("\n");
}

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}
//XXX

int id2iso_test_long_id2iso() {
    // var dec
    int found =1;  
    ibz_t temp,remainder,n;
    ibq_t ibq_norm;
    quat_alg_elem_t gen,gen_two,gen_check,quat_temp,gen_key;

    quat_left_ideal_t lideal_small,lideal_check,lideal_two;
    quat_left_ideal_t lideal_two_one,lideal_small_one;
    quat_order_t right_order;
    id2iso_compressed_long_two_isog_t zip;
    ec_isog_even_t isog_two_one;
    quat_alg_coord_t coeffs;
    ibz_mat_4x4_t reduced,gram;

    ec_basis_t basis_plus,basis_minus;
    ec_basis_t odd_basis[2];
    ec_curve_t curve = {0};
    ec_point_t kernel_dual = {0};
    ec_basis_t even_basis = {0};
 
    // var init
    ibq_init(&ibq_norm);
    ibz_init(&temp);ibz_init(&remainder);ibz_init(&n);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_two);
    quat_alg_elem_init(&gen_key);
    quat_alg_elem_init(&gen_check);
    quat_alg_elem_init(&quat_temp);
    quat_left_ideal_init(&lideal_small);
    quat_left_ideal_init(&lideal_check);
    quat_left_ideal_init(&lideal_two);
    quat_left_ideal_init(&lideal_two_one);
    quat_left_ideal_init(&lideal_small_one);
    quat_order_init(&right_order);
    quat_alg_coord_init(&coeffs);
    ibz_mat_4x4_init(&reduced);ibz_mat_4x4_init(&gram);
    // computation of lideal_small
    generate_random_prime(&n,1,KLPT_secret_key_prime_size);
    ibz_mul(&temp,&n,&TORSION_ODD);
    found = found && represent_integer(&gen_check,&temp,&QUATALG_PINFTY);
    assert(found);
    quat_lideal_create_from_primitive(&lideal_small,&gen_check,&n,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
    quat_alg_conj(&gen_check,&gen_check);
    quat_lideal_create_from_primitive(&lideal_check,&gen_check,&TORSION_ODD,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
    assert(quat_lideal_isom(&quat_temp,&lideal_small,&lideal_check,&QUATALG_PINFTY)); 

    // application of keygen_klpt to compute a generator of lideal_two 
    found = found && klpt_keygen_klpt(&gen,&lideal_small,&QUATALG_PINFTY);
    quat_alg_conj(&gen,&gen);
    quat_alg_elem_copy(&gen_key,&gen);

    

    // TODEBUG check that this is indeed an integer
    int length = KLPT_keygen_length/TORSION_PLUS_EVEN_POWER;
    assert(KLPT_keygen_length==TORSION_PLUS_EVEN_POWER*length); 

    // computation of the norm of lideal_two 
    quat_alg_norm(&ibq_norm,&gen,&QUATALG_PINFTY);
    ibq_to_ibz(&temp,&ibq_norm);
    ibz_div(&temp,&remainder,&temp,&lideal_small.norm);
    assert(0==ibz_cmp(&remainder,&ibz_const_zero));
    ibz_pow(&remainder,&ibz_const_two,KLPT_keygen_length);
    assert(0==ibz_cmp(&remainder,&temp));

    // computation of lideal_two
    quat_lideal_create_from_primitive(&lideal_two,&gen,&temp,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);

    // we compute a better generator
    int find_gen = 0;
    while (!find_gen) {
        ibz_rand_interval_minm_m(&coeffs[0],100);
        ibz_rand_interval_minm_m(&coeffs[1],100);
        ibz_rand_interval_minm_m(&coeffs[2],100);
        ibz_rand_interval_minm_m(&coeffs[3],100);
        ibz_mat_4x4_eval(&gen.coord,&lideal_two.lattice.basis,&coeffs);
        ibz_copy(&gen.denom,&lideal_two.lattice.denom);

        quat_alg_trace(&ibq_norm,&gen);
        ibq_to_ibz(&temp,&ibq_norm);
        find_gen = (0!=ibz_get(&temp)%2);
    }

    assert(quat_lideal_isom(&quat_temp,&lideal_two,&lideal_small,&QUATALG_PINFTY));
     assert(quat_lideal_isom(&quat_temp,&lideal_two,&lideal_check,&QUATALG_PINFTY));
    // computing the ideal of the first step 
    
    ibz_pow(&remainder,&ibz_const_two,TORSION_PLUS_EVEN_POWER);
    quat_lideal_create_from_primitive(&lideal_two_one,&gen,&remainder,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);


    // computing the right_order
    quat_lideal_right_order(&right_order,&lideal_two_one,&QUATALG_PINFTY);


    assert(quat_lattice_contains(&coeffs,&right_order,&gen,&QUATALG_PINFTY));

    // computing an ideal equivalent to lideal_two_one
    quat_lideal_reduce_basis(&reduced,&gram,&lideal_two_one,&QUATALG_PINFTY);
    found = found && klpt_lideal_equiv(&gen_two,&temp,&reduced,&gram,&lideal_two_one.norm,&lideal_two_one.lattice.denom,&QUATALG_PINFTY);
    quat_lideal_create_from_primitive(&lideal_small_one,&gen_two,&temp,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
    quat_alg_conj(&gen_two,&gen_two);
    assert(quat_lattice_contains(&coeffs,&lideal_two_one.lattice,&gen_two,&QUATALG_PINFTY));
    assert(quat_lideal_isom(&quat_temp,&lideal_two_one,&lideal_small_one,&QUATALG_PINFTY));

    quat_left_ideal_t ideal_test;
    quat_left_ideal_init(&ideal_test);
    quat_lideal_create_from_primitive(&ideal_test,&gen_two,&lideal_two_one.norm,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
    assert(quat_lideal_equals(&ideal_test,&lideal_two_one,&QUATALG_PINFTY));

    // casting gen into the right order of lideal_small_one
    quat_alg_elem_copy(&quat_temp,&gen_two);
    ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal_two_one.norm);
    ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal_small_one.norm);
    quat_alg_mul(&gen,&quat_temp,&gen,&QUATALG_PINFTY);
    quat_alg_conj(&quat_temp,&gen_two);
    quat_alg_mul(&gen,&gen,&quat_temp,&QUATALG_PINFTY);

    #ifndef NDEBUG 
        quat_lideal_right_order(&right_order,&lideal_small_one,&QUATALG_PINFTY);
        assert(quat_lattice_contains(&coeffs,&right_order,&gen,&QUATALG_PINFTY));
        int lideal_generator_ok1= quat_lideal_generator_coprime(&quat_temp,&lideal_small_one,&ibz_const_two,&QUATALG_PINFTY,0);
        assert(lideal_generator_ok1);
        quat_alg_mul(&quat_temp,&quat_temp,&gen,&QUATALG_PINFTY);
        ibz_pow(&temp,&ibz_const_two,TORSION_PLUS_EVEN_POWER*(length -1) );
        ibz_mul(&temp,&temp,&lideal_small_one.norm);
        quat_lideal_create_from_primitive(&ideal_test,&quat_temp,&temp,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
        assert(quat_lideal_isom(&quat_temp,&ideal_test,&lideal_check,&QUATALG_PINFTY));
        assert(quat_lideal_isom(&quat_temp,&ideal_test,&lideal_small,&QUATALG_PINFTY));
        if (!quat_lideal_isom(&quat_temp,&ideal_test,&lideal_check,&QUATALG_PINFTY)) {
            printf("the final check will fail ! \n");
        }
    #endif

    // copying the precomputed basis
    basis_plus = BASIS_ODD_PLUS;
    basis_minus = BASIS_ODD_MINUS;
    odd_basis[0] = basis_minus;
    odd_basis[1] = basis_plus;
    
    // copying the starting curve
    curve = CURVE_E0;

    // Computing the first isogeny 
    id2iso_ideal_to_isogeny_even(&isog_two_one,&lideal_two_one);

    // computation of the deterministic second point as the kernel dual
    ec_complete_basis_2(&even_basis,&curve,&isog_two_one.kernel);
    kernel_dual = even_basis.Q;

    // evaluation through the first isogeny
    ec_eval_even_basis(&curve,&isog_two_one,odd_basis,2);

    basis_minus = odd_basis[0]; basis_plus=odd_basis[1];
    ec_curve_t E;
    ec_eval_even(&E,&isog_two_one,&kernel_dual,1);


    id2iso_compressed_long_two_isog_init(&zip, length);

    found = found && id2iso_ideal_to_isogeny_two_long_power_of_2(&zip,&curve,&basis_minus,&basis_plus,&kernel_dual,&gen,length-1,&lideal_small_one,&lideal_two_one,&gen_two,&QUATALG_PINFTY);
    
    
    if ( found ) {
        ec_isog_odd_t isog_check;
        id2iso_ideal_to_isogeny_odd(&isog_check,&CURVE_E0,&BASIS_ODD_PLUS,&BASIS_ODD_MINUS,&lideal_check); 
        ec_curve_t new_curve=CURVE_E0;
        ec_eval_odd(&new_curve,&isog_check,&even_basis.P,1);

        // checking equality of j_inv
        fp2_t j1,j2;
        ec_j_inv(&j1,&new_curve);
        ec_j_inv(&j2,&curve);
        found = found && fp2_is_equal(&j1,&j2);
        if (!fp2_is_equal(&j1,&j2)) {
            printf(" the final check failed \n");
        }
    }
    else {
        printf("the long translation failed \n");
    }

  

    // var finalize    
    quat_alg_elem_finalize(&gen_key);
    quat_alg_elem_finalize(&gen_check);
    ibq_finalize(&ibq_norm);
    ibz_finalize(&temp);ibz_finalize(&remainder);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_two);
    quat_left_ideal_finalize(&lideal_two_one);
    quat_left_ideal_finalize(&lideal_small_one);
    quat_order_finalize(&right_order);
    quat_alg_coord_finalize(&coeffs);
    ibz_mat_4x4_finalize(&reduced);ibz_mat_4x4_finalize(&gram);
    id2iso_compressed_long_two_isog_finalize(&zip);

    return found;
}

int id2iso_test_id2iso() {
    int res = 1;
    printf("\n \nRunning id2iso tests for long id2iso \n \n");


     for (int i =0; i<3;i++) {
        res = res && id2iso_test_long_id2iso();
    }
     
    if (!res) {
        printf("ID2ISO unit test long id2iso failed\n");
    }

    return res;
}
