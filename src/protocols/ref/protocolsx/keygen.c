#include <protocols.h>


void secret_key_init(secret_key_t *sk) {
    quat_left_ideal_init(&sk->lideal_small);
    quat_left_ideal_init(&sk->lideal_two);
    quat_alg_elem_init(&sk->gen_two);
}

void secret_key_finalize(secret_key_t *sk) {
    quat_left_ideal_finalize(&sk->lideal_small);
    quat_left_ideal_finalize(&sk->lideal_two);
    quat_alg_elem_finalize(&sk->gen_two);
}


/**
 * @brief Computing a key pair
 *
 * @param pk : Output the public key 
 * @param sk : Output the secret key
 * @returns a bit indicating if the computation succeeded  
 * assumes that sk and pk have been initialized
    */
int protocols_keygen(public_key_t *pk, secret_key_t *sk) {

    // var dec
    int found = 1;  
    ibz_t temp,remainder;
    ibq_t ibq_norm;
    quat_alg_elem_t gen,gen_two;
    quat_alg_elem_t quat_temp;

    quat_left_ideal_t lideal_two_one,lideal_small_one;
    quat_order_t right_order;
    id2iso_compressed_long_two_isog_t zip;
    quat_alg_coord_t coeffs;
    ibz_mat_4x4_t reduced,gram;

    ec_basis_t basis_plus,basis_minus;
    ec_basis_t odd_basis[2];
    ec_curve_t curve = {0};
    ec_point_t kernel_dual = {0};
    ec_basis_t even_basis = {0};
    ec_isog_even_t isog_two_one;
 
 
    // var init
    ibq_init(&ibq_norm);
    ibz_init(&temp);ibz_init(&remainder);
    quat_alg_elem_init(&gen);quat_alg_elem_init(&quat_temp);
    quat_alg_elem_init(&gen_two);
    quat_left_ideal_init(&lideal_two_one);
    quat_left_ideal_init(&lideal_small_one);
    quat_order_init(&right_order);
    quat_alg_coord_init(&coeffs);
    ibz_mat_4x4_init(&reduced);ibz_mat_4x4_init(&gram);

    // computing the length of the walk
    int length = KLPT_keygen_length/TORSION_PLUS_EVEN_POWER; 
    assert(KLPT_keygen_length==TORSION_PLUS_EVEN_POWER*length);

    found = 0; 
    int cnt = 0;
    while (!found && cnt < SQISIGN_response_attempts) {
        cnt ++;

    // computation of lideal_small
    found = klpt_keygen_random_ideal(&sk->lideal_small,&STANDARD_EXTREMAL_ORDER,&QUATALG_PINFTY);
    if (!found) {
        continue;
    }
    

    // application of keygen_klpt to compute a generator of lideal_two 
    found = found && klpt_keygen_klpt(&gen,&sk->lideal_small,&QUATALG_PINFTY);
    quat_alg_normalize(&gen);
    assert(quat_lattice_contains(&coeffs,&sk->lideal_small.lattice,&gen,&QUATALG_PINFTY));
    #ifndef NDEBUG 
        quat_left_ideal_t ideal_test;
        quat_left_ideal_init(&ideal_test);
        quat_lideal_make_primitive_then_create(&ideal_test,&gen,&sk->lideal_small.norm,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
        assert(quat_lideal_equals(&ideal_test,&sk->lideal_small,&QUATALG_PINFTY));
        quat_left_ideal_finalize(&ideal_test);
    #endif

    // gen = conjugate(gen)
    quat_alg_conj(&gen,&gen);

    // sk.gen_two = gen
    quat_alg_elem_copy(&sk->gen_two,&gen);  

    

    // computation of the norm of lideal_two 
    quat_alg_norm(&ibq_norm,&gen,&QUATALG_PINFTY);
    ibq_to_ibz(&temp,&ibq_norm);
    ibz_div(&temp,&remainder,&temp,&sk->lideal_small.norm);
    assert(0==ibz_cmp(&remainder,&ibz_const_zero));

    // computation of lideal_two
    quat_lideal_make_primitive_then_create(&sk->lideal_two,&gen,&temp,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
    // for debug we check that the two ideals are equivalent
    assert(quat_lideal_isom(&quat_temp,&sk->lideal_two,&sk->lideal_small,&QUATALG_PINFTY));

    // we compute a better generator to make the translation 
    // for this, we look for a generator with odd trace
    int find_gen = 0;
    while (!find_gen) {
        ibz_rand_interval_minm_m(&coeffs[0],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[1],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[2],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[3],KLPT_equiv_bound_coeff);
        ibz_mat_4x4_eval(&gen.coord,&sk->lideal_two.lattice.basis,&coeffs);
        ibz_copy(&gen.denom,&sk->lideal_two.lattice.denom);

        quat_alg_trace(&ibq_norm,&gen);
        ibq_to_ibz(&temp,&ibq_norm);
        find_gen = (0!=ibz_get(&temp)%2);
    }

    // computing the ideal of the first step as lideal_two_one = O0 < gen,2^f >
    ibz_pow(&remainder,&ibz_const_two,TORSION_PLUS_EVEN_POWER);
    quat_lideal_create_from_primitive(&lideal_two_one,&gen,&remainder,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);

    // computing the right_order
    quat_lideal_right_order(&right_order,&lideal_two_one,&QUATALG_PINFTY);

    // for debug checking that gen is indeed contained in right_order
    assert(quat_lattice_contains(&coeffs,&right_order,&gen,&QUATALG_PINFTY));

    // computing an ideal equivalent to lideal_two_one
    quat_lideal_reduce_basis(&reduced,&gram,&lideal_two_one,&QUATALG_PINFTY);
    found = found && klpt_lideal_equiv(&gen_two,&temp,&reduced,&gram,&lideal_two_one.norm,&lideal_two_one.lattice.denom,&QUATALG_PINFTY);
    if (!found) {
        continue;
    }

    // lideal_small = O0 < gen_two , temp>
    quat_lideal_create_from_primitive(&lideal_small_one,&gen_two,&temp,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
    quat_alg_conj(&gen_two,&gen_two);
    assert(quat_lattice_contains(&coeffs,&lideal_two_one.lattice,&gen_two,&QUATALG_PINFTY));
    assert(quat_lideal_isom(&quat_temp,&lideal_two_one,&lideal_small_one,&QUATALG_PINFTY));

    // casting gen into the right order of lideal_small_one 
    // gen = gen_two * gen * gen_two^ {-1}
    quat_alg_elem_copy(&quat_temp,&gen_two);
    ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal_two_one.norm);
    ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal_small_one.norm);
    quat_alg_mul(&gen,&quat_temp,&gen,&QUATALG_PINFTY);
    quat_alg_conj(&quat_temp,&gen_two);
    quat_alg_mul(&gen,&gen,&quat_temp,&QUATALG_PINFTY);
    

    #ifndef NDEBUG 
        quat_left_ideal_init(&ideal_test);
        quat_lideal_right_order(&right_order,&lideal_small_one,&QUATALG_PINFTY);
        assert(quat_lattice_contains(&coeffs,&right_order,&gen,&QUATALG_PINFTY));
        int lideal_generator_ok = quat_lideal_generator_coprime(&quat_temp,&lideal_small_one,&ibz_const_two,&QUATALG_PINFTY,0);
        assert(lideal_generator_ok);
        quat_alg_mul(&quat_temp,&quat_temp,&gen,&QUATALG_PINFTY);
        ibz_pow(&temp,&ibz_const_two,TORSION_PLUS_EVEN_POWER*(length -1) );
        ibz_mul(&temp,&temp,&lideal_small_one.norm);
        quat_lideal_create_from_primitive(&ideal_test,&quat_temp,&temp,&STANDARD_EXTREMAL_ORDER.order,&QUATALG_PINFTY);
        assert(quat_lideal_isom(&quat_temp,&ideal_test,&sk->lideal_small,&QUATALG_PINFTY));
        quat_left_ideal_finalize(&ideal_test);
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
    ec_curve_t E = CURVE_E0;
    ec_eval_even(&E,&isog_two_one,&kernel_dual,1);

    // init of the zip
    id2iso_compressed_long_two_isog_init(&zip, length-1);


    // Computing the rest of the chain 
    found = found && id2iso_ideal_to_isogeny_two_long_power_of_2(&zip,&curve,&basis_minus,&basis_plus,&kernel_dual,&gen,length-1,&lideal_small_one,&lideal_two_one,&gen_two,&QUATALG_PINFTY);
    if (!found) {
        continue;
    }
    }

    sk->curve = curve;
    sk->kernel_dual = kernel_dual;
    sk->basis_minus = basis_minus;
    sk->basis_plus = basis_plus;

    // normalizing the curve
    ec_isom_t isom;
    ec_curve_normalize(&sk->curve,&isom,&sk->curve);
    ec_iso_eval(&sk->kernel_dual,&isom);
    ec_iso_eval(&sk->basis_minus.P,&isom);
    ec_iso_eval(&sk->basis_minus.Q,&isom);
    ec_iso_eval(&sk->basis_minus.PmQ,&isom);
    ec_iso_eval(&sk->basis_plus.P,&isom);
    ec_iso_eval(&sk->basis_plus.Q,&isom);
    ec_iso_eval(&sk->basis_plus.PmQ,&isom);

    // setting the public key
    pk->E = sk->curve;


    // var finalize
    ibq_finalize(&ibq_norm);
    ibz_finalize(&temp);ibz_finalize(&remainder);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_two);
    quat_left_ideal_finalize(&lideal_two_one);
    quat_left_ideal_finalize(&lideal_small_one);
    quat_order_finalize(&right_order);
    quat_alg_coord_finalize(&coeffs);
    ibz_mat_4x4_finalize(&reduced);ibz_mat_4x4_finalize(&gram);
    quat_alg_elem_finalize(&quat_temp);
    id2iso_compressed_long_two_isog_finalize(&zip);


    return found;
}
