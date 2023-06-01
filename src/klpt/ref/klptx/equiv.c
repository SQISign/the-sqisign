#include <klpt.h>
#include "tools.h"


/**
 * @brief Keygen random ideal
 * @param ideal : Output : random ideal
 * @param order : maximal extremal order
 * @param Bpoo the quaternion algebra
 * computes a keygen ideal
 */
int klpt_keygen_random_ideal(quat_left_ideal_t *ideal, const quat_p_extremal_maximal_order_t *order, const quat_alg_t *Bpoo) {

    int found;
    ibz_t temp;
    quat_alg_elem_t gamma;quat_alg_elem_t quat_temp;
    quat_alg_coord_t coeffs;
    ibz_t n;

    ibz_init(&temp);ibz_init(&n);
    quat_alg_elem_init(&gamma);quat_alg_elem_init(&quat_temp);
    quat_alg_coord_init(&coeffs);

    // we start by sampling the random norm n 
    generate_random_prime(&n,1,KLPT_secret_key_prime_size);

    // we start by sampling a random element of norm n 2^log(p)
    ibz_pow(&temp,&ibz_const_two,2*ibz_bitsize(&CHARACTERISTIC));
    ibz_mul(&temp,&temp,&n);
    found = represent_integer(&gamma,&temp,Bpoo);
    if (!found) {
        return 0;
    }
    // then we generate a random scalar
    ibz_rand_interval(&temp,&ibz_const_zero,&n);

    // if temp != 0 then we set gamma = gamma * (temp + i)
    if (0!=ibz_cmp(&temp,&ibz_const_zero)) {
        ibz_copy(&coeffs[0],&temp);
        ibz_copy(&coeffs[1],&ibz_const_one);
        ibz_copy(&coeffs[2],&ibz_const_zero);
        ibz_copy(&coeffs[3],&ibz_const_zero);
        order_elem_create(&quat_temp,order,&coeffs,Bpoo);
        quat_alg_mul(&gamma,&gamma,&quat_temp,Bpoo);
    }

    // we compute the output ideal O0 <gamma, n>
    quat_lideal_create_from_primitive(ideal,&gamma,&n,&order->order,Bpoo);


    ibz_finalize(&temp);ibz_finalize(&n);
    quat_alg_elem_finalize(&gamma);quat_alg_elem_finalize(&quat_temp);
    quat_alg_coord_finalize(&coeffs);
    return found;
}

/**
 * @brief Equivalent left ideal eichler randomized
 *
 * @param gen Output: generator of equiv  
 * @param n level of the eichler order
 * @param lideal quaternion ideal
 * @param Bpoo the quaternion algebra
 * finds an ideal randomized in the class of eichler orders of level n
 * assumes that n is prime
 */
void klpt_lideal_equiv_random_eichler(quat_alg_elem_t *gen, const ibz_t *n, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {
    

    // declaration of variables
    int ctr = 0;
    int check = 0;
    int found = 0;
    ibq_t ibq_norm;
    ibz_t temp;ibz_t disc;
    quat_alg_coord_t coeffs;
    quat_alg_elem_t beta,quat_temp,gamma;

    // var init
    quat_alg_coord_init(&coeffs);
    ibz_init(&temp);ibz_init(&disc);
    ibq_init(&ibq_norm);
    quat_alg_elem_init(&beta);
    quat_alg_elem_init(&gamma);
    quat_alg_elem_init(&quat_temp);

    // the loop : we are looking for our nice element
    // the number of tries is bounded
    while (!found && ctr < KLPT_equiv_num_iter) {
        ctr++;
        // we select our linear combination at random
        ibz_rand_interval_minm_m(&coeffs[0],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[1],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[2],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[3],KLPT_equiv_bound_coeff);

        // computation of the quaternion element
        // as beta = \sum_{0 <= i <= 3 } coeffs[i] * bas[i] 
        ibz_mat_4x4_eval(&beta.coord,&lideal->parent_order->basis,&coeffs);
        ibz_copy(&beta.denom,&lideal->parent_order->denom);

        
        // now we need to check that n is inert in the quadratic order ZZ[beta]
        // computing the discriminant n(beta) - 4*trace(beta)^2 
        quat_alg_norm(&ibq_norm,&beta,Bpoo);
        check = ibq_to_ibz(&disc,&ibq_norm);
        assert(&check);
        quat_alg_trace(&ibq_norm,&beta);
        check = ibq_to_ibz(&temp,&ibq_norm);
        assert(&check);
        ibz_mul(&temp,&temp,&temp);
        ibz_mul(&temp,&temp,&ibz_const_two);
        ibz_mul(&temp,&temp,&ibz_const_two);
        ibz_sub(&disc,&disc,&temp);
        ibz_neg(&disc,&disc);
        // testing that n is inert in ZZ[beta]
        ibz_mod(&temp,&disc,n);
        found = (0!=ibz_cmp(&temp,&ibz_const_zero)) && (-1 == ibz_legendre(&disc,n));
        
    }
    // we have found the nice element beta
    // now we need to find the nice element in lideal
    ctr = 0;
    found = 0;
    while (!found && ctr < KLPT_equiv_num_iter) {
        ctr++;
        // we select our linear combination at random
        ibz_rand_interval_minm_m(&coeffs[0],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[1],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[2],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[3],KLPT_equiv_bound_coeff);

        // computation of the quaternion element
        // as beta = \sum_{0 <= i <= 3 } coeffs[i] * bas[i] 
        ibz_mat_4x4_eval(&gamma.coord,&lideal->lattice.basis,&coeffs);
        ibz_copy(&gamma.denom,&lideal->lattice.denom);

        
        // now we need to check that n is coprime with norm of gamma
        quat_alg_norm(&ibq_norm,&gamma,Bpoo);
        ibq_to_ibz(&disc,&ibq_norm);
        
        ibz_gcd(&temp,&disc,n);
        // testing that the norm is coprime to 2*n
        found = (0==ibz_cmp(&temp,&ibz_const_one));
    }

    // now we can sample the random linear combination
    ibz_rand_interval(&temp,&ibz_const_zero,n);

    // if temp != 0 then we compute gamma as (temp + beta) * gamma  
    if (0!=ibz_cmp(&temp,&ibz_const_zero)) {
        quat_alg_scalar(&quat_temp,&temp,&ibz_const_one);
        quat_alg_add(&beta,&beta,&quat_temp);
        quat_alg_mul(&gamma,&beta,&gamma,Bpoo);
    } 

    // copying the result to gen
    quat_alg_conj(gen,&gamma);
    

    // freeing the variables
    quat_alg_coord_finalize(&coeffs);
    ibz_finalize(&temp);ibz_finalize(&disc);
    ibq_finalize(&ibq_norm);
    quat_alg_elem_finalize(&beta);
    quat_alg_elem_finalize(&gamma); 
    quat_alg_elem_finalize(&quat_temp);
} 

/**
 * @brief Equivalent left ideal
 *
 * @param gen Output: generator of equiv  
 * @param n Output : norm of the equivalent ideal
 * @param reduced the reduced basis for the ideal in input
 * @param gram the gram matrix of the reduced basis
 * @param lideal_norm the norm of the ideal represented by reduced basis
 * @param denom integer, the denominator of the ideal 
 * @param Bpoo the quaternion algebra
 * @return a bit indicating if the computation succeeded
 * Assumes that the basis is reduced
 * there is an equivalent ideal of norm n such that equiv_lideal = lideal * gen/Norm(lideal)
 * this search is randomized and may fail
 */
int klpt_lideal_equiv(quat_alg_elem_t *gen, ibz_t *n, const ibz_mat_4x4_t *reduced, const ibz_mat_4x4_t *gram, const ibz_t *lideal_norm, const ibz_t *denom, const quat_alg_t *Bpoo) {
    

    // declaration of variables
    int ctr = 0;
    int found = 0;
    ibz_t norm_v,remainder,quotient_norm_v,adjusted_norm;
    quat_alg_coord_t coeffs;

    // var init
    quat_alg_coord_init(&coeffs);
    ibz_init(&norm_v);
    ibz_init(&remainder);
    ibz_init(&quotient_norm_v);
    ibz_init(&adjusted_norm);

    // dividing by the denom to get the real norm from the gram matrix
    ibz_mul(&adjusted_norm,lideal_norm,denom);
    ibz_mul(&adjusted_norm,&adjusted_norm,denom);

    // the loop : we are looking for our nice element
    // the number of tries is bounded
    while (!found && ctr < KLPT_equiv_num_iter) {
        ctr++;
        // we select our linear combination at random
        ibz_rand_interval_minm_m(&coeffs[0],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[1],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[2],KLPT_equiv_bound_coeff);
        ibz_rand_interval_minm_m(&coeffs[3],KLPT_equiv_bound_coeff);

        //computation of the norm of the vector sampled
        quat_qf_eval(&norm_v,gram,&coeffs);

        // compute the norm of the equivalent ideal
        // can be improved by removing the power of two first and the odd part only if the trial division failed (this should always be called on an ideal of norm 2^x * N for some big prime N ) 
        ibz_div(&quotient_norm_v,&remainder,&norm_v,&adjusted_norm);

        // debug : check that the remainder is zero
	    assert(ibz_is_zero(&remainder));
        
        // pseudo-primality test
        if (ibz_probab_prime(&quotient_norm_v,KLPT_primality_num_iter)) { 
            // computes the generator using a matrix multiplication
            ibz_mat_4x4_eval(&gen->coord,reduced,&coeffs);

            ibz_copy(&gen->denom,denom);

            quat_alg_conj(gen,gen);
            ibz_copy(n,&quotient_norm_v);
            found = 1; break; 
        }

  
    }

    // freeing the variables
    ibz_finalize(&norm_v);
    ibz_finalize(&remainder);
    ibz_finalize(&quotient_norm_v);
    ibz_finalize(&adjusted_norm);
    quat_alg_coord_finalize(&coeffs);

    return found;    
} 
