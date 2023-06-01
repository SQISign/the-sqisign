#include <klpt.h>
#include "tools.h"


/**
 * @brief Deciding if a given vector is suitable for the signing klpt algorithm
 *
 * we test that gamma * mu is primitive in the end
 * @param mu Output: a quaternion algebra element
 * @param C the vector to be tested 
 * @param params parameters
 *  
 */
int condition_signing_klpt(quat_alg_elem_t *mu, const ibz_vec_2_t *C, const void* params) {
            // we start by casting params to the correct type
            const signing_klpt_param_t *signing_klpt_param = params;
            // var dec 
            int found;
            quat_alg_elem_t quat_temp;
            quat_alg_coord_t coeffs;
            ibq_t ibq_norm;
            ibz_t norm,rhs;
            ibz_t temp;
            ibz_t a,b;
            ibz_t ibz_q;

            // var init 
            quat_alg_elem_init(&quat_temp);
            quat_alg_coord_init(&coeffs);
            ibq_init(&ibq_norm);
            ibz_init(&norm);ibz_init(&rhs);
            ibz_init(&temp);
            ibz_init(&a);ibz_init(&b);
            ibz_init(&ibz_q);
            found = 0;


            // mu = j* ( C[0] + i * C[1])
            ibz_copy(&coeffs[0],&ibz_const_zero);
            ibz_copy(&coeffs[1],&ibz_const_zero);
            ibz_copy(&coeffs[2],&(*C)[0]);
            ibz_copy(&coeffs[3],&(*C)[1]);
            order_elem_create(mu,signing_klpt_param->order,&coeffs,&signing_klpt_param->Bpoo);
            


            // norm = n(mu) 
            quat_alg_norm(&ibq_norm,mu,&signing_klpt_param->Bpoo);
            found = ibq_to_ibz(&norm,&ibq_norm);
            if (!found) {
                assert(found);
            }

            // rhs = (target_norm - norm) / n²
            ibz_sub(&rhs,&signing_klpt_param->target_norm,&norm);
            ibz_mul(&temp,&signing_klpt_param->n,&signing_klpt_param->n);
            ibz_div(&rhs,&temp,&rhs,&temp);

            // for debug, we check that the remainder is zero
            assert(ibz_cmp(&temp,&ibz_const_zero)==0);

            found = found && ibz_cornacchia_extended(&a, &b,&rhs,SMALL_PRIMES_1MOD4,sizeof(SMALL_PRIMES_1MOD4)/sizeof(*SMALL_PRIMES_1MOD4), KLPT_primality_num_iter, &PROD_SMALL_PRIMES_3MOD4);

            if (found) {
                // computing the final value of mu
                // quat_temp = a + i*b
                ibz_copy(&coeffs[2],&ibz_const_zero);
                ibz_copy(&coeffs[3],&ibz_const_zero);
                ibz_copy(&coeffs[0],&a);
                ibz_copy(&coeffs[1],&b);
                ibz_mul(&coeffs[0],&coeffs[0],&signing_klpt_param->n);
                ibz_mul(&coeffs[1],&coeffs[1],&signing_klpt_param->n);
                order_elem_create(&quat_temp,signing_klpt_param->order,&coeffs,&signing_klpt_param->Bpoo);

                ibz_sub(&temp,&signing_klpt_param->target_norm,&norm);
                quat_alg_norm(&ibq_norm,&quat_temp,&signing_klpt_param->Bpoo);
                assert(ibq_to_ibz(&norm,&ibq_norm));
                assert(ibz_cmp(&norm,&temp)==0);
                
                // mu = mu +  quat_temp
                quat_alg_add(mu,mu,&quat_temp);

                #ifndef NDEBUG 
                    // for debug, we check that the norm is indeed the target norm 
                    quat_alg_norm(&ibq_norm,mu,&signing_klpt_param->Bpoo);
                    assert(ibq_to_ibz(&norm,&ibq_norm));
                    assert(ibz_cmp(&norm,&signing_klpt_param->target_norm)==0);
                #endif

                // mu = mu / 2
                ibz_mul(&mu->denom,&mu->denom,&ibz_const_two);
                // trying if mu is contained in signing_klpt_param-> order 
                // if not is is useless to continue
                found = quat_lattice_contains(&coeffs,&signing_klpt_param->order->order,mu,&signing_klpt_param->Bpoo);

        

                if (found) {
                    // mu is in the extremal order, we only need to verify that gamma * mu is primitive
                    quat_alg_mul(&quat_temp,&signing_klpt_param->gamma,mu,&signing_klpt_param->Bpoo);
                    found = quat_alg_is_primitive(&quat_temp,&signing_klpt_param->order->order,&signing_klpt_param->Bpoo); 

                    // verifying that gamma*mu*delta will also be primitive in the right order of ideal_start
                    quat_alg_mul(&quat_temp,&quat_temp,&signing_klpt_param->delta,&signing_klpt_param->Bpoo);
                    ibz_mul(&quat_temp.denom,&quat_temp.denom,&signing_klpt_param->equiv_n);
                    quat_alg_trace(&ibq_norm,&quat_temp);
                    if(!ibq_to_ibz(&norm,&ibq_norm)) {
                        assert(0);
                    }
                    found = found && ibz_get(&norm)%2 != 0;
                }
            }

            if (found) {
                quat_alg_norm(&ibq_norm,mu,&signing_klpt_param->Bpoo);
                assert(ibq_to_ibz(&norm,&ibq_norm));
                ibz_mul(&norm,&norm,&ibz_const_two);
                ibz_mul(&norm,&norm,&ibz_const_two);
                assert(ibz_cmp(&norm,&signing_klpt_param->target_norm)==0);
            }

            // var finalize 
            quat_alg_elem_finalize(&quat_temp);
            quat_alg_coord_finalize(&coeffs);
            ibq_finalize(&ibq_norm);
            ibz_finalize(&norm);ibz_finalize(&rhs);
            ibz_finalize(&temp);
            ibz_finalize(&a);ibz_finalize(&b);
            ibz_finalize(&ibz_q);

            return found;
}  



/**
 * @brief Deciding if a given vector is suitable for the signing klpt algorithm
 *
 * we test that gamma * mu is primitive in the end
 * @param mu Output: a quaternion algebra element
 * @param C the vector to be tested 
 * @param params parameters
 *  
 */
int condition_keygen_klpt(quat_alg_elem_t *mu, const ibz_vec_2_t *C, const void* params) {
            // we start by casting params to the correct type
            const signing_klpt_param_t *signing_klpt_param = params;
            // var dec 
            int found;
            quat_alg_elem_t quat_temp;
            quat_alg_coord_t coeffs;
            ibq_t ibq_norm;
            ibz_t norm,rhs;
            ibz_t temp;
            ibz_t a,b;
            ibz_t ibz_q;

            // var init 
            quat_alg_elem_init(&quat_temp);
            quat_alg_coord_init(&coeffs);
            ibq_init(&ibq_norm);
            ibz_init(&norm);ibz_init(&rhs);
            ibz_init(&temp);
            ibz_init(&a);ibz_init(&b);
            ibz_init(&ibz_q);
            found = 0;


            // mu = j* ( C[0] + i * C[1])
            ibz_copy(&coeffs[0],&ibz_const_zero);
            ibz_copy(&coeffs[1],&ibz_const_zero);
            ibz_copy(&coeffs[2],&(*C)[0]);
            ibz_copy(&coeffs[3],&(*C)[1]);
            order_elem_create(mu,signing_klpt_param->order,&coeffs,&signing_klpt_param->Bpoo);
            


            // norm = n(mu) 
            quat_alg_norm(&ibq_norm,mu,&signing_klpt_param->Bpoo);
            found = ibq_to_ibz(&norm,&ibq_norm);
            if (!found) {
                assert(found);
            }

            // rhs = (target_norm - norm) / n²
            ibz_sub(&rhs,&signing_klpt_param->target_norm,&norm);
            ibz_mul(&temp,&signing_klpt_param->n,&signing_klpt_param->n);
            ibz_div(&rhs,&temp,&rhs,&temp);

            // for debug, we check that the remainder is zero
            assert(ibz_cmp(&temp,&ibz_const_zero)==0);

            found = found && ibz_cornacchia_extended(&a, &b,&rhs,SMALL_PRIMES_1MOD4,sizeof(SMALL_PRIMES_1MOD4)/sizeof(*SMALL_PRIMES_1MOD4), KLPT_primality_num_iter, &PROD_SMALL_PRIMES_3MOD4);

            if (found) { 

                // computing the final value of mu
                // quat_temp = a + i*b
                ibz_copy(&coeffs[2],&ibz_const_zero);
                ibz_copy(&coeffs[3],&ibz_const_zero);
                ibz_copy(&coeffs[0],&a);
                ibz_copy(&coeffs[1],&b);
                ibz_mul(&coeffs[0],&coeffs[0],&signing_klpt_param->n);
                ibz_mul(&coeffs[1],&coeffs[1],&signing_klpt_param->n);
                order_elem_create(&quat_temp,signing_klpt_param->order,&coeffs,&signing_klpt_param->Bpoo);

                ibz_sub(&temp,&signing_klpt_param->target_norm,&norm);
                quat_alg_norm(&ibq_norm,&quat_temp,&signing_klpt_param->Bpoo);
                assert(ibq_to_ibz(&norm,&ibq_norm));
                assert(ibz_cmp(&norm,&temp)==0);
                
                // mu = mu +  quat_temp
                quat_alg_add(mu,mu,&quat_temp);

                // for debug, we check that the norm is indeed the target norm
                #ifndef NDEBUG
                    quat_alg_norm(&ibq_norm,mu,&signing_klpt_param->Bpoo);
                    assert(ibq_to_ibz(&norm,&ibq_norm));
                    assert(ibz_cmp(&norm,&signing_klpt_param->target_norm)==0);
                #endif

                // mu = mu / 2
                ibz_mul(&mu->denom,&mu->denom,&ibz_const_two);

                // trying if mu is contained in signing_klpt_param-> order 
                // if not is is useless to continue
                found = quat_lattice_contains(&coeffs,&signing_klpt_param->order->order,mu,&signing_klpt_param->Bpoo);

                if (found) {
                    // mu is in the extremal order, we only need to verify that gamma * mu is primitive
                    quat_alg_mul(&quat_temp,&signing_klpt_param->gamma,mu,&signing_klpt_param->Bpoo);
                    found = quat_alg_is_primitive(&quat_temp,&signing_klpt_param->order->order,&signing_klpt_param->Bpoo);
                }

            }


            // var finalize 
            quat_alg_elem_finalize(&quat_temp);
            quat_alg_coord_finalize(&coeffs);
            ibq_finalize(&ibq_norm);
            ibz_finalize(&norm);ibz_finalize(&rhs);
            ibz_finalize(&temp);
            ibz_finalize(&a);ibz_finalize(&b);
            ibz_finalize(&ibz_q);
            return found;
}  

/**
 * @brief Equivalent left ideal of power of two norm
 *
 * @param gen Output: generator of equiv  
 * @param lideal left O0 ideal
 * @param lideal_start left O0 ideal
 * @param delta quaternion algebra element
 * @param Bpoo the quaternion algebra
 *  
 * Let J = lideal_start
 * if K = intersect(lideal ,J) and I = conjugate (J) * K  
 * equiv_lideal is equivalent to the ideal I and we have the equality
 * equiv_lideal = I * gen/Norm(lideal) of norm n = 2^KLPT_signing_klpt_length
 * where conjugate(gen) lideal 
 * moreover we need that gen * delta is contained in the eichler order ZZ + J
 * Assumes that the ideals lideal and lideal_start has a "good" norm (ie that one needs not to apply lideal_equiv)
 * returns a bit indicating if the computation has succeeded
  */
int klpt_signing_klpt(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_left_ideal_t *lideal_start, const quat_alg_elem_t *delta, const quat_alg_t *Bpoo) {

    // var dec
    int found;
    int cnt;
    int e0,e1,center;
    ibz_t n_gamma,temp_rand,n_prod,n_mu,s1,s2;
    ibq_t ibq_n_quat_C,ibq_n_quat_C_start;
    ibz_t lambda,n_quat_C,n_quat_C_start;
    quat_alg_elem_t gamma,quat_one,mu;
    quat_alg_elem_t quat_C,quat_C_start;
    ibz_vec_2_t C,C_start;
    quat_alg_coord_t coeffs;
    signing_klpt_param_t param;

    // var init
    found = 0;
    cnt = 0;
    ibz_init(&n_gamma);
    ibz_init(&n_mu);
    ibz_init(&temp_rand);
    ibz_init(&n_prod);
    ibz_init(&lambda);
    ibz_init(&s1);ibz_init(&s2);
    ibz_init(&n_quat_C);ibz_init(&n_quat_C_start);
    ibq_init(&ibq_n_quat_C);ibq_init(&ibq_n_quat_C_start);
    ibz_vec_2_init(&C);
    ibz_vec_2_init(&C_start);
    quat_alg_elem_init(&mu);
    quat_alg_elem_init(&gamma);
    quat_alg_elem_init(&quat_one);
    quat_alg_elem_init(&quat_C);
    quat_alg_elem_init(&quat_C_start);
    quat_alg_coord_init(&coeffs);
    ibz_init(&param.n);
    ibz_init(&param.equiv_n);
    ibz_init(&param.target_norm);
    quat_alg_init_set(&param.Bpoo,&Bpoo->p);
    param.order = &STANDARD_EXTREMAL_ORDER;

    // quaternion element for 1 
    quat_alg_scalar(&quat_one,&ibz_const_one,&ibz_const_one);

    // initializing the center of the sample interval as log(p) - log( n(lideal) ) + KLPT_gamma_exponent_center_shift
    center = ibz_bitsize(&Bpoo->p) - ibz_bitsize(&lideal->norm) + KLPT_gamma_exponent_center_shift ;
    
    // the main loop, iterate until a solution is found or the bound is reached  
    while (!found && cnt < KLPT_signing_num_gamma_trial) {
        // incrementing the counter
        cnt++;
        // we start by choosing the exponent of gamma at random inside some interval (depending on the choice of parameters this interval might actually be of size 1)
        // sampling the exponent
        ibz_rand_interval_i(&temp_rand,center - KLPT_gamma_exponent_interval_size, center + KLPT_gamma_exponent_interval_size); 

        // computation of the norm
        ibz_pow(&n_gamma,&ibz_const_two,ibz_get(&temp_rand));
        ibz_mul(&n_gamma,&n_gamma,&lideal->norm);

        // computing the value of gamma
        found = represent_integer(&gamma,&n_gamma,Bpoo); 

        if (!found ) {
            assert(found);
        }
        else {
            // computing the norm of the remaining part
            ibz_pow(&temp_rand,&ibz_const_two,KLPT_signing_klpt_length);
            ibz_mul(&temp_rand,&lideal->norm,&temp_rand);
            ibz_div(&n_mu,&temp_rand,&temp_rand,&n_gamma);

            // computing the first linear combination mod n(lideal) 

            found = solve_combi_eichler(&C,&STANDARD_EXTREMAL_ORDER,&gamma,&quat_one,lideal,Bpoo,1); 
            if (!found) {
                assert(found);
            }
            // computing the element quat_C = j* (C[0] + i C[1]); 
            ibz_copy(&coeffs[0],&ibz_const_zero);
            ibz_copy(&coeffs[1],&ibz_const_zero);
            ibz_copy(&coeffs[2],&C[0]);
            ibz_copy(&coeffs[3],&C[1]);
            order_elem_create(&quat_C,&STANDARD_EXTREMAL_ORDER,&coeffs,Bpoo);

            // computing the norm n_quat_C of quat_C 
            quat_alg_norm(&ibq_n_quat_C,&quat_C,Bpoo);
            found = ibq_to_ibz(&n_quat_C,&ibq_n_quat_C);
            if (!found) {
                assert(found);
            }

            // check that the quadratic reduosity condition is verified 
            // n_mu / n_quat_C must be square mod lideal.norm
            ibz_invmod(&temp_rand,&n_quat_C,&lideal->norm); 
            ibz_mul(&temp_rand,&n_mu,&temp_rand);
            found = ibz_sqrt_mod_p(&s1,&temp_rand,&lideal->norm);

            if (found) {
                // computing the second linear combination mod n(lideal_start)
                found  = solve_combi_eichler(&C_start,&STANDARD_EXTREMAL_ORDER,&gamma,delta,lideal_start,Bpoo,0); 

                
                // computing the element quat_C = j* (C[0] + i C[1]); 
                ibz_copy(&coeffs[2],&C_start[0]);
                ibz_copy(&coeffs[3],&C_start[1]);
                order_elem_create(&quat_C_start,&STANDARD_EXTREMAL_ORDER,&coeffs,Bpoo);

                // computing the norm n_quat_C of quat_C 
                quat_alg_norm(&ibq_n_quat_C_start,&quat_C_start,Bpoo);
                found = ibq_to_ibz(&n_quat_C_start,&ibq_n_quat_C_start);
                if (!found) {
                    assert(found);
                }

                // checking the quadratic reduosity condition
                // n_mu / n_quat_C_start must be square mod lideal_start.norm
                ibz_invmod(&temp_rand,&n_quat_C_start,&lideal_start->norm); 
                ibz_mul(&temp_rand,&n_mu,&temp_rand);
                found = ibz_sqrt_mod_p(&s2,&temp_rand,&lideal_start->norm);
                if (found) {
                    // computationn of the product 
                    ibz_mul(&n_prod,&lideal->norm,&lideal_start->norm);

                    // computation of the value of lambda
                    ibz_crt(&lambda,&s1,&s2,&lideal->norm,&lideal_start->norm);

                    // multiply lambda by 2 and the final norm by 4
                    ibz_mul(&temp_rand,&n_mu,&ibz_const_two);
                    ibz_mul(&n_mu,&temp_rand,&ibz_const_two);
                    ibz_mul(&lambda,&lambda,&ibz_const_two); 

                    // replace C by CRT of C and C_start 
                    ibz_crt(&C[0],&C[0],&C_start[0],&lideal->norm,&lideal_start->norm);
                    ibz_crt(&C[1],&C[1],&C_start[1],&lideal->norm,&lideal_start->norm);

                    // computing the element quat_C = j* (C[0] + i C[1]); 
                    ibz_copy(&coeffs[0],&ibz_const_zero);
                    ibz_copy(&coeffs[1],&ibz_const_zero);
                    ibz_copy(&coeffs[2],&C[0]);
                    ibz_copy(&coeffs[3],&C[1]);
                    order_elem_create(&quat_C,&STANDARD_EXTREMAL_ORDER,&coeffs,Bpoo);

                    // computing the norm n_quat_C of quat_C 
                    quat_alg_norm(&ibq_n_quat_C,&quat_C,Bpoo);
                    found = ibq_to_ibz(&n_quat_C,&ibq_n_quat_C);
                    if (!found) {
                        assert(found);
                    }

                    // preparation of the strong approx computation
                    // setting param 
                    ibz_copy(&param.target_norm,&n_mu);
                    ibz_copy(&param.n,&n_prod);
                    ibz_copy(&param.equiv_n,&lideal->norm);
                    param.gamma = gamma;
                    param.delta = *delta;

                    // computing the strong approximation
                    found = strong_approximation(&mu,&n_mu,&STANDARD_EXTREMAL_ORDER,&C,&n_quat_C,&lambda,&n_prod,KLPT_signing_number_strong_approx,condition_signing_klpt,&param,Bpoo); 
                    
                    if (found) {
                        

                        

                        quat_alg_norm(&ibq_n_quat_C,&mu,Bpoo);
                        assert( ibq_to_ibz(&n_quat_C,&ibq_n_quat_C));
                        ibz_mul(&n_quat_C,&n_quat_C,&ibz_const_two);
                        ibz_mul(&n_quat_C,&n_quat_C,&ibz_const_two);
                        assert(0==ibz_cmp(&n_quat_C,&n_mu));

                        // computation of the final output gen = gamma * mu
                        quat_alg_mul(gen,&gamma,&mu,Bpoo);

                    } 
                    

                }
            }
        }
    
    }

    // var finalize
    ibz_finalize(&n_gamma);
    ibz_finalize(&temp_rand);
    ibz_finalize(&n_prod);
    ibz_finalize(&lambda);
    ibz_finalize(&s1);ibz_finalize(&s2);
    ibz_finalize(&n_quat_C);ibz_finalize(&n_quat_C_start);
    ibq_finalize(&ibq_n_quat_C);ibq_finalize(&ibq_n_quat_C_start);
    ibz_vec_2_finalize(&C);
    ibz_vec_2_finalize(&C_start);
    quat_alg_elem_finalize(&gamma);
    quat_alg_elem_finalize(&quat_one);
    quat_alg_elem_finalize(&mu);
    ibz_finalize(&n_mu);
    quat_alg_elem_finalize(&quat_C);
    quat_alg_elem_finalize(&quat_C_start);
    quat_alg_coord_finalize(&coeffs);
    ibz_finalize(&param.target_norm);
    ibz_finalize(&param.equiv_n);
    ibz_finalize(&param.n);
    quat_alg_finalize(&param.Bpoo);

    return found;


}


/**
 * @brief Equivalent left O0-ideal of power of two norm
 *
 * @param gen Output: generator of equiv  
 * @param lideal left O0 ideal
 * @param delta quaternion algebra element
 * @param Bpoo the quaternion algebra
 *  
 * 
 * equiv_lideal is equivalent to the ideal lideal and we have the equality
 * equiv_lideal = lideal * gen/Norm(lideal) of norm n = 2^KLPT_keygen_length
 * where conjugate(gen) is in lideal 
 * Assumes that the ideals lideal has a "good" norm (ie that one needs not to apply lideal_equiv)
 * returns a bit indicating if the computation has succeeded
 */
int klpt_keygen_klpt(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo) {
    // var dec
    int found;
    int cnt;
    int e0,e1,center;
    ibz_t n_gamma,temp_rand,n_prod,n_mu;
    ibq_t ibq_n_quat_C;
    ibz_t lambda,n_quat_C;
    quat_alg_elem_t gamma,quat_one,mu;
    quat_alg_elem_t quat_C;
    ibz_vec_2_t C;
    quat_alg_coord_t coeffs;
    signing_klpt_param_t param;

    // var init
    found = 0;
    cnt = 0;
    ibz_init(&n_gamma);
    ibz_init(&n_mu);
    ibz_init(&temp_rand);
    ibz_init(&n_prod);
    ibz_init(&lambda);
    ibz_init(&n_quat_C);
    ibq_init(&ibq_n_quat_C);
    ibz_vec_2_init(&C);
    quat_alg_elem_init(&mu);
    quat_alg_elem_init(&gamma);
    quat_alg_elem_init(&quat_one);
    quat_alg_elem_init(&quat_C);
    quat_alg_coord_init(&coeffs);
    ibz_init(&param.n);
    ibz_init(&param.target_norm);
    quat_alg_init_set(&param.Bpoo,&Bpoo->p);
    param.order = &STANDARD_EXTREMAL_ORDER;

    // quaternion element for 1 
    quat_alg_scalar(&quat_one,&ibz_const_one,&ibz_const_one);

    // initializing the center of the sample interval as log(p) - log( n(lideal) ) + KLPT_gamma_exponent_center_shift
    center = ibz_bitsize(&Bpoo->p) - ibz_bitsize(&lideal->norm) + KLPT_gamma_exponent_center_shift;

    
    // the main loop, iterate until a solution is found or the bound is reached  
    while (!found && cnt < KLPT_keygen_num_gamma_trial) {
        // incrementing the counter
        cnt++;
        // we start by choosing the exponent of gamma at random inside some interval (depending on the choice of parameters this interval might actually be of size 1)
        // sampling the exponent
        ibz_rand_interval_i(&temp_rand,center - KLPT_gamma_exponent_interval_size, center + KLPT_gamma_exponent_interval_size); 

        // computation of the norm
        ibz_pow(&n_gamma,&ibz_const_two,ibz_get(&temp_rand));
        ibz_mul(&n_gamma,&n_gamma,&lideal->norm);

        // computing the value of gamma
        found = represent_integer(&gamma,&n_gamma,Bpoo); 

        if (found) {
            // computing the norm of the remaining part
            ibz_pow(&temp_rand,&ibz_const_two,KLPT_keygen_length);
            ibz_mul(&temp_rand,&lideal->norm,&temp_rand);
            ibz_div(&n_mu,&temp_rand,&temp_rand,&n_gamma);

            // computing the linear combination mod n(lideal) 

            found = solve_combi_eichler(&C,&STANDARD_EXTREMAL_ORDER,&gamma,&quat_one,lideal,Bpoo,1); 
            if (!found) {
                assert(found);
            }
            // computing the element quat_C = j* (C[0] + i C[1]); 
            ibz_copy(&coeffs[0],&ibz_const_zero);
            ibz_copy(&coeffs[1],&ibz_const_zero);
            ibz_copy(&coeffs[2],&C[0]);
            ibz_copy(&coeffs[3],&C[1]);
            order_elem_create(&quat_C,&STANDARD_EXTREMAL_ORDER,&coeffs,Bpoo);

            // computing the norm n_quat_C of quat_C 
            quat_alg_norm(&ibq_n_quat_C,&quat_C,Bpoo);
            found = ibq_to_ibz(&n_quat_C,&ibq_n_quat_C);
            if (!found) {
                assert(found);
            }

            // check that the quadratic reduosity condition is verified 
            // n_mu / n_quat_C must be square mod lideal.norm
            ibz_invmod(&temp_rand,&n_quat_C,&lideal->norm); 
            ibz_mul(&temp_rand,&n_mu,&temp_rand);
            found = ibz_sqrt_mod_p(&lambda,&temp_rand,&lideal->norm);
            if (found) {


                // multiply lambda by 2 and the final norm by 4
                ibz_mul(&temp_rand,&n_mu,&ibz_const_two);
                ibz_mul(&n_mu,&temp_rand,&ibz_const_two);
                ibz_mul(&lambda,&lambda,&ibz_const_two);

                // preparation of the strong approx computation
                // setting param 
                ibz_copy(&param.target_norm,&n_mu);
                ibz_copy(&param.n,&lideal->norm);
                param.gamma = gamma;

                // computing the strong approximation
                found = strong_approximation(&mu,&n_mu,&STANDARD_EXTREMAL_ORDER,&C,&n_quat_C,&lambda,&lideal->norm,KLPT_signing_number_strong_approx,condition_keygen_klpt,&param,Bpoo); 
                    
                if (found) {
                        
                    quat_alg_norm(&ibq_n_quat_C,&mu,Bpoo);
                    assert( ibq_to_ibz(&n_quat_C,&ibq_n_quat_C));
                    ibz_mul(&n_quat_C,&n_quat_C,&ibz_const_two);
                    ibz_mul(&n_quat_C,&n_quat_C,&ibz_const_two);
                    assert(0==ibz_cmp(&n_quat_C,&n_mu));

                    // computation of the final output gen = gamma * mu
                    quat_alg_mul(gen,&gamma,&mu,Bpoo);
                } 

            }
        }
    }

    // var finalize
    ibz_finalize(&n_gamma);
    ibz_finalize(&temp_rand);
    ibz_finalize(&n_prod);
    ibz_finalize(&lambda);
    ibz_finalize(&n_quat_C);
    ibq_finalize(&ibq_n_quat_C);
    ibz_vec_2_finalize(&C);
    quat_alg_elem_finalize(&gamma);
    quat_alg_elem_finalize(&quat_one);
    quat_alg_elem_finalize(&mu);
    ibz_finalize(&n_mu);
    quat_alg_elem_finalize(&quat_C);
    quat_alg_coord_finalize(&coeffs);

    ibz_finalize(&param.n);
    ibz_finalize(&param.target_norm);
    quat_alg_finalize(&param.Bpoo);

    return found;

}
