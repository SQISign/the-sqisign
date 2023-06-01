#include <klpt.h>
#include "tools.h"

/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief the eichler norm equation implementation 
 */


/**
 * @brief Deciding if a given vector is suitable for the eichler norm algorithm
 *
 * @param mu Output: a quaternion algebra element
 * @param C the vector to be tested 
 * @param params parameters
 *  
 */
int condition_eichler( quat_alg_elem_t *mu, const ibz_vec_2_t *C, const void* params) {
            
            // we start by casting myparams to the correct type
            eichler_norm_param_t *eichler_norm_param = (eichler_norm_param_t*)params;
            // var dec 
            int found = 0;
            quat_alg_elem_t quat_temp;
            quat_alg_coord_t coeffs;
            ibz_t norm,rhs;
            ibz_t temp;
            ibq_t ibq_norm;
            ibz_t a,b;
            ibz_t ibz_q;
            quat_left_ideal_t id1,id2;

            // var init 
            quat_alg_elem_init(&quat_temp);
            quat_alg_coord_init(&coeffs);
            ibz_init(&norm);ibz_init(&rhs);
            ibz_init(&temp);
            ibq_init(&ibq_norm);
            ibz_init(&a);ibz_init(&b);
            ibz_init(&ibz_q);
            quat_left_ideal_init(&id1);
            quat_left_ideal_init(&id2);

            ibz_set(&ibz_q,eichler_norm_param->order.q);

            // mu = j* ( C[0] + i * C[1])
            ibz_copy(&coeffs[0],&ibz_const_zero);
            ibz_copy(&coeffs[1],&ibz_const_zero);
            ibz_copy(&coeffs[2],&(*C)[0]);
            ibz_copy(&coeffs[3],&(*C)[1]);


            order_elem_create(mu,&eichler_norm_param->order,&coeffs,&eichler_norm_param->Bpoo);
            
            // for debug check that mu is in the correct order
            assert(quat_lattice_contains(&coeffs,&eichler_norm_param->right_order,mu,&eichler_norm_param->Bpoo));

            // norm = n(mu) 
            quat_alg_norm(&ibq_norm,mu,&eichler_norm_param->Bpoo);
            if (!ibq_to_ibz(&norm,&ibq_norm)) {
                assert(0);
            }
            
            #ifndef NDEBUG 
                ibz_mul(&rhs,&(*C)[0],&(*C)[0]);
                ibz_set(&temp,eichler_norm_param->order.q);
                ibz_mul(&temp,&temp,&(*C)[1]);
                ibz_mul(&temp,&temp,&(*C)[1]);
                ibz_add(&temp,&temp,&rhs);
                ibz_mul(&temp,&temp,&eichler_norm_param->Bpoo.p);
                assert(0==ibz_cmp(&temp,&norm));
            #endif

            // rhs = (target_norm - norm) / n² 
            ibz_sub(&rhs,&eichler_norm_param->target_norm,&norm);
            ibz_mul(&temp,&eichler_norm_param->n,&eichler_norm_param->n);
            ibz_div(&rhs,&temp,&rhs,&temp);

            

            // for debug check that the remainder is zero
            assert(ibz_cmp(&temp,&ibz_const_zero)==0);

            // covers the case where we adjusted the norm of the target by 4, to help obtaining a primitive element
            // will always happen when eichler_nom_param.order.q = 1 mod 4
            if (ibz_get(&eichler_norm_param->target_norm)%2==0) {
                
                // when eichler_nom_param.order.q = 3 mod 4, we make some adjustements to help obtaining a primitive element and 
                // make possible the resolution of the binary quadratic equation below 
                // in this first case we will find an element of the form 4 * M for some prime M 
                if (eichler_norm_param->order.q%8 == 3 && ibz_get(&rhs)%4 == 0) {
                    ibz_set(&temp,4);
                    ibz_div(&norm,&temp,&rhs,&temp);
                    assert(0==ibz_cmp(&temp,&ibz_const_zero));
                } 
                // in this first case we will find an element of the form 8 * M for some prime M 
                else if ((eichler_norm_param->order.q%8 == 7 && ibz_get(&rhs)%8 == 0)) {
                    ibz_set(&temp,8);
                    ibz_div(&norm,&temp,&rhs,&temp);
                    assert(0==ibz_cmp(&temp,&ibz_const_zero));
                }
                else {
                    ibz_copy(&norm,&rhs);
                }
                assert(ibz_cmp(&norm,&ibz_const_zero)>=0);

                // TODUPDATE we can probably filter out some cases depending on the reduosiy 
                if ((eichler_norm_param->order.q == 1 && ibz_cornacchia_extended(&a, &b,&norm, SMALL_PRIMES_1MOD4 , sizeof(SMALL_PRIMES_1MOD4)/sizeof(*SMALL_PRIMES_1MOD4), KLPT_primality_num_iter, &PROD_SMALL_PRIMES_3MOD4)) 
                    || (eichler_norm_param->order.q%4 == 1 && ibz_probab_prime(&norm,KLPT_primality_num_iter) && ibz_cornacchia_prime(&a,&b, &ibz_q, &norm))
                    || (eichler_norm_param->order.q%8 == 3 && ibz_get(&rhs)%4 == 0 && ibz_probab_prime(&norm,KLPT_primality_num_iter) && ibz_cornacchia_special_prime(&a,&b,&ibz_q,&norm,2))
                    || (eichler_norm_param->order.q%8 == 7 && ibz_get(&rhs)%8 ==0 && ibz_probab_prime(&norm,KLPT_primality_num_iter) && ibz_cornacchia_special_prime(&a,&b,&ibz_q,&norm,3))

                ){                   

                    ibz_copy(&coeffs[2],&ibz_const_zero);
                    ibz_copy(&coeffs[3],&ibz_const_zero);
                    ibz_copy(&coeffs[0],&a);
                    ibz_copy(&coeffs[1],&b);
                    ibz_mul(&coeffs[0],&coeffs[0],&eichler_norm_param->n);
                    ibz_mul(&coeffs[1],&coeffs[1],&eichler_norm_param->n);
                    
                    // computing the final value of mu
                    // quat_temp = n * (a + i*b)
                    order_elem_create(&quat_temp,&eichler_norm_param->order,&coeffs,&eichler_norm_param->Bpoo);
            
                    // mu = mu +  quat_temp
                    quat_alg_add(mu,mu,&quat_temp);


                    #ifndef NDEBUG 
                        // for debug we check check that the norm is indeed the target norm
                        quat_alg_norm(&ibq_norm,mu,&eichler_norm_param->Bpoo);
                        ibq_to_ibz(&norm,&ibq_norm);
                        assert(ibz_cmp(&norm,&eichler_norm_param->target_norm)==0);
                        assert(quat_lattice_contains(&coeffs,&eichler_norm_param->order.order,mu,&eichler_norm_param->Bpoo));
                    #endif
                    
                    // mu = mu / 2 and checking that the result in still in the order
                    ibz_mul(&mu->denom,&mu->denom,&ibz_const_two);
                    found = quat_lattice_contains(&coeffs,&eichler_norm_param->order.order,mu,&eichler_norm_param->Bpoo);

                }
            }

            // this second case is only possible when q = 3 mod 4 
            // here, we have not adjusted the target_norm, so it is odd  
            else {
                // same as before, we differentiate several cases to enable the computation and boost the probability of finding a primitive element
                // in particular, we will try to find a,b of the right norm such that a+ib/2 is in the maximal extremal order
                if (eichler_norm_param->order.q%8 == 3 && ibz_get(&rhs)%2 != 0) {
                    ibz_copy(&norm,&rhs);
                }
                else if ((eichler_norm_param->order.q%8 == 3 && ibz_get(&rhs)%4 == 0)) {
                    ibz_set(&temp,4);
                    ibz_div(&norm,&temp,&rhs,&temp);
                    assert(0==ibz_cmp(&temp,&ibz_const_zero));
                }
                else if ((eichler_norm_param->order.q%8 == 7 && ibz_get(&rhs)%4 == 2)) {
                    ibz_set(&temp,2);
                    ibz_div(&norm,&temp,&rhs,&temp);
                    assert(0==ibz_cmp(&temp,&ibz_const_zero));
                }

                // TODUPDATE we can probably filter out some cases depending on the reduosiy 
                if ((eichler_norm_param->order.q%8 == 3 && (ibz_get(&rhs)%2 != 0) && ibz_probab_prime(&norm,KLPT_primality_num_iter) && ibz_cornacchia_special_prime(&a,&b,&ibz_q,&norm,2))
                    || (eichler_norm_param->order.q%8 == 3 && ibz_get(&rhs)%4 ==0 && ibz_probab_prime(&norm,KLPT_primality_num_iter) && ibz_cornacchia_special_prime(&a,&b,&ibz_q,&norm,4))
                    || (eichler_norm_param->order.q%8 == 7 && ibz_get(&rhs)%4 ==2 && ibz_probab_prime(&norm,KLPT_primality_num_iter) && ibz_cornacchia_special_prime(&a,&b,&ibz_q,&norm,3))

                ) {
                    ibz_copy(&coeffs[2],&ibz_const_zero);
                    ibz_copy(&coeffs[3],&ibz_const_zero);
                    ibz_copy(&coeffs[0],&a);
                    ibz_copy(&coeffs[1],&b);
                    ibz_mul(&coeffs[0],&coeffs[0],&eichler_norm_param->n);
                    ibz_mul(&coeffs[1],&coeffs[1],&eichler_norm_param->n);
                    
                    // computing the final value of mu
                    // quat_temp = n * (a + i*b) / 2 
                    order_elem_create(&quat_temp,&eichler_norm_param->order,&coeffs,&eichler_norm_param->Bpoo);
                    ibz_mul(&quat_temp.denom,&quat_temp.denom,&ibz_const_two);


                    // debug testing that quat_temp is in the correct order
                    assert(quat_lattice_contains(&coeffs,&eichler_norm_param->order.order,&quat_temp,&eichler_norm_param->Bpoo));

                    // mu = mu +  quat_temp
                    quat_alg_add(mu,mu,&quat_temp);

                    // debug testing that mu is in the correct order
                    assert(quat_lattice_contains(&coeffs,&eichler_norm_param->order.order,mu,&eichler_norm_param->Bpoo));

                    #ifndef NDEBUG
                        // debug testing the norm
                        quat_alg_norm(&ibq_norm,mu,&eichler_norm_param->Bpoo);
                        ibq_to_ibz(&norm,&ibq_norm);
                        assert(ibz_cmp(&norm,&eichler_norm_param->target_norm)==0);
                    #endif

                    found = 1;

                }

            }
            if (found) {
                    
                    // for debug we check that the element is primitive and contained in the correct order 
                    assert(quat_alg_is_primitive(mu,&eichler_norm_param->order.order,&eichler_norm_param->Bpoo));
                    assert(quat_lattice_contains(&coeffs,&eichler_norm_param->right_order,mu,&eichler_norm_param->Bpoo));
                    assert(quat_lattice_contains(&coeffs,&eichler_norm_param->order.order,mu,&eichler_norm_param->Bpoo));
                    assert(quat_lattice_contains(&coeffs,&eichler_norm_param->right_order,&eichler_norm_param->gen_constraint,&eichler_norm_param->Bpoo));
                    assert(quat_alg_is_primitive(&eichler_norm_param->gen_constraint,&eichler_norm_param->right_order,&eichler_norm_param->Bpoo));

                    // mu is in the extremal order, we only to verify the final condition
                    // computing the first ideal to compare id1 = order < gen_constraint , 2 > 
                    quat_lideal_make_primitive_then_create(&id1,&eichler_norm_param->gen_constraint,&ibz_const_two,&eichler_norm_param->right_order,&eichler_norm_param->Bpoo);

                    // computing the generator of the second ideal quat_temp = gen_constraint * conjugate( mu ) 
                    quat_alg_conj(&quat_temp,mu);
                    quat_alg_mul(&quat_temp,&eichler_norm_param->gen_constraint,&quat_temp,&eichler_norm_param->Bpoo);
                    quat_alg_normalize(&quat_temp);

                    // computing id2 = order < quat_temp ,2>
                    quat_lideal_make_primitive_then_create(&id2,&quat_temp,&ibz_const_two,&eichler_norm_param->right_order,&eichler_norm_param->Bpoo);

                    // testing equality 
                    found = !quat_lideal_equals(&id1,&id2,&eichler_norm_param->Bpoo);

            }
            

            
            // var finalize 
            quat_alg_elem_finalize(&quat_temp);
            quat_alg_coord_finalize(&coeffs);
            ibz_finalize(&norm);ibz_finalize(&rhs);
            ibz_finalize(&temp);
            ibq_finalize(&ibq_norm);    
            ibz_finalize(&a);ibz_finalize(&b);
            ibz_finalize(&ibz_q);
            quat_left_ideal_finalize(&id1);
            quat_left_ideal_finalize(&id2);

            return found;
}  


/**
 * @brief Finding an endomorphism of smooth norm inside a fixed eichler order
 *
 * @param beta Output: the linear combination 
 * @param n_beta Output: the norm of bet<a
 * @param order the special extremal order used for the computation 
 * @param gen the generator of the order-ideal
 * @param n the norm of the ideal 
 * @param gen_constraint generator of the constraint (belong to the right order of the ideal generated by gen), may be modified during the computation (division by a small scalar)
 * @param Bpoo the quaternion algebra
 *  
 * Let J = order < gen,n> and O = OR(J) and ideal_constraint = O < gen_constraint,2>
 * Find beta of norm defined n_quat_C in OR(J) \ (ZZ + ideal_constraint) (by design the norm of beta is odd)
 * gen_constraint is assumed to be primitive!
 * assumes n is prime 
 * returns a bit indicating if the computation succeeded
 */
int eichler_special_norm_fixed(quat_alg_elem_t *beta,ibz_t *n_beta, const quat_p_extremal_maximal_order_t *order, const quat_alg_elem_t *gen, const ibz_t *n, quat_alg_elem_t *gen_constraint, const quat_alg_t *Bpoo) {

    // var dec 
    int cnt;
    ibz_t nj; // norm of j
    ibq_t temp_ibq;
    ibz_t temp,lambda;
    ibz_t n_quat_C;  
    int size_list;
    quat_left_ideal_t lideal;
    quat_alg_coord_t coeffs;
    quat_alg_elem_t quat_1,quat_C,scal,quat_temp,mu;
    ibz_vec_2_t C; // the linear combination
    int is_divisible = 0;
    int is_n_quat_C_square;
    int log_margin;
    int found;
    int abort = 0;
    eichler_norm_param_t param;
    
    // var init
    found = 1; 
    cnt = 0;
    ibz_init(&nj);
    ibz_init(&n_quat_C);
    ibz_t n_mu_list[KLPT_eichler_number_mu_norm]; 

    ibz_init(&temp);
    ibz_init(&lambda);
    ibq_init(&temp_ibq);
    quat_left_ideal_init(&lideal);
    quat_alg_elem_init(&quat_1);
    quat_alg_elem_init(&quat_C);
    quat_alg_elem_init(&scal);
    quat_alg_elem_init(&quat_temp);
    quat_alg_elem_init(&mu);
    quat_alg_coord_init(&coeffs);
    ibz_vec_2_init(&C);
    for (int ind=0;ind<KLPT_eichler_number_mu_norm;ind++){
        ibz_init(&n_mu_list[ind]);
    }
    
    // init param
    quat_alg_init_set(&param.Bpoo,&Bpoo->p);
    ibz_init(&param.n);
    ibz_init(&param.target_norm);
    quat_alg_elem_init(&param.gen_constraint);
    param.order = *order;
    quat_order_init(&param.right_order);    

    //the quaternion element 1
    quat_alg_scalar(&quat_1,&ibz_const_one,&ibz_const_one);
    
    // computation the ideal lideal
    quat_lideal_make_primitive_then_create(&lideal,gen,n,&order->order,Bpoo);

    // computation of the right order of lideal
    quat_lideal_right_order(&param.right_order,&lideal,Bpoo);
   

    // removing the potential scalar factors of gen_constraint
    // TODUPDATE optimize this to avoid doing it when unnecessary (using the trace ?)
    quat_alg_make_primitive(&coeffs,&temp,gen_constraint,&param.right_order,Bpoo);
    ibz_mul(&gen_constraint->denom,&gen_constraint->denom,&temp);
    quat_alg_normalize(gen_constraint);

    // copying to eichler norm param
    quat_alg_elem_copy(&param.gen_constraint,gen_constraint);

    // for debug, we check that gen_constraint is in the correct right order 
    assert(quat_lattice_contains(&coeffs,&param.right_order,gen_constraint,Bpoo));
    assert(quat_alg_is_primitive(gen_constraint,&param.right_order,Bpoo));
    
    // computation of the linear combination
    found = found && solve_combi_eichler(&C,order,&quat_1,&quat_1,&lideal,Bpoo,is_divisible);
    assert(found);

    // compute quat_C = j * (C[0] + i * C[1])
    ibz_copy(&coeffs[0],&ibz_const_zero);
    ibz_copy(&coeffs[1],&ibz_const_zero);
    ibz_copy(&coeffs[2],&C[0]);
    ibz_copy(&coeffs[3],&C[1]);
    order_elem_create(&quat_C,order,&coeffs,Bpoo);
    assert(quat_lattice_contains(&coeffs,&param.right_order,&quat_C,Bpoo));

    // computation of the norm of beta 
    quat_alg_norm(&temp_ibq,&quat_C,Bpoo);
    if (!ibq_to_ibz(&n_quat_C,&temp_ibq)) {
      assert(0);
    }
    
    
    ibz_mul(&lambda,&C[0],&C[0]);
    ibz_set(&temp,order->q);
    ibz_mul(&temp,&temp,&C[1]);
    ibz_mul(&temp,&temp,&C[1]);
    ibz_add(&temp,&temp,&lambda);
    ibz_mul(&temp,&temp,&Bpoo->p);
    assert(0==ibz_cmp(&temp,&n_quat_C));

    // computation of the legendre symbol
    is_n_quat_C_square = ibz_legendre(&n_quat_C,n);

    //computation of the log_margin = log (p) + 3 log(n) + KLPT_strong_approx_log_margin
    log_margin = ibz_bitsize(&Bpoo->p) + 3* ibz_bitsize(n) + KLPT_eichler_strong_approx_log_margin;

    // computation of a list of possible values for n_mu
    size_list = norm_list_computation(n_mu_list,KLPT_eichler_number_mu_norm,log_margin,n,is_n_quat_C_square);

    // if no suitable norm has been found we can abort directly 
    if (!size_list) {
        abort = 1;
    }

    if(!abort){
        found = 0;

        // strong_approximation loop
        while (!found && cnt < size_list) {

            // computation of the value of lambda
            ibz_invmod(&temp,&n_quat_C,n);
            ibz_mul(&temp,&temp,&n_mu_list[cnt]);
            ibz_mod(&temp,&temp,n);
            found=ibz_sqrt_mod_p(&lambda,&temp,n);
            assert(found);

            // if q = 1 mod 4 multiply lambda by 2 and the final norm by 4
            if (order->q%4 == 1) {
            // if (1) {
                ibz_mul(&temp,&n_mu_list[cnt],&ibz_const_two);
                ibz_mul(&temp,&temp,&ibz_const_two);
                ibz_mul(&lambda,&lambda,&ibz_const_two); 
            }
            else {
                ibz_copy(&temp,&n_mu_list[cnt]);
            }
            // preparation of the strong approx computation
            // setting param
            ibz_copy(&param.target_norm,&temp);
            ibz_copy(&param.n,n);


            // actual strong_approximation
            found = strong_approximation(&mu,&temp,order,&C,&n_quat_C,&lambda,n,KLPT_eichler_number_strong_approx,condition_eichler,&param,Bpoo);      
            if (!found && order->q%4==3) {
                // if q =3 mod 4, we try again by multiplying the norm by 4
                ibz_mul(&temp,&n_mu_list[cnt],&ibz_const_two);
                ibz_mul(&temp,&temp,&ibz_const_two);
                ibz_mul(&lambda,&lambda,&ibz_const_two); 
                ibz_copy(&param.target_norm,&temp);
                found = strong_approximation(&mu,&temp,order,&C,&n_quat_C,&lambda,n,KLPT_eichler_number_strong_approx,condition_eichler,&param,Bpoo);      

            }

            if (found) {
                // computing the final beta = mu
                ibz_copy(&beta->coord[0],&mu.coord[0]);
                ibz_copy(&beta->coord[1],&mu.coord[1]);
                ibz_copy(&beta->coord[2],&mu.coord[2]);
                ibz_copy(&beta->coord[3],&mu.coord[3]);
                ibz_copy(&beta->denom,&mu.denom);
                ibz_copy(n_beta,&n_mu_list[cnt]);
            }
            else {
                // incrementing the counter
                cnt++;
            }
        }
    } else {
        found = 0;
    }

    // var finalize
    ibz_finalize(&nj);
    ibz_finalize(&n_quat_C);
    quat_left_ideal_finalize(&lideal);
    quat_alg_elem_finalize(&quat_1);
    quat_alg_elem_finalize(&quat_C);
    quat_alg_elem_finalize(&scal);
    quat_alg_elem_finalize(&quat_temp);
    ibz_finalize(&temp);
    ibz_finalize(&lambda);
    ibq_finalize(&temp_ibq);
    quat_alg_elem_finalize(&mu);
    quat_alg_coord_finalize(&coeffs);
    ibz_vec_2_finalize(&C);
    for (int ind=0;ind<KLPT_eichler_number_mu_norm;ind++){
        ibz_finalize(&n_mu_list[ind]);
    }
    quat_alg_finalize(&param.Bpoo);
    ibz_finalize(&param.n);
    ibz_finalize(&param.target_norm);
    quat_order_finalize(&param.right_order);
    quat_alg_elem_finalize(&param.gen_constraint);

    
    return found;

}


/**
 * @brief Finding an endomorphism of smooth norm inside an eichler order
 *
 * @param beta Output: the endomorphism 
 * @param n_beta Output: the norm of beta
 * @param gen Output: the generator of equivalent ideal 
 * @param n Output: the norm of the equivalent ideal
 * @param lideal left O0 ideal
 * @param gen_constraint generator of a left OR(lideal)-ideal, is contained inside OR(lideal) may be modified during the computation (division by some scalar)
 * @param Bpoo the quaternion algebra
 *  
 * Let I = lideal and K = O_R(I) < gen_constraint,2> 
 * Find beta of norm n_beta in OR(J) \ O (by design the norm of beta is odd)
 * where J = I * gen /Norm(I) of norm n, a prime number  
 * and O is the embedding of ZZ + K inside OR(J) (induced by the isomorphism from OR(I) to OR(J))
 * returns a bit indicating if the computation succeeded
 * in most cases beta will be contained in the order ZZ + J, but in some extreme cases (where we use another special extremal order) it might not be true. This does not change anything for the applications of eichler_special_norm. 
 */
int klpt_eichler_special_norm(quat_alg_elem_t *beta, ibz_t *n_beta, quat_alg_elem_t *gen, ibz_t *n, quat_left_ideal_t *lideal, quat_alg_elem_t *gen_constraint, const quat_alg_t *Bpoo) {

    // var declaration
    int bitsize;
    ibz_t mod2;
    int found,found_equiv;
    int cnt;
    quat_alg_elem_t gen_lideal;
    ibz_t n_lideal;
    quat_left_ideal_t connecting_lideal;
    quat_order_t right_order;
    quat_alg_elem_t new_gen_constraint,quat_temp;
    const quat_p_extremal_maximal_order_t *alternate_order = NULL; 
    ibz_mat_4x4_t reduced,gram;
    quat_alg_coord_t coeffs;
    
    // var init
    found = 0; found_equiv = 0;
    cnt = 0;
    ibz_init(&mod2);
    ibz_init(&n_lideal);
    quat_alg_elem_init(&gen_lideal);
    quat_alg_elem_init(&new_gen_constraint);
    quat_alg_elem_init(&quat_temp);
    quat_left_ideal_init(&connecting_lideal);
    quat_order_init(&right_order);
    ibz_mat_4x4_init(&reduced);
    ibz_mat_4x4_init(&gram);
    quat_alg_coord_init(&coeffs);

    // we start by deciding if we can use J = I 
    // for that we try to see if the value of n(I) is small and odd
    bitsize = ibz_bitsize(&lideal->norm);
    if (ibz_get(&lideal->norm)%2 == 1 && bitsize < KLPT_eichler_smallnorm_bitsize ) { 
        // in that case, we should have enough room to find an output with J = I
        // compute the generator and copy the norm 
        int lideal_generator_ok = quat_lideal_generator(gen,lideal,Bpoo,0);
        assert(lideal_generator_ok);
        ibz_copy(n,&lideal->norm);
        found = eichler_special_norm_fixed(beta,n_beta,&STANDARD_EXTREMAL_ORDER,gen,n,gen_constraint,Bpoo);
        
        // we adjust the value of gen to match the requirement J = I * gen /Norm(I) of norm n with J = I
        quat_alg_scalar(gen,n,&ibz_const_one);
        found_equiv = 1;

        // if it failed we used alternate orders
        if (!found ) {

            cnt = 0;
            // we compute the right order of lideal
            quat_lideal_right_order(&right_order,lideal,Bpoo);

            while ( !found && cnt < KLPT_eichler_num_equiv_ideal*NUM_ALTERNATE_EXTREMAL_ORDERS ) {

                if (cnt%KLPT_eichler_num_equiv_ideal == 0) {
                    // we initialize the alternate order and the other values
                    int index = cnt/KLPT_eichler_num_equiv_ideal;

                    alternate_order = &ALTERNATE_EXTREMAL_ORDERS[index];

                    // fist we need to compute connecting_lideal = connecting_ideal ( right_order, alternate_order )
                    quat_connecting_ideal(&connecting_lideal,&ALTERNATE_EXTREMAL_ORDERS[index].order,&right_order,Bpoo);

                    // we reduce the basis and compute its gram matrix
                    // basis reduction step
                    quat_lideal_reduce_basis(&reduced,&gram,&connecting_lideal,Bpoo);

                }
                cnt++;
            
                // we compute an equivalent ideal
                found = klpt_lideal_equiv(&gen_lideal,&n_lideal,&reduced,&gram,&connecting_lideal.norm,&connecting_lideal.lattice.denom,Bpoo);
                if (!found) {
                    continue;
                }

                // filtering if the bitsize of n is too big
                if (ibz_bitsize(&CHARACTERISTIC)+ 3*ibz_bitsize(&n_lideal) > 2*ibz_bitsize(&TORSION_ODD)) {
                    found = 0;
                    continue;

                }
                // we adapt the value of gen_constraint to be contained inside the right order of equiv_lideal
                // new_gen_constraint =   ( conj(gen_lideal) / n(gen_lideal) ) *  new_gen_constraint * gen_lideal
                quat_alg_conj(&quat_temp,&gen_lideal);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,&n_lideal);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,&connecting_lideal.norm);
                quat_alg_mul(&new_gen_constraint,&quat_temp,gen_constraint,Bpoo);
                quat_alg_mul(&new_gen_constraint,&new_gen_constraint,&gen_lideal,Bpoo);

                // we try to solve the eichler norm equation 
                found = eichler_special_norm_fixed(beta,n_beta,alternate_order,&gen_lideal,&n_lideal,&new_gen_constraint,Bpoo);

                assert(!found || quat_lattice_contains(&coeffs,&alternate_order->order,beta,Bpoo));

            }


            if (found ) {
                // we apply the isomorphism to embed beta into the correct order
                
                // first, we send beta to right_order
                // beta = gen_lideal * beta * conj (gen_lideal) / n( gen_lideal )  
                quat_alg_mul(beta,beta,&quat_temp,Bpoo);
                quat_alg_mul(beta,&gen_lideal,beta,Bpoo);
                assert(quat_lattice_contains(&coeffs,&right_order,beta,Bpoo));

                // second, we send beta to the right order of the equivalent ideal generated by gen
                // beta = conj ( gen ) / n ( gen ) * beta * gen 
                quat_alg_conj(&quat_temp,gen);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,n);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal->norm);
                quat_alg_mul(beta,&quat_temp,beta,Bpoo);
                quat_alg_mul(beta,beta,gen,Bpoo);

            }
        }
    }
    else {
        // in all other cases we try to run the eichler order norm equation for the standard extremal order first
        // before trying alternates order

        // we compute a reduced basis and its gram matrix
        quat_lideal_reduce_basis(&reduced,&gram,lideal,Bpoo);

        // we will try different ideals equivalent to lideal 
        while ( !found && cnt < KLPT_eichler_num_equiv_ideal ) {
           
            cnt++;
            // we start by computing an equivalent ideal.
            // equiv_lideal =   lideal * gen/Norm(lideal)
            // gen is contained in equiv_lideal and conj(lideal)

            found_equiv = klpt_lideal_equiv(gen,n,&reduced,&gram,&lideal->norm,&lideal->lattice.denom,Bpoo);
            if (!found_equiv) {
                found = 0;
                continue;
            }
            
            // filtering if the bitsize of n is too big
            if (ibz_bitsize(&CHARACTERISTIC)+ 3*ibz_bitsize(n) > 2*ibz_bitsize(&TORSION_ODD)) { 
                found = 0;
                continue;
            }

            // we adapt the value of gen_constraint to be contained inside the right order of equiv_lideal
            // new_gen_constraint =   ( conj(gen) / n(gen) ) *  gen_constraint * gen
            quat_alg_conj(&new_gen_constraint,gen);
            ibz_mul(&new_gen_constraint.denom,&new_gen_constraint.denom,n);
            ibz_mul(&new_gen_constraint.denom,&new_gen_constraint.denom,&lideal->norm);
            quat_alg_mul(&new_gen_constraint,&new_gen_constraint,gen_constraint,Bpoo);
            quat_alg_mul(&new_gen_constraint,&new_gen_constraint,gen,Bpoo);


            // we try to solve the eichler norm equation 
            found = eichler_special_norm_fixed(beta,n_beta,&STANDARD_EXTREMAL_ORDER,gen,n,&new_gen_constraint,Bpoo);

        }

        // then we try the same with but with an alternate special extremal orders 
        // we will need to find a connecting ideal and then try the same procedure with ideals equivalent to the connecting ideal
        if (!found && found_equiv) {

            cnt = 0;
            // we compute the right order of lideal
            quat_lideal_right_order(&right_order,lideal,Bpoo);

            while ( !found && cnt < KLPT_eichler_num_equiv_ideal*NUM_ALTERNATE_EXTREMAL_ORDERS ) {

                if (cnt%KLPT_eichler_num_equiv_ideal == 0) {
                    // we initialize the alternate order and the other values
                    int index = cnt/KLPT_eichler_num_equiv_ideal;

                    alternate_order = &ALTERNATE_EXTREMAL_ORDERS[index];

                    // fist we need to compute connecting_lideal = connecting_ideal ( right_order, alternate_order )
                    quat_connecting_ideal(&connecting_lideal,&ALTERNATE_EXTREMAL_ORDERS[index].order,&right_order,Bpoo);

                    // we reduce the basis and compute its gram matrix
                    // basis reduction step
                    quat_lideal_reduce_basis(&reduced,&gram,&connecting_lideal,Bpoo);

                }
                cnt++;
            
                // we compute an equivalent ideal
                found = klpt_lideal_equiv(&gen_lideal,&n_lideal,&reduced,&gram,&connecting_lideal.norm,&connecting_lideal.lattice.denom,Bpoo);
                if (!found) {
                    continue;
                }

                // filtering if the bitsize of n is too big
                if ( ibz_bitsize(&CHARACTERISTIC)+ 3*ibz_bitsize(&n_lideal) > 2*ibz_bitsize(&TORSION_ODD)) {
                    found = 0;
                    continue;

                }

                // we adapt the value of gen_constraint to be contained inside the right order of equiv_lideal
                // new_gen_constraint =   ( conj(gen_lideal) / n(gen_lideal) ) *  new_gen_constraint * gen_lideal
                quat_alg_conj(&quat_temp,&gen_lideal);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,&n_lideal);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,&connecting_lideal.norm);
                quat_alg_mul(&new_gen_constraint,&quat_temp,gen_constraint,Bpoo);
                quat_alg_mul(&new_gen_constraint,&new_gen_constraint,&gen_lideal,Bpoo);

                // we try to solve the eichler norm equation 
                found = eichler_special_norm_fixed(beta,n_beta,alternate_order,&gen_lideal,&n_lideal,&new_gen_constraint,Bpoo);

                assert(!found || quat_lattice_contains(&coeffs,&alternate_order->order,beta,Bpoo));

            }


            if (found  ) {
                // we apply the isomorphism to embed beta into the correct order
                
                // first, we send beta to right_order
                // beta = gen_lideal * beta * conj (gen_lideal) / n( gen_lideal )  
                quat_alg_mul(beta,beta,&quat_temp,Bpoo);
                quat_alg_mul(beta,&gen_lideal,beta,Bpoo);
                assert(quat_lattice_contains(&coeffs,&right_order,beta,Bpoo));

                // second, we send beta to the right order of the equivalent ideal generated by gen
                // beta = conj ( gen ) / n ( gen ) * beta * gen 
                quat_alg_conj(&quat_temp,gen);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,n);
                ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal->norm);
                quat_alg_mul(beta,&quat_temp,beta,Bpoo);
                quat_alg_mul(beta,beta,gen,Bpoo);

            }

        }
    }


    // var finalization
    ibz_finalize(&mod2);
    ibz_finalize(&n_lideal);
    quat_alg_elem_finalize(&gen_lideal);
    quat_alg_elem_finalize(&quat_temp);
    quat_alg_elem_finalize(&new_gen_constraint);
    quat_left_ideal_finalize(&connecting_lideal);
    quat_order_finalize(&right_order);
    ibz_mat_4x4_finalize(&reduced);
    ibz_mat_4x4_finalize(&gram);
    quat_alg_coord_finalize(&coeffs);
    return found & found_equiv;
} 
