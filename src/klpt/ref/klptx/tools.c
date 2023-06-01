#include<quaternion.h>
#include <klpt.h>
#include "tools.h"



/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief the eichler norm equation implementation 
 */


void quat_alg_elem_copy(quat_alg_elem_t *copy,const quat_alg_elem_t *copied) {
    ibz_copy(&copy->denom,&copied->denom);
    ibz_copy(&copy->coord[0],&copied->coord[0]);
    ibz_copy(&copy->coord[1],&copied->coord[1]);
    ibz_copy(&copy->coord[2],&copied->coord[2]);
    ibz_copy(&copy->coord[3],&copied->coord[3]);
}
void quat_left_ideal_copy(quat_left_ideal_t *copy,const quat_left_ideal_t *copied) {
    copy->parent_order = copied->parent_order;
    ibz_copy(&copy->norm,&copied->norm);
    ibz_copy(&copy->lattice.denom,&copied->lattice.denom);
    for (int i=0;i<4;i++) {
        for (int j=0;j<4;j++) {
                ibz_copy(&copy->lattice.basis[i][j],&copied->lattice.basis[i][j]); 
        }
    }
}


/**
 * @brief Create an element of a extremal maximal order from its coefficients
 *
 * @param elem Output: the quaternion element
 * @param order the order
 * @param coeffs the vector of 4 ibz coefficients
 * @param Bpoo quaternion algebra
 *  
 * elem = x + i*y + j*z + j*i*t  
 * where coeffs = [x,y,z,t] and i = order.i, j = order.j
 *  
   */
void order_elem_create(quat_alg_elem_t *elem, const quat_p_extremal_maximal_order_t *order, const quat_alg_coord_t *coeffs, const quat_alg_t *Bpoo) {

    // var dec 
    quat_alg_elem_t quat_temp;

    // var init
    quat_alg_elem_init(&quat_temp);

    // elem = x 
    quat_alg_scalar(elem,&(*coeffs)[0],&ibz_const_one);

    // quat_temp = i*y
    quat_alg_scalar(&quat_temp,&((*coeffs)[1]),&ibz_const_one);
    quat_alg_mul(&quat_temp,&order->i,&quat_temp,Bpoo);
    
    // elem = x + i*y
    quat_alg_add(elem,elem,&quat_temp);
    
    // quat_temp = z * j
    quat_alg_scalar(&quat_temp,&(*coeffs)[2],&ibz_const_one);
    quat_alg_mul(&quat_temp,&order->j,&quat_temp,Bpoo);

    // elem = x + i* + z*j
    quat_alg_add(elem,elem,&quat_temp);
    
    // quat_temp = t * j * i
    quat_alg_scalar(&quat_temp,&(*coeffs)[3],&ibz_const_one);
    quat_alg_mul(&quat_temp,&order->j,&quat_temp,Bpoo);
    quat_alg_mul(&quat_temp,&quat_temp,&order->i,Bpoo);

    // elem =  x + i*y + j*z + j*i*t  
    quat_alg_add(elem,elem,&quat_temp);

    quat_alg_elem_finalize(&quat_temp);
}   

/**
 * @brief Finding the good linear combination
 *
 * @param C Output: the linear combination
 * @param beta the endomorphism 
 * @param order the order
 * @param n the norm
 * @param e the exponent
 * @param gen_start element of order, generate a left O-ideal of norm n = 2^exp
 * @param gen_end element of order, generate a left O-ideal of norm n = 2^exp
 * @param Bpoo quaternion algebra
 *  
 * lideal_start = order  < gen_start,n>
 * lideal_end = order  < gen_end,n>
 * This algorithm finds two integers C[1], C[2] such that the pushforward of the ideal lideal_start by the endomorphism C[1] + C[2] beta is equal to lideal_end
 * beta is an endomorphism of order, it has been chosen to ensure that a solution exists
 * Returns a bit indicating if the computation succeeded 
   */
int klpt_find_linear_comb(ibz_vec_2_t *C,const quat_alg_elem_t *beta, const quat_order_t *order, const ibz_t *n, const unsigned short exp, const quat_alg_elem_t *gen_start, const quat_alg_elem_t *gen_end,const quat_alg_t *Bpoo) {

    // var dec
    int found;
    quat_left_ideal_t lideal_end;
    quat_alg_elem_t temp,gamma,bas2,bas3;
    quat_alg_coord_t v0,v1,v2,v3;
    ibz_mat_4x4_t system;
     
    // var init
    found = 0;
    quat_alg_coord_init(&v1);quat_alg_coord_init(&v0);
    quat_alg_coord_init(&v2);quat_alg_coord_init(&v3);
    quat_alg_elem_init(&gamma);
    quat_alg_elem_init(&temp);
    quat_alg_elem_init(&bas2);
    quat_alg_elem_init(&bas3);
    quat_left_ideal_init(&lideal_end);
    ibz_mat_4x4_init(&system);

    assert(quat_lattice_contains(&v0,order,gen_end,Bpoo));

    //computing lideal_start from gen_start
    quat_lideal_create_from_primitive(&lideal_end,gen_end,n,order,Bpoo);


    // compute gamma = gen_start * \overline { beta } 
    quat_alg_conj(&temp,beta);
    quat_alg_mul(&gamma,gen_start,&temp,Bpoo);
    quat_alg_normalize(&gamma);
    

    // computing the system to solve
    // possibly solving the system over the basis of Bpoo should be enough ?? (up to rescaling to allow for integer solutions) 

    // converting the basis in actual quaternion element bas2,bas3
    ibz_copy(&bas2.denom,&lideal_end.lattice.denom);
    ibz_copy(&bas3.denom,&lideal_end.lattice.denom);
    ibz_copy(&bas2.coord[0],&lideal_end.lattice.basis[0][2]);
    ibz_copy(&bas2.coord[1],&lideal_end.lattice.basis[1][2]);
    ibz_copy(&bas2.coord[2],&lideal_end.lattice.basis[2][2]);
    ibz_copy(&bas2.coord[3],&lideal_end.lattice.basis[3][2]);
    ibz_copy(&bas3.coord[0],&lideal_end.lattice.basis[0][3]);
    ibz_copy(&bas3.coord[1],&lideal_end.lattice.basis[1][3]);
    ibz_copy(&bas3.coord[2],&lideal_end.lattice.basis[2][3]);
    ibz_copy(&bas3.coord[3],&lideal_end.lattice.basis[3][3]);

    

    // setting the four vectors needed for the system
    found = quat_lattice_contains(&v0,order,gen_start,Bpoo);
    assert(found);
    found = found && quat_lattice_contains(&v1,order,&gamma,Bpoo);
    assert(found);
    found = found && quat_lattice_contains(&v2,order,&bas2,Bpoo);
    assert(found);
    found = found && quat_lattice_contains(&v3,order,&bas3,Bpoo);
    assert(found);

    // creation of the actual system
    ibz_copy(&system[0][0],&v0[0]);
    ibz_copy(&system[0][1],&v1[0]);
    ibz_copy(&system[0][2],&v2[0]);
    ibz_copy(&system[0][3],&v3[0]);
    ibz_copy(&system[1][0],&v0[1]);
    ibz_copy(&system[1][1],&v1[1]);
    ibz_copy(&system[1][2],&v2[1]);
    ibz_copy(&system[1][3],&v3[1]);
    ibz_copy(&system[2][0],&v0[2]);
    ibz_copy(&system[2][1],&v1[2]);
    ibz_copy(&system[2][2],&v2[2]);
    ibz_copy(&system[2][3],&v3[2]);
    ibz_copy(&system[3][0],&v0[3]);
    ibz_copy(&system[3][1],&v1[3]);
    ibz_copy(&system[3][2],&v2[3]);
    ibz_copy(&system[3][3],&v3[3]);
    // ibz_copy(&system[0][0],&gen_start->coord[0]);
    // ibz_copy(&system[0][1],&gamma.coord[0]);
    // ibz_copy(&system[0][2],&bas2.coord[0]);
    // ibz_copy(&system[0][3],&bas3.coord[0]);
    // ibz_copy(&system[1][0],&gen_start->coord[1]);
    // ibz_copy(&system[1][1],&gamma.coord[1]);
    // ibz_copy(&system[1][2],&bas2.coord[1]);
    // ibz_copy(&system[1][3],&bas3.coord[1]);
    // ibz_copy(&system[2][0],&gen_start->coord[2]);
    // ibz_copy(&system[2][1],&gamma.coord[2]);
    // ibz_copy(&system[2][2],&bas2.coord[2]);
    // ibz_copy(&system[2][3],&bas3.coord[2]);
    // ibz_copy(&system[3][0],&gen_start->coord[3]);
    // ibz_copy(&system[3][1],&gamma.coord[3]);
    // ibz_copy(&system[3][2],&bas2.coord[3]);
    // ibz_copy(&system[3][3],&bas3.coord[3]);

    // finding the kernel of the system 
    found = ibz_4x4_right_ker_mod_power_of_2(&v0,&system,exp);

    // computation of the result
    ibz_copy(&(*C)[0],&v0[0]);
    ibz_copy(&(*C)[1],&v0[1]);
    // var fin
    quat_alg_coord_finalize(&v1);quat_alg_coord_finalize(&v0);
    quat_alg_coord_finalize(&v2);quat_alg_coord_finalize(&v3);
    quat_alg_elem_finalize(&gamma);
    quat_alg_elem_finalize(&temp);
    quat_alg_elem_finalize(&bas2);
    quat_alg_elem_finalize(&bas3);
    quat_left_ideal_finalize(&lideal_end);
    ibz_mat_4x4_finalize(&system);

    return found;
} 

/**
 * @brief Finding the good linear span(j,i*j) combination
 *
 * @param C Output: two integers modulo n  
 * @param order special extremal order 
 * @param delta quaternion element
 * @param gamma quaternion element
 * @param lideal generator of left order ideal
 * @param Bpoo quaternion algebra
 * @param is_divisible intean indicating if Norm(gamma) is divisible by n
 *  
 * This algorithm finds two integers C[0], C[1] such that gamma * j * (C[0] + i C[1]]) * delta is in the eichler order ZZ + lideal  
 * if gamma have norm divisible by n(lideal) (indicated by the intean is_divisible), then the result will be contained in J 
 * Failure should not be possible
 * Assumes n is prime !
 */
int solve_combi_eichler(ibz_vec_2_t *C, const quat_p_extremal_maximal_order_t *order, const quat_alg_elem_t *gamma, const quat_alg_elem_t *delta, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo, int is_divisible){


    //var dec  
    int found;
    quat_alg_elem_t gamma_delta, gamma_j,gamma_j_delta,gamma_ji_delta,denom1,denom2,denom3,bas2,bas3;

    //var init
    found = 0;
    quat_alg_elem_init(&gamma_j);
    quat_alg_elem_init(&gamma_delta);
    quat_alg_elem_init(&gamma_j_delta);
    quat_alg_elem_init(&gamma_ji_delta);
    quat_alg_elem_init(&denom1);
    quat_alg_elem_init(&denom2);
    quat_alg_elem_init(&denom3);
    quat_alg_elem_init(&bas3);
    quat_alg_elem_init(&bas2);

    // computation of the system to solve 

    //computation of the algebra elements we need
    // we multiply everything by gamma and delta
    quat_alg_mul(&gamma_j,gamma,&order->j,Bpoo);
    quat_alg_mul(&gamma_ji_delta,&gamma_j,&order->i,Bpoo);
    quat_alg_mul(&gamma_j_delta,&gamma_j,delta,Bpoo);
    quat_alg_mul(&gamma_ji_delta,&gamma_ji_delta,delta,Bpoo);
    quat_alg_normalize(&gamma_ji_delta);
    quat_alg_normalize(&gamma_j_delta);

    // we create the two elements that are interesting to us
    ibz_copy(&bas2.denom,&ibz_const_one);
    ibz_copy(&bas3.denom,&ibz_const_one);
    ibz_copy(&bas2.coord[0],&lideal->lattice.basis[0][2]);
    ibz_copy(&bas2.coord[1],&lideal->lattice.basis[1][2]);
    ibz_copy(&bas2.coord[2],&lideal->lattice.basis[2][2]);
    ibz_copy(&bas2.coord[3],&lideal->lattice.basis[3][2]);
    ibz_copy(&bas3.coord[0],&lideal->lattice.basis[0][3]);
    ibz_copy(&bas3.coord[1],&lideal->lattice.basis[1][3]);
    ibz_copy(&bas3.coord[2],&lideal->lattice.basis[2][3]);
    ibz_copy(&bas3.coord[3],&lideal->lattice.basis[3][3]);

    // we compute the denominators
    quat_alg_scalar(&denom1,&gamma_ji_delta.denom,&ibz_const_one); 
    quat_alg_scalar(&denom2,&lideal->lattice.denom,&ibz_const_one);
    quat_alg_scalar(&denom3,&gamma_j_delta.denom,&ibz_const_one);

    // performing the necessary multiplication by crossed denom
    quat_alg_mul(&gamma_j_delta,&gamma_j_delta,&denom2,Bpoo); quat_alg_mul(&gamma_j_delta,&gamma_j_delta,&denom1,Bpoo);
    quat_alg_mul(&gamma_ji_delta,&gamma_ji_delta,&denom2,Bpoo);quat_alg_mul(&gamma_ji_delta,&gamma_ji_delta,&denom3,Bpoo);
    quat_alg_mul(&bas2,&bas2,&denom1,Bpoo);quat_alg_mul(&bas2,&bas2,&denom3,Bpoo);
    quat_alg_mul(&bas3,&bas3,&denom1,Bpoo);quat_alg_mul(&bas3,&bas3,&denom3,Bpoo);

    if (is_divisible) {
        // declaration and init of the system and kernel
        ibz_mat_4x4_t system;
        ibz_mat_4x4_init(&system);
        quat_alg_coord_t ker;
        quat_alg_coord_init(&ker); 

        // we fill the system
        ibz_copy(&system[0][0],&gamma_j_delta.coord[0]);
        ibz_copy(&system[0][1],&gamma_ji_delta.coord[0]);
        ibz_copy(&system[0][2],&bas2.coord[0]);
        ibz_copy(&system[0][3],&bas3.coord[0]);
        ibz_copy(&system[1][0],&gamma_j_delta.coord[1]);
        ibz_copy(&system[1][1],&gamma_ji_delta.coord[1]);
        ibz_copy(&system[1][2],&bas2.coord[1]);
        ibz_copy(&system[1][3],&bas3.coord[1]);
        ibz_copy(&system[2][0],&gamma_j_delta.coord[2]);
        ibz_copy(&system[2][1],&gamma_ji_delta.coord[2]);
        ibz_copy(&system[2][2],&bas2.coord[2]);
        ibz_copy(&system[2][3],&bas3.coord[2]);
        ibz_copy(&system[3][0],&gamma_j_delta.coord[3]);
        ibz_copy(&system[3][1],&gamma_ji_delta.coord[3]);
        ibz_copy(&system[3][2],&bas2.coord[3]);
        ibz_copy(&system[3][3],&bas3.coord[3]);
    
        // ibz_mat_4x4_print(&system);

        // resolution of the system
        found = ibz_4x4_right_ker_mod_prime(&ker,&system,&lideal->norm);
        assert(found);

        // computation of the result
        ibz_copy(&(*C)[0],&ker[0]);
        ibz_copy(&(*C)[1],&ker[1]);

        // finalization of the kernel and system
        ibz_mat_4x4_finalize(&system);
        quat_alg_coord_finalize(&ker);

    }
    else {
        // somme missing operations
        // gamma_delta = 1
        quat_alg_scalar(&gamma_delta,&ibz_const_one,&ibz_const_one);

        // multiply by the denom
        quat_alg_mul(&gamma_delta,&gamma_delta,&denom1,Bpoo);
        quat_alg_mul(&gamma_delta,&gamma_delta,&denom2,Bpoo);
        quat_alg_mul(&gamma_delta,&gamma_delta,&denom3,Bpoo);

        // declaration and init of the system and 
        ibz_mat_4x5_t system;
        ibz_mat_4x5_init(&system);
        ibz_vec_5_t ker;
        ibz_vec_5_init(&ker);

        // we fill the system
        ibz_copy(&system[0][0],&gamma_j_delta.coord[0]);
        ibz_copy(&system[0][1],&gamma_ji_delta.coord[0]);
        ibz_copy(&system[0][2],&gamma_delta.coord[0]);
        ibz_copy(&system[0][3],&bas2.coord[0]);
        ibz_copy(&system[0][4],&bas3.coord[0]);
        ibz_copy(&system[1][0],&gamma_j_delta.coord[1]);
        ibz_copy(&system[1][1],&gamma_ji_delta.coord[1]);
        ibz_copy(&system[1][2],&gamma_delta.coord[1]);
        ibz_copy(&system[1][3],&bas2.coord[1]);
        ibz_copy(&system[1][4],&bas3.coord[1]);
        ibz_copy(&system[2][0],&gamma_j_delta.coord[2]);
        ibz_copy(&system[2][1],&gamma_ji_delta.coord[2]);
        ibz_copy(&system[2][2],&gamma_delta.coord[2]);
        ibz_copy(&system[2][3],&bas2.coord[2]);
        ibz_copy(&system[2][4],&bas3.coord[2]);
        ibz_copy(&system[3][0],&gamma_j_delta.coord[3]);
        ibz_copy(&system[3][1],&gamma_ji_delta.coord[3]);
        ibz_copy(&system[3][2],&gamma_delta.coord[3]);
        ibz_copy(&system[3][3],&bas2.coord[3]);
        ibz_copy(&system[3][4],&bas3.coord[3]);

        // resolution of the system 
        found = ibz_4x5_right_ker_mod_prime(&ker,&system,&lideal->norm);
        assert(found);
        // computation of the result
        ibz_copy(&(*C)[0],&ker[0]);
        ibz_copy(&(*C)[1],&ker[1]);

        // finalization of the kernel and system
        ibz_mat_4x5_finalize(&system);
        ibz_vec_5_finalize(&ker);


    }

    // var finalize 
    quat_alg_elem_finalize(&gamma_j);
    quat_alg_elem_finalize(&gamma_delta);
    quat_alg_elem_finalize(&gamma_j_delta);
    quat_alg_elem_finalize(&gamma_ji_delta);
    quat_alg_elem_finalize(&denom1);
    quat_alg_elem_finalize(&denom2);
    quat_alg_elem_finalize(&denom3);
    quat_alg_elem_finalize(&bas3);
    quat_alg_elem_finalize(&bas2);
    return found;

}


/**
 * @brief Representing an integer by the quadratic norm form of a maximal extremal order 
 *
 * @param gamma Output: a quaternion element    
 * @param n_gamma Outut: norm of gamma (also part of the input, it is the target norm a multiple of the final norm)
 * @param Bpoo the quaternion algebra
 * @return 1 if the computation succeeded
 *  
 * This algorithm finds a primitive quaternion element gamma of n_gamma inside the standard maximal extremal order
 */
int represent_integer(quat_alg_elem_t *gamma, ibz_t *n_gamma, const quat_alg_t *Bpoo){
  
    // var dec 
    int found;
    int cnt;
    ibz_t cornacchia_target;
    ibz_t adjusted_n_gamma;
    ibz_t bound,sq_bound,temp;
    quat_alg_coord_t coeffs; // coeffs = [x,y,z,t]
    quat_alg_elem_t quat_temp;
    int prime_list_length;short prime_list[1]; ibz_t prod_bad_primes; // TODO adapt
    prime_list[0] = 5; prime_list_length = 1;  
    ibz_init(&prod_bad_primes);ibz_copy(&prod_bad_primes,&ibz_const_one); 
    
    
    // var init
    found = 0;
    cnt = 0;
    ibz_init(&bound);ibz_init(&temp);ibz_init(&sq_bound);
    quat_alg_coord_init(&coeffs);
    quat_alg_elem_init(&quat_temp);
    ibz_init(&adjusted_n_gamma);
    ibz_init(&cornacchia_target);
    
    //adjusting the norm of gamma (multiplied by 4 to find a solution in the full maximal order)
    ibz_mul(&adjusted_n_gamma,n_gamma,&ibz_const_two);
    ibz_mul(&adjusted_n_gamma,&adjusted_n_gamma,&ibz_const_two);

    // computation of the first bound = sqrt (adjust_n_gamma / p )  
    ibz_div(&sq_bound,&bound,&adjusted_n_gamma,&Bpoo->p);
    ibz_sqrt_floor(&bound,&sq_bound);
    // entering the main loop
    while (!found && cnt < KLPT_repres_num_gamma_trial) { 
        cnt ++;
        // we start by sampling the first coordinate 
        ibz_rand_interval(&coeffs[2],&ibz_const_one,&bound); 
        // then, we sample the second coordinate 
        // computing the second bound in temp as sqrt( (adjust_n_gamma - coeffs[2]²)/p )
        ibz_mul(&cornacchia_target,&coeffs[2],&coeffs[2]);
        ibz_sub(&temp,&sq_bound,&cornacchia_target);
        ibz_sqrt_floor(&temp,&temp);

        if (ibz_cmp(&temp,&ibz_const_zero) ==0) {
            continue;
        } 
        // sampling the second value 
        ibz_rand_interval(&coeffs[3],&ibz_const_one,&temp);

        // compute cornacchia_target = n_gamma - p * (z² + t²)
        ibz_mul(&temp,&coeffs[3],&coeffs[3]);
        ibz_add(&cornacchia_target,&cornacchia_target,&temp);
        ibz_mul(&cornacchia_target,&cornacchia_target,&Bpoo->p);
        ibz_sub(&cornacchia_target,&adjusted_n_gamma,&cornacchia_target);
        // applying cornacchia extended

        found = ibz_cornacchia_extended(&coeffs[0],&coeffs[1],&cornacchia_target,prime_list,prime_list_length, KLPT_primality_num_iter, &prod_bad_primes);

        // check that we can divide by two at least once
        // we must have x = t mod  2 and y = z mod 2
        found = found && (ibz_get(&coeffs[0])%2 == ibz_get(&coeffs[3])%2) && (ibz_get(&coeffs[1])%2 == ibz_get(&coeffs[2])%2);

    }
    if (found) {
        // translate x,y,z,t into the quaternion element gamma
        order_elem_create(gamma,&STANDARD_EXTREMAL_ORDER,&coeffs,Bpoo);
        // making gamma primitive
        // coeffs contains the coefficients of primitivized gamma in the basis of order
        quat_alg_make_primitive(&coeffs,&temp,gamma,&STANDARD_EXTREMAL_ORDER.order,Bpoo);

        // new gamma 
        ibz_mat_4x4_eval(&coeffs,&STANDARD_EXTREMAL_ORDER.order.basis,&coeffs);
        ibz_copy(&gamma->coord[0],&coeffs[0]);
        ibz_copy(&gamma->coord[1],&coeffs[1]);
        ibz_copy(&gamma->coord[2],&coeffs[2]);
        ibz_copy(&gamma->coord[3],&coeffs[3]);
        ibz_copy(&gamma->denom,&STANDARD_EXTREMAL_ORDER.order.denom);

        // adjust the norm of gamma by dividing by the scalar temp²
        ibz_mul(&temp,&temp,&temp);
        ibz_div(n_gamma,&bound,&adjusted_n_gamma,&temp);

    }    
    // var finalize
    ibz_finalize(&bound);ibz_finalize(&temp);ibz_finalize(&sq_bound);
    quat_alg_coord_finalize(&coeffs);
    quat_alg_elem_finalize(&quat_temp);
    ibz_finalize(&adjusted_n_gamma);
    ibz_finalize(&cornacchia_target);
    ibz_finalize(&prod_bad_primes);

    return found;
}

/**
 * @brief Computing the strong approximation
 *
 * @param mu Output: a quaternion element    
 * @param n_mu 4 times target norm of mu
 * @param order special extremal order 
 * @param C coefficients to be lifted
 * @param n_quat_C norm of the quaternion element corresponding to C
 * @param lambda an integer used in the computation
 * @param n modulo for the lift
 * @param max_tries int, bound on the number of tries before we abort
 * @param condition a filter function returning whether a vector should be output or not
 * @param params extra parameters passed to `condition`. May be NULL.
 * @param Bpoo the quaternion algebra
 *  
 * This algorithm finds a quaternion element mu of norm equal to n_mu/4 such that 
 * mu = lambda* (C[0] + i * C[1])*j   mod   ( n * order) 
 * the value of lambda has been computed in a way to make this computation possible
 * the value of mu is computed by the function condition (which depends on some additional params) and must satisfy a set of constraints defined by condition
 * the number of trials is max_tries
 * Assumes n is either a prime or a product of big prime ! this means the probability to encounter non-invertible elements is considered negligible
 */
int strong_approximation(quat_alg_elem_t *mu, const ibz_t *n_mu, const quat_p_extremal_maximal_order_t *order, const ibz_vec_2_t *C, const ibz_t *n_quat_C, const ibz_t *lambda, const ibz_t *n, int max_tries,
int (*condition)(quat_alg_elem_t*, const ibz_vec_2_t*, const void*), const void *params, const quat_alg_t *Bpoo) {
    
    // var dec
    ibz_t ibz_q;
    ibz_t p_lambda_2;
    ibz_t coeff_c0,coeff_c1;
    ibz_t cst_term,temp;
    ibz_t bound;
    ibz_mat_2x2_t lat;
    ibz_vec_2_t target;
    int found;

    // var init 
    ibz_init(&ibz_q); ibz_init(&temp);ibz_init(&bound);
    ibz_init(&p_lambda_2);
    ibz_init(&coeff_c0);ibz_init(&coeff_c1);
    ibz_init(&cst_term);
    ibz_mat_2x2_init(&lat);
    ibz_vec_2_init(&target);

    // computation of various ibz
    ibz_set(&ibz_q,order->q);

    // p_lambda_2 = p * lambda * 2
    ibz_mul(&p_lambda_2,lambda,&Bpoo->p);
    ibz_mul(&p_lambda_2,&p_lambda_2,&ibz_const_two);

    //coeff_c0 = p_lambda_2 * C[0]
    ibz_mul(&coeff_c0,&p_lambda_2,&(*C)[0]);
    ibz_mod(&coeff_c0,&coeff_c0,n);
    //coeff_c1 = 1/(q * p_lambda_2 * C[1]) mod n 
    ibz_mul(&coeff_c1,&p_lambda_2,&(*C)[1]);
    ibz_mul(&coeff_c1,&coeff_c1,&ibz_q);
    ibz_mod(&coeff_c1,&coeff_c1,n);
    ibz_invmod(&coeff_c1,&coeff_c1,n);


    // cst_term = (n_mu-lambda^2 * n_quat_c)/n
    ibz_mul(&temp,lambda,lambda);
    ibz_mul(&temp,&temp,n_quat_C);
    ibz_sub(&temp,n_mu,&temp);
    ibz_div(&cst_term,&temp,&temp,n);
    ibz_mod(&cst_term,&cst_term,n);

    // for debug check that the remainder is zero
    assert(ibz_cmp(&temp,&ibz_const_zero)==0);

    // computation of the lattice 
    // first_column = n * ( 1, - coeff_c0 * coeff_c1 mod n)
    // we start by inverting coeff_c1 
    ibz_copy(&lat[0][0],n);
    ibz_copy(&lat[1][0],&coeff_c1);
    ibz_neg(&lat[1][0],&lat[1][0]);
    ibz_mul(&lat[1][0],&lat[1][0],&coeff_c0);
    ibz_mod(&lat[1][0],&lat[1][0],n);
    ibz_mul(&lat[1][0],&lat[1][0],n);
    
    // second_colum = ( 0, n² ) 
    ibz_copy(&lat[0][1],&ibz_const_zero);
    ibz_copy(&lat[1][1],n);
    ibz_mul(&lat[1][1],&lat[1][1],n);

    // computation of the target vector = (lambda * C[0], lambda * C[1] + n * cst_term * coeff_c1 )
    ibz_mul(&target[0],lambda,&(*C)[0]);
    ibz_mul(&target[1],lambda,&(*C)[1]);
    ibz_mul(&temp,&cst_term,&coeff_c1);
    ibz_mod(&temp,&temp,n);
    ibz_mul(&temp,&temp,n);
    ibz_add(&target[1],&target[1],&temp);

    // computation of the norm log_bound
    ibz_div(&bound,&temp,n_mu,&Bpoo->p);
    int dist_bound = ibz_bitsize(&bound) - 1;
    // to avoid that dist_bound is too big
    if  (dist_bound > 3*ibz_bitsize(n) + 10) {
        dist_bound= 3*ibz_bitsize(n) + 10;
    }

    // enumeration of small vectors until a suitable one is found
    found = quat_2x2_lattice_enumerate_cvp_filter(mu,&lat,&target,order->q,dist_bound,condition,params,max_tries);

    // var finalize
    ibz_finalize(&ibz_q);ibz_finalize(&temp);ibz_finalize(&bound);
    ibz_finalize(&p_lambda_2);
    ibz_finalize(&coeff_c0);ibz_finalize(&coeff_c1);
    ibz_finalize(&cst_term);
    ibz_mat_2x2_finalize(&lat);
    ibz_vec_2_finalize(&target);
    
    return found;
}


/**
 * @brief Computing a list of good norm
 *
 * @param list Output: a list of ibz
 * @param list_size int, size of the list
 * @param log_bound int, logarithmic bound 
 * @param n ibz, the modulo
 * @param target_reduosity int, a legendre symbol
 * @returns the number of integers above the bound that were able to set
 *  
 * This algorithm finds a list of size list_size  
 * of ibz dividing TORSION_ODD^2 
 * such that the log of these ibz is bigger than log_bound   
 * and the legendre symbol of these integers mod n is equal to target_reduosity
 * assumes that TORSION_ODD is divided by a big power of 3, so we look for T² and T² divided by small powers of 3 (possibly divided by another factor to satisfy the reduosity condition). 
 */
int norm_list_computation(ibz_t *list, int list_size, int log_bound, const ibz_t *n, int target_reduosity) {
    
    // var dec
    ibz_t temp,remainder; 
    ibz_t elli;
    int is_smaller;
    int cnt = 0;
    int num_found=0;
    int found = 0;
    int log_freedom =0;
    ibz_t const *const T = &TORSION_ODD;
    size_t const num_fac = sizeof(TORSION_ODD_PRIMES) / sizeof(*TORSION_ODD_PRIMES);
    uint64_t const *const fac = TORSION_ODD_PRIMES;


    // var init
    ibz_init(&temp);ibz_init(&remainder);
    ibz_init(&elli);

    assert(fac[0]==3);  // by construction of p

    // temp = T²
    ibz_mul(&temp, T, T);

    // log_freedom is the bitsize of temp
    log_freedom = ibz_bitsize(&temp);
    // indicates if we have a chance to find an integer bigger than the bound
    is_smaller = (log_bound < log_freedom);
    // elli = fac[0]; 
    ibz_set(&elli,fac[0]);
    int small_reduosity = ibz_legendre(&elli,n);

    // in that case we look for a number that is a square
    if (target_reduosity == 1  ) {
        
        while (is_smaller && num_found < list_size ) {
            // copy the value to the output list 
            ibz_copy(&(list[num_found]),&temp);
            num_found++;

            // we set up to compute another value 
            if (small_reduosity == 1) {
                // we can afford to remove only fac[0] since it is a square
                ibz_set(&remainder,fac[0]);
            }
            else {
                // we need to remove fac[0]² because fac[0] is not a square
                ibz_set(&remainder,fac[0]*fac[0]);
            }
            ibz_div(&temp,&remainder,&temp,&remainder);
            log_freedom = ibz_bitsize(&temp);
            is_smaller = (log_bound < log_freedom);
        }
        
    }
    else if (is_smaller) {
        
        // first, we look for a factor that is not a square that we can remove from T² to obtain an integer with the right reduosity
        while (!found && num_found < num_fac ) {
            
            ibz_set(&elli,fac[num_found]);
            if (ibz_legendre(&elli,n) == -1 ) {
                ibz_div(&temp,&remainder,&temp,&elli);
                assert(ibz_cmp(&remainder,&ibz_const_zero)==0);
                log_freedom = ibz_bitsize(&temp);
                is_smaller = (log_bound < log_freedom);
                found = 1;
                
            }
            num_found++;
        }
        if (!found) {
            // no suitable factor was found, we cannot find any suitable integer
            return 0;
        }
        // now we actually try to find integers 
        num_found =0;
        while (is_smaller && num_found < list_size) {
            // copy the good value to the output list 
            ibz_copy(&list[num_found],&temp);            
            num_found++;

            // we set up to compute another value 
            // we need to remove fac[0]²
            if (small_reduosity == 1) {
                // we can afford to remove only fac[0] since it is a square
                ibz_set(&remainder,fac[0]);
            }
            else {
                // we need to remove fac[0]² because fac[0] is not a square
                ibz_set(&remainder,fac[0]*fac[0]);
            }
            ibz_div(&temp,&remainder,&temp,&remainder);
            assert(ibz_cmp(&remainder,&ibz_const_zero)==0);
            log_freedom = ibz_bitsize(&temp);
            is_smaller = (log_bound < log_freedom);
        }
    }

    

    // var finalize
    ibz_finalize(&temp);ibz_finalize(&remainder);
    ibz_finalize(&elli);

    return num_found;

}
