/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief The norm equation algorithms 
 */

#ifndef KLPT_H
#define KLPT_H

#include <intbig.h>
#include <quaternion.h>
#include <klpt_constants.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <stdio.h>


/*************************** Functions *****************************/

/** @defgroup klpt_klpt Functions and types for KLPT
 * @{
*/

/** @defgroup klpt_extremals Init for extremal orders
 * @{
*/
static inline int generate_random_prime(ibz_t *p , int is3mod4 ,int bitsize) {
    int found = 0;
    ibz_t two_pow,two_powp;

    ibz_init(&two_pow);ibz_init(&two_powp);
    ibz_pow(&two_pow,&ibz_const_two,bitsize);
    ibz_pow(&two_powp,&ibz_const_two,bitsize+1);

    int cnt = 0;
    while (!found && cnt < KLPT_random_prime_attempts*bitsize) {
        cnt ++ ;
        ibz_rand_interval(p,&two_pow,&two_powp);
        
        found = ibz_probab_prime(p,30) && (!is3mod4 || (ibz_get(p)%4 == 3)) ;
    }
    ibz_finalize(&two_pow);ibz_finalize(&two_powp); 
    return found;

}


void quat_alg_elem_copy(quat_alg_elem_t *copy,const quat_alg_elem_t *copied);
void quat_left_ideal_copy(quat_left_ideal_t *copy,const quat_left_ideal_t *copied);


/**
 * @brief Representing an integer by the quadratic norm form of a maximal extremal order 
 *
 * @param gamma Output: a quaternion element    
 * @param n_gamma Output : target norm of gamma (it is also an input, the final value will be a divisor of the initial value)
 * @param Bpoo the quaternion algebra
 *  
 * This algorithm finds a primitive quaternion element gamma of n_gamma inside the standard maximal extremal order
 * Failure is possible
 */
int represent_integer(quat_alg_elem_t *gamma, ibz_t *n_gamma, const quat_alg_t *Bpoo);

/** @}
*/

/** @defgroup klpt_equiv Finding equivalent left ideals
 * @{
*/

/**
 * @brief Keygen random ideal
 * @param ideal : Output : random ideal
 * @param order : maximal extremal order
 * @param Bpoo the quaternion algebra
 * computes a keygen ideal
 */
int klpt_keygen_random_ideal(quat_left_ideal_t *ideal, const quat_p_extremal_maximal_order_t *order, const quat_alg_t *Bpoo);

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
void klpt_lideal_equiv_random_eichler(quat_alg_elem_t *gen, const ibz_t *n, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo);

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
int klpt_lideal_equiv(quat_alg_elem_t *gen, ibz_t *n, const ibz_mat_4x4_t *reduced, const ibz_mat_4x4_t *gram, const ibz_t *lideal_norm, const ibz_t *denom, const quat_alg_t *Bpoo); 


/** @}
*/


/** @defgroup klpt_klpt_general Finding equivalent left ideals of power of 2 norm
 * @{
*/
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
 * equiv_lideal = I * gen/Norm(lideal) of norm n = 2^KLPT_signing_length
 * where conjugate(gen) lideal 
 * moreover we need that gen * delta is contained in the eichler order ZZ + J
 * Assumes that the ideals lideal and lideal_start has a "good" norm (ie that one needs not to apply lideal_equiv)
 * returns a bit indicating if the computation has succeeded
 */
int klpt_signing_klpt(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_left_ideal_t *lideal_start, const quat_alg_elem_t *delta, const quat_alg_t *Bpoo);

/**
 * @brief Equivalent left O0-ideal of power of two norm
 *
 * @param gen Output: generator of equiv  
 * @param lideal left O0 ideal
 * @param Bpoo the quaternion algebra
 *  
 * 
 * equiv_lideal is equivalent to the ideal lideal and we have the equality
 * equiv_lideal = lideal * gen/Norm(lideal) of norm n = 2^KLPT_keygen_length
 * where conjugate(gen) is in lideal 
 * Assumes that the ideals lideal has a "good" norm (ie that one needs not to apply lideal_equiv)
 * returns a bit indicating if the computation has succeeded
 */
int klpt_keygen_klpt(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo);


/** @}
*/


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
int klpt_eichler_special_norm(quat_alg_elem_t *beta, ibz_t *n_beta, quat_alg_elem_t *gen, ibz_t *n, quat_left_ideal_t *lideal, quat_alg_elem_t *gen_constraint, const quat_alg_t *Bpoo); 


/** @}
*/

/** @}
*/



/** @ingroup klpt_klpt
 * @defgroup klpt_tools Various tools for KLPT and norm equations
 * @{
*/

/**
 * @brief Finding the good linear combination
 *
 * @param C Output: the linear combination
 * @param beta the endomorphism 
 * @param order the order
 * @param n the norm
 * @param exp the exponent
 * @param gen_start element of order, generate a left O-ideal of norm n = 2^exp
 * @param gen_end element of order, generate a left O-ideal of norm n = 2^exp
 * @param Bpoo quaternion algebra
 *  
 * lideal_start = order  <gen_start,n>
 * lideal_end = order  <gen_end,n>
 * This algorithm finds two integers C[1], C[2] such that the pushforward of the ideal lideal_start by the endomorphism C[1] + C[2] beta is equal to lideal_end
 * beta is an endomorphism of order
 * Returns a bit indicating if the computation succeeded 
   */
int klpt_find_linear_comb(ibz_vec_2_t *C,const quat_alg_elem_t *beta, const quat_order_t *order, const ibz_t *n, const unsigned short exp, const quat_alg_elem_t *gen_start, const quat_alg_elem_t *gen_end,const quat_alg_t *Bpoo);


/** @}
*/


#endif
