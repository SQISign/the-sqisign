#include <intbig.h>
#include <quaternion.h>
#include <klpt.h>

/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief the klpt_tools header 
 */

#ifndef KLPT_TOOLS_H
#define KLPT_TOOLS_H

/** @internal
 * @ingroup klpt_klpt
 * @defgroup klpt_internal KLPT module internal functions and types
 * @{
 */

/** @internal
 * @defgroup klpt_param_t Types for klpt parameters
 * @{
*/

/** @brief Type for klpt parameter * 
 * @typedef signing_klpt_param_t
 *
 * @struct signing_klpt_param
 * 
 * the structure used in klpt for the call to strong approx
*/
typedef struct signing_klpt_param {
    

    quat_alg_elem_t gamma;
    quat_alg_elem_t delta;
    ibz_t target_norm;
    quat_alg_t Bpoo;
    const quat_p_extremal_maximal_order_t *order;
    ibz_t n;
    ibz_t equiv_n;


} signing_klpt_param_t;


/** @brief Type for quaternion algebras
 * 
 * @typedef eichler_norm_param_t
 *
 * @struct eichler_norm_param
 * 
 * the structure used in eichler_norm_param for the call to strong approx 
 */
typedef struct eichler_norm_param {

    quat_alg_elem_t gen_constraint;
    
    
    ibz_t target_norm;
    quat_alg_t Bpoo;
    quat_p_extremal_maximal_order_t order;
    ibz_t n;
    quat_order_t right_order;

} eichler_norm_param_t;


/**
 * @}
*/

/** @}
*/


/** @internal
 * @ingroup klpt_tools
 * @{
*/

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
void order_elem_create(quat_alg_elem_t *elem, const quat_p_extremal_maximal_order_t *order, const quat_alg_coord_t *coeffs, const quat_alg_t *Bpoo);

/**
 * @brief Finding the good linear span(j,i*j) combination
 *
 * @param C Output: two integers modulo n  
 * @param order special extremal order 
 * @param delta quaternion element
 * @param gamma quaternion element
 * @param lideal left order-ideal
 * @param Bpoo quaternion algebra
 * @param is_divisible boolean indicating if Norm(gamma) is divisible by n
 *  
 * This algorithm finds two integers C[0], C[1] such that gamma * j * (C[0] + i C[1]]) * delta is in the eichler order ZZ + lideal  
 * if gamma have norm divisible by n(lideal) (indicated by the boolean is_divisible), then the result will be contained in J 
 * Failure should not be possible
 * Assumes n is prime !
 */
int solve_combi_eichler(ibz_vec_2_t *C, const quat_p_extremal_maximal_order_t *order, const quat_alg_elem_t *gamma, const quat_alg_elem_t *delta, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo, int is_divisible);


int ibz_cornacchia_special_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p, const int exp_adjust);




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
int (*condition)(quat_alg_elem_t*, const ibz_vec_2_t*, const void*), const void *params, const quat_alg_t *Bpoo);


/**
 * @brief Computing a list of good norm
 *
 * @param list Output: a list of ibz
 * @param list_size int, size of the list
 * @param log_bound int, logarithmic bound 
 * @param n ibz, the modulo
 * @param target_reduosity int, a legendre symbol
 * @returns a bit indicating if at least one factor was found
 * 
 * This algorithm finds a list of size list_size  
 * of ibz dividing T^2 
 * such that the log of these ibz is smaller than log_bound and their legendre symbol   
 * and the legendre symbol of these ibz mod n is equal to target_reduosity
 */
int norm_list_computation(ibz_t *list, int list_size, int log_bound, const ibz_t *n, int target_reduosity);

/** @}
*/


#endif
