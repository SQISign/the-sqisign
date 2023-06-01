/** @file
 *
 * @authors Sina Schaeffler
 *
 * @brief Declarations of tests of quaternion algebra operations
 */

#ifndef QUATERNION_TESTS_H
#define QUATERNION_TESTS_H

#include <quaternion.h>
#include <stdio.h>
#include "../internal.h"

/** @internal 
 * @ingroup quat_quat
 * @defgroup quat_tests Quaternion module test functions
 * @{
 */

/** @brief Test initializers and finalizers for quaternion algebra types
 * 
 * Test initializers and finalizers for the following types: 
 * 
 * quat_alg_t
 * 
 * quat_alg_elem_t
 * 
 * quat_alg_coord_t
 * 
 * ibz_vec_4_t
 * 
 * ibz_vec_5_t
 * 
 * ibz_mat_2x2_t
 * 
 * ibz_mat_4x4_t
 * 
 * ibz_mat_4x5_t
 * 
 * ibz_mat_4x8_t
 * 
 * quat_lattice_t
 * 
 * quat_order_t
 * 
 * quat_left_ideal_t
 */
int quat_test_finit();

/** @brief Test integer, quadratic form and matrix functions for dimension 4 from the quaternion module
 * 
 * Runs unit tests for the following functions
 * 
 * int ibz_4x5_right_ker_mod_prime(ibz_vec_5_t *ker, const ibz_mat_4x5_t *mat, const ibz_t *p);
 * 
 * int ibz_4x4_right_ker_mod_prime(ibz_vec_4_t *ker, const ibz_mat_4x4_t *mat, const ibz_t *p);
 * 
 * int ibz_4x4_right_ker_mod_power_of_2(ibz_vec_4_t *ker, const ibz_mat_4x4_t *mat, unsigned short exp);
 * 
 * void ibz_mat_4x4_eval(quat_alg_coord_t  *res, const ibz_mat_4x4_t *mat, const quat_alg_coord_t *vec);
 * 
 * void quat_qf_eval(ibz_t *res, const ibz_mat_4x4_t *qf, const quat_alg_coord_t *coord);
 * 
 * void ibz_content(ibz_t *content, const quat_alg_coord_t *v);
 */
int quat_test_dim4();

/** @brief Test integer, lattice and matrix functions for dimension 2 from the quaternion module
 * 
 * Runs unit tests for the following functions
 * 
 * void ibz_mat_2x2_eval(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec);
 * 
 * void ibz_2x2_mul_mod(ibz_mat_2x2_t *prod, const ibz_mat_2x2_t *mat_a, const ibz_mat_2x2_t *mat_b, const ibz_t *m);
 * 
 * int ibz_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);
 * 
 * int quat_2x2_lattice_enumerate_cvp_filter(quat_alg_elem_t *res, const ibz_mat_2x2_t *lat_basis, const ibz_vec_2_t *target,unsigned int qf, unsigned int dist_bound, int (*condition)(quat_alg_elem_t* , const ibz_vec_2_t*, const void*), const void* params, unsigned int max_tries);
 */
int quat_test_dim2();

/** @brief Test integer functions
 * 
 * Runs unit tests for the following functions
 * 
 * int ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p);
 * 
 * int ibz_cornacchia_special_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p, const int exp_adjust);
 * 
 * int ibz_cornacchia_extended(ibz_t *x, ibz_t *y, const ibz_t *n, const short *prime_list, const int prime_list_length, short primality_test_iterations, const ibz_t *bad_primes_prod); 
 */
int quat_test_integers();

/** @brief Test operations on quaternion algebra elements
 * 
 * Runs unit tests for the following functions
 * 
 * void quat_alg_add(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
 * 
 * void quat_alg_sub(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
 * 
 * void quat_alg_mul(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b, const quat_alg_t *alg);
 * 
 * void quat_alg_norm(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_t *alg);
 * 
 * void quat_alg_trace(quat_alg_elem_t *res, const quat_alg_elem_t *a);
 * 
 * void quat_alg_scalar(quat_alg_elem_t *elem, const ibz_t *numerator, const ibz_t *denominator);
 * 
 * void quat_alg_conj(quat_alg_elem_t *conj, const quat_alg_elem_t *x);
 * 
 * void quat_alg_make_primitive(quat_alg_coord_t *primitive_x, ibz_t *content, const quat_alg_elem_t *x, const quat_order_t *order, const quat_alg_t *alg){
 * 
 * int quat_alg_is_primitive(const quat_alg_elem_t *x, const quat_order_t *order, const quat_alg_t *alg):
 * 
 * void quat_alg_normalize(quat_alg_elem_t *x);
 * 
 * int quat_alg_elem_is_zero(const quat_alg_elem_t *x);
 * 
 * int quat_alg_coord_is_zero(const quat_alg_coord_t *x);
 */
int quat_test_algebra();

/** @brief Test operations on lattices
 * 
 * Runs unit tests for the following functions
 * 
 * void quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);
 * 
 * void quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);
 * 
 * void quat_lattice_mul(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2, const quat_alg_t *alg);
 * 
 * int  quat_lattice_contains(quat_alg_coord_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t *x, const quat_alg_t *alg);
 * 
 * void quat_lattice_hnf(quat_lattice_t *lat);
 */
int quat_test_lattice();

/** @brief Test operations on left ideals and their creation
 * 
 * Runs unit tests for the following functions
 * 
 * void quat_lideal_create_principal(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const quat_order_t *order, const quat_alg_t *alg);
 * 
 * void quat_lideal_create_from_primitive(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const ibz_t *N, const quat_order_t *order, const quat_alg_t *alg); 
 * 
 * void quat_lideal_make_primitive_then_create(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const ibz_t *N, const quat_order_t *order, const quat_alg_t *alg); 
 * 
 * void quat_lideal_random_2e(quat_left_ideal_t *lideal, const quat_order_t *order, const quat_alg_t *alg, int64_t e); 
 * 
 * int quat_lideal_generator(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const quat_alg_t *alg, int bound);
 * 
 * int quat_lideal_generator_coprime(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const ibz_t *n, const quat_alg_t *alg, int bound);
 * 
 * int quat_lideal_mul(quat_left_ideal_t *product, const quat_left_ideal_t *lideal, const quat_alg_elem_t *alpha, const quat_alg_t *alg, int bound);
 * 
 * void quat_lideal_add(quat_left_ideal_t *sum, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
 * 
 * void quat_lideal_inter(quat_left_ideal_t *intersection, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
 * 
 * int quat_lideal_equals(const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
 * 
 * int quat_lideal_isom(quat_alg_elem_t *iso, const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t *alg);
 * 
 * void quat_lideal_right_order(quat_order_t *order, const quat_left_ideal_t *lideal, const quat_alg_t *alg);
 * 
 * void quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced, ibz_mat_4x4_t *gram, const quat_left_ideal_t *lideal, const quat_alg_t *alg); //replaces lideal_lll
 * 
 * void quat_connecting_ideal(quat_left_ideal_t *connecting_ideal, const quat_order_t *O1, const quat_order_t *O2, const quat_alg_t *alg);
*/
int quat_test_lideal();

/** @brief test matkermod and Howell form
 */
int quat_test_matkermod();

/** @brief Test with randomization for complex functions where this is possible
 * 
 */
int quat_test_with_randomization();


/** @}
 */

#endif
