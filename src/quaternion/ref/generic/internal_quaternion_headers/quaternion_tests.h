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
#include "internal.h"

/** @internal
 * @ingroup quat_helpers
 * @defgroup quat_tests Quaternion module test functions
 * @{
 */

/** @internal
 * @ingroup quat_tests
 * @defgroup quat_test_inputs Quaternion module random test input generation
 * @{
 */

/**
 * @brief Generates list of random ideals of a special extremal order
 *
 * @param ideals Output: Array of iterations left ideals of the given special extremal order and of
 * norms of size norm_bitsize
 * @param norm_bitsize Bitsize of the norms of the outut ideals
 * @param iterations Number of ideals to sample. Most be smaller than ideals is long
 * @param params quat_represent_integer_params_t parameters for the left order of the ideals to be
 * sampled.
 * @return 0 if success, 1 if failure
 */
int quat_test_input_random_ideal_generation(quat_left_ideal_t *ideals,
                                            int norm_bitsize,
                                            int iterations,
                                            const quat_represent_integer_params_t *params);

/**
 * @brief Generates list of random ideals of a special extremal order
 *
 * @param lattices Output: Array of iterations left ideals of the given special extremal order and
 * of norms of size norm_bitsize, given only by their lattices
 * @param norms Output: Array which will contain the norms of the sampled ideals, in the same order
 * as their lattices are. Can be NULL, in which case no norm is output. Otherwise, it must be an
 * array of length at least iterations
 * @param norm_bitsize Bitsize of the norms of the outut ideals
 * @param iterations Number of ideals to sample. Most be smaller than lattices is long
 * @param params quat_represent_integer_params_t parameters for the left order of the ideals to be
 * sampled.
 * @return 0 if success, 1 if failure
 */
int quat_test_input_random_ideal_lattice_generation(quat_lattice_t *lattices,
                                                    ibz_t *norms,
                                                    int norm_bitsize,
                                                    int iterations,
                                                    const quat_represent_integer_params_t *params);

/**
 * @brief Generates list of random lattices
 *
 * @param lattices Output: Array of iterations lattices
 * @param bitsize Bitsize of the coefficients of a random basis of the lattices
 * @param iterations Number of lattices to sample.Most be smaller than lattices is long
 * @param in_hnf If not 0, the lattices are put in HNF before being outputs. Their coefficients in
 * this basis might be larger than bitsize.
 * @return 0 if success, 1 if failure
 */
int quat_test_input_random_lattice_generation(quat_lattice_t *lattices, int bitsize, int iterations, int in_hnf);
/**
 * @}
 */

/** @brief Test for integer functions
 *
 * void ibz_init(ibz_t *x);
 *
 * void ibz_finalize(ibz_t *x);
 *
 * void ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b);
 *
 * void ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b);
 *
 * void ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b);
 *
 * void ibz_neg(ibz_t *neg, const ibz_t *a);
 *
 * void ibz_abs(ibz_t *abs, const ibz_t *a);
 *
 * void ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b);
 *
 * void ibz_div_2exp(ibz_t *quotient, const ibz_t *a, uint32_t exp);
 *
 * void ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b);
 *
 * unsigned long int ibz_mod_ui(const mpz_t *n, unsigned long int d);
 *
 * int ibz_divides(const ibz_t *a, const ibz_t *b);
 *
 * void ibz_pow(ibz_t *pow, const ibz_t *x, uint32_t e);
 *
 * void ibz_pow_mod(ibz_t *pow, const ibz_t *x, const ibz_t *e, const ibz_t *m);
 *
 * int ibz_cmp(const ibz_t *a, const ibz_t *b);
 *
 * int ibz_is_zero(const ibz_t *x);
 *
 * int ibz_is_one(const ibz_t *x);
 *
 * int ibz_cmp_int32(const ibz_t *x, int32_t y);
 *
 * int ibz_is_even(const ibz_t *x);
 *
 * int ibz_is_odd(const ibz_t *x);
 *
 * void ibz_set(ibz_t *i, int32_t x);
 *
 * void ibz_copy(ibz_t *target, const ibz_t *value);
 *
 * void ibz_swap(ibz_t *a, ibz_t *b);
 *
 * void ibz_copy_digits(ibz_t *target, const digit_t *dig, int dig_len);
 *
 * void ibz_to_digits(digit_t *target, const ibz_t *ibz);
 *
 * int32_t ibz_get(const ibz_t *i);
 *
 * int ibz_rand_interval(ibz_t *rand, const ibz_t *a, const ibz_t *b);
 *
 * int ibz_rand_interval_minm_m(ibz_t *rand, int32_t m);
 *
 * int ibz_bitsize(const ibz_t *a);
 *
 * void ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b);
 *
 * int ibz_invmod(ibz_t *inv, const ibz_t *a, const ibz_t *mod);
 *
 * void ibz_sqrt_floor(ibz_t *sqrt, const ibz_t *a);
 */
int ibz_test_intbig(void);

/** @brief Test for implementations of GMP functions missing from the mini-GMP API
 *
 * int mpz_legendre(const mpz_t a, const mpz_t p);
 *
 * double mpz_get_d_2exp(signed long int *exp, const mpz_t op);
 */
int mini_gmp_test(void);

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
 * ibz_vec_2_t
 *
 * ibz_vec_4_t
 *
 * ibz_mat_2x2_t
 *
 * ibz_mat_4x4_t
 *
 * quat_lattice_t
 *
 * quat_lattice_t
 *
 * quat_left_ideal_t
 */
int quat_test_finit(void);

/** @brief Test integer, quadratic form and matrix functions for dimension 4 from the quaternion
 * module
 *
 * Runs unit tests for the following functions
 *
 * void ibz_mat_4x4_eval(quat_alg_coord_t  *res, const ibz_mat_4x4_t *mat, const quat_alg_coord_t
 * *vec);
 *
 * void quat_qf_eval(ibz_t *res, const ibz_mat_4x4_t *qf, const quat_alg_coord_t *coord);
 *
 * void ibz_vec_4_content(ibz_t *content, const quat_alg_coord_t *v);
 */
int quat_test_dim4(void);

/** @brief Test integer, lattice and matrix functions for dimension 2 from the quaternion module
 *
 * Runs unit tests for the following functions
 *
 * void ibz_mat_2x2_copy(ibz_vec_2_t *copy, const ibz_mat_2x2_t *copied);
 *
 * void ibz_mat_2x2_eval(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec);
 *
 * void ibz_mat_4x4_eval_t(ibz_vec_4_t *res, const ibz_vec_4_t *vec, const ibz_mat_4x4_t *mat);
 *
 * void ibz_2x2_mul_mod(ibz_mat_2x2_t *prod, const ibz_mat_2x2_t *mat_a, const ibz_mat_2x2_t *mat_b,
 * const ibz_t *m);
 *
 * int ibz_mat_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);
 */
int quat_test_dim2(void);

/** @brief Test integer functions
 *
 * Runs unit tests for the following functions
 *
 * int ibz_generate_random_prime(ibz_t *p, int is3mod4, int bitsize);
 *
 * int ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p);
 */
int quat_test_integers(void);

/** @brief Test operations on quaternion algebra elements
 *
 * Runs unit tests for the following functions
 *
 * void quat_alg_add(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
 *
 * void quat_alg_sub(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
 *
 * void quat_alg_mul(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b, const
 * quat_alg_t *alg);
 *
 * void quat_alg_norm(ibz_t *res_num, ibz_t *res_denom, const quat_alg_elem_t *a, const quat_alg_t
 * *alg);
 *
 * void quat_alg_scalar(quat_alg_elem_t *elem, const ibz_t *numerator, const ibz_t *denominator);
 *
 * void quat_alg_conj(quat_alg_elem_t *conj, const quat_alg_elem_t *x);
 *
 * void quat_alg_make_primitive(quat_alg_coord_t *primitive_x, ibz_t *content, const quat_alg_elem_t
 * *x, const quat_lattice_t *order, const quat_alg_t *alg){
 *
 * void quat_alg_normalize(quat_alg_elem_t *x);
 *
 * int quat_alg_elem_is_zero(const quat_alg_elem_t *x);
 *
 * int quat_alg_coord_is_zero(const quat_alg_coord_t *x);
 *
 * void quat_alg_elem_copy(quat_alg_elem_t *copy, const quat_alg_elem_t *copied);
 */
int quat_test_algebra(void);

/** @brief Test operations on lattices
 *
 * Runs unit tests for the following functions
 *
 * void quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
 * *lat2);
 *
 * void quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
 * *lat2);
 *
 * void quat_lattice_conjugate_without_hnf(quat_lattice_t *conj, const quat_lattice_t *lat);
 *
 * void quat_lattice_alg_elem_mul(quat_lattice_t *prod, const quat_lattice_t *lat, const
 * quat_alg_elem_t *elem, const quat_alg_t *alg);
 *
 * void quat_lattice_mul(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t
 * *lat2, const quat_alg_t *alg);
 *
 * int  quat_lattice_contains(quat_alg_coord_t *coord, const quat_lattice_t *lat, const
 * quat_alg_elem_t *x, const quat_alg_t *alg);
 *
 * void quat_lattice_hnf(quat_lattice_t *lat);
 */
int quat_test_lattice(void);

/** @brief Test for lattice reduction and functions based on it
 *
 * int quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced, ibz_mat_4x4_t *gram, const quat_left_ideal_t
 * *lideal, const quat_alg_t *alg);
 *
 * int quat_lideal_lideal_mul_reduced(quat_left_ideal_t *prod, ibz_mat_4x4_t *gram, const
 * quat_left_ideal_t *lideal1,const quat_left_ideal_t *lideal2, const quat_alg_t *alg);
 *
 * int quat_lideal_prime_norm_reduced_equivalent(quat_left_ideal_t *lideal, const quat_alg_t *alg,
 * const int primality_num_iter, const int equiv_bound_coeff);
 */
int quat_test_lll(void);

/** @brief Test operations on left ideals and their creation
 *
 * Runs unit tests for the following functions
 *
 * void quat_lideal_copy(quat_left_ideal_t *copy, const quat_left_ideal_t *copied)
 *
 * void quat_lideal_create_principal(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const
 * quat_lattice_t *order, const quat_alg_t *alg);
 *
 * void quat_lideal_create(quat_left_ideal_t *lideal, const quat_alg_elem_t *x, const
 * ibz_t *N, const quat_lattice_t *order, const quat_alg_t *alg);
 *
 * int quat_lideal_generator(quat_alg_elem_t *gen, const quat_left_ideal_t *lideal, const
 * quat_alg_t);
 *
 * void quat_lideal_add(quat_left_ideal_t *sum, const quat_left_ideal_t *I1, const quat_left_ideal_t
 * *I2, const quat_alg_t *alg);
 *
 * void quat_lideal_inter(quat_left_ideal_t *intersection, const quat_left_ideal_t *I1, const
 * quat_left_ideal_t *I2, const quat_alg_t *alg);
 *
 * int quat_lideal_equals(const quat_left_ideal_t *I1, const quat_left_ideal_t *I2, const quat_alg_t
 * *alg);
 *
 */
int quat_test_lideal(void);

/** @brief Test operations on represent integer
 *
 * int quat_sampling_random_ideal_O0_given_norm(quat_left_ideal_t *lideal,const ibz_t *norm,int
 * is_prime,const quat_represent_integer_params_t *params,int prime_sampling_attempts);
 *
 * void quat_change_to_O0_basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el);
 *
 * int quat_represent_integer(quat_alg_elem_t *gamma, const ibz_t *n_gamma,  int non_diag, const
 * quat_represent_integer_params_t *params);
 */
int quat_test_normeq(void);

/** @brief Test functions for sampling lattice points of bounded norm
 *
 * int quat_lattice_sample_from_ball(ibz_vec_4_t *x, const quat_lattice_t *lattice, const quat_alg_t
 * *alg, const ibz_t *radius);
 */
int quat_test_lat_ball(void);

/** @brief Test for hnf core and verification
 *
 */
int quat_test_hnf(void);

/** @brief Test with randomization for complex functions where this is possible
 *
 */
int quat_test_with_randomization(void);

/** @}
 */

#endif
