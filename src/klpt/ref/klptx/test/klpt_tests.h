#ifndef KLPT_TESTS_H
#define KLPT_TESTS_H

// #include <klpt.h>
#include <quaternion.h>
#include <intbig.h>
#include <stdio.h>
#include "../tools.h" 
#include <rng.h>
#include <quaternion_data.h>

/** @internal 
 * @ingroup klpt_klpt
 * @defgroup klpt_tests KLPT module test functions
 * @{
 */

/** @brief Test for the tools of the KLPT module, covers the following functions :
 * 
 * int represent_integer(quat_alg_elem_t *gamma, ibz_t *n_gamma, const quat_alg_t *Bpoo);
 * int solve_combi_eichler(ibz_vec_2_t *C, const quat_p_extremal_maximal_order_t *order, const quat_alg_elem_t *gamma, const quat_alg_elem_t *delta, const quat_left_ideal_t *lideal, const quat_alg_t *Bpoo, int is_divisible);
 * int klpt_find_linear_comb(ibz_vec_2_t *C,const quat_alg_elem_t *beta, const quat_order_t *order, const ibz_t *n, const unsigned short exp, const quat_alg_elem_t *gen_start, const quat_alg_elem_t *gen_end,const quat_alg_t *Bpoo);
 * 
 */
int klpt_test_tools();

/** @brief Test for the equiv of the KLPT module, covers the following functions :
 * 
 * int klpt_lideal_equiv(quat_alg_elem_t *gen, ibz_t *n, const quat_lattice_t *reduced, const ibz_mat_4x4_t *gram, const ibz_t *lideal_norm, const quat_alg_t *Bpoo) 
 * 
 */
int klpt_test_equiv();

int klpt_test_klpt();

int klpt_test_eichler();

/** @}
 */

#endif
