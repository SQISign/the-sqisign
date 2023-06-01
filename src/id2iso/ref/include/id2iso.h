/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief The id2iso algorithms 
 */

#ifndef ID2ISO_H
#define ID2ISO_H

#include <intbig.h>
#include <quaternion.h>
#include <klpt.h>
#include <ec.h>

/** @defgroup id2iso_id2iso Ideal to isogeny conversion
 * @{
*/

/** @defgroup id2iso_iso_types Types for isogenies needed for the id2iso  
 * @{
*/

/** @brief Type for long chain of two isogenies
 * 
 * @typedef id2iso_long_two_isog
 * 
 * Represented as a vector ec_isog_even_t
*/
typedef struct id2iso_long_two_isog {
    unsigned short length; ///< the number of smaller two isogeny chains
    ec_isog_even_t *chain; ///< the chain of two isogeny
} id2iso_long_two_isog_t;

/** @brief Type for compressed long chain of two isogenies
 * 
 * @typedef id2iso_compressed_long_two_isog 
 * 
 * Represented as a vector of intbig
*/
typedef struct id2iso_compressed_long_two_isog {
    unsigned short length; ///< the number of smaller two isogeny chains
    ibz_t *zip_chain; ///< the chain of two isogeny, compressed
    unsigned char bit_first_step ; ///< the bit for the first step
} id2iso_compressed_long_two_isog_t;

/** @}
*/

/*************************** Functions *****************************/


/** @defgroup id2iso Constructors and Destructors
 * @{
*/
void id2iso_long_two_isog_init(id2iso_long_two_isog_t *isog, const size_t length);
void id2iso_long_two_isog_finalize(id2iso_long_two_isog_t *isog);

void id2iso_compressed_long_two_isog_init(id2iso_compressed_long_two_isog_t *zip, const size_t length);
void id2iso_compressed_long_two_isog_finalize(id2iso_compressed_long_two_isog_t *zip);

/** @}
*/

/** @defgroup id2iso_others Other functions needed for id2iso
 * @{
*/

/**
 * @brief Translating an ideal of norm a big power of 2 into the corresponding isogeny
 *
 * @param isog_zip Output : compression of the output isogeny  
 * @param basis_minus : odd torsion basis (in the end, this will be the basis pushed through the output isogeny)
 * @param basis_plus : odd torsion basis (in the end, this will be the basis pushed through the output isogeny)
 * @param domain : the starting curve (in the end, this will be the codomain of the output isogeny)
 * @param kernel_dual : the dual of the kernel of the last step of isog_start_two (in the end this will be the kernel of the dual of the last step of the output isogeny) 
 * @param gen_input : quaternion element, element of a maximal order O, generator of the O-ideal to be translated 
 * @param length : the length of the chain to be translated
 * @param lideal_start_small : a small ideal equivalent to lideal_start_two of right order equal to O
 * @param lideal_start_two : O0-ideal of norm a power of 2 equivalent to lideal_start_small, corresponding to an isogeny isog_start_two
 * @param gen_two element of O0, generator of lideal_start_two
 * @param Bpoo : the quaternion algebra
 * @returns a bit indicating if the computation succeeded
 *  /!\ the composition of isog_start_two and isog might be backtracking
 * lideal_start_two = O0 < gen_two , 2^*> 
 * lideal_start_small = O0 < conj(gen_two), * >
 * lideal_start_small = lideal_lideal_start_two * conj(gen_two) / 2^*  
 * The ideal to be translated is equal to O < gen_input, 2^e> where O = OR(lideal_start_small)  
 * 
 * assumes that the ideal given in input has exactly norm 2^e where e = length * f (where f = TORSION_PLUS_EVEN_POWER)
 * when used for compressing an isogeny, we assume that the curve given in input is normalized
 */
int id2iso_ideal_to_isogeny_two_long_power_of_2(id2iso_compressed_long_two_isog_t *isog_zip, ec_curve_t *domain, ec_basis_t *basis_minus, ec_basis_t *basis_plus, ec_point_t *kernel_dual, const quat_alg_elem_t *gen_input, const int length, const quat_left_ideal_t *lideal_start_small, const quat_left_ideal_t *lideal_start_two, const quat_alg_elem_t *gen_two, const quat_alg_t *Bpoo);

/**
 * @brief Translating an ideal of odd norm dividing p²-1 into the corresponding isogeny
 *
 * @param isog Output : the output isogeny 
 * @param basis_minus : a basis of ec points
 * @param basis_plus : a basis of ec points
 * @param domain : an elliptic curve
 * @param lideal_input : O0-ideal corresponding to the ideal to be translated 
 *  
 * compute  the isogeny starting from domain corresponding to ideal_input 
 * the coefficients extracted from the ideal are to be applied to basis_minus and basis_plus to compute the kernel of the isogeny. 
 *
 */

void id2iso_ideal_to_isogeny_odd(ec_isog_odd_t *isog, const ec_curve_t *domain, const ec_basis_t *basis_plus,const ec_basis_t *basis_minus, const quat_left_ideal_t *lideal_input);

/**
 * @brief Translating an ideal of norm a power of two dividing p²-1 into the corresponding isogeny
 *
 * @param isog Output : the output isogeny 
 * @param lideal_input : O0-ideal corresponding to the ideal to be translated 
 *  
 */

void id2iso_ideal_to_isogeny_even(ec_isog_even_t *isog, const quat_left_ideal_t *lideal_input);


/**
 * @brief Translating a kernel on the curve E0, represented as two vectors with respect to the precomputed 2^f- and 3^e-torsion bases, into the corresponding O0-ideal
 *
 * @param lideal Output : the output O0-ideal
 * @param vec2 : length-2 vector giving the 2-power part of the kernel with respect to the precomputed TORSION_PLUS_2POWER basis
 * @param vec3 : length-2 vector giving the 3-power part of the kernel with respect to the precomputed TORSION_PLUS_3POWER basis
 *
 */
void id2iso_kernel_dlogs_to_ideal(quat_left_ideal_t *lideal, const ibz_vec_2_t *vec2, const ibz_vec_2_t *vec3);

/** @}
 */
/** @}
 */



#endif
