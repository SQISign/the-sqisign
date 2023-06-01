/** @file
 * 
 * @authors Antonin Leroux
 * 
 * @brief The protocols
 */

#ifndef PROTOCOLS_H
#define PROTOCOLS_H

#include <intbig.h>
#include <quaternion.h>
#include <klpt_constants.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <rng.h>
#include <endomorphism_action.h>
#include <encoded_sizes.h>
#include <fp_constants.h>

#include <stdio.h>
#include <klpt.h>
#include <id2iso.h>

/** @defgroup protocols_protocols SQIsign protocols
 * @{
*/
/** @defgroup protocols_t Types for sqisign protocols
 * @{
*/

/** @brief Type for the signature
 * 
 * @typedef signature_t
 * 
 * @struct signature
 * 
*/


typedef struct signature {
    id2iso_compressed_long_two_isog_t zip;  /// the compressed isogeny 
    digit_t r[NWORDS_FIELD];          /// first scalar encoding the challenge
    struct {
        unsigned char select23; // this&1 is bit2, this&2 is bit3
        digit_t scalar2[NWORDS_FIELD];
        digit_t scalar3[NWORDS_FIELD];
    } s;                                    /// second scalar encoding the challenge
} signature_t;

/** @brief Type for the public keys
 * 
 * @typedef public_key_t
 * 
 * @struct public_key
 * 
*/
typedef struct public_key {
	ec_curve_t E; /// the normalized A coefficient of the Montgomery curve
} public_key_t;

/** @brief Type for the secret keys
 * 
 * @typedef secret_key_t
 * 
 * @struct secret_key
 * 
*/
typedef struct secret_key {
	ec_curve_t curve; /// the j-invariant
    quat_left_ideal_t lideal_small; /// ideal of prime norm
    quat_left_ideal_t lideal_two; /// ideal of norm a power of 2
    quat_alg_elem_t gen_two; /// generator of lideal_two
    ec_point_t kernel_dual; /// kernel of the dual of the last step of the phi_two  
    ec_basis_t basis_plus; /// basis pushed through phi_two
    ec_basis_t basis_minus; /// basis pushed through phi_two
} secret_key_t;


/** @}
*/


/*************************** Functions *****************************/


void secret_key_init(secret_key_t *sk);
void secret_key_finalize(secret_key_t *sk);

void signature_init(signature_t *sig);
void signature_finalize(signature_t *sig);

/** @defgroup signature The signature protocol
 * @{
*/

/**
 * @brief Computing a key pair
 *
 * @param pk : Output the public key 
 * @param sk : Output the secret key
 * @returns a bit indicating if the computation succeeded  
    */
int protocols_keygen(public_key_t *pk, secret_key_t *sk);


/**
 * @brief Computing a commitment isogeny from the starting curve
 *
 * @param ideal Output: the left O0-ideal defining the commitment isogeny
 * @param E1 Output: the codomain curve, normalized à la ec_curve_normalize()
 * @param basis Input and output: a basis on E0 which will be pushed through the commitment i
sogeny
 */
void protocols_commit(quat_left_ideal_t *ideal, ec_curve_t *E1, ec_basis_t *basis);

/**
 * @brief Hashing a commitment curve and a message to a challenge
 *
 * @param scalars Output: scalars defining the challenge kernel in terms of a deterministic t
orsion basis of the commitment curve
 * @param curve: the commitment curve
 * @param message: byte string to be signed
 * @param length: length of the message
 */
void hash_to_challenge(ibz_vec_2_t *scalars, const ec_curve_t *curve, const unsigned char *message, size_t length);


/**
 * @brief Computing a challenge isogeny from a commitment curve
 *
 * @param ideal Output: the left O0-ideal defining the challenge isogeny (pulled back through
 the commitment)
 * @param sig Output: signature structure in which the members r and s will be filled by this
 function
 * @param E1: commitment curve
 * @param pushedbasis6: image of the precomputed fixed torsion basis #BASIS_CHALLENGE under t
he commitment isogeny
 * @param hash: message hash according to hash_to_challenge()
 * @param out_E2: optional output (can be NULL) to store the value of the challenge curve for
 testing and debugging
 */
void protocols_challenge(quat_left_ideal_t *ideal, signature_t *sig, const ec_curve_t *E1, const ec_basis_t *pushedbasis6, const ibz_vec_2_t *hash, ec_curve_t *out_E2);


/**
 * @brief Computing a signature
 *
 * @param sig Output: the signature
 * @param pk the public key 
 * @param sk the secret key
 * @param m the message
 * @param l the message length
 * @returns a bit indicating if the computation succeeded  
    */
int protocols_sign(signature_t *sig,const public_key_t *pk, const secret_key_t *sk,const unsigned char* m, size_t l);


/** @}
*/


/** @defgroup verification The verification protocol
 * @{
*/


/**
 * @brief Recovering the challenge curve from a signature
 *
 * @param E2 Output: the challenge curve, will be normalized à la ec_curve_normalize()
 * @param dual Output: point defining the first step of the dual of the response isogeny
 * @param sig the signature
 * @param pk the public key
 */
void protocols_verif_unpack_chall(ec_curve_t *E2, ec_point_t *dual, const signature_t *sig, const public_key_t *pk);

/**
 * @brief Verifying the challenge from a signature
 *
 * @param sig the signature
 * @param E2 the challenge curve
 * @param dual point defining the first step of the dual of the response isogeny
 * @param m the message
 * @param l the message length
 * @returns a single bit indicating whether the given data forms correct signature data
 */
int protocols_verif_from_chall(const signature_t *sig, ec_curve_t const *E2, const ec_point_t *dual, const unsigned char* m, size_t l);


/**
 * @brief Verifying a signature
 *
 * @param sig: the signature
 * @param pk the public key 
 * @param m the message
 * @param l length of the message
 * @returns a bit indicating if the verification succeeded  
    */
int protocols_verif(const signature_t *sig,const public_key_t *pk,const unsigned char* m, size_t l);


/** @}
*/


/*************************** Encoding *****************************/


/** @defgroup encoding Encoding and decoding functions
 * @{
*/

/**
 * @brief Encodes a public key as a byte array
 *
 * @param enc : Output the encoded public key
 * @param pk : public key
    */
void public_key_encode(unsigned char *enc, const public_key_t* pk);

/**
 * @brief Encodes a secret key as a byte array
 *
 * @param enc : Output the encoded secret key (including public key)
 * @param sk : secret key
 * @param pk : public key
  */
void secret_key_encode(unsigned char *enc, const secret_key_t* sk, const public_key_t* pk);

/**
 * @brief Encodes a signature as a byte array
 *
 * @param enc : Output the encoded signature
 * @param sig : signature
    */
void signature_encode(unsigned char* enc, const signature_t* sig);

/**
 * @brief Encodes a public key as a byte array
 *
 * @param pk : Output the decoded public key
 * @param enc : encoded public key
    */
void public_key_decode(public_key_t* pk, const unsigned char *enc);

/**
 * @brief Encodes a secret key (and public key) as a byte array
 *
 * @param sk : Output the decoded secret key
 * @param pk : Output the decoded public key contained in enc
 * @param enc : encoded secret key
    */
void secret_key_decode(secret_key_t* sk, public_key_t* pk, const unsigned char *enc);

/**
 * @brief Encodes a signature as a byte array
 *
 * @param sig : Output the decoded signature
 * @param enc : encoded signature
    */
void signature_decode(signature_t* sig, const unsigned char* enc);

/** @}
*/
/** @}
*/

#endif
