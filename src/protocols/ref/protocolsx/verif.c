#include <protocols.h>
#include <inttypes.h>

static inline void copy_point(ec_point_t* A, const ec_point_t* B){
    fp2_copy(&A->x, &B->x);
    fp2_copy(&A->z, &B->z);
}

//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.im[i]);
    printf("\n");
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}
//XXX

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}


void protocols_verif_unpack_chall(ec_curve_t *E2, ec_point_t *dual, const signature_t *sig, const public_key_t *pk)
{
    // Generate 2^f-basis
    ec_basis_t B2;
    ec_curve_to_basis_2(&B2, &pk->E);
    if(sig->zip.bit_first_step){
        // Swap P and Q
        ec_point_t tmp;
        copy_point(&tmp, &B2.P);
        copy_point(&B2.P, &B2.Q);
        copy_point(&B2.Q, &tmp);
    }

    // 2^f-isogeny chains
    ec_point_t K;
    digit_t scalar[NWORDS_FIELD];
    ec_curve_t E;
    ec_isog_even_t isog2;
    fp2_copy(&E.A, &pk->E.A);
    fp2_copy(&E.C, &pk->E.C);
    isog2.length = POWER_OF_2;
    for(size_t i = 0; i < sig->zip.length-1; i++){        
        ibz_to_digit_array(scalar, &sig->zip.zip_chain[i]);
        ec_ladder3pt(&K, scalar, &B2.Q, &B2.P, &B2.PmQ, &E);
        copy_point(&isog2.kernel, &K);
        fp2_copy(&isog2.curve.A, &E.A);
        fp2_copy(&isog2.curve.C, &E.C);
        ec_eval_even(&E, &isog2, &B2.P, 1);
        ec_complete_basis_2(&B2, &E, &B2.P);
    }   
    ibz_to_digit_array(scalar, &sig->zip.zip_chain[sig->zip.length-1]);
    ec_ladder3pt(&K, scalar, &B2.Q, &B2.P, &B2.PmQ, &E);
    copy_point(&isog2.kernel, &K);
    fp2_copy(&isog2.curve.A, &E.A);
    fp2_copy(&isog2.curve.C, &E.C);
    ec_eval_even(&E, &isog2, &B2.P, 1);

    copy_point(dual,&B2.P);
    // Normalize curve
    ec_isom_t isom;
    ec_curve_normalize(E2, &isom, &E);
    ec_iso_eval(dual,&isom);

}

int protocols_verif_from_chall(const signature_t *sig, ec_curve_t const *E2, const ec_point_t *dual, const unsigned char* m, size_t l)
{   

    // Generate 2^f and 3^g bases
    ec_basis_t B2, B3, B6;
    ec_curve_to_basis_6(&B6, E2);
    copy_point(&B3.P, &B6.P);
    copy_point(&B3.Q, &B6.Q);
    copy_point(&B3.PmQ, &B6.PmQ);
    for(size_t i = 0; i < POWER_OF_2; i++){
        ec_dbl(&B3.P, E2, &B3.P);
        ec_dbl(&B3.Q, E2, &B3.Q);
        ec_dbl(&B3.PmQ, E2, &B3.PmQ);
    }
    copy_point(&B2.P, &B6.P);
    copy_point(&B2.Q, &B6.Q);
    copy_point(&B2.PmQ, &B6.PmQ);
    for(size_t i = 0; i < POWER_OF_3; i++){
        ec_point_t tmp;
        ec_dbl(&tmp, E2, &B2.P);
        ec_add(&B2.P, &tmp, &B2.P, &B2.P);
        ec_dbl(&tmp, E2, &B2.Q);
        ec_add(&B2.Q, &tmp, &B2.Q, &B2.Q);
        ec_dbl(&tmp, E2, &B2.PmQ);
        ec_add(&B2.PmQ, &tmp, &B2.PmQ, &B2.PmQ);
    }

    bool const bit2 = sig->s.select23 & 1;
    bool const bit3 = sig->s.select23 & 2;

    if(!bit2){
        // Swap P2 and Q2
        ec_point_t tmp;
        copy_point(&tmp, &B2.P);
        copy_point(&B2.P, &B2.Q);
        copy_point(&B2.Q, &tmp);
    }
    if(!bit3){
        // Swap P3 and Q3
        ec_point_t tmp;
        copy_point(&tmp, &B3.P);
        copy_point(&B3.P, &B3.Q);
        copy_point(&B3.Q, &tmp);
    }

    // Kernels for the dual of the challenge isogeny
    ec_point_t K2, K3;
    ec_ladder3pt(&K2, sig->s.scalar2, &B2.P, &B2.Q, &B2.PmQ, E2);
    ec_ladder3pt(&K3, sig->s.scalar3, &B3.P, &B3.Q, &B3.PmQ, E2);

    // Evaluating the dual isogeny
    ec_isog_odd_t isog3={0};
    ec_point_t push_points[2];
    if (bit2 == bit3) {
        push_points[0] = bit2 ? B6.Q : B6.P;
    }
    else {
        digit_t scalars[2][NWORDS_FIELD];
        ibz_to_digit_array(scalars[bit3], &TORSION_PLUS_2POWER);
        ibz_to_digit_array(scalars[bit2], &TORSION_PLUS_3POWER);
        ec_biscalar_mul(&push_points[0], E2, scalars[0], scalars[1], &B6);
    }
    copy_point(&push_points[1], &K3);
    ec_isog_even_t isog2;
    isog2.length = POWER_OF_2;
    copy_point(&isog2.kernel, &K2);


    // entering the cyclicity test zone
    ec_point_t test1,test2;
    test1 = *dual;
    test2 = K2;
    for (int i =0 ; i < TORSION_PLUS_EVEN_POWER-1; i++) {
        ec_dbl(&test1,E2,&test1);
        ec_dbl(&test2,E2,&test2);
    }
    assert(!fp2_is_zero(&test1.z));
    assert(!fp2_is_zero(&test2.z));
    if (ec_is_equal(&test1,&test2)) {
        return 0;
    }      

    fp2_copy(&isog2.curve.A, &E2->A);
    fp2_copy(&isog2.curve.C, &E2->C);
    ec_curve_t E;
    ec_eval_even(&E, &isog2, push_points, 2);
    isog3.degree[0] = POWER_OF_3;
    copy_point(&isog3.ker_plus, &push_points[1]);
    fp2_copy(&isog3.curve.A, &E.A);
    fp2_copy(&isog3.curve.C, &E.C);
    ec_eval_odd(&E, &isog3, push_points, 1);
assert(!fp2_is_zero(&E.C));

    // Normalize curve
    ec_curve_t E1;
    ec_isom_t isom;
    ec_curve_normalize(&E1, &isom, &E);
    ec_iso_eval(&push_points[0], &isom);
assert(!fp2_is_zero(&E1.C));

    // Generate 2^f*3^g basis
    ec_curve_to_basis_6(&B6, &E1);

    // Recover original challenge kernel point from the hash
    ec_point_t K;
    ibz_vec_2_t vec;
    ibz_vec_2_init(&vec);
    hash_to_challenge(&vec, &E1, m, l);
    digit_t scalars[2][NWORDS_FIELD];
    ibz_to_digit_array(scalars[0], &vec[0]);
    ibz_to_digit_array(scalars[1], &vec[1]);
    ec_biscalar_mul(&K, &E1, scalars[0], scalars[1], &B6);
    ibz_vec_2_finalize(&vec);
assert(!fp2_is_zero(&K.z));

    // Multiply Q by r
    ec_mul(&push_points[0], &E1, sig->r, &push_points[0]);
assert(!fp2_is_zero(&push_points[0].z));

    return ec_is_equal(&K, &push_points[0]);
}


/** @defgroup verification The verification protocol
 * @{
*/

/**
 * @brief Verifying a signature
 *
 * @param sig: the signature
 * @param pk the public key 
 * @param m the message
 * @returns a bit indicating if the verification succeeded  
    */
int protocols_verif(const signature_t *sig, const public_key_t *pk, const unsigned char* m, size_t l)
{
    //TODO general input checks?

    ec_curve_t E2;
    ec_point_t dual;
    protocols_verif_unpack_chall(&E2, &dual, sig, pk);

    return protocols_verif_from_chall(sig, &E2, &dual, m, l);
}

