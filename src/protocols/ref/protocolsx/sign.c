#include <protocols.h>


void signature_init(signature_t *sig)
{
    id2iso_compressed_long_two_isog_init(&sig->zip, SQISIGN_signing_length);
}

void signature_finalize(signature_t *sig)
{
    id2iso_compressed_long_two_isog_finalize(&sig->zip);
}


void protocols_commit(quat_left_ideal_t *ideal, ec_curve_t *E1, ec_basis_t *basis)
{
    ibz_vec_2_t vec;
    ibz_vec_2_init(&vec);

    ibz_t tmp;
    ibz_init(&tmp);

    // sample the random scalars
    do {
        ibz_rand_interval(&vec[0], &ibz_const_one, &DEGREE_COMMITMENT);
        ibz_rand_interval(&vec[1], &ibz_const_one, &DEGREE_COMMITMENT);
        ibz_gcd(&tmp, &vec[0], &vec[1]);
        ibz_gcd(&tmp, &tmp, &DEGREE_COMMITMENT);
    } while (!ibz_is_one(&tmp));

    // deduce the commitment ideal
    {
        quat_alg_elem_t gen;
        quat_alg_elem_init(&gen);

        ibz_set(&gen.denom, 1);
        for (unsigned i = 0; i < 4; ++i) {
            ibz_mul(&gen.coord[i], &COMMITMENT_IDEAL_DISTORTION_ENDO.coord[i], &vec[1]);
            if (i)
                ibz_neg(&gen.coord[i], &gen.coord[i]);  // conjugate
        }
        ibz_add(&gen.coord[0], &gen.coord[0], &vec[0]);
        // now gen = a + b*dual(theta) where vec=(a,b) and theta is the distortion map

        quat_alg_mul(&gen, &COMMITMENT_IDEAL_UNDISTORTED_GEN, &gen, &QUATALG_PINFTY);

        #ifndef NDEBUG
            quat_alg_coord_t coeffs;
            ibz_t temp,remainder;
            quat_alg_coord_init(&coeffs);
            ibz_init(&temp);ibz_init(&remainder);
            quat_alg_make_primitive(&coeffs,&temp,&gen,&MAXORD_O0,&QUATALG_PINFTY);
            ibz_mul(&gen.denom,&gen.denom,&temp);
            quat_alg_normalize(&gen);
            assert(quat_alg_is_primitive(&gen,&MAXORD_O0,&QUATALG_PINFTY));
            ibq_t ibq_norm;
            ibq_init(&ibq_norm);
            quat_alg_norm(&ibq_norm,&gen,&QUATALG_PINFTY);
            assert(ibq_to_ibz(&temp,&ibq_norm));
            ibq_finalize(&ibq_norm);
            ibz_gcd(&temp,&temp,&DEGREE_COMMITMENT);
            ibz_div(&temp,&remainder,&temp,&DEGREE_COMMITMENT);
            assert(0==ibz_cmp(&remainder,&ibz_const_zero));
            ibz_finalize(&temp);ibz_finalize(&remainder);
            quat_alg_coord_finalize(&coeffs);
        #endif


        quat_lideal_make_primitive_then_create(ideal, &gen, &DEGREE_COMMITMENT, &MAXORD_O0, &QUATALG_PINFTY);
        assert(!ibz_cmp(&ideal->norm, &DEGREE_COMMITMENT));

        quat_alg_elem_finalize(&gen);
    }

    // deduce the commitment isogeny
    ec_isog_odd_t isogeny;
    {
        isogeny.curve = CURVE_E0;
        for (size_t i = 0; i < sizeof(TORSION_ODD_PRIMES)/sizeof(*TORSION_ODD_PRIMES); ++i)
            isogeny.degree[i] = DEGREE_COMMITMENT_POWERS[i];

        digit_t scalars_plus[2][NWORDS_FIELD], scalars_minus[2][NWORDS_FIELD];
        ibz_mod(&tmp, &vec[0], &DEGREE_COMMITMENT_PLUS);
        ibz_to_digit_array(scalars_plus[0], &tmp);
        ibz_mod(&tmp, &vec[1], &DEGREE_COMMITMENT_PLUS);
        ibz_to_digit_array(scalars_plus[1], &tmp);
        ibz_mod(&tmp, &vec[0], &DEGREE_COMMITMENT_MINUS);
        ibz_to_digit_array(scalars_minus[0], &tmp);
        ibz_mod(&tmp, &vec[1], &DEGREE_COMMITMENT_MINUS);
        ibz_to_digit_array(scalars_minus[1], &tmp);

        ec_biscalar_mul(&isogeny.ker_plus, &CURVE_E0, scalars_plus[0], scalars_plus[1], &BASIS_COMMITMENT_PLUS);
        ec_biscalar_mul(&isogeny.ker_minus, &CURVE_E0, scalars_minus[0], scalars_minus[1], &BASIS_COMMITMENT_MINUS);
    }

    ec_curve_t E1comp;
    ec_eval_odd_basis(&E1comp, &isogeny, basis, 1);

    // normalize E1 with the pushed basis
    {
        ec_isom_t isom;
        ec_curve_normalize(E1, &isom, &E1comp);
        ec_iso_eval(&basis->P, &isom);
        ec_iso_eval(&basis->Q, &isom);
        ec_iso_eval(&basis->PmQ, &isom);
    }

    ibz_finalize(&tmp);
    ibz_vec_2_finalize(&vec);
}

void protocols_challenge(quat_left_ideal_t *ideal, signature_t *sig, const ec_curve_t *E1, const ec_basis_t *pushedbasis6, const ibz_vec_2_t *hash, ec_curve_t *out_E2)
{
    // compute deterministic and pushed 2*3*-torsion bases on E1
    ec_basis_t E1basis6, E1basis2, E1basis3;
    ec_basis_t pushedbasis2, pushedbasis3;
    ec_curve_to_basis_6(&E1basis6, E1);
    {
        digit_t scalar[NWORDS_FIELD];
        ibz_to_digit_array(scalar, &TORSION_PLUS_3POWER);
        ec_mul(&E1basis2.P, E1, scalar, &E1basis6.P);
        ec_mul(&E1basis2.Q, E1, scalar, &E1basis6.Q);
        ec_mul(&E1basis2.PmQ, E1, scalar, &E1basis6.PmQ);
        ec_mul(&pushedbasis2.P, E1, scalar, &pushedbasis6->P);
        ec_mul(&pushedbasis2.Q, E1, scalar, &pushedbasis6->Q);
        ec_mul(&pushedbasis2.PmQ, E1, scalar, &pushedbasis6->PmQ);
        ibz_to_digit_array(scalar, &TORSION_PLUS_2POWER);
        ec_mul(&E1basis3.P, E1, scalar, &E1basis6.P);
        ec_mul(&E1basis3.Q, E1, scalar, &E1basis6.Q);
        ec_mul(&E1basis3.PmQ, E1, scalar, &E1basis6.PmQ);
        ec_mul(&pushedbasis3.P, E1, scalar, &pushedbasis6->P);
        ec_mul(&pushedbasis3.Q, E1, scalar, &pushedbasis6->Q);
        ec_mul(&pushedbasis3.PmQ, E1, scalar, &pushedbasis6->PmQ);
    }

    // compute the kernel of the challenge isogeny
    ec_point_t ker;
    {
        digit_t scalars[2][NWORDS_FIELD];
        ibz_to_digit_array(scalars[0], &(*hash)[0]);
        ibz_to_digit_array(scalars[1], &(*hash)[1]);
        ec_biscalar_mul(&ker, E1, scalars[0], scalars[1], &E1basis6);
    }
    ec_point_t ker2, ker3;
    {
        digit_t scalar[NWORDS_FIELD];
        ibz_to_digit_array(scalar, &TORSION_PLUS_3POWER);
        ec_mul(&ker2, E1, scalar, &ker);
        ibz_to_digit_array(scalar, &TORSION_PLUS_2POWER);
        ec_mul(&ker3, E1, scalar, &ker);
    }

    // compute ideal corresponding to the challenge isogeny, pulled back to O0
    {
        // compute the logarithm of the kernel with respect to the pushed basis
        ibz_vec_2_t vec2, vec3;
        ibz_vec_2_init(&vec2);
        ibz_vec_2_init(&vec3);
        digit_t log[2][NWORDS_FIELD];

        ec_dlog_2(log[0], log[1], &pushedbasis2, &ker2, E1);

        ibz_copy_digit_array(&vec2[0], log[0]);
        ibz_copy_digit_array(&vec2[1], log[1]);

        ec_dlog_3(log[0], log[1], &pushedbasis3, &ker3, E1);

        ibz_copy_digit_array(&vec3[0], log[0]);
        ibz_copy_digit_array(&vec3[1], log[1]);

        // now compute the ideal
        id2iso_kernel_dlogs_to_ideal(ideal, &vec2, &vec3);

        ibz_vec_2_finalize(&vec2);
        ibz_vec_2_finalize(&vec3);
        assert(!ibz_cmp(&ideal->norm, &DEGREE_CHALLENGE));
    }

    // compute the isogeny, evaluate at suitable point to find kernel of dual
    ec_curve_t Emid, E2comp, E2;
    ec_point_t pts[2];
    ec_isom_t E2isom;
    {
        ec_isog_even_t isog2 = {.curve = *E1, .length = TORSION_PLUS_EVEN_POWER, .kernel = ker2};

        // find a point independent to the kernel to push it through
        {
            bool const bit2 = !ibz_divides(&(*hash)[0], &ibz_const_two);
            bool const bit3 = !ibz_divides(&(*hash)[0], &ibz_const_three);
            if (bit2 == bit3) {
                pts[0] = bit2 ? E1basis6.Q : E1basis6.P;
            }
            else {
                digit_t scalars[2][NWORDS_FIELD];
                ibz_to_digit_array(scalars[bit3], &TORSION_PLUS_2POWER);
                ibz_to_digit_array(scalars[bit2], &TORSION_PLUS_3POWER);
                ec_biscalar_mul(&pts[0], E1, scalars[0], scalars[1], &E1basis6);
            }
        }
        pts[1] = ker;

        ec_eval_even(&Emid, &isog2, pts, 2);
        assert(!fp2_is_zero(&Emid.C));

        assert(TORSION_ODD_PRIMES[0] == 3);
        ec_isog_odd_t isog3 = {.curve = Emid, .degree = {TORSION_ODD_POWERS[0]}, .ker_plus = pts[1]};

        ec_eval_odd(&E2comp, &isog3, pts, 1);
        assert(!fp2_is_zero(&E2comp.C));

        ec_curve_normalize(&E2, &E2isom, &E2comp);
        if (out_E2) *out_E2 = E2;

        ec_iso_eval(&pts[0], &E2isom);
    }
    // now the kernel of the dual is in pts[0]

    // decompose kernel of dual over deterministic 2*3*-torsion basis on E2
    ibz_vec_2_t vec2, vec3;
    ibz_vec_2_init(&vec2);
    ibz_vec_2_init(&vec3);
    ec_basis_t E2basis6;
    {
        ec_basis_t E2basis2, E2basis3;
        ec_curve_to_basis_6(&E2basis6, &E2);

        // set up the 2-power and 3-power bases as images of the 2*3*-torsion one
        {
            digit_t scalar[NWORDS_FIELD];
            ibz_to_digit_array(scalar, &TORSION_PLUS_3POWER);
            ec_mul(&E2basis2.P, &E2, scalar, &E2basis6.P);
            ec_mul(&E2basis2.Q, &E2, scalar, &E2basis6.Q);
            ec_mul(&E2basis2.PmQ, &E2, scalar, &E2basis6.PmQ);
            ibz_to_digit_array(scalar, &TORSION_PLUS_2POWER);
            ec_mul(&E2basis3.P, &E2, scalar, &E2basis6.P);
            ec_mul(&E2basis3.Q, &E2, scalar, &E2basis6.Q);
            ec_mul(&E2basis3.PmQ, &E2, scalar, &E2basis6.PmQ);
        }

        // compute the 2-power and 3-power parts of the kernel of the dual
        {
            digit_t scalar[NWORDS_FIELD];
            ibz_to_digit_array(scalar, &TORSION_PLUS_2POWER);
            ec_mul(&pts[1], &E2, scalar, &pts[0]);
            ibz_to_digit_array(scalar, &TORSION_PLUS_3POWER);
            ec_mul(&pts[0], &E2, scalar, &pts[0]);
        }

        // now compute the logarithms
        {
            digit_t scalars[2][NWORDS_FIELD];

            ec_dlog_2(scalars[0], scalars[1], &E2basis2, &pts[0], &E2);
            ibz_copy_digit_array(&vec2[0], scalars[0]);
            ibz_copy_digit_array(&vec2[1], scalars[1]);
//ibz_vec_2_print(&vec2);
#ifndef NDEBUG
{
    ec_point_t tmp;
    ec_biscalar_mul(&tmp, &E2, scalars[0], scalars[1], &E2basis2);
    assert(ec_is_equal(&tmp, &pts[0]));
}
#endif

            ec_dlog_3(scalars[0], scalars[1], &E2basis3, &pts[1], &E2);
            ibz_copy_digit_array(&vec3[0], scalars[0]);
            ibz_copy_digit_array(&vec3[1], scalars[1]);
//ibz_vec_2_print(&vec3);
#ifndef NDEBUG
{
    ec_point_t tmp;
    ec_biscalar_mul(&tmp, &E2, scalars[0], scalars[1], &E2basis3);
    assert(ec_is_equal(&tmp, &pts[1]));
}
#endif
        }
    }

    // encode this projective vector in the signature field s
    bool const bit2 = !ibz_divides(&vec2[0], &ibz_const_two);
    bool const bit3 = !ibz_divides(&vec3[0], &ibz_const_three);
    {
        assert(!ibz_divides(&vec2[!bit2], &ibz_const_two));
        assert(!ibz_divides(&vec3[!bit3], &ibz_const_three));

        {
            ibz_t inv;
            ibz_init(&inv);

            ibz_invmod(&inv, &vec2[!bit2], &TORSION_PLUS_2POWER);
            ibz_mul(&vec2[bit2], &vec2[bit2], &inv);
            ibz_mod(&vec2[bit2], &vec2[bit2], &TORSION_PLUS_2POWER);
            ibz_set(&vec2[!bit2], 1);

            ibz_invmod(&inv, &vec3[!bit3], &TORSION_PLUS_3POWER);
            ibz_mul(&vec3[bit3], &vec3[bit3], &inv);
            ibz_mod(&vec3[bit3], &vec3[bit3], &TORSION_PLUS_3POWER);
            ibz_set(&vec3[!bit3], 1);

            ibz_finalize(&inv);
        }
//printf("bit2=%u bit3=%u\n", (unsigned) bit2, (unsigned) bit3);
//ibz_vec_2_print(&vec2);
//ibz_vec_2_print(&vec3);

        sig->s.select23 = bit2 | ((unsigned) bit3 << 1);
        ibz_to_digit_array(sig->s.scalar2, &vec2[bit2]);
        ibz_to_digit_array(sig->s.scalar3, &vec3[bit3]);
    }

    // find another independent point R deterministically
    ec_point_t R;
    if (bit2 == bit3) {
        R = bit2 ? E2basis6.Q : E2basis6.P;
    }
    else {
        digit_t scalars[2][NWORDS_FIELD];
        ibz_to_digit_array(scalars[bit3], &TORSION_PLUS_2POWER);
        ibz_to_digit_array(scalars[bit2], &TORSION_PLUS_3POWER);
        ec_biscalar_mul(&R, &E2, scalars[0], scalars[1], &E2basis6);
    }

    // evaluate the dual of the challenge isogeny
    {
        {
            ec_isom_t E2isom_inv = E2isom;
            ec_iso_inv(&E2isom_inv);
            ec_iso_eval(&pts[0], &E2isom_inv);
            ec_iso_eval(&pts[1], &E2isom_inv);
            ec_iso_eval(&R, &E2isom_inv);
        }

        assert(TORSION_ODD_PRIMES[0] == 3);
        ec_isog_odd_t isog3dual = {.curve = E2comp, .degree = {TORSION_ODD_POWERS[0]}, .ker_plus = pts[1]};
        ec_curve_t Emidcomp;
        pts[1] = R;
        ec_eval_odd(&Emidcomp, &isog3dual, pts, 2);
        assert(!fp2_is_zero(&Emidcomp.C));
#ifndef NDEBUG
{
    fp2_t j1, j2;
    ec_j_inv(&j1, &Emid);
    ec_j_inv(&j2, &Emidcomp);
    assert(fp2_is_equal(&j1, &j2));
}
#endif

        R = pts[1];

        ec_isog_even_t isog2dual = {.curve = Emidcomp, .length = TORSION_PLUS_EVEN_POWER, .kernel = pts[0]};
        ec_curve_t E1comp;
        ec_eval_even(&E1comp, &isog2dual, &R, 1);
        assert(!fp2_is_zero(&E1comp.C));
#ifndef NDEBUG
{
    fp2_t j1, j2;
    ec_j_inv(&j1, E1);
    ec_j_inv(&j2, &E1comp);
    assert(fp2_is_equal(&j1, &j2));
}
#endif

        {
            ec_isom_t E1isom;
            ec_isomorphism(&E1isom, &E1comp, E1);
            ec_iso_eval(&R, &E1isom);
        }
    }

    ibz_t r;
    ibz_init(&r);
    //TODO this should just be a 1-dimensional logarithm in the 2*3*-torsion
    {
        ec_point_t R2, R3;
        {
            digit_t scalar[NWORDS_FIELD];
            ibz_to_digit_array(scalar, &TORSION_PLUS_3POWER);
            ec_mul(&R2, E1, scalar, &R);
            ibz_to_digit_array(scalar, &TORSION_PLUS_2POWER);
            ec_mul(&R3, E1, scalar, &R);
        }

        ibz_vec_2_t log2, log3;
        ibz_vec_2_init(&log2);
        ibz_vec_2_init(&log3);
        {
            digit_t scalars[2][NWORDS_FIELD] = {0};
            ec_dlog_2(scalars[0], scalars[1], &E1basis2, &R2, E1);
#ifndef NDEBUG
{
    ec_point_t tmp;
    ec_biscalar_mul(&tmp, E1, scalars[0], scalars[1], &E1basis2);
    assert(ec_is_equal(&tmp, &R2));
}
#endif
            ibz_copy_digit_array(&log2[0], scalars[0]);
            ibz_copy_digit_array(&log2[1], scalars[1]);
            memset(scalars, 0, sizeof(scalars));
            ec_dlog_3(scalars[0], scalars[1], &E1basis3, &R3, E1);
#ifndef NDEBUG
{
    ec_point_t tmp;
    ec_biscalar_mul(&tmp, E1, scalars[0], scalars[1], &E1basis3);
    assert(ec_is_equal(&tmp, &R3));
}
#endif
            ibz_copy_digit_array(&log3[0], scalars[0]);
            ibz_copy_digit_array(&log3[1], scalars[1]);
#ifndef NDEBUG
{
ibz_t lhs, rhs;
ibz_init(&lhs);
ibz_init(&rhs);
// check log2 is a multiple of hash modulo 2*
ibz_mul(&lhs, &log2[0], &(*hash)[1]);
ibz_mul(&rhs, &log2[1], &(*hash)[0]);
ibz_mod(&lhs, &lhs, &TORSION_PLUS_2POWER);
ibz_mod(&rhs, &rhs, &TORSION_PLUS_2POWER);
assert(!ibz_cmp(&lhs, &rhs));
// check log3 is a multiple of hash modulo 3*
ibz_mul(&lhs, &log3[0], &(*hash)[1]);
ibz_mul(&rhs, &log3[1], &(*hash)[0]);
ibz_mod(&lhs, &lhs, &TORSION_PLUS_3POWER);
ibz_mod(&rhs, &rhs, &TORSION_PLUS_3POWER);
assert(!ibz_cmp(&lhs, &rhs));
ibz_finalize(&lhs);
ibz_finalize(&rhs);
}
#endif

            {
                ibz_t r2, r3;
                ibz_init(&r2);
                ibz_init(&r3);
                bool bit = ibz_divides(&log2[0], &ibz_const_two);
                assert(!ibz_divides(&log2[bit], &ibz_const_two));
                ibz_invmod(&r2, &log2[bit], &TORSION_PLUS_2POWER);
                ibz_mul(&r2, &(*hash)[bit], &r2);
                ibz_mod(&r2, &r2, &TORSION_PLUS_2POWER);
//ibz_printf("r2 = %#Zx\n", &r2);
#ifndef NDEBUG
{
    digit_t scalar[NWORDS_FIELD];
    ibz_to_digit_array(scalar, &r2);
    ec_point_t T1, T2;
    ec_mul(&T1, E1, scalar, &R);
    ibz_to_digit_array(scalar, &TORSION_PLUS_3POWER);
    ec_mul(&T1, E1, scalar, &T1);
    ec_mul(&T2, E1, scalar, &ker);
    assert(ec_is_equal(&T1, &T2));
}
#endif
                bit = ibz_divides(&log3[0], &ibz_const_three);
                assert(!ibz_divides(&log3[bit], &ibz_const_three));
                ibz_invmod(&r3, &log3[bit], &TORSION_PLUS_3POWER);
                ibz_mul(&r3, &(*hash)[bit], &r3);
                ibz_mod(&r3, &r3, &TORSION_PLUS_3POWER);
//ibz_printf("r3 = %#Zx\n", &r3);
#ifndef NDEBUG
{
    digit_t scalar[NWORDS_FIELD];
    ibz_to_digit_array(scalar, &r3);
    ec_point_t T1, T2;
    ec_mul(&T1, E1, scalar, &R);
    ibz_to_digit_array(scalar, &TORSION_PLUS_2POWER);
    ec_mul(&T1, E1, scalar, &T1);
    ec_mul(&T2, E1, scalar, &ker);
    assert(ec_is_equal(&T1, &T2));
}
#endif
                ibz_crt(&r, &r2, &r3, &TORSION_PLUS_2POWER, &TORSION_PLUS_3POWER);
                // check if we mixed up the signs and correct it in that case
                {
                    digit_t scalar[NWORDS_FIELD];
                    ibz_to_digit_array(scalar, &r);
                    ec_point_t T;
                    ec_mul(&T, E1, scalar, &R);
                    if (!ec_is_equal(&T, &ker)) {
                        ibz_sub(&r3, &TORSION_PLUS_3POWER, &r3);
                        ibz_crt(&r, &r2, &r3, &TORSION_PLUS_2POWER, &TORSION_PLUS_3POWER);
                    }
                }
                ibz_finalize(&r2);
                ibz_finalize(&r3);
            }
        }

        ibz_vec_2_finalize(&log2);
        ibz_vec_2_finalize(&log3);
        ibz_vec_2_finalize(&vec2);
        ibz_vec_2_finalize(&vec3);
    }
#ifndef NDEBUG
{
    digit_t scalar[NWORDS_FIELD];
    ibz_to_digit_array(scalar, &r);
    ec_point_t T;
    ec_mul(&T, E1, scalar, &R);
    assert(ec_is_equal(&T, &ker));
assert(!ibz_divides(&r, &ibz_const_two));
assert(!ibz_divides(&r, &ibz_const_three));
}
#endif

    ibz_to_digit_array(sig->r, &r);

    ibz_finalize(&r);
}



int protocols_sign(signature_t *sig,const public_key_t *pk, const secret_key_t *sk, const unsigned char* m, size_t l) {

    // var dec
    quat_order_t right_order_key;
    quat_left_ideal_t ideal_commit,ideal_challenge;
    quat_alg_elem_t gen,gen_small,gen_answer,delta;
    ibq_t ibq_norm;
    ibz_t norm,temp,remainder;
    ibz_mat_4x4_t reduced,gram;
    quat_alg_coord_t coeffs;

    quat_left_ideal_t ideal_comchall,ideal_seccomchall;
    quat_left_ideal_t ideal_eichler_rand;
    quat_left_ideal_t ideal_pullback,ideal_input_klpt;
    quat_alg_elem_t gen_challenge;

    
    // var init
    ibq_init(&ibq_norm);
    ibz_init(&norm);ibz_init(&temp);
    ibz_init(&remainder);
    quat_order_init(&right_order_key);
    quat_left_ideal_init(&ideal_commit);
    quat_left_ideal_init(&ideal_challenge);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_small);
    quat_alg_elem_init(&delta);
    ibz_mat_4x4_init(&reduced);
    ibz_mat_4x4_init(&gram);
    quat_alg_coord_init(&coeffs);
    quat_left_ideal_init(&ideal_comchall);
    quat_left_ideal_init(&ideal_seccomchall);
    quat_left_ideal_init(&ideal_eichler_rand);
    quat_left_ideal_init(&ideal_pullback);
    quat_left_ideal_init(&ideal_input_klpt);
    quat_alg_elem_init(&gen_challenge);


    ec_curve_t curve = sk->curve;
    ec_basis_t basis_plus = sk->basis_plus;
    ec_basis_t basis_minus = sk->basis_minus;
    ec_point_t kernel_dual = sk->kernel_dual;

    // Commitment (computes ideal_commit and pushes the 2*3*-torsion basis)
    ec_curve_t E1;
    ec_basis_t E1basis = BASIS_CHALLENGE;
    {
        protocols_commit(&ideal_commit, &E1, &E1basis);

        assert(!fp2_is_zero(&E1.C));
        assert(!fp2_is_zero(&E1basis.P.z) && !fp2_is_zero(&E1basis.Q.z) && !fp2_is_zero(&E1basis.PmQ.z));
    }
  

    // Challenge (computes ideal_challenge and sig->r, sig->s)
    ec_curve_t E2_for_testing;  //XXX
    {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);

        hash_to_challenge(&vec, &E1, m, l);

        protocols_challenge(&ideal_challenge, sig, &E1, &E1basis, &vec, &E2_for_testing);

        ibz_vec_2_finalize(&vec);
    }
    assert(ibz_get(&ideal_commit.norm)%2!=0);
        
    // ideal_comchall is the intersectin of ideal_commit and ideal_challenge
    quat_lideal_inter(&ideal_comchall,&ideal_challenge,&ideal_commit,&QUATALG_PINFTY);

    #ifndef NDEBUG 
        // debug the norm
        ibz_mul(&temp,&ideal_commit.norm,&ideal_challenge.norm);
        assert(0==ibz_cmp(&temp,&ideal_comchall.norm));
    #endif

    int lideal_generator_ok = quat_lideal_generator(&gen,&ideal_comchall,&QUATALG_PINFTY,0);
    assert(lideal_generator_ok);
    // computing the generator of ideal seccomchall as gen = sk.gen_two * gen
    quat_alg_elem_copy(&gen_challenge,&gen);
    quat_alg_norm(&ibq_norm,&gen_challenge,&QUATALG_PINFTY);
    ibq_to_ibz(&norm,&ibq_norm);

    // computing a good norm element for sk.lideal_small
    lideal_generator_ok = quat_lideal_generator_coprime(&gen_small,&sk->lideal_small,&norm,&QUATALG_PINFTY,0);
    assert(lideal_generator_ok);
    quat_alg_conj(&gen_small,&gen_small);
    quat_alg_normalize(&gen_small);

    // gen = gen_small * gen_challenge
    quat_alg_mul(&gen,&gen_small,&gen_challenge,&QUATALG_PINFTY);
    ibz_mul(&temp,&ideal_comchall.norm,&sk->lideal_small.norm);
    quat_alg_normalize(&gen);

    // computing of the right order of lideal_small
    quat_lideal_right_order(&right_order_key,&sk->lideal_small,&QUATALG_PINFTY);
    assert(quat_lattice_contains(&coeffs,&right_order_key,&gen,&QUATALG_PINFTY));

    // creation of ideal_seccomchall
    quat_lideal_make_primitive_then_create(&ideal_seccomchall,&gen,&temp,&right_order_key,&QUATALG_PINFTY);

    // copying a generator for later
    quat_alg_elem_copy(&gen_challenge,&gen);


    #ifndef NDEBUG 
        // debug the norm
        ibz_mul(&temp,&ideal_comchall.norm,&sk->lideal_small.norm);
        assert(0==ibz_cmp(&temp,&ideal_seccomchall.norm));
    #endif

    int found = 0;

    for (unsigned attempt = 0; !found && attempt < SQISIGN_response_attempts; ++attempt) { 
        // TODUPDATE add a random ideal of norm 2^SQISIGN_random_length (right now SQISIGN_random_length )

        // now we are computing a small ideal equivalent to ideal_seccomchall 
        //TODUPDATE replace with another randomization to ensure a good distribution in the eichler order class set   
        quat_lideal_reduce_basis(&reduced,&gram,&ideal_seccomchall,&QUATALG_PINFTY);
        found = klpt_lideal_equiv(&gen,&temp,&reduced,&gram,&ideal_seccomchall.norm,&ideal_seccomchall.lattice.denom,&QUATALG_PINFTY);
        if (!found) {
            continue;
        }
        // quat_alg_elem_copy(&gen_challenge,&gen);
        quat_alg_conj(&gen_challenge,&gen);

        // computing ideal_eichler_rand = right_order_key < gen, temp >
        quat_lideal_create_from_primitive(&ideal_eichler_rand,&gen,&temp,&right_order_key,&QUATALG_PINFTY);

        assert(quat_lideal_isom(&delta,&ideal_eichler_rand,&ideal_seccomchall,&QUATALG_PINFTY));

        // computing the pullback through sk->lideal_small 
        quat_alg_mul(&gen,&gen_small,&gen,&QUATALG_PINFTY);
        quat_alg_normalize(&gen);   
        quat_lideal_create_from_primitive(&ideal_pullback,&gen,&ideal_eichler_rand.norm,&MAXORD_O0,&QUATALG_PINFTY);

        // computing the final input ideal as equivalent to ideal_pullback
        quat_lideal_reduce_basis(&reduced,&gram,&ideal_pullback,&QUATALG_PINFTY);
        found = found && klpt_lideal_equiv(&gen,&temp,&reduced,&gram,&ideal_pullback.norm,&ideal_pullback.lattice.denom,&QUATALG_PINFTY);
        // ruling out failure, or extreme values of gen
        if (!found || 0==ibz_cmp(&gen.coord[0],&ibz_const_zero) || 0==ibz_cmp(&gen.coord[1],&ibz_const_zero)) {
            found = 0;
            continue;
            
        }

        // creating this ideal
        quat_lideal_create_from_primitive(&ideal_input_klpt,&gen,&temp,&MAXORD_O0,&QUATALG_PINFTY);
        assert(quat_lideal_isom(&delta,&ideal_input_klpt,&ideal_pullback,&QUATALG_PINFTY));

        // applying klpt 
        quat_alg_conj(&delta,&gen);
        #ifndef NDEBUG 
            assert(quat_lattice_contains(&coeffs,&ideal_pullback.lattice,&delta,&QUATALG_PINFTY));
        #endif 
    

        // applying signing klpt
        found = klpt_signing_klpt(&gen,&ideal_input_klpt,&sk->lideal_small,&delta,&QUATALG_PINFTY);
        if (!found) {
            continue;

        }
        assert(quat_lattice_contains(&coeffs,&ideal_input_klpt.lattice,&gen,&QUATALG_PINFTY));

        // gen = gen *delta /norm(ideal_input_klpt)
        quat_alg_mul(&gen,&gen,&delta,&QUATALG_PINFTY);
        ibz_mul(&gen.denom,&gen.denom,&ideal_input_klpt.norm);
        quat_alg_normalize(&gen);

        // for debug we check that gen is contained in the right order and ideal
        assert(quat_lattice_contains(&coeffs,&ideal_pullback.lattice,&gen,&QUATALG_PINFTY));
        assert(quat_lattice_contains(&coeffs,&right_order_key,&gen,&QUATALG_PINFTY));
    
        // gen = conjugate(gen) 
        // gen should be the generator of the ideal to be translated
        quat_alg_conj(&gen,&gen); 


        #ifndef NDEBUG 
            quat_left_ideal_t ideal_signing_test;
            quat_left_ideal_init(&ideal_signing_test);
            ibz_pow(&temp,&ibz_const_two,KLPT_signing_klpt_length);
            quat_lideal_create_from_primitive(&ideal_signing_test,&gen,&temp,&right_order_key,&QUATALG_PINFTY);
            assert(quat_lideal_isom(&delta,&ideal_signing_test,&ideal_seccomchall,&QUATALG_PINFTY));
            assert(quat_lideal_isom(&delta,&ideal_eichler_rand,&ideal_seccomchall,&QUATALG_PINFTY));
            assert(quat_lideal_isom(&delta,&ideal_eichler_rand,&ideal_seccomchall,&QUATALG_PINFTY));
            assert(0==ibz_cmp(&temp,&ideal_signing_test.norm));
            quat_alg_conj(&delta,&gen);
            quat_lideal_create_from_primitive(&ideal_signing_test,&delta,&ideal_eichler_rand.norm,&right_order_key,&QUATALG_PINFTY);
            assert(quat_lideal_equals(&ideal_signing_test,&ideal_eichler_rand,&QUATALG_PINFTY));
            quat_left_ideal_finalize(&ideal_signing_test);
        #endif
        
        // checking cyclicity
        quat_alg_conj(&delta,&gen);
        assert(quat_lattice_contains(&coeffs,&ideal_eichler_rand.lattice,&delta,&QUATALG_PINFTY));
        // quat_alg_elem_copy(&delta,&gen);
        quat_alg_mul(&delta,&delta,&gen_challenge,&QUATALG_PINFTY);
        ibz_mul(&delta.denom,&delta.denom,&ideal_eichler_rand.norm);
        quat_alg_normalize(&delta);
        assert(quat_lattice_contains(&coeffs,&right_order_key,&delta,&QUATALG_PINFTY));
        found = found && quat_alg_is_primitive(&delta,&right_order_key,&QUATALG_PINFTY);
        if (!found) { 
            continue;
        }

        assert(SQISIGN_signing_length == KLPT_signing_klpt_length/TORSION_PLUS_EVEN_POWER);

        found = found && id2iso_ideal_to_isogeny_two_long_power_of_2(&sig->zip,&curve,&basis_minus,&basis_plus,&kernel_dual,&gen,SQISIGN_signing_length,&sk->lideal_small,&sk->lideal_two,&sk->gen_two,&QUATALG_PINFTY);
    }

    // if (!found)
    //     return found;

    // quick check to see if we got the correct curve in the end
#ifndef NDEBUG
    if (found) {
        fp2_t jchall, jresp;
        ec_j_inv(&jchall, &E2_for_testing);
        ec_j_inv(&jresp, &curve);
        assert(fp2_is_equal(&jchall, &jresp));
    } else {
        printf("signature failed \n");
    }
#endif

    // var finalize 
    ibq_finalize(&ibq_norm);
    ibz_finalize(&norm);ibz_finalize(&temp);
    ibz_finalize(&remainder);
    quat_order_finalize(&right_order_key);
    quat_left_ideal_finalize(&ideal_commit);
    quat_left_ideal_finalize(&ideal_challenge);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_small);
    quat_alg_elem_finalize(&delta);
    ibz_mat_4x4_finalize(&reduced);
    ibz_mat_4x4_finalize(&gram);
    quat_alg_coord_finalize(&coeffs);
    quat_left_ideal_finalize(&ideal_comchall);
    quat_left_ideal_finalize(&ideal_seccomchall);
    quat_left_ideal_finalize(&ideal_eichler_rand);
    quat_left_ideal_finalize(&ideal_pullback);
    quat_left_ideal_finalize(&ideal_input_klpt);
    quat_alg_elem_finalize(&gen_challenge);

    return found;
}

