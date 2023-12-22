#include <quaternion.h>
#include <ec.h>
#include <endomorphism_action.h>
#include <id2iso.h>
#include <inttypes.h>
#include <locale.h> 
#include <bench.h>

static __inline__ uint64_t rdtsc(void)
{
    return (uint64_t) cpucycles();
}

static int test_point_order_twof(ec_point_t *P, ec_curve_t *E) {
    ec_point_t test = *P;
    assert(!fp2_is_zero(&test.z));
    for (int i = 0;i<TORSION_PLUS_EVEN_POWER-1;i++) {
        ec_dbl(&test,E,&test);
    }
    assert(!fp2_is_zero(&test.z));
    ec_dbl(&test,E,&test);
    return (fp2_is_zero(&test.z));
}


//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.im[i]);
    printf("\n");
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}

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
//XXX

void id2iso_long_two_isog_init(id2iso_long_two_isog_t *isog, const size_t length)
{
    isog->length = length;
    isog->chain = malloc(length * sizeof(*isog->chain));
}

void id2iso_long_two_isog_finalize(id2iso_long_two_isog_t *isog)
{
    free(isog->chain);
}

void id2iso_compressed_long_two_isog_init(id2iso_compressed_long_two_isog_t *zip, const size_t length)
{
    zip->length = length;
    zip->zip_chain = malloc(length * sizeof(*zip->zip_chain));
    for (size_t i = 0; i < length; ++i)
        ibz_init(&zip->zip_chain[i]);
}

void id2iso_compressed_long_two_isog_finalize(id2iso_compressed_long_two_isog_t *zip)
{
    for (size_t i = 0; i < zip->length; ++i)
        ibz_finalize(&zip->zip_chain[i]);
    free(zip->zip_chain);
}

/**
 * @brief evaluation of an endomorphism beta of norm dividing T² on the 2^f torsion
 *
 * @param action_matrix Output : the matrix of the action 
 * @param ideal_beta_L : a left O0-ideal for the left part of the endomorphism
 * @param ideal_beta_R : a left 00-ideal for the right part of the endomorphism
 * @param trace : the trace of the endomorphism mod 2^f
 * @param Bpoo : the quaternion algebra
 * @param domain : the starting curve
 * @param basis_minus : odd torsion basis
 * @param basis_plus : odd torsion basis
 * @param two_basis : a basis of the 2^f torsion
 * @param ibz_two_f : the value of 2^f
 * @returns a bit indicating if the computation succeeded
 * the matrix given in output correspond to the action of the endomorphism beta on the 2^f torsion  
 * the trace is used to lift the sign ambiguity
 */
void endomorphism_evaluation(ibz_mat_2x2_t *action_matrix, const quat_left_ideal_t *ideal_beta_L, const quat_left_ideal_t *ideal_beta_R, const ibz_t *trace, const quat_alg_t *Bpoo, const ec_curve_t *domain, const ec_basis_t *basis_minus, const ec_basis_t *basis_plus, const ec_basis_t *two_basis,const ibz_t *ibz_two_f) {

    // var dec
    ec_isog_odd_t isog_beta_L,isog_beta_R; 
    ec_basis_t basis_R;ec_point_t basis_L[2]; 
    ec_curve_t codomain_L,codomain_R; 
    ec_isom_t isom;
    digit_t x1[NWORDS_ORDER] = {0},x2[NWORDS_ORDER] = {0},x3[NWORDS_ORDER] = {0},x4[NWORDS_ORDER] = {0};  
    ibz_t i1,i2,i3,i4,trace_temp;
    ibz_init(&i1);ibz_init(&i2);ibz_init(&i3);ibz_init(&i4);ibz_init(&trace_temp);

    // translate the ideals ideal_beta_R[ind] and ideal_beta_L[ind] to isogenies 
    id2iso_ideal_to_isogeny_odd(&isog_beta_L,domain,basis_plus,basis_minus,ideal_beta_L);
    id2iso_ideal_to_isogeny_odd(&isog_beta_R,domain,basis_plus,basis_minus,ideal_beta_R);

    // init of points and basis
    basis_R = *two_basis;
    basis_L[0] = two_basis->P;
    basis_L[1] = two_basis->Q; 

    codomain_R = *domain;
    codomain_L = *domain;
    assert(test_point_order_twof(&basis_R.Q,&codomain_R));
    assert(test_point_order_twof(&basis_R.PmQ,&codomain_R));
    assert(test_point_order_twof(&basis_R.P,&codomain_R));
    assert(test_point_order_twof(&basis_L[0],&codomain_L));

    // evaluation of the 2^f basis 
    // evaluating the right isogeny on the full basis
    ec_eval_odd_basis(&codomain_R,&isog_beta_R,&basis_R,1);

    // evaluation the left isogeny on the incomplete basis
    ec_eval_odd(&codomain_L,&isog_beta_L,basis_L,2);

    // checking that the points still have the correct order
    assert(test_point_order_twof(&basis_L[0],&codomain_L));

    // for debug, we check that the two codomains are isomorphic 
    #ifndef NDEBUG 
        fp2_t j_R,j_L;
        ec_j_inv(&j_R, &codomain_R);
        ec_j_inv(&j_L, &codomain_L);
        assert(fp2_is_equal(&j_R,&j_L));
    #endif 

    // computing the isomorphism
    ec_isomorphism(&isom, &codomain_L, &codomain_R);

    // applying the isomorphism
    ec_iso_eval(&basis_L[0],&isom);
    ec_iso_eval(&basis_L[1],&isom);

    assert(test_point_order_twof(&basis_L[1],&codomain_R));
    assert(test_point_order_twof(&basis_R.P,&codomain_R));
    assert(test_point_order_twof(&basis_R.Q,&codomain_R));
    assert(test_point_order_twof(&basis_R.PmQ,&codomain_R));
    assert(ec_is_on_curve(&codomain_R,&basis_R.P));
    assert(ec_is_on_curve(&codomain_R,&basis_R.Q));
    assert(ec_is_on_curve(&codomain_R,&basis_R.PmQ));
    assert(ec_is_on_curve(&codomain_R,&basis_L[0]));

    // DLP step
    // we do the DLP and get the result up to sign 
    // first DLP
    assert(test_point_order_twof(&basis_L[0],&codomain_R));
    ec_dlog_2(x1,x2,&basis_R,&basis_L[0],&codomain_R);

    // copying the digits
    ibz_copy_digit_array(&i2,x2);
    ibz_copy_digit_array(&i1,x1);

    // second DLP
    ec_dlog_2(x3,x4,&basis_R,&basis_L[1],&codomain_R);

    //copyting the digits
    ibz_copy_digit_array(&i3,x3);
    ibz_copy_digit_array(&i4,x4);

    #ifndef NDEBUG
        ec_point_t test;
        ec_biscalar_mul(&test,&codomain_R,x1,x2,&basis_R);
        assert(ec_is_equal(&test,&basis_L[0]));
        ec_biscalar_mul(&test,&codomain_R,x3,x4,&basis_R);
        assert(ec_is_equal(&test,&basis_L[1]));
    #endif


    ibz_copy_digit_array(&i1,x1);
    ibz_copy_digit_array(&i2,x2);
    ibz_copy_digit_array(&i3,x3);
    ibz_copy_digit_array(&i4,x4);

    // multiplication by the norm of ideal_beta_R
    ibz_mul(&i1,&i1,&ideal_beta_R->norm);
    ibz_mul(&i2,&i2,&ideal_beta_R->norm);
    ibz_mul(&i3,&i3,&ideal_beta_R->norm);
    ibz_mul(&i4,&i4,&ideal_beta_R->norm);

    // computation mod 2^f
    ibz_mod(&i1,&i1,ibz_two_f);
    ibz_mod(&i2,&i2,ibz_two_f);
    ibz_mod(&i3,&i3,ibz_two_f);
    ibz_mod(&i4,&i4,ibz_two_f);

    // compute the trace 
    ibz_add(&trace_temp,&i1,&i4);
    ibz_mod(&trace_temp,&trace_temp,ibz_two_f);

    if (ibz_cmp(&trace_temp,trace)!=0){
        ibz_neg(&i1,&i1);
        ibz_neg(&i2,&i2);
        ibz_mod(&i1,&i1,ibz_two_f);
        ibz_mod(&i2,&i2,ibz_two_f);
        ibz_add(&trace_temp,&i1,&i4);
        ibz_mod(&trace_temp,&trace_temp,ibz_two_f);
        if (ibz_cmp(&trace_temp,trace)!=0) {
            ibz_neg(&i4,&i4);
            ibz_neg(&i3,&i3);
            ibz_mod(&i4,&i4,ibz_two_f);
            ibz_mod(&i3,&i3,ibz_two_f);
            ibz_add(&trace_temp,&i1,&i4);
            ibz_mod(&trace_temp,&trace_temp,ibz_two_f);
            if (ibz_cmp(trace,&trace_temp)!=0) {
                ibz_neg(&i1,&i1);
                ibz_neg(&i2,&i2);
                ibz_mod(&i1,&i1,ibz_two_f);
                ibz_mod(&i2,&i2,ibz_two_f);
                ibz_add(&trace_temp,&i1,&i4);
                ibz_mod(&trace_temp,&trace_temp,ibz_two_f);
            }
        }
    }
    // for debug checking final equality
    assert(0==ibz_cmp(&trace_temp,trace));

    // we have the correct values, now we can copy the output
    ibz_copy(&(*action_matrix)[0][0],&i1);
    ibz_copy(&(*action_matrix)[0][1],&i2);
    ibz_copy(&(*action_matrix)[1][0],&i3);
    ibz_copy(&(*action_matrix)[1][1],&i4);


    // var finalize
    ibz_finalize(&i1);ibz_finalize(&i2);ibz_finalize(&i3);ibz_finalize(&i4);
    ibz_finalize(&trace_temp);

}


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
 * when used for compressing an isogeny, we assume that the curve given in input is normalized!!
 */
int id2iso_ideal_to_isogeny_two_long_power_of_2(id2iso_compressed_long_two_isog_t *isog_zip, ec_curve_t *domain, ec_basis_t *basis_minus, ec_basis_t *basis_plus, ec_point_t *kernel_dual, const quat_alg_elem_t *gen_input, const int length, const quat_left_ideal_t *lideal_start_small, const quat_left_ideal_t *lideal_start_two, const quat_alg_elem_t *gen_two, const quat_alg_t *Bpoo) {
    
    // var dec 
    int found;
    quat_alg_elem_t beta;
    quat_alg_elem_t gen,gen_constraint;
    ibz_t n_beta,n;
    quat_left_ideal_t lideal_small;
    quat_left_ideal_t ideal_beta_L[length],ideal_beta_R[length];
    ibq_t ibq_trace;
    ibz_t beta_trace[length];
    ibz_vec_2_t linear_comb[length];
    quat_alg_elem_t gen_input_temp,quat_temp;
    quat_alg_elem_t gen_two_temp;
    quat_left_ideal_t lideal_equiv;quat_left_ideal_init(&lideal_equiv); quat_order_t right_order;quat_order_init(&right_order); // TODECIDE we might avoid this by solving solve linear combination without using the right order 
    quat_alg_coord_t coeffs;

    int const f = TORSION_PLUS_EVEN_POWER;
    ibz_t ibz_two_f,temp,remainder;
    ec_basis_t temp_odd_basis[2]; 
    ec_curve_t temp_domain; ec_point_t temp_kernel_dual; ec_basis_t temp_two_basis;
    ec_isog_even_t isog;
    ec_point_t pushed_points[7];
    ibz_mat_2x2_t action_matrix;
    digit_t digit_a[NWORDS_ORDER] = {0}; digit_t digit_b[NWORDS_ORDER] = {0};
    ibz_t ibz_a,ibz_b;

    // var init
    found = 0;
    quat_alg_elem_init(&beta);
    quat_alg_elem_init(&gen);quat_alg_elem_init(&gen_constraint);
    quat_alg_elem_init(&gen_input_temp);quat_alg_elem_init(&quat_temp);
    quat_alg_elem_init(&gen_two_temp);
    ibz_init(&n_beta); ibz_init(&n);
    ibz_init(&temp); ibz_init(&ibz_two_f);
    ibz_init(&remainder);
    ibq_init(&ibq_trace);
    ibz_init(&ibz_a);ibz_init(&ibz_b);
    ibz_mat_2x2_init(&action_matrix);
    quat_left_ideal_init(&lideal_small);
    quat_alg_coord_init(&coeffs);
    for (int ind=0; ind<length;ind++) {
        quat_left_ideal_init(&ideal_beta_L[ind]);
        quat_left_ideal_init(&ideal_beta_R[ind]);
        ibz_init(&beta_trace[ind]);
        ibz_vec_2_init(&linear_comb[ind]);
    }

    //ibz_two_f = 2^f
    ibz_pow(&ibz_two_f,&ibz_const_two,f);
    
    // first, we perform all operations over the quaternions to verify that no failure occurs before doing the costly ec operations    
    
    // initializing the quaternion elements and ideals we are going to use 
    // here is a list of all the important elements :
    // integer n,
    // quaternion element gen
    // quaternion element gen_input_temp
    // quaternion element gen_two_temp
    // ideal lideal_small  
 
    // Before iteration ind: 
    //  - lideal_small is a small O0-ideal (in the sense that it has somewhat small norm for ideals of the same class) whose right order is the endomorphism ring of the domain of the (ind+1)-th isogeny to be translated 
    //    in particular at ind = 0, it is simply lideal_start_small   
    //  - gen is the generator of an ideal of norm n equivalent to lideal_start_small. We call it lideal_equiv even though this ideal is never explicitly computed, and it is equal to lideal_equiv = lideal_small * gen / norm    
    //    (lideal_small) ( = O0 < gen, n(gen)/n(lideal_small) > when ind is bigger than 0). conj(gen) is contained in lideal_small.  
    //    At ind = 0, we have lideal_equiv = lideal_small  
    //  - gen_input_temp is an element of the right order of lideal_equiv generating the ideal of norm 2^f corresponding to the (ind+1)-th isogeny to be translated.
    //  - gen_two_temp is an element of the right order of lideal_equiv such that the principal ideal O0 < gen_two_temp > = lideal_two * conj(lideal_equiv) where lideal_two corresponds to the isogeny of degree 2^*  connecting E0 to the domain of the (ind+1)-th isogeny to be translated

    // At the beginning of iteration ind:
    //  - lideal_small is set to be lideal_equiv times an ideal of norm 2^f corresponding to the (ind+1)-th isogeny to be translated
    // During iteration ind, all the other values will be updated  

    // The main goal of iteration ind is to find an endomorphism beta contained inside the right order of lideal_equiv of step ind+1 (which is isomorphic to the endomorphism ring of the codomain of the (ind+1)-th isogeny to be translated) 
    // this endomorphism beta will be used to compute the two ideals ideal_beta_L[ind],ideal_beta_R[ind] and the trace trace_beta[ind] which are required to actually compute the (ind+1)-th isogeny to be translated. 
    // Two coefficients linear_comb[ind] will also be necessary, it is also computed in the loop below.  

    // gen_constraint  = conjugate( gen_two ) 
    quat_alg_conj(&gen_constraint,gen_two);

    // gen = norm(lideal_start_small)
    quat_alg_scalar(&gen,&lideal_start_small->norm,&ibz_const_one);

    // n = norm(lideal_small) T
    ibz_copy(&n,&lideal_start_small->norm);

    // gen_input_temp = gen_input
    quat_alg_elem_copy(&gen_input_temp,gen_input);


    // lideal_small = lideal_start_small 
    quat_left_ideal_copy(&lideal_small,lideal_start_small);
    
    // gen_two_temp = gen_two
    quat_alg_elem_copy(&gen_two_temp,gen_two);

    // for debug we check that all values are contained in the correct orders and lattices
    assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen_two_temp,Bpoo));
    assert(quat_lattice_contains(&coeffs,&lideal_start_two->lattice,&gen_two_temp,Bpoo));
    assert(!quat_lattice_contains(&coeffs,&lideal_start_small->lattice,&gen_two_temp,Bpoo));
    quat_alg_conj(&quat_temp,&gen_two_temp);
    assert(quat_lattice_contains(&coeffs,&lideal_start_small->lattice,&quat_temp,Bpoo));

    // timing
    setlocale(LC_NUMERIC, "");
    uint64_t t0, t1;

    // we will loop through and try to find an endomorphism beta for each
    for (int ind=0; ind<length;ind++) {
        if (ind > 0) {
            // we update the various ideals and quaternion elements
        
            // if we set lideal_equiv = ideal_small * gen / norm(lideal_small) = O0 < gen, n(gen)/n(lideal_small) >
            // we have set (at the end of the previous iteration) gen_input_temp = conj ( gen ) / norm(gen ) * gen_input_temp * gen 
            // and lideal_small = lideal_equiv * O_R(lideal_equiv) < gen_input_temp , 2^f >  = O0 <  make_primitive (gen * gen_input_temp) , n(gen) * 2^f /n (lideal_small) >

            // update lideal_small
            // quat_temp = gen * gen_input_temp
            quat_alg_mul(&quat_temp,&gen,&gen_input_temp, Bpoo);
            quat_alg_normalize(&quat_temp);

            assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&quat_temp,Bpoo));
            assert(quat_alg_is_primitive(&gen_input_temp,&right_order,Bpoo));
            assert(quat_lattice_contains(&coeffs,&lideal_equiv.lattice,&quat_temp,Bpoo));

            // temp = n * 2^f
            ibz_mul(&temp,&n,&ibz_two_f);
            quat_lideal_make_primitive_then_create(&lideal_small,&quat_temp,&temp,&STANDARD_EXTREMAL_ORDER.order,Bpoo);

        }
        // copying n in temp
        ibz_copy(&temp,&n);

        // trying to solve the norm equation algorithm 
        found = klpt_eichler_special_norm(&beta,&n_beta,&gen,&n,&lideal_small,&gen_constraint, Bpoo); 
        if (!found) {
            break;
        }

        // the solution has been found, we proceed to the remaining computations 
        // this includes computation of an updated values of gen_two_temp and gen_input_temp,
        // of the ideals ideal_beta_L[ind] and ideal_beta_R[ind],
        // of the trace of beta as beta_trace[ind], and of the linear_combination linear_comb[ind]  
        
         #ifndef NDEBUG
            quat_alg_conj(&quat_temp,&gen);
            quat_lattice_contains(&coeffs,&lideal_small.lattice,&quat_temp,Bpoo);
         #endif 



        // computation of lideal_equiv as lideal_small* gen / norm(lideal_small)
        quat_alg_elem_copy(&quat_temp,&gen);
        ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal_small.norm);
        int lideal_mul_ok = quat_lideal_mul(&lideal_equiv,&lideal_small,&quat_temp,Bpoo,0);
        assert(lideal_mul_ok);
        

        // computing the right order of lideal_equiv 
        quat_lideal_right_order(&right_order,&lideal_equiv,Bpoo);
        assert(quat_lattice_contains(&coeffs,&right_order,&beta,Bpoo));

        if (ind > 0) {
            // we can adjust the scalar of gen_input_temp, because it was adjusted for gen_constraint during the execution of eichler norm fixed
            ibz_copy(&gen_input_temp.denom,&gen_constraint.denom);
        }
        

        // update gen_input_temp = conj ( gen ) / norm(gen ) * gen_input_temp * gen
        quat_alg_normalize(&gen);
        quat_alg_conj(&quat_temp,&gen);
        assert(quat_lattice_contains(&coeffs,&lideal_small.lattice,&quat_temp,Bpoo));
        ibz_mul(&quat_temp.denom,&quat_temp.denom,&n);
        ibz_mul(&quat_temp.denom,&quat_temp.denom,&lideal_small.norm);
        quat_alg_mul(&gen_input_temp,&quat_temp,&gen_input_temp,Bpoo);
        quat_alg_mul(&gen_input_temp,&gen_input_temp,&gen,Bpoo);
        quat_alg_normalize(&gen_input_temp);

        // update gen_constraint
        if (ind > 0) {
            // gen_constraint is simply the conjugate of gen_input_temp
            quat_alg_conj(&gen_constraint,&gen_input_temp);
            
        }
        else {
            // gen_constraint = conj (gen) /norm(gen) * gen_constraint * gen
            quat_alg_mul(&gen_constraint,&quat_temp,&gen_constraint,Bpoo);
            quat_alg_mul(&gen_constraint,&gen_constraint,&gen,Bpoo);
            quat_alg_normalize(&gen_constraint);
        }
        // for debug we check that our updated element are contained in the correct orders
        assert(quat_lattice_contains(&coeffs,&right_order,&gen_constraint,Bpoo));
        assert(quat_lattice_contains(&coeffs,&right_order,&gen_input_temp,Bpoo));
        // assert(quat_alg_is_primitive(&gen_input_temp,&right_order,Bpoo));
        assert(quat_lattice_contains(&coeffs,&right_order,&beta,Bpoo));
        assert(quat_alg_is_primitive(&beta,&right_order,Bpoo));
        
        // for debug we check that gen_two_temp is indeed contained inside order.order
        assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen_two_temp,Bpoo));
        
        // update gen_two_temp =   conj( gen ) * gen_two_temp / (lideal_equiv) 
        quat_alg_conj(&quat_temp,&gen);
        quat_alg_mul(&gen_two_temp,&quat_temp,&gen_two_temp,Bpoo);
        ibz_mul(&gen_two_temp.denom,&gen_two_temp.denom,&temp);
        quat_alg_normalize(&gen_two_temp);
        assert(quat_lattice_contains(&coeffs,&STANDARD_EXTREMAL_ORDER.order,&gen_two_temp,Bpoo));

        // computation of the two ideals
        // ideal_beta_L = O0 < beta * gen_two_temp, N>, ideal_beta_R = O0 < conj(beta) * gen_two_temp, N'> with N N' = n_beta

        // quat_temp = beta * gen_two_temp
        quat_alg_mul(&quat_temp,&beta,&gen_two_temp,&QUATALG_PINFTY);

        // computation of the norm of this first ideal
        ibz_gcd(&temp,&n_beta,&TORSION_ODD);       
        // creation of ideal_beta_L
        quat_lideal_make_primitive_then_create(&ideal_beta_L[ind],&quat_temp,&temp,&STANDARD_EXTREMAL_ORDER.order,Bpoo);
        assert(0==ibz_cmp(&temp,&ideal_beta_L[ind].norm));

        // computation of the generator of ideal_beta_R 
        // as conj(beta) * gen_two_temp 
        quat_alg_conj(&quat_temp,&beta);
        quat_alg_mul(&quat_temp,&quat_temp,&gen_two_temp, Bpoo);
        quat_alg_normalize(&quat_temp);

        // computation of the norm of ideal_beta_R 
        ibz_div(&temp,&remainder,&n_beta,&ideal_beta_L[ind].norm);

        // creation of ideal_beta_R 
        quat_lideal_make_primitive_then_create(&ideal_beta_R[ind],&quat_temp,&temp,&STANDARD_EXTREMAL_ORDER.order,Bpoo);
        assert(0==ibz_cmp(&temp,&ideal_beta_R[ind].norm)); 

        #ifndef NDEBUG 
            ibz_mul(&temp,&ideal_beta_R[ind].norm,&ideal_beta_L[ind].norm);
            assert(0==ibz_cmp(&temp,&n_beta));
        #endif

        // computation of the trace of beta
        quat_alg_trace(&ibq_trace,&beta);
        if (!ibq_to_ibz(&beta_trace[ind],&ibq_trace)) {
            assert(0);
        }
        else {
            // adjusting the trace mod 2^f
            ibz_mod(&beta_trace[ind],&beta_trace[ind],&ibz_two_f);
        } 
        
        // removing scalar factors if needed 
        // TODUPDATE : it is only needed if trace is even so we may check that to avoid doing that operation every time
        quat_alg_make_primitive(&coeffs,&temp,&gen_input_temp,&right_order,Bpoo);
        ibz_mul(&gen_input_temp.denom,&gen_input_temp.denom,&temp);
        quat_alg_normalize(&gen_input_temp);

        // computation of the linear_combination such that linear_comb[ind][0] + beta * linear_comb[ind][1] sends the kernel of the dual of the 2^f isogeny of step ind-1 to the kernel of the 2^f isogeny of step ind
        found = klpt_find_linear_comb(&linear_comb[ind],&beta,&right_order,&ibz_two_f,f,&gen_constraint,&gen_input_temp,Bpoo);
        if (!found) {
            break;
            // TODEBUG
        }

        int lideal_generator_ok;
        lideal_generator_ok = quat_lideal_generator_coprime(&gen,&lideal_equiv,&ibz_const_two,Bpoo,0);
        assert(lideal_generator_ok);
        quat_alg_normalize(&gen);

        // ensure that we still have the correct value of gen_constraint
        quat_alg_conj(&gen_constraint,&gen_input_temp);
        

    }
    // the quaternion computation has succeeded, we can proceed to perform all elliptic curve and isogeny operations 
    if (found) {
        // init all ec values needed for this computation

        // domain 
        temp_domain = *domain;

        // kernel dual
        temp_kernel_dual = *kernel_dual;

        // odd torsion points 
        temp_odd_basis[0] = *basis_minus;
        temp_odd_basis[1] = *basis_plus;

        for (int ind =0;ind<length;ind ++) {
            
            // completing the 2^f basis 
            assert(test_point_order_twof(&temp_kernel_dual,&temp_domain));
            ec_complete_basis_2(&temp_two_basis,&temp_domain,&temp_kernel_dual);
            assert(ec_is_equal(&temp_kernel_dual,&temp_two_basis.P));
            assert(test_point_order_twof(&temp_two_basis.P,&temp_domain));
            assert(test_point_order_twof(&temp_two_basis.Q,&temp_domain));
            assert(test_point_order_twof(&temp_two_basis.PmQ,&temp_domain));
            
            // for each step, we must evaluate the endomorphism, compute the next isogny, push the torsion and compute the compression

            // evaluation of the endomorphism 
            endomorphism_evaluation(&action_matrix,&ideal_beta_L[ind],&ideal_beta_R[ind],&beta_trace[ind],Bpoo,&temp_domain,&temp_odd_basis[0],&temp_odd_basis[1], &temp_two_basis,&ibz_two_f);    

            // computation of the linear combination of points to get the kernel 
            // linear_comb[ind][0] = linear_comb[ind][1] * action_matrix[0][0] + linear_comb[ind][0]
            ibz_mul(&temp,&linear_comb[ind][1],&action_matrix[0][0]);
            ibz_add(&linear_comb[ind][0],&temp,&linear_comb[ind][0]);
            ibz_mod(&linear_comb[ind][0],&linear_comb[ind][0],&ibz_two_f);

            // linear_comb[ind][1] = linear_comb[ind][1] * action_matrix[0][1]     
            ibz_mul(&linear_comb[ind][1],&linear_comb[ind][1],&action_matrix[0][1]);
            // ibz_neg(&linear_comb[ind][1],&linear_comb[ind][1]); 
            ibz_mod(&linear_comb[ind][1],&linear_comb[ind][1],&ibz_two_f);
            

            // computation of the 2^f isogeny and evaluation of the points 
            // in that case we can rescale the basis and compute the zipped value directly
            if (ind > 0) {                
                assert(0!=ibz_get(&linear_comb[ind][1])%2);
                ibz_invmod(&linear_comb[ind][1],&linear_comb[ind][1],&ibz_two_f);
                ibz_mul(&linear_comb[ind][0],&linear_comb[ind][0],&linear_comb[ind][1]);
                ibz_mod(&linear_comb[ind][0],&linear_comb[ind][0],&ibz_two_f);
                ibz_copy(&linear_comb[ind][1],&ibz_const_one);
                ibz_copy(&isog_zip->zip_chain[ind],&linear_comb[ind][0]);
            }

            // computation of the kernel
            // changing the ibzs to digits
            ibz_to_digit_array(digit_a,&linear_comb[ind][0]);
            ibz_to_digit_array(digit_b,&linear_comb[ind][1]);

            isog.kernel = temp_two_basis.P;

            ec_biscalar_mul(&isog.kernel,&temp_domain,digit_a,digit_b,&temp_two_basis);
            assert(test_point_order_twof(&isog.kernel,&temp_domain));

            // setting the isogeny
            isog.curve = temp_domain;
            isog.length = f;

            // compression for the first step 
            // the first step is different for the compression
            if (ind == 0) {
                // first we compute a new deterministic basis of the curve 
                ec_curve_to_basis_2(&temp_two_basis,&temp_domain);

                // then we perform a dlp to express isog.kernel in this basis
                ec_dlog_2(digit_a,digit_b,&temp_two_basis,&isog.kernel,&temp_domain);

                // translate the digit_t as ibz_t
                ibz_copy_digit_array(&ibz_a,digit_a);
                ibz_copy_digit_array(&ibz_b,digit_b);
                // testing the sign

        
                #ifndef NDEBUG 
                    ec_point_t test;
                    ec_biscalar_mul(&test,&temp_domain,digit_a,digit_b,&temp_two_basis);
                    assert(ec_is_equal(&test,&isog.kernel));
                #endif

                // then we encode the isogeny as a scalar and a bit
                // set the bit_first_step as 0 when we can get a generator as temp_two_basis.P + s temp_two_basis.Q
                // if not, then this means we can write as Q + s P
                isog_zip->bit_first_step = ibz_get(&ibz_a)%2!=0; 
                if (isog_zip->bit_first_step) {
                    // compute the scalar s as ibz_b / ibz_a mod 2^f 
                    if(!ibz_invmod(&ibz_a,&ibz_a,&ibz_two_f)) {
                        assert(0);
                    }
                    ibz_mul(&ibz_b,&ibz_b,&ibz_a);
                    ibz_mod(&ibz_b,&ibz_b,&ibz_two_f);
                    ibz_copy(&isog_zip->zip_chain[ind],&ibz_b);
                    // and we can set the kernel dual to be temp_two_basis.Q
                    temp_kernel_dual = temp_two_basis.Q;

                    ibz_to_digit_array(digit_b,&ibz_b);
                    ibz_copy_digit_array(&ibz_b,digit_b);
                    // ec_ladder3pt(&isog.kernel,digit_b,&temp_two_basis.P,&temp_two_basis.Q,&temp_two_basis.PmQ,&temp_domain);
                }
                else {
                    // compute the scalar s as ibz_a / ibz_b mod 2^f 
                    if (!ibz_invmod(&ibz_b,&ibz_b,&ibz_two_f)) {
                        assert(0);
                    }
                    ibz_mul(&ibz_a,&ibz_b,&ibz_a);
                    ibz_mod(&ibz_a,&ibz_a,&ibz_two_f);
                    ibz_copy(&isog_zip->zip_chain[ind],&ibz_a);
                    // and we can set the kernel dual to be temp_two_basis.P
                    temp_kernel_dual = temp_two_basis.P;
                    
                    ibz_to_digit_array(digit_a,&ibz_a);
                    ibz_copy_digit_array(&ibz_a,digit_a);
                    // ec_ladder3pt(&isog.kernel,digit_a,&temp_two_basis.Q,&temp_two_basis.P,&temp_two_basis.PmQ,&temp_domain);
                } 

            }
            
            assert(test_point_order_twof(&temp_kernel_dual,&temp_domain));
            // evaluating the isogeny
            // TODECIDE : we may use a non-zero eval for this one
            pushed_points[0] = temp_odd_basis[0].P;
            pushed_points[1] = temp_odd_basis[0].Q;
            pushed_points[2] = temp_odd_basis[0].PmQ;
            pushed_points[3] = temp_odd_basis[1].P;
            pushed_points[4] = temp_odd_basis[1].Q;
            pushed_points[5] = temp_odd_basis[1].PmQ;
            pushed_points[6] = temp_kernel_dual;
            ec_eval_even(&temp_domain,&isog,pushed_points,7);
            temp_odd_basis[0].P = pushed_points[0];
            temp_odd_basis[0].Q = pushed_points[1];
            temp_odd_basis[0].PmQ = pushed_points[2];
            temp_odd_basis[1].P = pushed_points[3];
            temp_odd_basis[1].Q = pushed_points[4];
            temp_odd_basis[1].PmQ = pushed_points[5];
            temp_kernel_dual = pushed_points[6];
            assert(test_point_order_twof(&temp_kernel_dual,&temp_domain)); 
            
        }

        *domain = temp_domain;
        *kernel_dual = temp_kernel_dual;
        *basis_minus = temp_odd_basis[0];
        *basis_plus = temp_odd_basis[1];

       
    }

     // var finalize 
    quat_alg_elem_finalize(&beta);
    quat_alg_elem_finalize(&gen);quat_alg_elem_finalize(&gen_constraint);
    quat_alg_elem_finalize(&gen_input_temp);quat_alg_elem_finalize(&quat_temp);
    quat_alg_elem_finalize(&gen_two_temp);
    ibz_finalize(&n_beta); ibz_finalize(&n);
    ibz_finalize(&temp); ibz_finalize(&ibz_two_f);
    ibz_finalize(&remainder);
    ibz_finalize(&ibz_a);ibz_finalize(&ibz_b);
    ibq_finalize(&ibq_trace);
    ibz_mat_2x2_finalize(&action_matrix);
    quat_left_ideal_finalize(&lideal_small);
    quat_alg_coord_finalize(&coeffs);
    quat_left_ideal_finalize(&lideal_equiv); quat_order_finalize(&right_order);
    for (int ind=0; ind<length;ind++) {
        quat_left_ideal_finalize(&ideal_beta_L[ind]);
        quat_left_ideal_finalize(&ideal_beta_R[ind]);
        ibz_finalize(&beta_trace[ind]);
        ibz_vec_2_finalize(&linear_comb[ind]);
    }

    return found;

   

}





//TODO this should probably be in the quaternion module
static void from_1ijk_to_O0basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el)
{
    ibz_copy(&(*vec)[2], &el->coord[2]);
    ibz_copy(&(*vec)[3], &el->coord[3]);
    ibz_sub(&(*vec)[0], &el->coord[0], &el->coord[3]);
    ibz_sub(&(*vec)[1], &el->coord[1], &el->coord[2]);
    if (ibz_is_one(&el->denom)) {
        ibz_add(&(*vec)[2], &(*vec)[2], &(*vec)[2]);
        ibz_add(&(*vec)[3], &(*vec)[3], &(*vec)[3]);
    } else {
        assert(!ibz_cmp(&el->denom, &ibz_const_two));
        ibz_div_2exp(&(*vec)[0], &(*vec)[0], 1);
        ibz_div_2exp(&(*vec)[1], &(*vec)[1], 1);
    }
}

void id2iso_ideal_to_kernel_dlogs_even(ibz_vec_2_t *vec, const quat_left_ideal_t *lideal)
{
    ibz_t tmp;
    ibz_init(&tmp);

    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    // construct the matrix of the dual of alpha on the 2^f-torsion
    {
        quat_alg_elem_t alpha;
        quat_alg_elem_init(&alpha);

        int lideal_generator_ok;
        lideal_generator_ok =  quat_lideal_generator(&alpha, lideal, &QUATALG_PINFTY,0);
        assert(lideal_generator_ok);
        quat_alg_conj(&alpha, &alpha);

        ibz_vec_4_t coeffs;
        ibz_vec_4_init(&coeffs);
        from_1ijk_to_O0basis(&coeffs, &alpha);

        for (unsigned i = 0; i < 2; ++i) {
            ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
            for (unsigned j = 0; j < 2; ++j) {
                ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
            }
        }

        ibz_vec_4_finalize(&coeffs);
        quat_alg_elem_finalize(&alpha);
    }

    // find the kernel of alpha modulo the norm of the ideal
    {
        ibz_t const *const norm = &lideal->norm;

        ibz_mod(&(*vec)[0], &mat[0][0], norm);
        ibz_mod(&(*vec)[1], &mat[1][0], norm);
        ibz_gcd(&tmp, &(*vec)[0], &(*vec)[1]);
        if (!(ibz_get(&tmp) & 1)) {
            ibz_mod(&(*vec)[0], &mat[0][1], norm);
            ibz_mod(&(*vec)[1], &mat[1][1], norm);
        }
#ifndef NDEBUG
        ibz_gcd(&tmp, &(*vec)[0], norm);
        ibz_gcd(&tmp, &(*vec)[1], &tmp);
        assert(!ibz_cmp(&tmp, &ibz_const_one));
#endif
    }

    ibz_mat_2x2_finalize(&mat);
    ibz_finalize(&tmp);
}



void id2iso_ideal_to_isogeny_even(ec_isog_even_t *isog, const quat_left_ideal_t *lideal_input)
{
    // compute length
    isog->length = 0;
    ibz_t norm;
    ibz_init(&norm);
    ibz_copy(&norm, &lideal_input->norm);
    while (!ibz_is_one(&norm)) {
        assert(!ibz_is_zero(&norm) && !(ibz_get(&norm) & 1));
        ibz_div_2exp(&norm, &norm, 1);
        ++isog->length;
    }
    ibz_finalize(&norm);
    assert(isog->length <= TORSION_PLUS_EVEN_POWER);

    digit_t scalars[2][NWORDS_FIELD];
    {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);

        id2iso_ideal_to_kernel_dlogs_even(&vec, lideal_input);

        // multiply out unnecessary cofactor from 2^f-torsion basis
        for (size_t i = isog->length; i < TORSION_PLUS_EVEN_POWER; ++i) {
            ibz_add(&vec[0], &vec[0], &vec[0]);
            ibz_add(&vec[1], &vec[1], &vec[1]);
        }

        ibz_to_digit_array(scalars[0], &vec[0]);
        ibz_to_digit_array(scalars[1], &vec[1]);

        ibz_vec_2_finalize(&vec);
    }

    isog->curve = CURVE_E0;
    ec_biscalar_mul(&isog->kernel, &isog->curve, scalars[0], scalars[1], &BASIS_EVEN);

}



void id2iso_ideal_to_kernel_dlogs_odd(ibz_vec_2_t *vec, ec_degree_odd_t *deg, const quat_left_ideal_t *lideal)
{
    ibz_t tmp;
    ibz_init(&tmp);

    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);

    // construct the matrix of the dual of alpha on the T-torsion
    {
        quat_alg_elem_t alpha;
        quat_alg_elem_init(&alpha);

        int lideal_generator_ok;
        lideal_generator_ok = quat_lideal_generator(&alpha, lideal, &QUATALG_PINFTY,0);
        assert(lideal_generator_ok);
        assert(ibz_divides(&ibz_const_two, &alpha.denom));  // denominator is invertible mod T, ignore

        for (unsigned i = 0; i < 2; ++i) {
            ibz_add(&mat[i][i], &mat[i][i], &alpha.coord[0]);
            for (unsigned j = 0; j < 2; ++j) {
                ibz_mul(&tmp, &ACTION_I[i][j], &alpha.coord[1]);
                ibz_sub(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_J[i][j], &alpha.coord[2]);
                ibz_sub(&mat[i][j], &mat[i][j], &tmp);
                ibz_mul(&tmp, &ACTION_K[i][j], &alpha.coord[3]);
                ibz_sub(&mat[i][j], &mat[i][j], &tmp);
//                ibz_mod(&mat[i][j], &mat[i][j], &TORSION_ODD);
            }
        }

        quat_alg_elem_finalize(&alpha);
    }

    // determine prime powers in the norm of the ideal
    size_t numpp = 0;
    #define NUMPP (sizeof(TORSION_ODD_PRIMEPOWERS) / sizeof(*TORSION_ODD_PRIMEPOWERS))
    ibz_t pps[NUMPP];
    ibz_vec_2_t vs[NUMPP];
    {
        ibz_t const *const norm = &lideal->norm;

        for (size_t i = 0; i < NUMPP; ++i) {
            ibz_gcd(&tmp, norm, &TORSION_ODD_PRIMEPOWERS[i]);
            (*deg)[i] = 0;
            if (!ibz_is_one(&tmp)) {
                ibz_init(&pps[numpp]);
                ibz_copy(&pps[numpp], &tmp);
                ibz_vec_2_init(&vs[numpp]);
                ++numpp;

                // compute valuation
                ibz_t l, r;
                ibz_init(&l);
                ibz_init(&r);
                ibz_set(&l, TORSION_ODD_PRIMES[i]);
                do {
                    ++(*deg)[i];
                    ibz_div(&tmp, &r, &tmp, &l);
                    assert(!ibz_is_zero(&tmp) && ibz_is_zero(&r));
                } while (!ibz_is_one(&tmp));
                ibz_finalize(&r);
                ibz_finalize(&l);
            }
        }
    }
    #undef NUMPP

    // find the kernel of alpha modulo each prime power
    {
        for (size_t i = 0; i < numpp; ++i) {
            ibz_mod(&vs[i][0], &mat[0][0], &pps[i]);
            ibz_mod(&vs[i][1], &mat[1][0], &pps[i]);
            ibz_gcd(&tmp, &vs[i][0], &pps[i]);
            ibz_gcd(&tmp, &vs[i][1], &tmp);
            if (ibz_cmp(&tmp, &ibz_const_one)) {
                ibz_mod(&vs[i][0], &mat[0][1], &pps[i]);
                ibz_mod(&vs[i][1], &mat[1][1], &pps[i]);
            }
#ifndef NDEBUG
            ibz_gcd(&tmp, &vs[i][0], &pps[i]);
            ibz_gcd(&tmp, &vs[i][1], &tmp);
            assert(!ibz_cmp(&tmp, &ibz_const_one));
#endif
        }
    }

    // now CRT them together
    {
        //TODO use a product tree instead
        ibz_t mod;
        ibz_init(&mod);
        ibz_set(&mod, 1);
        ibz_set(&(*vec)[0], 0);
        ibz_set(&(*vec)[1], 0);
        for (size_t i = 0; i < numpp; ++i) {
            //TODO use vector CRT
            ibz_crt(&(*vec)[0], &(*vec)[0], &vs[i][0], &mod, &pps[i]);
            ibz_crt(&(*vec)[1], &(*vec)[1], &vs[i][1], &mod, &pps[i]);
            //TODO optionally return lcm from CRT and use it
            ibz_mul(&mod, &mod, &pps[i]);
            ibz_finalize(&pps[i]);
        }
        ibz_finalize(&mod);
    }

    for (size_t i = 0; i < numpp; ++i)
        ibz_vec_2_finalize(&vs[i]);

    ibz_mat_2x2_finalize(&mat);

    ibz_finalize(&tmp);
}




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

void id2iso_ideal_to_isogeny_odd(ec_isog_odd_t *isog, const ec_curve_t *domain, const ec_basis_t *basis_plus,const ec_basis_t *basis_minus, const quat_left_ideal_t *lideal_input)
{
    digit_t scalars_plus[2][NWORDS_FIELD], scalars_minus[2][NWORDS_FIELD];
    {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);

        id2iso_ideal_to_kernel_dlogs_odd(&vec, &isog->degree, lideal_input);

        ibz_t tmp;
        ibz_init(&tmp);

        // multiply out unnecessary cofactor from T-torsion basis
        assert(sizeof(isog->degree)/sizeof(*isog->degree)
                == sizeof(TORSION_ODD_PRIMEPOWERS)/sizeof(*TORSION_ODD_PRIMEPOWERS));
        for (size_t i = 0; i < sizeof(isog->degree)/sizeof(*isog->degree); ++i) {
            assert(isog->degree[i] <= TORSION_ODD_POWERS[i]);
            if (isog->degree[i] == TORSION_ODD_POWERS[i])
                continue;
            ibz_set(&tmp, TORSION_ODD_PRIMES[i]);
            ibz_pow(&tmp, &tmp, TORSION_ODD_POWERS[i] - isog->degree[i]);
            ibz_mul(&vec[0], &vec[0], &tmp);
            ibz_mul(&vec[1], &vec[1], &tmp);
        }

        ibz_mod(&tmp, &vec[0], &TORSION_ODD_PLUS);
        ibz_to_digit_array(scalars_plus[0], &tmp);
        ibz_mod(&tmp, &vec[1], &TORSION_ODD_PLUS);
        ibz_to_digit_array(scalars_plus[1], &tmp);
        ibz_mod(&tmp, &vec[0], &TORSION_ODD_MINUS);
        ibz_to_digit_array(scalars_minus[0], &tmp);
        ibz_mod(&tmp, &vec[1], &TORSION_ODD_MINUS);
        ibz_to_digit_array(scalars_minus[1], &tmp);

        ibz_finalize(&tmp);

        ibz_vec_2_finalize(&vec);
    }

    isog->curve = *domain;
    ec_biscalar_mul(&isog->ker_plus, domain, scalars_plus[0], scalars_plus[1], basis_plus);
    ec_biscalar_mul(&isog->ker_minus, domain, scalars_minus[0], scalars_minus[1], basis_minus);
}


void id2iso_kernel_dlogs_to_ideal(quat_left_ideal_t *lideal, const ibz_vec_2_t *vec2, const ibz_vec_2_t *vec3)
{
    ibz_vec_2_t ker;
    ibz_vec_2_init(&ker);
    ibz_crt(&ker[0], &(*vec2)[0], &(*vec3)[0], &TORSION_PLUS_2POWER, &TORSION_PLUS_3POWER);
    ibz_crt(&ker[1], &(*vec2)[1], &(*vec3)[1], &TORSION_PLUS_2POWER, &TORSION_PLUS_3POWER);

    // algorithm: apply endomorphisms 1 and j+(1+k)/2 to the kernel point,
    // the result should form a basis of the respective torsion subgroup.
    // then apply i to the kernel point and decompose over said basis.
    // hence we have an equation a*P + b*[j+(1+k)/2]P == [i]P, which will
    // easily reveal an endomorphism that kills P.

    ibz_vec_2_t vec;
    ibz_vec_2_init(&vec);

    {
        ibz_mat_2x2_t mat;
        ibz_mat_2x2_init(&mat);

        ibz_copy(&mat[0][0], &ker[0]);
        ibz_copy(&mat[1][0], &ker[1]);

        ibz_mat_2x2_eval(&vec, &ACTION_J, &ker);
        ibz_copy(&mat[0][1], &vec[0]);
        ibz_copy(&mat[1][1], &vec[1]);
        ibz_mat_2x2_eval(&vec, &ACTION_GEN4, &ker);
        ibz_add(&mat[0][1], &mat[0][1], &vec[0]);
        ibz_add(&mat[1][1], &mat[1][1], &vec[1]);
        ibz_mod(&mat[0][1], &mat[0][1], &TORSION_PLUS_23POWER);
        ibz_mod(&mat[1][1], &mat[1][1], &TORSION_PLUS_23POWER);

        ibz_mat_2x2_t inv;
        ibz_mat_2x2_init(&inv);
        {
            int inv_ok = ibz_2x2_inv_mod(&inv, &mat, &TORSION_PLUS_23POWER);
            assert(inv_ok);
        }
        ibz_mat_2x2_finalize(&mat);

        ibz_mat_2x2_eval(&vec, &ACTION_I, &ker);
        ibz_mat_2x2_eval(&vec, &inv, &vec);

        ibz_mat_2x2_finalize(&inv);
    }

    ibz_vec_2_finalize(&ker);

    // final result: a - i + b*(j+(1+k)/2)
    quat_alg_elem_t gen;
    quat_alg_elem_init(&gen);
    ibz_set(&gen.denom, 2);
    ibz_add(&gen.coord[0], &vec[0], &vec[0]);
    ibz_set(&gen.coord[1], -2);
    ibz_add(&gen.coord[2], &vec[1], &vec[1]);
    ibz_copy(&gen.coord[3], &vec[1]);
    ibz_add(&gen.coord[0], &gen.coord[0], &vec[1]);
    ibz_vec_2_finalize(&vec);

    quat_lideal_create_from_primitive(lideal, &gen, &TORSION_PLUS_23POWER, &MAXORD_O0, &QUATALG_PINFTY);

    assert(0 == ibz_cmp(&lideal->norm, &TORSION_PLUS_23POWER));

    quat_alg_elem_finalize(&gen);
}


