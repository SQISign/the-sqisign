#include "quaternion_tests.h"
#include <rng.h>

// int ibz_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);
int quat_test_randomized_ibz_mat_2x2_inv_mod()
{
    int res = 0;
    ibz_t m, det, gcd;
    ibz_mat_2x2_t a, inv, id, prod;
    int rand[4];
    int64_t rand_m;
    ibz_init(&m);
    ibz_init(&det);
    ibz_init(&gcd);
    ibz_mat_2x2_init(&a);
    ibz_mat_2x2_init(&inv);
    ibz_mat_2x2_init(&prod);
    ibz_mat_2x2_init(&id);
    ibz_mat_2x2_set(&id, 1, 0, 0, 1);


    for (int iter = 0; iter < 100; iter++){
        // generate random matrix and modulo, with modulo larger than 2
        int randret = randombytes((unsigned char*)rand, 4*sizeof(int));
        if (randret != 0)
            return 1;
        ibz_mat_2x2_set(&a, rand[0],rand[1],rand[2],rand[3]);
        rand_m = 0;
        while(rand_m < 2) {
            int randret = randombytes((unsigned char*)&rand_m, sizeof(int64_t));
            if (randret != 0)
                return 1;
        }
        ibz_set(&m,rand_m);

        // compute det
        ibz_mat_2x2_det_from_ibz(&det,&(a[0][0]),&(a[0][1]),&(a[1][0]),&(a[1][1]));
        // is it prime to mod
        ibz_gcd(&gcd,&det,&m);
        if(ibz_is_one(&gcd)){
            // matrix should be invertible mod m
            if (ibz_2x2_inv_mod(&inv, &a, &m))
            {
                // ibz_2x2_mul_mod(&prod,&a,&inv, &m);
                ibz_2x2_mul_mod(&prod, &inv, &a, &m);
                res = res || ibz_cmp(&(prod[0][0]), &(id[0][0]));
                res = res || ibz_cmp(&(prod[0][1]), &(id[0][1]));
                res = res || ibz_cmp(&(prod[1][0]), &(id[1][0]));
                res = res || ibz_cmp(&(prod[1][1]), &(id[1][1]));
            }
            else
            {
                res = 1;
            }
        } else {
            res = res || ibz_2x2_inv_mod(&inv, &a, &m);
        }
    }

    if (res != 0)
    {
        printf("Quaternion unit test with randomization for ibz_mat_2x2_inv_mod failed\n");
    }
    ibz_mat_2x2_finalize(&a);
    ibz_mat_2x2_finalize(&inv);
    ibz_mat_2x2_finalize(&prod);
    ibz_mat_2x2_finalize(&id);
    ibz_finalize(&m);
    ibz_finalize(&det);
    ibz_finalize(&gcd);
    return (res);
}

// int quat_2x2_lattice_enumerate_cvp_filter(quat_alg_elem_t *res, const ibz_mat_2x2_t *lat_basis, const ibz_vec_2_t *target,unsigned int qf, unsigned int dist_bound, int (*condition)(quat_alg_elem_t* , const ibz_vec_2_t*, const void*), const void* params, unsigned int max_tries);
int quat_test_randomized_2x2_lattice_enumerate_cvp_filter()
{
    // test only wether a result is returned which verifies the conditions, if true is returned
    // cannot test wether false is retured correctly or not
    int res = 0;
    ibz_t p, bound, norm_q, norm;
    quat_alg_elem_t cvp_res;
    ibz_mat_2x2_t basis;
    ibz_vec_2_t target, diff, found;
    ibz_mat_2x2_init(&basis);
    quat_alg_elem_init(&cvp_res);
    ibz_vec_2_init(&target);
    ibz_vec_2_init(&diff);
    ibz_vec_2_init(&found);
    ibz_init(&p);
    ibz_init(&norm_q);
    ibz_init(&norm);
    ibz_init(&bound);
    unsigned int q;
    unsigned int dist_bound;
    void *params = (void *)&p;
    unsigned int max_tries;
    int rand[9];

    // fix p since only used inside condition which is internal, unused by any external function
    ibz_set(&p, 3);


    for (int iter = 0; iter < 20; iter++){
        // generate random matrix  with non-0 det
        ibz_set(&bound,0);
        while(ibz_is_zero(&bound)){
            int randret = randombytes((unsigned char*)rand, 9*sizeof(int));
            if (randret != 0)
                return 1;
            ibz_mat_2x2_set(&basis, rand[0],rand[1],rand[2],rand[3]);
            // check det is not 0
            ibz_mat_2x2_det_from_ibz(&bound,&(basis[0][0]),&(basis[0][1]),&(basis[1][0]),&(basis[1][1]));
        }
        // set target
        ibz_vec_2_set(&target, rand[4], rand[5]);
        //set other params randomly in reasonable ranges
        q = ((unsigned int)rand[6]) % 1024;
        if(q == 0){
            q = 1;
        }
        dist_bound = ((unsigned int)rand[7]) % 128;
        max_tries = ((unsigned int)rand[8]) % 16384;
        if (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries)){
            // used cvp_res and condition to exfiltrate standard coords of distance to target from close vector found
            ibz_copy(&(diff[0]), &(cvp_res.coord[0]));
            ibz_copy(&(diff[1]), &(cvp_res.coord[1]));
            // compute close vector in lattice from target and diff
            ibz_sub(&(found[0]), &(target[0]), &(diff[0]));
            ibz_sub(&(found[1]), &(target[1]), &(diff[1]));
            // norm bound on diff ok?
            ibz_set(&norm_q, q);
            ibz_set(&bound, 2);
            ibz_pow(&bound, &bound, dist_bound);
            quat_dim2_lattice_norm(&norm, &(diff[0]), &(diff[1]), &norm_q);
            res = res || !(ibz_cmp(&norm, &bound) < 0);
            // Condition ok
            res = res || (!quat_dim2_lattice_test_cvp_condition(&cvp_res, &diff, params));
            // Is in lattice
            res = res || (!quat_dim2_lattice_contains(&basis, &(found[0]), &(found[1])));
        }
    }

    if (res != 0)
    {
        printf("Quaternion unit test with randomization for 2x2_lattice_enumerate_cvp_filter failed\n");
    }
    ibz_mat_2x2_finalize(&basis);
    quat_alg_elem_finalize(&cvp_res);
    ibz_vec_2_finalize(&target);
    ibz_vec_2_finalize(&diff);
    ibz_vec_2_finalize(&found);
    ibz_finalize(&p);
    ibz_finalize(&norm_q);
    ibz_finalize(&norm);
    ibz_finalize(&bound);
    return res;
}

//int ibz_4x4_inv_with_det_as_denom(ibz_mat_4x4_t *inv, ibz_t *det, const ibz_mat_4x4_t mat);
int quat_test_randomized_ibz_mat_4x4_inv_with_det_as_denom(){
    int res = 0;
    ibz_t det;
    ibz_mat_4x4_t mat, inv;
    ibz_mat_4x4_t det_id, prod;
    ibz_mat_4x4_init(&det_id);
    ibz_mat_4x4_init(&prod);
    ibz_init(&det);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&inv);

    for (int r = 0; r < 10; r++) {
        do {
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    ibz_rand_interval_i(&mat[i][j], -(1 << 20), 1 << 20);
        } while (!ibz_mat_4x4_inv_with_det_as_denom(&inv, &det, &mat));
        ibz_mat_4x4_identity(&det_id);
        ibz_mat_4x4_scalar_mul(&det_id,&det,&det_id);
        ibz_mat_4x4_mul(&prod, &inv, &mat);
        res = res || !ibz_mat_4x4_equal(&det_id,&prod);
        ibz_mat_4x4_mul(&prod, &mat, &inv);
        res = res || !ibz_mat_4x4_equal(&det_id,&prod);
    }
    
    if (res != 0){
        printf("Quaternion unit test with randomization for ibz_mat_4x4_inv_with_det_as_denom failed\n");
    }
    ibz_mat_4x4_finalize(&det_id);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&inv);
    ibz_finalize(&det);
    return res;
}

//void ibz_mat_4x8_hnf_core(ibz_mat_4x4_t *hnf, const ibz_mat_4x8_t *generators);
int quat_test_randomized_ibz_mat_4x8_hnf_core(){
    // only partial test, since lattice equality cannot be tested without hnf, so only one inclusion is checked
    int res = 0;
    quat_alg_elem_t vec;
    quat_lattice_t lat;
    ibz_mat_4x8_t mat;
    ibz_mat_4x4_t hnf,cmp;
    quat_lattice_init(&lat);
    quat_alg_elem_init(&vec);
    ibz_mat_4x8_init(&mat);
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&cmp);
    int64_t rand[8][4];
    int det_non_0;

    for (int iter = 0; iter < 100; iter++){
        int randret = randombytes((unsigned char*)rand, 8*4*sizeof(int64_t));
        if (randret != 0)
            return 1;

        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 8; j++){
                ibz_set(&(mat[i][j]),rand[j][i]);
            }
        }
        ibz_mat_4x8_hnf_core(&hnf,&mat);
        res = res || (!ibz_mat_4x4_is_hnf(&hnf));
        // also should test that they generate the same lattice. However, can only check one inclusion efficiently (and this only if full rank), so do so
        det_non_0 = 1;
        for(int i = 0; i <4; i ++){
            det_non_0 = det_non_0 && ibz_is_zero(&(hnf[i][i]));
        }
        if(det_non_0){
            ibz_mat_4x4_copy(&(lat.basis),&hnf);
            ibz_set(&(lat.denom),1);
            for(int i = 0; i <8; i ++){
                quat_alg_elem_copy_ibz(&vec,&(lat.denom), &(mat[i][0]), &(mat[i][1]), &(mat[i][2]), &(mat[i][3]));
                res = res || !quat_lattice_contains_without_alg(NULL,&lat,&vec);
            }
        }
    }

    if (res != 0){
        printf("Quaternion unit test with randomization for ibz_mat_4x8_hnf_core failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&vec);
    ibz_mat_4x8_finalize(&mat);
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&cmp);
    return(res);
}

//int ibz_4x5_right_ker_mod_prime(ibz_vec_5_t *ker, const ibz_mat_4x5_t *mat, const ibz_t *p);
int quat_test_randomized_ibz_4x5_right_ker_mod_prime(){
    // check that if the function returns 1, the vector is in the kernel
    int res = 0;
    ibz_t prime;
    ibz_mat_4x5_t mat;
    ibz_vec_5_t ker;
    ibz_t  prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_mat_4x5_init(&mat);
    ibz_vec_5_init(&ker);
    ibz_init(&prime);
    int64_t rand[5][4];
    int64_t rand_p;


    for (int iter = 0; iter < 100; iter++){
        // generate random matrix and modulo, with modulo larger than 2
        int randret = randombytes((unsigned char*)rand, 4*5*sizeof(int64_t));
        if (randret != 0)
            return 1;
        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 5; j++){
                ibz_set(&(mat[i][j]),rand[j][i]);
            }
        }
        rand_p = 0;
        while(rand_p ==0) {
            int randret = randombytes((unsigned char*)&rand_p, sizeof(int64_t));
            if (randret != 0)
                return 1;
            ibz_set(&prime,rand_p);
            if(!ibz_probab_prime(&prime,20))
                rand_p = 0;
        }
        if (ibz_4x5_right_ker_mod_prime(&ker,&mat,&prime)){
            for(int i = 0; i < 4; i++){
                ibz_set(&sum,0);
                for(int j = 0; j < 5; j++){
                    ibz_mul(&prod,&(mat[i][j]),&(ker[j]));
                    ibz_add(&sum,&sum,&prod);
                    ibz_mod(&sum,&sum,&prime);
                }
                res = res || !ibz_is_zero(&sum);
            }
        }
    }

    if (res != 0){
        printf("Quaternion unit test with randomization for ibz_4x5_right_ker_mod_prime failed\n");
    }
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    ibz_vec_5_finalize(&ker);
    ibz_mat_4x5_finalize(&mat);
    ibz_finalize(&prime);
    return(res);
}


//int ibz_4x4_right_ker_mod_prime(ibz_vec_4_t *ker, const ibz_mat_4x4_t *mat, const ibz_t *p);
int quat_test_randomized_ibz_4x4_right_ker_mod_prime(){
    int res = 0;
    ibz_t prime;
    ibz_mat_4x4_t mat, rank3, rank2;
    ibz_vec_4_t ker;
    ibz_t  prod, sum,det;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_init(&det);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&rank3);
    ibz_mat_4x4_init(&rank2);
    ibz_vec_4_init(&ker);
    ibz_init(&prime);
    int64_t rand[4][4];
    int64_t rand_p;
    ibz_mat_4x4_identity(&rank3);
    ibz_set(&(rank3[3][3]),0);
    ibz_mat_4x4_identity(&rank2);
    ibz_set(&(rank2[3][3]),0);
    ibz_set(&(rank2[2][2]),0);

    // test case where always 1 solution
    for (int iter = 0; iter < 100; iter++){
        rand_p = 0;
        while(rand_p ==0) {
            int randret = randombytes((unsigned char*)&rand_p, sizeof(int64_t));
            if (randret != 0)
                return 1;
            ibz_set(&prime,rand_p);
            if(!ibz_probab_prime(&prime,20))
                rand_p = 0;
        }
        // generate random invertible matrix
        ibz_set(&det,0);
        while(ibz_is_zero(&det)){
            int randret = randombytes((unsigned char*)rand, 4*4*sizeof(int64_t));
            if (randret != 0)
                return 1;
            for (int i = 0; i < 4; i++){
                for (int j = 0; j < 4; j++){
                    ibz_set(&(mat[i][j]),rand[j][i]);
                }
            }
            ibz_mat_4x4_inv_with_det_as_denom(NULL,&det,&mat);
            ibz_mod(&det,&det,&prime);
        }
        // rank 4 matrix does not work
        res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);
        // multiply with rank3 to get random rank 3 matrix
        ibz_mat_4x4_mul(&mat,&mat,&rank3);

        if (ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime)){
            for(int i = 0; i < 4; i++){
                ibz_set(&sum,0);
                for(int j = 0; j < 4; j++){
                    ibz_mul(&prod,&(mat[i][j]),&(ker[j]));
                    ibz_add(&sum,&sum,&prod);
                    ibz_mod(&sum,&sum,&prime);
                }
                res = res || !ibz_is_zero(&sum);
            }
        } else {
            res = 1;
        }
        // multiply with rank2 to get matrix of too small rank
        ibz_mat_4x4_mul(&mat,&mat,&rank2);
        res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);
        ibz_mat_4x4_mul(&mat,&mat,&rank2);
        res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);
        ibz_mat_4x4_mul(&mat,&mat,&rank2);
        res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);
    }

    if (res != 0){
        printf("Quaternion unit test with randomization for ibz_4x4_right_ker_mod_prime failed\n");
    }
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    ibz_finalize(&det);
    ibz_vec_4_finalize(&ker);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&rank3);
    ibz_mat_4x4_finalize(&rank2);
    ibz_finalize(&prime);
    return(res);
}

//int ibz_4x4_right_ker_mod_power_of_2(ibz_vec_4_t *ker, const ibz_mat_4x4_t *mat, unsigned short exp);
int quat_test_randomized_ibz_4x4_right_ker_mod_power_of_2(){
    // this only tests that if a vector is returned, it fullfills the conditions
    int res = 0;
    int zero = 1;
    short exp;
    int64_t rand[4][4];
    unsigned char rand_exp;
    ibz_mat_4x4_t mat, rank3;
    ibz_vec_4_t ker;
    ibz_vec_4_t prod;
    ibz_t q, r, two, det;
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&rank3);
    ibz_vec_4_init(&ker);
    ibz_vec_4_init(&prod);
    ibz_init(&q);
    ibz_init(&r);
    ibz_init(&det);
    ibz_init(&two);
    ibz_set(&two,2);
    ibz_mat_4x4_identity(&rank3);
    ibz_set(&(rank3[3][3]),0);

    for (int iter = 0; iter < 100; iter++){
        rand_exp = 0;
        while(rand_exp <=0) {
            int randret = randombytes((unsigned char*)&rand_exp, sizeof(unsigned char));
            if (randret != 0)
                return 1;
            exp = (short)rand_exp;
        }
        // generate random invertible matrix
        ibz_set(&det,0);
        while(ibz_is_zero(&det)){
            int randret = randombytes((unsigned char*)rand, 4*4*sizeof(int64_t));
            if (randret != 0)
                return 1;
            for (int i = 0; i < 4; i++){
                for (int j = 0; j < 4; j++){
                    ibz_set(&(mat[i][j]),rand[j][i]);
                }
            }
            ibz_mat_4x4_inv_with_det_as_denom(NULL,&det,&mat);
            ibz_div(&q,&det,&det,&two);
        }


        // rank 4 matrix does not work
        res = res || ibz_4x4_right_ker_mod_power_of_2(&ker,&mat,exp);
        // multiply with rank3 to get random rank 3 matrix
        ibz_mat_4x4_mul(&mat,&mat,&rank3);

        if (ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, exp)){
            zero = 1;
            ibz_mat_4x4_eval(&prod,&mat,&ker);
            for (int i = 0; i < 4; i++){
                ibz_div(&q,&r,&(ker[i]),&two);
                zero = zero && ibz_is_zero(&r);
                ibz_pow(&q,&two,exp);
                ibz_mod(&(prod[i]),&(prod[i]),&q);
            }
            res |= !quat_alg_coord_is_zero(&prod);
            res |= zero;
        }
    }

    
    if (res != 0){
        printf("Quaternion unit test with randomization for ibz_4x4_right_ker_mod_power_of_2 failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&rank3);
    ibz_vec_4_finalize(&ker);
    ibz_vec_4_finalize(&prod);
    ibz_finalize(&q);
    ibz_finalize(&r);
    ibz_finalize(&two);
    ibz_finalize(&det);
    return(res);
}


//int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int precision);
int quat_test_randomized_lattice_lll(){
    int res = 0;
    quat_lattice_t lat, test;
    ibz_mat_4x4_t red;
    ibz_t num, denom, q, det;
    ibq_t coeff;
    int64_t rand[4][4];
    uint32_t rand_q, rand_prec, rand_denom;
    ibz_init(&num);
    ibz_init(&denom);
    ibz_init(&q);
    ibz_init(&det);
    ibq_init(&coeff);
    ibz_mat_4x4_init(&red);
    quat_lattice_init(&lat);
    quat_lattice_init(&test);
    ibz_set(&num,3);
    ibz_set(&denom,4);
    ibq_set(&coeff,&num,&denom);

    for (int iter = 0; iter < 20; iter++){
        rand_denom = 0;
        while(rand_denom ==0) {
            int randret = randombytes((unsigned char*)&rand_denom, sizeof(uint32_t));
            if (randret != 0)
                return 1;
        }
        int randret = randombytes((unsigned char*)&rand_q, sizeof(uint32_t));
        if (randret != 0)
            return 1;
        // generate random invertible matrix
        ibz_set(&det,0);
        while(ibz_is_zero(&det)){
            int randret = randombytes((unsigned char*)rand, 4*4*sizeof(int64_t));
            if (randret != 0)
                return 1;
            for (int i = 0; i < 4; i++){
                for (int j = 0; j < 4; j++){
                    ibz_set(&(lat.basis[i][j]),rand[j][i]);
                }
            }
            ibz_mat_4x4_inv_with_det_as_denom(NULL,&det,&(lat.basis));
        }

        // set lattice
        ibz_set(&lat.denom, rand_denom);
        quat_lattice_hnf(&lat);
        // set other parameter
        ibz_set(&q,rand_q % 1024);
        //reduce
        res = res || quat_lattice_lll(&red,&lat,&q,0);
        // test lll reduced
        res = res || !quat_dim4_lll_verify(&red,&coeff,&q);
        // test lattice equality
        ibz_copy(&(test.denom),&(lat.denom));
        ibz_mat_4x4_copy(&(test.basis),&(red));
        quat_lattice_hnf(&test);
        res = res || !quat_lattice_equal(&test,&lat);
    }
    
    if (res != 0){
        printf("Quaternion unit test with randomization for lattice_lll failed\n");
    }
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibz_finalize(&q);
    ibz_finalize(&det);
    ibq_finalize(&coeff);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&test);
    return(res);
}

//int quat_lattice_contains_without_alg(quat_alg_coord_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t *x);
int quat_test_randomized_lattice_contains_without_alg(){
    // only tests the case where the element is in the lattice
    int res = 0;
    ibz_t det;
    quat_alg_elem_t x, cmp;
    quat_alg_coord_t coord;
    quat_lattice_t lat;
    int64_t rand_denom;
    int64_t rand[5][4];
    ibz_init(&det);
    quat_alg_coord_init(&coord);
    quat_alg_elem_init(&cmp);
    quat_alg_elem_init(&x);
    quat_lattice_init(&lat);

    for (int iter = 0; iter < 10; iter++){
        rand_denom = 0;
        while(rand_denom ==0) {
            int randret = randombytes((unsigned char*)&rand_denom, sizeof(int64_t));
            if (randret != 0)
                return 1;
        }
        // generate random invertible matrix
        ibz_set(&det,0);
        while(ibz_is_zero(&det)){
            int randret = randombytes((unsigned char*)rand, 5*4*sizeof(int64_t));
            if (randret != 0)
                return 1;
            for (int i = 0; i < 4; i++){
                for (int j = 0; j < 4; j++){
                    ibz_set(&(lat.basis[i][j]),rand[j][i]);
                }
                ibz_set(&(coord[i]),rand[4][i]);
            }
            ibz_mat_4x4_inv_with_det_as_denom(NULL,&det,&(lat.basis));
        }

        ibz_set(&(lat.denom),rand_denom);
        quat_lattice_hnf(&lat);
        ibz_mat_4x4_eval(&(x.coord),&(lat.basis),&coord);
        ibz_copy(&(x.denom),&(lat.denom));
        ibz_vec_4_set(&coord,1,0,1,0);
        if(quat_lattice_contains_without_alg(&coord,&lat,&x)){
            ibz_mat_4x4_eval(&(cmp.coord),&(lat.basis),&coord);
            ibz_copy(&(cmp.denom),&(lat.denom));
            quat_alg_sub(&cmp,&x,&cmp);
            res = res || !quat_alg_elem_is_zero(&cmp);
        } else {
            res = 1;
        }

    }

    if (res != 0){
        printf("Quaternion unit test with randomization for lattice_contains_without_alg failed\n");
    }
    ibz_finalize(&det);
    quat_alg_coord_finalize(&coord);
    quat_alg_elem_finalize(&x);
    quat_alg_elem_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return(res);
}


//int ibz_cornacchia_extended(ibz_t *x, ibz_t *y, const ibz_t *n, const short *prime_list, const int prime_list_length, short primality_test_iterations, const ibz_t *bad_primes_prod); 
int quat_test_randomized_ibz_cornacchia_extended(){
    int res = 0;
    ibz_t x,y,n, prod,c_res,bad, p;
    int counter = 1;
    int counter_p = 3;
    int counter_good = 1;
    short rand_exps[100];
    int64_t rand_fact;
    short primes[100]; 	 	 
    short good_primes[100]; 	 	 	 
    int primes_length = 100;
    short iterations = 20;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&n);
    ibz_init(&p);
    ibz_init(&prod);
    ibz_init(&c_res);
    ibz_init(&bad);
    ibz_set(&bad,1);
    good_primes[0] = 2;
    primes[0] = 2;
    while(counter <100){
        ibz_set(&p,counter_p);
        if(ibz_probab_prime(&p,20)){
            primes[counter]=(short) counter_p;
            if(counter_p%4 == 3){
                ibz_mul(&bad,&bad,&p);
            } else {
                good_primes[counter_good]=(short) counter_p;
                counter_good++;
            }
            counter++;
        }
        counter_p++;
    }

    for (int iter = 0; iter < 1000; iter++){
        ibz_set(&n,1);
        int randret = randombytes((unsigned char*)&rand_exps, 100*sizeof(short));
        for(int i = 0; i  < counter_good; i++){
            ibz_set(&p,good_primes[i]);
            ibz_pow(&p,&p,((int)((unsigned short)rand_exps[i])%16));
            ibz_mul(&n,&n,&p);
        }
        // sample large prime as last factor
        rand_fact = 0;
        while(rand_fact == 0) {
            int randret = randombytes((unsigned char*)&rand_fact, sizeof(int64_t));
            if (randret != 0)
                return 1;
            if(rand_fact < 0) 
                rand_fact = - rand_fact;
            ibz_set(&p,rand_fact);
            if(!ibz_probab_prime(&p,20))
                rand_fact = 0;
        }
        ibz_mul(&n,&n,&p);

        ibz_set(&prod,4);
        ibz_mod(&prod,&p,&prod);
        if((ibz_probab_prime(&p,100) || (ibz_is_one(&p))) && (ibz_is_one(&prod))){
            res = res || (!ibz_cornacchia_extended(&x,&y,&n,primes,100,100,&bad));
            if(res){
                ibz_mul(&c_res,&x,&x);
                ibz_mul(&prod,&y,&y);
                ibz_add(&c_res,&c_res,&prod);
                res = res || ibz_cmp(&n,&c_res);
            }
            //(this test depends on the primality test)
            // now it is not a prime factor any more
            ibz_set(&p,rand_fact);
            ibz_mul(&n,&n,&p);
            res = res || ibz_cornacchia_extended(&x,&y,&n,primes,100,1000,&bad);
            // multiply with random bad primes
            ibz_div(&n,&p,&n,&p);
            int randret = randombytes((unsigned char*)&rand_exps, 100*sizeof(short));
            if (randret != 0)
                return 1;
            for(int i = 0; i  < 100; i++){
                ibz_set(&p,primes[i]);
                ibz_pow(&p,&p,((int)((unsigned short)rand_exps[i])%8));
                ibz_mul(&n,&n,&p);
            }
            ibz_gcd(&p,&n,&bad);
            if(ibz_is_one(&p)){
                res = res || !ibz_cornacchia_extended(&x,&y,&n,primes,100,100,&bad);
                if(res){
                    ibz_mul(&c_res,&x,&x);
                    ibz_mul(&prod,&y,&y);
                    ibz_add(&c_res,&c_res,&prod);
                    res = res || ibz_cmp(&n,&c_res);
                }
            } else {
                res = res || ibz_cornacchia_extended(&x,&y,&n,primes,100,100,&bad);
            } 
        } else {
            res = res || ibz_cornacchia_extended(&x,&y,&n,primes,100,1000,&bad);
            for(int i = 0; i  < 100; i++){
                ibz_set(&p,primes[i]);
                ibz_pow(&p,&p,((int)((unsigned short)rand_exps[i])%8));
                ibz_mul(&n,&n,&p);
            }
            ibz_gcd(&p,&n,&bad);
            res = res || ibz_cornacchia_extended(&x,&y,&n,primes,100,100,&bad);
        } 
    }

    if (res != 0){
        printf("Quaternion unit test with randomization for ibz_cornacchia_extended failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&n);
    ibz_finalize(&p);
    ibz_finalize(&prod);
    ibz_finalize(&c_res);
    ibz_finalize(&bad);
    return res;
}


// run all previous tests
int quat_test_with_randomization(){
    int res = 0;
    printf("\nRunning randomized tests from quaternion module\n");
    res = res | quat_test_randomized_ibz_mat_2x2_inv_mod();
    res = res | quat_test_randomized_2x2_lattice_enumerate_cvp_filter();
    res = res | quat_test_randomized_ibz_mat_4x4_inv_with_det_as_denom();
    res = res | quat_test_randomized_ibz_mat_4x8_hnf_core();
    res = res | quat_test_randomized_ibz_4x5_right_ker_mod_prime();
    res = res | quat_test_randomized_ibz_4x4_right_ker_mod_prime();
    res = res | quat_test_randomized_ibz_4x4_right_ker_mod_power_of_2();
    res = res | quat_test_randomized_lattice_lll();
    res = res | quat_test_randomized_lattice_contains_without_alg();
    res = res | quat_test_randomized_ibz_cornacchia_extended();
    return(res);
}
