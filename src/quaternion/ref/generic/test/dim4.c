#include "quaternion_tests.h"

// internal helper function
//void ibz_mat_4x4_mul(ibz_mat_4x4_t *res, const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b)
int quat_test_dim4_ibz_mat_4x4_mul(){
    int res = 0;
    ibz_mat_4x4_t a, b, prod, cmp;
    ibz_mat_4x4_init(&a);
    ibz_mat_4x4_init(&b);
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&cmp);
    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            ibz_set(&(a[i][j]), 0);
            ibz_set(&(b[i][j]), 0);
            ibz_set(&(cmp[i][j]), 0);
        }
    }

    ibz_set(&(a[0][0]), 1);
    ibz_set(&(a[0][1]), 2);
    ibz_set(&(a[0][2]), 1);
    ibz_set(&(a[1][1]), 1);
    ibz_set(&(a[1][2]), 3);
    ibz_set(&(a[2][2]), 1);
    ibz_set(&(a[2][3]), 4);
    ibz_set(&(a[3][3]), 1);
    ibz_set(&(b[0][0]), -1);
    ibz_set(&(b[1][0]), 1);
    ibz_set(&(b[1][1]), -2);
    ibz_set(&(b[2][0]), 1);
    ibz_set(&(b[2][2]), 3);
    ibz_set(&(b[3][0]), 1);
    ibz_set(&(b[3][1]), 5);
    ibz_set(&(b[3][2]), -1);
    ibz_set(&(b[3][3]), 2);
    ibz_set(&(cmp[0][0]), 4);
    ibz_set(&(cmp[0][1]), -4);
    ibz_set(&(cmp[0][2]), 1);
    ibz_set(&(cmp[1][0]), 10);
    ibz_set(&(cmp[1][1]), -2);
    ibz_set(&(cmp[1][2]), 3);
    ibz_set(&(cmp[2][0]), 7);
    ibz_set(&(cmp[2][1]), 20);
    ibz_set(&(cmp[2][2]), -3);
    ibz_set(&(cmp[2][3]), 8);
    ibz_set(&(cmp[3][0]), 1);
    ibz_set(&(cmp[3][1]), 5);
    ibz_set(&(cmp[3][2]), -1);
    ibz_set(&(cmp[3][3]), 2);
    ibz_mat_4x4_mul(&prod, &a, &b);
    res = res || ibz_mat_4x4_equal(&cmp,&prod);

    ibz_mat_4x4_mul(&b, &a, &b);
    res = res || ibz_mat_4x4_equal(&b,&cmp);

    ibz_set(&(b[0][0]), -1);
    ibz_set(&(b[1][0]), 1);
    ibz_set(&(b[1][1]), -2);
    ibz_set(&(b[2][0]), 1);
    ibz_set(&(b[2][2]), 3);
    ibz_set(&(b[3][0]), 1);
    ibz_set(&(b[3][1]), 5);
    ibz_set(&(b[3][2]), -1);
    ibz_set(&(b[3][3]), 2);
    ibz_mat_4x4_mul(&a, &a, &b);
    res = res || ibz_mat_4x4_equal(&a,&cmp);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_mul failed\n");
    }
    ibz_mat_4x4_finalize(&a);
    ibz_mat_4x4_finalize(&b);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&cmp);
    return res;
}

// helper functions for lattices
//void ibz_vec_4_set(ibz_vec_4_t *vec, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3)
int quat_test_dim4_ibz_vec_4_set(){
    int res = 0;
    ibz_vec_4_t a;
    ibz_vec_4_init(&a);
    ibz_vec_4_set(&a,1,2,3,4);
    res = res || !(1 == ibz_get(&(a[0])));
    res = res || !(2 == ibz_get(&(a[1])));
    res = res || !(3 == ibz_get(&(a[2])));
    res = res || !(4 == ibz_get(&(a[3])));
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_vec_4_set failed\n");
    }
    ibz_vec_4_finalize(&a);
    return(res);
}


//void ibz_vec_4_copy(ibz_vec_4_t *new, const ibz_vec_4_t  *vec);
int quat_test_dim4_ibz_vec_4_copy(){
    int res = 0;
    ibz_vec_4_t a,b;
    ibz_vec_4_init(&a);
    ibz_vec_4_init(&b);
    ibz_vec_4_set(&a,1,2,3,4);
    ibz_vec_4_copy(&b,&a);
    res = res || !(1 == ibz_get(&(b[0])));
    res = res || !(2 == ibz_get(&(b[1])));
    res = res || !(3 == ibz_get(&(b[2])));
    res = res || !(4 == ibz_get(&(b[3])));
    ibz_vec_4_copy(&a,&a);
    res = res || !(1 == ibz_get(&(a[0])));
    res = res || !(2 == ibz_get(&(a[1])));
    res = res || !(3 == ibz_get(&(a[2])));
    res = res || !(4 == ibz_get(&(a[3])));

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_vec_4_copy failed\n");
    }
    ibz_vec_4_finalize(&a);
    ibz_vec_4_finalize(&b);
    return(res);
}


//void ibz_vec_4_negate(ibz_vec_4_t *neg, const ibz_vec_4_t  *vec);
int quat_test_dim4_ibz_vec_4_negate(){
    int res = 0;
    ibz_vec_4_t a,b;
    ibz_vec_4_init(&a);
    ibz_vec_4_init(&b);
    ibz_vec_4_set(&a,1,2,3,4);
    ibz_vec_4_negate(&b,&a);
    res = res || !(-1 == ibz_get(&(b[0])));
    res = res || !(-2 == ibz_get(&(b[1])));
    res = res || !(-3 == ibz_get(&(b[2])));
    res = res || !(-4 == ibz_get(&(b[3])));
    ibz_vec_4_negate(&a,&a);
    res = res || !(-1 == ibz_get(&(a[0])));
    res = res || !(-2 == ibz_get(&(a[1])));
    res = res || !(-3 == ibz_get(&(a[2])));
    res = res || !(-4 == ibz_get(&(a[3])));
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_vec_4_negate failed\n");
    }
    ibz_vec_4_finalize(&a);
    ibz_vec_4_finalize(&b);
    return(res);
}

//void ibz_vec_4_linear_combination(ibz_vec_4_t *lc, const ibz_t *coeff_a, const ibz_vec_4_t  *vec_a, const ibz_t *coeff_b, const ibz_vec_4_t *vec_b){
int quat_test_dim4_ibz_vec_4_linear_combination(){
    int res = 0;
    ibz_vec_4_t a, b, lc, cmp;
    ibz_t ca, cb;
    ibz_init(&ca);
    ibz_init(&cb);
    ibz_vec_4_init(&a);
    ibz_vec_4_init(&b);
    ibz_vec_4_init(&lc);
    ibz_vec_4_init(&cmp);
    ibz_vec_4_set(&a,1,2,3,4);
    ibz_vec_4_set(&b,-2,1,3,-3);
    ibz_set(&ca,2);
    ibz_set(&cb,-1);
    ibz_vec_4_set(&cmp,4,3,3,11);
    ibz_vec_4_linear_combination(&lc,&ca,&a,&cb,&b);
    for(int i = 0; i < 4; i++){
        res = res || ibz_cmp(&(lc[i]),&(cmp[i]));
    }
    ibz_vec_4_set(&cmp,1,2,3,4);
    ibz_vec_4_linear_combination(&a,&ca,&a,&cb,&a);
    for(int i = 0; i < 4; i++){
        res = res || ibz_cmp(&(a[i]),&(cmp[i]));
    }
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_vec_4_linear_combination failed\n");
    }
    ibz_finalize(&ca);
    ibz_finalize(&cb);
    ibz_vec_4_finalize(&a);
    ibz_vec_4_finalize(&b);
    ibz_vec_4_finalize(&lc);
    ibz_vec_4_finalize(&cmp);
    return(res);
}

//int ibz_vec_4_scalar_div(ibz_vec_4_t *quot, const ibz_t *scalar, const ibz_vec_4_t *vec);
int quat_test_dim4_ibz_vec_4_scalar_div(){
    int res = 0;
    int s;
    ibz_t scalar;
    ibz_vec_4_t quot,vec;
    ibz_vec_4_init(&vec);
    ibz_vec_4_init(&quot);
    ibz_init(&scalar);

    s = 5;
    ibz_set(&scalar,s);
    for(int i = 0; i < 4; i++){
        ibz_set(&(vec[i]),(i)*s);
    }
    res = res || !ibz_vec_4_scalar_div(&quot,&scalar,&vec);
    for(int i = 0; i < 4; i++){
        res = res || (ibz_get(&(quot[i])) !=i);
    }

    res = res || ibz_vec_4_scalar_div(&vec,&scalar,&vec);
    for(int i = 0; i < 4; i++){
        res = res || (ibz_get(&(vec[i])) !=i);
    }

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_vec_4_scalar_div failed\n");
    }
    ibz_vec_4_finalize(&vec);
    ibz_vec_4_finalize(&quot);
    ibz_finalize(&scalar);
    return(res);
}

//void ibz_mat_4x4_copy(ibz_mat_4x4_t *new, const ibz_mat_4x4_t *mat);
int quat_test_dim4_ibz_mat_4x4_copy(){
    int res = 0;
    ibz_mat_4x4_t mat, new;
    ibz_mat_4x4_init(&new);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_zero(&mat);
    ibz_set(&(mat[0][0]),1);
    ibz_set(&(mat[0][1]),2);
    ibz_set(&(mat[0][2]),-7);
    ibz_set(&(mat[0][3]), 77);
    ibz_set(&(mat[2][0]),13);
    ibz_set(&(mat[1][1]),20);
    ibz_set(&(mat[3][2]),-77);
    ibz_set(&(mat[3][3]), 7);
    ibz_mat_4x4_copy(&new,&mat);
    res = res || !ibz_mat_4x4_equal(&new,&mat);
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_copy failed\n");
    }
    ibz_mat_4x4_finalize(&new);
    ibz_mat_4x4_finalize(&mat);
    return(res);
}

//void ibz_mat_4x4_negate(ibz_mat_4x4_t *neg, const ibz_mat_4x4_t *mat);
int quat_test_dim4_ibz_mat_4x4_negate(){
    int res = 0;
    ibz_mat_4x4_t mat, neg,cmp;
    ibz_mat_4x4_init(&neg);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_zero(&cmp);
    ibz_mat_4x4_zero(&mat);
    ibz_set(&(mat[0][0]),1);
    ibz_set(&(cmp[0][0]),-1);
    ibz_set(&(mat[0][1]),2);
    ibz_set(&(cmp[0][1]),-2);
    ibz_set(&(mat[0][2]),-7);
    ibz_set(&(cmp[0][2]),7);
    ibz_set(&(mat[0][3]), 77);
    ibz_set(&(cmp[0][3]), -77);
    ibz_set(&(mat[2][0]),13);
    ibz_set(&(cmp[2][0]),-13);
    ibz_set(&(mat[1][1]),20);
    ibz_set(&(cmp[1][1]),-20);
    ibz_set(&(mat[3][2]),-77);
    ibz_set(&(cmp[3][2]),77);
    ibz_set(&(mat[3][3]), 7);
    ibz_set(&(cmp[3][3]), -7);
    ibz_mat_4x4_negate(&neg,&mat);
    res = res || !ibz_mat_4x4_equal(&neg,&cmp);
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_negate failed\n");
    }
    ibz_mat_4x4_finalize(&neg);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&cmp);
    return(res);
}

//void ibz_mat_4x4_transpose(ibz_mat_4x4_t *transposed, const ibz_mat_4x4_t *mat)
int quat_test_dim4_ibz_mat_4x4_transpose(){
    int res = 0;
    ibz_mat_4x4_t mat, transposed, cmp;
    ibz_mat_4x4_init(&transposed);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_zero(&mat);
    ibz_mat_4x4_zero(&cmp);
    ibz_set(&(mat[0][0]),1);
    ibz_set(&(cmp[0][0]),1);
    ibz_set(&(mat[0][1]),2);
    ibz_set(&(cmp[1][0]),2);
    ibz_set(&(mat[0][2]),-7);
    ibz_set(&(cmp[2][0]),-7);
    ibz_set(&(mat[0][3]), 77);
    ibz_set(&(cmp[3][0]), 77);
    ibz_set(&(mat[2][0]),13);
    ibz_set(&(cmp[0][2]),13);
    ibz_set(&(mat[1][1]),20);
    ibz_set(&(cmp[1][1]),20);
    ibz_set(&(mat[3][2]),-77);
    ibz_set(&(cmp[2][3]),-77);
    ibz_set(&(mat[3][3]), 7);
    ibz_set(&(cmp[3][3]), 7);
    ibz_mat_4x4_transpose(&transposed,&mat);
    res = res || !ibz_mat_4x4_equal(&transposed,&cmp);
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_transpose failed\n");
    }
    ibz_mat_4x4_finalize(&transposed);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&cmp);
    return(res);
}

//void ibz_mat_4x4_zero(ibz_mat_4x4_t *zero);
int quat_test_dim4_ibz_mat_4x4_zero(){
    int res = 0;
    ibz_mat_4x4_t mat, cmp;
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_init(&mat);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(cmp[i][j]),0);
        }
    }
    ibz_mat_4x4_zero(&mat);
    res = res || !ibz_mat_4x4_equal(&cmp,&mat);
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_zero failed\n");
    }
    ibz_mat_4x4_finalize(&cmp);
    ibz_mat_4x4_finalize(&mat);
    return(res);
}

//void ibz_mat_4x4_identity(ibz_mat_4x4_t *id);
int quat_test_dim4_ibz_mat_4x4_identity(){
    int res = 0;
    ibz_mat_4x4_t mat, cmp;
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_init(&mat);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(cmp[i][j]),0);
        }
        ibz_set(&(cmp[i][i]),1);
    }
    ibz_mat_4x4_identity(&mat);
    res = res || !ibz_mat_4x4_equal(&cmp,&mat);
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_identity failed\n");
    }
    ibz_mat_4x4_finalize(&cmp);
    ibz_mat_4x4_finalize(&mat);
    return(res);
}

//int ibz_mat_4x4_is_identity(const ibz_mat_4x4_t *mat);
int quat_test_dim4_ibz_mat_4x4_is_identity(){
    int res = 0;
    ibz_mat_4x4_t mat;
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_identity(&mat);
    res = res || !ibz_mat_4x4_is_identity(&mat);
    ibz_set(&(mat[0][1]),1);
    res = res || ibz_mat_4x4_is_identity(&mat);
    ibz_set(&(mat[0][1]),0);
    ibz_set(&(mat[3][3]),0);
    res = res || ibz_mat_4x4_is_identity(&mat);
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_is_identity failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    return(res);
}

//int ibz_mat_4x4_equal(const ibz_mat_4x4_t *mat1, const ibz_mat_4x4_t *mat2);
int quat_test_dim4_ibz_mat_4x4_equal(){
    int res = 0;
    ibz_mat_4x4_t a,b;
    ibz_mat_4x4_init(&a);
    ibz_mat_4x4_init(&b);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(a[i][j]),i+j);
            ibz_set(&(b[i][j]),i+j);
        }
    }
    res = res || (!ibz_mat_4x4_equal(&a,&b));

    ibz_set(&(b[2][2]),2);
    res = res  || ibz_mat_4x4_equal(&a,&b);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_equal failed\n");
    }
    ibz_mat_4x4_finalize(&a);
    ibz_mat_4x4_finalize(&b);
    return(res);
}

//void ibz_mat_4x4_scalar_mul(ibz_mat_4x4_t *prod, const ibz_t *scalar, const ibz_mat_4x4_t *mat);
int quat_test_dim4_ibz_mat_4x4_scalar_mul(){
    int res = 0;
    int s;
    ibz_t scalar;
    ibz_mat_4x4_t prod,mat,cmp;
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_init(&prod);
    ibz_init(&scalar);

    s = 5;
    ibz_set(&scalar,s);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),i+j);
            ibz_set(&(cmp[i][j]),(i+j)*s);
        }
    }
    ibz_mat_4x4_scalar_mul(&prod,&scalar,&mat);
    res = res || (!ibz_mat_4x4_equal(&prod,&cmp));


    ibz_mat_4x4_scalar_mul(&mat,&scalar,&mat);
    res = res || (!ibz_mat_4x4_equal(&mat,&cmp));

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_scalar_mul failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&cmp);
    ibz_mat_4x4_finalize(&prod);
    ibz_finalize(&scalar);
    return(res);
}

// void ibz_mat_4x4_gcd(ibz_t *gcd, const ibz_mat_4x4_t *mat);
int quat_test_dim4_ibz_mat_4x4_gcd(){
    int res = 0;
    int d;
    ibz_t cmp, gcd;
    ibz_mat_4x4_t mat;
    ibz_mat_4x4_init(&mat);
    ibz_init(&cmp);
    ibz_init(&gcd);
    
    d = 2;
    ibz_set(&cmp,d);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),d*i*j);
        }
    }
    ibz_mat_4x4_gcd(&gcd,&mat);
    res = res || ibz_cmp(&gcd,&cmp);

    d = 21;
    ibz_set(&cmp,d);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),d*i*j);
        }
    }
    ibz_mat_4x4_gcd(&gcd,&mat);
    res = res || ibz_cmp(&gcd,&cmp);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_gcd failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_finalize(&cmp);
    ibz_finalize(&gcd);
    return(res);
}

//int ibz_mat_4x4_scalar_div(ibz_mat_4x4_t *quot, const ibz_t *scalar, const ibz_mat_4x4_t *mat);
int quat_test_dim4_ibz_mat_4x4_scalar_div(){
    int res = 0;
    int s;
    ibz_t scalar;
    ibz_mat_4x4_t quot,mat,cmp;
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&cmp);
    ibz_mat_4x4_init(&quot);
    ibz_init(&scalar);

    s = 5;
    ibz_set(&scalar,s);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),(i+j)*s);
            ibz_set(&(cmp[i][j]),(i+j));
        }
    }
    ibz_mat_4x4_scalar_div(&quot,&scalar,&mat);
    res = res || (!ibz_mat_4x4_equal(&quot,&cmp));

    ibz_mat_4x4_scalar_div(&mat,&scalar,&mat);
    res = res || (!ibz_mat_4x4_equal(&mat,&cmp));

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_scalar_div failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&cmp);
    ibz_mat_4x4_finalize(&quot);
    ibz_finalize(&scalar);
    return(res);
}


//int ibz_mat_4x4_is_hnf(const ibz_mat_4x4_t *mat);
int quat_test_dim4_is_hnf(){
    int res = 0;
    ibz_mat_4x4_t mat;
    ibz_mat_4x4_init(&mat);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    res = res || (!ibz_mat_4x4_is_hnf(&mat));
    ibz_set(&(mat[0][0]),7);
    ibz_set(&(mat[0][1]),6);
    ibz_set(&(mat[0][2]),5);
    ibz_set(&(mat[0][3]),4);
    ibz_set(&(mat[1][1]),6);
    ibz_set(&(mat[1][2]),5);
    ibz_set(&(mat[1][3]),4);
    ibz_set(&(mat[2][2]),5);
    ibz_set(&(mat[2][3]),4);
    ibz_set(&(mat[3][3]),4);
    res = res || (!ibz_mat_4x4_is_hnf(&mat));

    ibz_set(&(mat[0][0]),7);
    ibz_set(&(mat[0][1]),0);
    ibz_set(&(mat[0][2]),5);
    ibz_set(&(mat[0][3]),4);
    ibz_set(&(mat[1][1]),0);
    ibz_set(&(mat[1][2]),0);
    ibz_set(&(mat[1][3]),0);
    ibz_set(&(mat[2][2]),5);
    ibz_set(&(mat[2][3]),4);
    ibz_set(&(mat[3][3]),4);
    res = res || (!ibz_mat_4x4_is_hnf(&mat));

    // negative tests
    ibz_set(&(mat[0][0]),7);
    ibz_set(&(mat[0][1]),0);
    ibz_set(&(mat[0][2]),5);
    ibz_set(&(mat[0][3]),4);
    ibz_set(&(mat[1][1]),1);
    ibz_set(&(mat[1][2]),5);
    ibz_set(&(mat[1][3]),9);
    ibz_set(&(mat[2][2]),5);
    ibz_set(&(mat[2][3]),4);
    ibz_set(&(mat[3][3]),4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));

    ibz_set(&(mat[0][0]),7);
    ibz_set(&(mat[0][1]),0);
    ibz_set(&(mat[0][2]),5);
    ibz_set(&(mat[0][3]),4);
    ibz_set(&(mat[1][1]),1);
    ibz_set(&(mat[1][2]),-5);
    ibz_set(&(mat[1][3]),1);
    ibz_set(&(mat[2][2]),5);
    ibz_set(&(mat[2][3]),4);
    ibz_set(&(mat[3][3]),4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));


    ibz_set(&(mat[0][0]),7);
    ibz_set(&(mat[0][1]),0);
    ibz_set(&(mat[0][2]),5);
    ibz_set(&(mat[0][3]),4);
    ibz_set(&(mat[1][0]),2);
    ibz_set(&(mat[1][1]),3);
    ibz_set(&(mat[1][2]),1);
    ibz_set(&(mat[1][3]),1);
    ibz_set(&(mat[2][2]),5);
    ibz_set(&(mat[2][3]),4);
    ibz_set(&(mat[3][3]),4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));


    ibz_set(&(mat[0][0]),7);
    ibz_set(&(mat[0][1]),0);
    ibz_set(&(mat[0][2]),5);
    ibz_set(&(mat[0][3]),4);
    ibz_set(&(mat[1][0]),2);
    ibz_set(&(mat[1][1]),3);
    ibz_set(&(mat[1][2]),-1);
    ibz_set(&(mat[1][3]),7);
    ibz_set(&(mat[2][2]),0);
    ibz_set(&(mat[2][3]),0);
    ibz_set(&(mat[3][3]),4);
    res = res || (ibz_mat_4x4_is_hnf(&mat));

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_is_hnf failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    return(res);
}


//void ibz_mat_4x8_hnf_core(ibz_mat_4x4_t *hnf, const ibz_mat_4x8_t *generators);
int quat_test_dim4_ibz_mat_4x8_hnf_core(){
    int res = 0;
    ibz_mat_4x8_t mat;
    ibz_mat_4x4_t hnf,cmp;
    ibz_mat_4x8_init(&mat);
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&cmp);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 8; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    ibz_mat_4x8_hnf_core(&hnf,&mat);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));
    // also should test that they generate the same lattice. Since HNF is unique, copute the HNF for test vectors might be ok

    ibz_set(&(mat[0][2]),2);
    ibz_set(&(mat[1][3]),3);
    ibz_set(&(mat[0][4]),4);
    ibz_set(&(mat[3][2]),5);
    ibz_set(&(mat[3][7]),6);
    ibz_set(&(mat[1][7]),7);
    ibz_set(&(mat[1][3]),8);
    ibz_set(&(mat[1][1]),9);
    ibz_set(&(mat[0][6]),10);
    ibz_set(&(mat[0][5]),11);
    ibz_set(&(mat[0][0]),12);
    ibz_mat_4x8_hnf_core(&hnf,&mat);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));

    ibz_set(&(mat[2][5]),1);
    ibz_set(&(mat[2][0]),2);
    ibz_mat_4x8_hnf_core(&hnf,&mat);
    res = res || (!ibz_mat_4x4_is_hnf(&hnf));

    // test equality of result to a known hnf
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 8; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    ibz_set(&(mat[0][0]),4);
    ibz_set(&(mat[0][2]),3);
    ibz_set(&(mat[0][4]),1);
    ibz_set(&(mat[0][7]),-1);
    ibz_set(&(mat[1][1]),5);
    ibz_set(&(mat[1][5]),-2);
    ibz_set(&(mat[2][2]),3);
    ibz_set(&(mat[2][6]),1);
    ibz_set(&(mat[2][5]),1);
    ibz_set(&(mat[3][3]),7);
    ibz_set(&(mat[3][7]),-3);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            ibz_set(&(cmp[i][j]),0);
        }
        ibz_set(&(cmp[i][i]),1);
    }
    ibz_mat_4x8_hnf_core(&hnf,&mat);
    res = res || (!ibz_mat_4x4_equal(&cmp,&hnf));

    // test known hnf encountered in
    // https://github.com/SQISign/sqisign-nist/issues/38#issuecomment-1554585079
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 8; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    ibz_set(&(mat[0][4]),438);
    ibz_set(&(mat[1][4]),400);
    ibz_set(&(mat[2][4]),156);
    ibz_set(&(mat[3][4]),-2);
    ibz_set(&(mat[0][5]),-400);
    ibz_set(&(mat[1][5]),438);
    ibz_set(&(mat[2][5]),2);
    ibz_set(&(mat[3][5]),156);
    ibz_set(&(mat[0][6]),-28826);
    ibz_set(&(mat[1][6]),-148);
    ibz_set(&(mat[2][6]),220);
    ibz_set(&(mat[3][6]),-122);
    ibz_set(&(mat[0][7]),586);
    ibz_set(&(mat[1][7]),-28426);
    ibz_set(&(mat[2][7]),278);
    ibz_set(&(mat[3][7]),218);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            ibz_set(&(cmp[i][j]),0);
        }
    }
    ibz_set(&(cmp[0][0]),2321156);
    ibz_set(&(cmp[1][1]),2321156);
    ibz_set(&(cmp[0][2]),620252);
    ibz_set(&(cmp[1][2]),365058);
    ibz_set(&(cmp[2][2]),2);
    ibz_set(&(cmp[0][3]),1956098);
    ibz_set(&(cmp[1][3]),620252);
    ibz_set(&(cmp[3][3]),2);
    ibz_mat_4x8_hnf_core(&hnf,&mat);
    res = res || (!ibz_mat_4x4_equal(&cmp,&hnf));

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x8_hnf_core failed\n");
    }
    ibz_mat_4x8_finalize(&mat);
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&cmp);
    return(res);
}

//void ibz_mat_4x4_hnf_mod(ibz_mat_4x4_t *hnf, const ibz_mat_4x4_t *mat, const ibz_t *mod);
int quat_test_dim4_ibz_mat_4x4_hnf_mod(){
    int res = 0;
    ibz_mat_4x4_t hnf, mat, cmp;
    ibz_t mod;
    ibz_mat_4x4_init(&hnf);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&cmp);
    ibz_init(&mod);
    
    ibz_set(&mod,3);
    ibz_set(&(mat[0][0]), 4);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 3);
    ibz_set(&(mat[0][3]), 0);
    ibz_set(&(mat[1][0]), 0);
    ibz_set(&(mat[1][1]), 5);
    ibz_set(&(mat[1][2]), 0);
    ibz_set(&(mat[1][3]), 0);
    ibz_set(&(mat[2][0]), 3);
    ibz_set(&(mat[2][1]), -1);
    ibz_set(&(mat[2][2]), 3);
    ibz_set(&(mat[2][3]), -5);
    ibz_set(&(mat[3][0]), 7);
    ibz_set(&(mat[3][1]), 0);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 11);

    ibz_set(&(cmp[0][0]), 3);
    ibz_set(&(cmp[0][1]), 2);
    ibz_set(&(cmp[0][2]), 1);
    ibz_set(&(cmp[0][3]), 1);
    ibz_set(&(cmp[1][0]), 0);
    ibz_set(&(cmp[1][1]), 1);
    ibz_set(&(cmp[1][2]), 0);
    ibz_set(&(cmp[1][3]), 0);
    ibz_set(&(cmp[2][0]), 0);
    ibz_set(&(cmp[2][1]), 0);
    ibz_set(&(cmp[2][2]), 1);
    ibz_set(&(cmp[2][3]), 0);
    ibz_set(&(cmp[3][0]), 0);
    ibz_set(&(cmp[3][1]), 0);
    ibz_set(&(cmp[3][2]), 0);
    ibz_set(&(cmp[3][3]), 1);

    ibz_mat_4x4_hnf_mod(&hnf, &mat, &mod);
    res = res || !ibz_mat_4x4_equal(&hnf,&cmp);

    ibz_set(&mod,100);
    ibz_set(&(mat[0][0]), 12);
    ibz_set(&(mat[0][1]), -111);
    ibz_set(&(mat[0][2]), 12);
    ibz_set(&(mat[0][3]), -345);
    ibz_set(&(mat[1][0]), 1134);
    ibz_set(&(mat[1][1]), -2);
    ibz_set(&(mat[1][2]), 56);
    ibz_set(&(mat[1][3]), 72);
    ibz_set(&(mat[2][0]), 8);
    ibz_set(&(mat[2][1]), 231);
    ibz_set(&(mat[2][2]), -21);
    ibz_set(&(mat[2][3]), 22);
    ibz_set(&(mat[3][0]), 30);
    ibz_set(&(mat[3][1]), 34);
    ibz_set(&(mat[3][2]), 5);
    ibz_set(&(mat[3][3]), -33);

    ibz_set(&(cmp[0][0]), 4);
    ibz_set(&(cmp[0][1]), 2);
    ibz_set(&(cmp[0][2]), 3);
    ibz_set(&(cmp[0][3]), 3);
    ibz_set(&(cmp[1][0]), 0);
    ibz_set(&(cmp[1][1]), 2);
    ibz_set(&(cmp[1][2]), 0);
    ibz_set(&(cmp[1][3]), 0);
    ibz_set(&(cmp[2][0]), 0);
    ibz_set(&(cmp[2][1]), 0);
    ibz_set(&(cmp[2][2]), 1);
    ibz_set(&(cmp[2][3]), 0);
    ibz_set(&(cmp[3][0]), 0);
    ibz_set(&(cmp[3][1]), 0);
    ibz_set(&(cmp[3][2]), 0);
    ibz_set(&(cmp[3][3]), 1);

    ibz_mat_4x4_hnf_mod(&hnf, &mat, &mod);
    res = res || !ibz_mat_4x4_equal(&hnf,&cmp);

    ibz_set(&mod,3);
    ibz_set(&(cmp[0][0]), 3);
    ibz_set(&(cmp[0][1]), 0);
    ibz_set(&(cmp[0][2]), 0);
    ibz_set(&(cmp[0][3]), 0);
    ibz_set(&(cmp[1][0]), 0);
    ibz_set(&(cmp[1][1]), 3);
    ibz_set(&(cmp[1][2]), 0);
    ibz_set(&(cmp[1][3]), 1);
    ibz_set(&(cmp[2][0]), 0);
    ibz_set(&(cmp[2][1]), 0);
    ibz_set(&(cmp[2][2]), 1);
    ibz_set(&(cmp[2][3]), 0);
    ibz_set(&(cmp[3][0]), 0);
    ibz_set(&(cmp[3][1]), 0);
    ibz_set(&(cmp[3][2]), 0);
    ibz_set(&(cmp[3][3]), 1);

    ibz_mat_4x4_hnf_mod(&mat, &mat, &mod);
    res = res || !ibz_mat_4x4_equal(&mat,&cmp);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_hnf_mod failed\n");
    }
    ibz_mat_4x4_finalize(&hnf);
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&cmp);
    ibz_finalize(&mod);
    return(res);
}


// test for lll verification
//void ibq_vec_4_copy_ibz(ibq_t (*vec)[4], const ibz_t *coeff0, const ibz_t *coeff1,const ibz_t *coeff2,const ibz_t *coeff3);
int quat_test_dim4_ibq_vec_4_copy_ibz(){
    int res = 0;
    ibq_t vec[4];
    ibz_vec_4_t vec_z;
    ibz_vec_4_init(&vec_z);
    ibq_init(&(vec[0]));
    ibq_init(&(vec[1]));
    ibq_init(&(vec[2]));
    ibq_init(&(vec[3]));
    ibz_vec_4_set(&vec_z,2,3,4,5);
    ibq_vec_4_copy_ibz(&vec,&(vec_z[0]),&(vec_z[1]),&(vec_z[2]),&(vec_z[3]));
    for(int i = 0; i <4; i++){
        ibq_to_ibz(&(vec_z[i]),&(vec[i]));
        res = res || !(ibz_get(&(vec_z[i]))==i+2);
    }

    if (res != 0){
        printf("Quaternion unit test dim4_ibq_vec_4_copy_ibz failed\n");
    }
    ibz_vec_4_finalize(&vec_z);
    ibq_finalize(&(vec[0]));
    ibq_finalize(&(vec[1]));
    ibq_finalize(&(vec[2]));
    ibq_finalize(&(vec[3]));
    return(res);
}

//void quat_dim4_lll_bilinear(ibq_t *b, const ibq_t (*vec0)[4], const ibq_t (*vec1)[4], const ibz_t *q);
int quat_test_dim4_lll_bilinear(){
    int res = 0;
    ibz_vec_4_t init_helper;
    ibq_t vec0[4];
    ibq_t vec1[4];
    ibz_t q;
    ibq_t cmp, b;
    ibz_vec_4_init(&init_helper);
    ibq_init(&cmp);
    ibq_init(&b);
    ibz_init(&q);
    for(int i = 0; i <4; i++){
        ibq_init(&(vec0[i]));
        ibq_init(&(vec1[i]));
    }
    ibz_vec_4_set(&init_helper,1,2,3,4);
    ibq_vec_4_copy_ibz(&vec0,&(init_helper[0]),&(init_helper[1]),&(init_helper[2]),&(init_helper[3]));
    ibz_vec_4_set(&init_helper,9,-8,7,-6);
    ibq_vec_4_copy_ibz(&vec1,&(init_helper[0]),&(init_helper[1]),&(init_helper[2]),&(init_helper[3]));
    for(int i = 0; i <4; i++){
        ibq_inv(&(vec0[i]),&(vec0[i]));
    }
    ibz_set(&q,3);
    ibz_vec_4_set(&init_helper,15,2,0,0);
    ibq_set(&cmp,&(init_helper[0]),&(init_helper[1]));
    quat_dim4_lll_bilinear(&b,&vec0,&vec1,&q);
    res = res || (ibq_cmp(&b,&cmp));

    if (res != 0){
        printf("Quaternion unit test quat_dim4_lll_bilinear failed\n");
    }
    ibq_finalize(&cmp);
    ibq_finalize(&b);
    ibz_finalize(&q);
    ibz_vec_4_finalize(&init_helper);
    for(int i = 0; i <4; i++){
        ibq_finalize(&(vec0[i]));
        ibq_finalize(&(vec1[i]));
    }
    return(res);
}

//void quat_dim4_gram_schmidt_transposed_with_ibq(ibq_t (*orthogonalised_transposed)[4][4], const ibz_mat_4x4_t *mat, const ibz_t *q);
int quat_test_dim4_gram_schmidt_transposed_with_ibq(){
    int res = 0;
    int zero;
    ibq_t ot[4][4];
    ibq_t cmp[4][4];
    ibz_mat_4x4_t mat;
    ibz_t q, num, denom;
    ibq_t b;
    ibz_init(&q);
    ibz_init(&num);
    ibz_init(&denom);
    ibq_init(&b);
    ibz_mat_4x4_init(&mat);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibq_init(&(cmp[i][j]));
            ibq_init(&(ot[i][j]));
        }
    }


    ibz_mat_4x4_zero(&mat);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]), i*i+(j+5)*j-2+(i==j));
        }
    }
    ibz_set(&q,3);
    quat_dim4_gram_schmidt_transposed_with_ibq(&ot,&mat,&q);
    // test orthogonality
    for(int i = 0; i < 4; i++){
        for(int j = i+1; j < 4; j++){
            quat_dim4_lll_bilinear(&b,&(ot[i]),&(ot[j]),&q);
            res = res || !ibq_is_zero(&b);
        }
    }
    // test first vector is identical to mat
    for(int i = 0; i < 4; i++){
        ibq_to_ibz(&q,&(ot[0][i]));
        res = res || ibz_cmp(&q,&(mat[i][0]));
    }
    // test no zero vector
    for(int i = 0; i < 4; i++){
        zero = 1;
        for(int j = 0; j < 4; j++){
            zero = zero && ibq_is_zero(&(ot[i][j]));
        }
        res = res || zero;
    }

    ibz_set(&(mat[0][0]), 1);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 1);
    ibz_set(&(mat[0][3]), 0);
    ibz_set(&(mat[1][0]), 0);
    ibz_set(&(mat[1][1]), 1);
    ibz_set(&(mat[1][2]), 0);
    ibz_set(&(mat[1][3]), 1);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), 0);
    ibz_set(&(mat[2][2]), 2);
    ibz_set(&(mat[2][3]), 0);
    ibz_set(&(mat[3][0]), 0);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 2);
    ibz_set(&denom,1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibq_set(&(cmp[i][j]),&(mat[j][i]),&denom);
        }
    }
    ibz_set(&denom,3);
    ibz_set(&num,-2);
    ibq_set(&(cmp[2][0]),&num,&denom);
    ibq_set(&(cmp[3][1]),&num,&denom);
    ibz_set(&num,1);
    ibq_set(&(cmp[2][2]),&num,&denom);
    ibq_set(&(cmp[3][3]),&num,&denom);
    ibz_set(&q,2);
    quat_dim4_gram_schmidt_transposed_with_ibq(&ot,&mat,&q);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res = res || ibq_cmp(&(cmp[i][j]),&(ot[i][j]));
        }
    }

    ibz_set(&(mat[0][0]), 1);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 1);
    ibz_set(&(mat[0][3]), 0);
    ibz_set(&(mat[1][0]), 0);
    ibz_set(&(mat[1][1]), 1);
    ibz_set(&(mat[1][2]), 0);
    ibz_set(&(mat[1][3]), 1);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), 0);
    ibz_set(&(mat[2][2]), 2);
    ibz_set(&(mat[2][3]), 1);
    ibz_set(&(mat[3][0]), 0);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 2);
    ibz_set(&denom,1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibq_set(&(cmp[i][j]),&(mat[j][i]),&denom);
        }
    }
    ibz_set(&denom,3);
    ibz_set(&num,-2);
    ibq_set(&(cmp[2][0]),&num,&denom);
    ibq_set(&(cmp[3][1]),&num,&denom);
    ibz_set(&num,1);
    ibq_set(&(cmp[2][2]),&num,&denom);
    ibq_set(&(cmp[3][3]),&num,&denom);
    ibz_set(&num,0);
    ibq_set(&(cmp[3][0]),&num,&denom);
    ibq_set(&(cmp[3][2]),&num,&denom);
    ibz_set(&q,2);
    quat_dim4_gram_schmidt_transposed_with_ibq(&ot,&mat,&q);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res = res || ibq_cmp(&(cmp[i][j]),&(ot[i][j]));
        }
    }

    if (res != 0){
        printf("Quaternion unit test dim4_gram_schmidt_transposed_with_ibq failed\n");
    }
    ibz_finalize(&q);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibq_finalize(&b);
    ibz_mat_4x4_finalize(&mat);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibq_finalize(&(ot[i][j]));
            ibq_finalize(&(cmp[i][j]));
        }
    }
    return(res);
}



//int quat_dim4_lll_verify(const ibz_mat_4x4_t *mat, const ibq_t *coeff, const ibz_t *q);
int quat_test_dim4_lll_verify(){
    int res = 0;
    ibz_mat_4x4_t mat;
    ibz_t q, coeff_num, coeff_denom;
    ibq_t coeff;
    ibz_mat_4x4_init(&mat);
    ibz_init(&q);
    ibq_init(&coeff);
    ibz_init(&coeff_num);
    ibz_init(&coeff_denom);

    ibz_mat_4x4_identity(&mat);
    ibz_set(&q,1);
    ibz_set(&coeff_num,1);
    ibz_set(&coeff_denom,1);
    ibq_set(&coeff,&coeff_num,&coeff_denom);
    res = res || !quat_dim4_lll_verify(&mat,&coeff,&q);

    // non reduced
    ibz_set(&q,1);
    ibz_set(&coeff_num,99);
    ibz_set(&coeff_denom,100);
    ibq_set(&coeff,&coeff_num,&coeff_denom);
    ibz_mat_4x4_identity(&mat);
    ibz_set(&(mat[0][0]),50);
    res = res || quat_dim4_lll_verify(&mat,&coeff,&q);

    // non reduced
    ibz_set(&q,1);
    ibz_set(&coeff_num,99);
    ibz_set(&coeff_denom,100);
    ibq_set(&coeff,&coeff_num,&coeff_denom);
    ibz_mat_4x4_identity(&mat);
    ibz_set(&(mat[0][1]),4);
    res = res || quat_dim4_lll_verify(&mat,&coeff,&q);

    // reduced
    ibz_set(&q,1);
    ibz_set(&coeff_num,99);
    ibz_set(&coeff_denom,100);
    ibq_set(&coeff,&coeff_num,&coeff_denom);
    ibz_set(&(mat[0][0]), 0);
    ibz_set(&(mat[0][1]), 2);
    ibz_set(&(mat[0][2]), 1);
    ibz_set(&(mat[0][3]), -7);
    ibz_set(&(mat[1][0]), 2);
    ibz_set(&(mat[1][1]), -1);
    ibz_set(&(mat[1][2]), -1);
    ibz_set(&(mat[1][3]), -6);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), -2);
    ibz_set(&(mat[2][2]), 4);
    ibz_set(&(mat[2][3]), 1);
    ibz_set(&(mat[3][0]), 1);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 13);
    res = res || !quat_dim4_lll_verify(&mat,&coeff,&q);

    // reduced: non-1 norm
    ibz_set(&q,3);
    ibz_set(&coeff_num,99);
    ibz_set(&coeff_denom,100);
    ibq_set(&coeff,&coeff_num,&coeff_denom);
    ibz_set(&(mat[0][0]), 0);
    ibz_set(&(mat[0][1]), 2);
    ibz_set(&(mat[0][2]), 3);
    ibz_set(&(mat[0][3]), -14);
    ibz_set(&(mat[1][0]), 2);
    ibz_set(&(mat[1][1]), -1);
    ibz_set(&(mat[1][2]), -4);
    ibz_set(&(mat[1][3]), -8);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), -2);
    ibz_set(&(mat[2][2]), 1);
    ibz_set(&(mat[2][3]), 0);
    ibz_set(&(mat[3][0]), 1);
    ibz_set(&(mat[3][1]), 1);
    ibz_set(&(mat[3][2]), 0);
    ibz_set(&(mat[3][3]), 7);
    res = res || !quat_dim4_lll_verify(&mat,&coeff,&q);

    // reduced: non-1 norm
    ibz_set(&q,103);
    ibz_set(&coeff_num,99);
    ibz_set(&coeff_denom,100);
    ibq_set(&coeff,&coeff_num,&coeff_denom);
    ibz_set(&(mat[0][0]), 3);
    ibz_set(&(mat[0][1]), 0);
    ibz_set(&(mat[0][2]), 90);
    ibz_set(&(mat[0][3]), -86);
    ibz_set(&(mat[1][0]), 11);
    ibz_set(&(mat[1][1]), 15);
    ibz_set(&(mat[1][2]), 12);
    ibz_set(&(mat[1][3]), 50);
    ibz_set(&(mat[2][0]), 1);
    ibz_set(&(mat[2][1]), -2);
    ibz_set(&(mat[2][2]), 0);
    ibz_set(&(mat[2][3]), 3);
    ibz_set(&(mat[3][0]), -1);
    ibz_set(&(mat[3][1]), 0);
    ibz_set(&(mat[3][2]), 5);
    ibz_set(&(mat[3][3]), 5);
    res = res || !quat_dim4_lll_verify(&mat,&coeff,&q);

    if (res != 0){
        printf("Quaternion unit test quat_dim4_lll_verify failed\n");
    }
    ibz_finalize(&q);
    ibq_finalize(&coeff);
    ibz_finalize(&coeff_num);
    ibz_finalize(&coeff_denom);
    ibz_mat_4x4_finalize(&mat);
    return(res);
}

//void ibz_inv_dim4_make_coeff_pmp(ibz_t *coeff, const ibz_t *a1, const ibz_t *a2, const ibz_t *b1, const ibz_t *b2, const ibz_t *c1, const ibz_t *c2);
int quat_test_dim4_ibz_inv_dim4_make_coeff_pmp(){
    int res = 0;
    ibz_t coeff,cmp,a1,a2,b1,b2, c1, c2;
    ibz_init(&a1);
    ibz_init(&a2);
    ibz_init(&b1);
    ibz_init(&b2);
    ibz_init(&c1);
    ibz_init(&c2);
    ibz_init(&coeff);
    ibz_init(&cmp);

    ibz_set(&a1,0);
    ibz_set(&a2,3);
    ibz_set(&b1,-1);
    ibz_set(&b2,0);
    ibz_set(&c1,-1);
    ibz_set(&c2,0);
    ibz_set(&cmp,0);
    ibz_inv_dim4_make_coeff_pmp(&coeff,&a1,&a2,&b1,&b2, &c1,&c2);
    res = res || ibz_cmp(&cmp,&coeff);

    ibz_set(&a1,2);
    ibz_set(&a2,3);
    ibz_set(&b1,-1);
    ibz_set(&b2,1);
    ibz_set(&c1,-4);
    ibz_set(&c2,2);
    ibz_set(&cmp,-1);
    ibz_inv_dim4_make_coeff_pmp(&coeff,&a1,&a2,&b1,&b2, &c1,&c2);
    res = res || ibz_cmp(&cmp,&coeff);

    ibz_set(&a1,2);
    ibz_set(&cmp,4);
    ibz_inv_dim4_make_coeff_pmp(&a1,&a1,&a1,&a1,&a1, &a1,&a1);
    res = res || ibz_cmp(&cmp,&a1);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_inv_dim4_make_coeff_pmp failed\n");
    }
    ibz_finalize(&a1);
    ibz_finalize(&a2);
    ibz_finalize(&b1);
    ibz_finalize(&b2);
    ibz_finalize(&c1);
    ibz_finalize(&c2);
    ibz_finalize(&cmp);
    ibz_finalize(&coeff);
    return res;
}

//void ibz_inv_make_coeff_mpm(ibz_t *coeff, const ibz_t *a1, const ibz_t *a2, const ibz_t *b1, const ibz_t *b2, const ibz_t *c1, const ibz_t *c2);
int quat_test_dim4_ibz_inv_dim4_make_coeff_mpm(){
    int res = 0;
    ibz_t coeff,cmp,a1,a2,b1,b2, c1, c2;
    ibz_init(&a1);
    ibz_init(&a2);
    ibz_init(&b1);
    ibz_init(&b2);
    ibz_init(&c1);
    ibz_init(&c2);
    ibz_init(&coeff);
    ibz_init(&cmp);

    ibz_set(&a1,0);
    ibz_set(&a2,3);
    ibz_set(&b1,-1);
    ibz_set(&b2,0);
    ibz_set(&c1,-1);
    ibz_set(&c2,0);
    ibz_set(&cmp,0);
   ibz_inv_dim4_make_coeff_mpm(&coeff,&a1,&a2,&b1,&b2, &c1,&c2);
    res = res || ibz_cmp(&cmp,&coeff);

    ibz_set(&a1,2);
    ibz_set(&a2,3);
    ibz_set(&b1,-1);
    ibz_set(&b2,1);
    ibz_set(&c1,-4);
    ibz_set(&c2,2);
    ibz_set(&cmp,1);
   ibz_inv_dim4_make_coeff_mpm(&coeff,&a1,&a2,&b1,&b2, &c1,&c2);
    res = res || ibz_cmp(&cmp,&coeff);


    ibz_set(&a1,2);
    ibz_set(&cmp,-4);
   ibz_inv_dim4_make_coeff_mpm(&a1,&a1,&a1,&a1,&a1, &a1,&a1);
    res = res || ibz_cmp(&cmp,&a1);


    if (res != 0){
        printf("Quaternion unit test dim4_ibz_inv_dim4_make_coeff_mpm failed\n");
    }
    ibz_finalize(&a1);
    ibz_finalize(&a2);
    ibz_finalize(&b1);
    ibz_finalize(&b2);
    ibz_finalize(&c1);
    ibz_finalize(&c2);
    ibz_finalize(&cmp);
    ibz_finalize(&coeff);
    return res;
}

// returns 1 if inverse is valid, 0 otherwise
int quat_test_dim4_validate_mat_4x4_rational_inv_if_exists(const ibz_mat_4x4_t *mat, const ibz_t *det, const ibz_mat_4x4_t *inv){
    int res = 1;
    ibz_mat_4x4_t det_id, prod;
    ibz_mat_4x4_init(&det_id);
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_identity(&det_id);
    ibz_mat_4x4_scalar_mul(&det_id,det,&det_id);
    ibz_mat_4x4_mul(&prod, inv, mat);
    res = res && ibz_mat_4x4_equal(&det_id,&prod);
    ibz_mat_4x4_mul(&prod, mat, inv);
    res = res && ibz_mat_4x4_equal(&det_id,&prod);
    ibz_mat_4x4_finalize(&det_id);
    ibz_mat_4x4_finalize(&prod);
    return res;
}

//int ibz_4x4_inv_with_det_as_denom(ibz_mat_4x4_t *inv, ibz_t *det, const ibz_mat_4x4_t mat);
int quat_test_dim4_ibz_mat_4x4_inv_with_det_as_denom(){
    int res = 0;
    ibz_t det;
    ibz_mat_4x4_t mat, inv;
    ibz_init(&det);
    ibz_mat_4x4_init(&mat);
    ibz_mat_4x4_init(&inv);

    ibz_mat_4x4_zero(&mat);
    res = res || ibz_mat_4x4_inv_with_det_as_denom(&inv,&det,&mat);
    res = res || !ibz_is_zero(&det);
    ibz_mat_4x4_identity(&mat);
    if(ibz_mat_4x4_inv_with_det_as_denom(&inv,&det,&mat)){
        res = res || !ibz_is_one(&det);
        res = res || !quat_test_dim4_validate_mat_4x4_rational_inv_if_exists(&inv,&det,&mat);
    } else {
        res = 1;
    }
    ibz_set(&(mat[0][0]),2);
    ibz_set(&(mat[0][1]),-17);
    ibz_set(&(mat[0][2]),3);
    ibz_set(&(mat[0][3]),5);
    ibz_set(&(mat[1][1]),-2);
    ibz_set(&(mat[1][2]),3);
    ibz_set(&(mat[1][3]),2);
    ibz_set(&(mat[2][2]),-3);
    ibz_set(&(mat[2][3]),0);
    ibz_set(&(mat[3][3]),1);
    if(ibz_mat_4x4_inv_with_det_as_denom(&inv,&det,&mat)){
        res = res || (ibz_get(&det)!=12);
        res = res || !quat_test_dim4_validate_mat_4x4_rational_inv_if_exists(&inv,&det,&mat);
    } else {
        res = 1;
    }
    ibz_set(&(mat[3][0]),1);
    ibz_set(&(mat[3][1]),8);
    ibz_set(&(mat[3][2]),-9);
    ibz_set(&(mat[2][0]),3);
    ibz_set(&(mat[2][1]),0);
    ibz_set(&(mat[1][0]),4);
    if(ibz_mat_4x4_inv_with_det_as_denom(&inv,&det,&mat)){
        res = res || (ibz_get(&det)!=-1503);
        res = res || !quat_test_dim4_validate_mat_4x4_rational_inv_if_exists(&inv,&det,&mat);
    } else {
        res = 1;
    }
    
    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_inv_with_det_as_denom failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_mat_4x4_finalize(&inv);
    ibz_finalize(&det);
    return res;
}

// larger matrix modular kernel

// returns 1 if product is 0 mod p, 0 otherwise
int quat_test_dim4_4x5_vec_5_mul_is_zero_mod(const ibz_mat_4x5_t *mat, const ibz_vec_5_t *vec, ibz_t *p){
    int res = 1;
    ibz_t  prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    for(int i = 0; i < 4; i++){
        ibz_set(&sum,0);
        for(int j = 0; j < 5; j++){
            ibz_mul(&prod,&((*mat)[i][j]),&((*vec)[j]));
            ibz_add(&sum,&sum,&prod);
            ibz_mod(&sum,&sum,p);
        }
        res = res && ibz_is_zero(&sum);
    }
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    return(res);
}

//int ibz_4x5_right_ker_mod_prime(ibz_vec_5_t *ker, const ibz_mat_4x5_t *mat, const ibz_t *p);
int quat_test_dim4_ibz_4x5_right_ker_mod_prime(){
    int res = 0;
    ibz_t prime;
    ibz_mat_4x5_t mat;
    ibz_vec_5_t ker;
    ibz_mat_4x5_init(&mat);
    ibz_vec_5_init(&ker);
    ibz_init(&prime);

    // kernel has dimension 1
    ibz_set(&prime,5);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 5; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    ibz_set(&(mat[0][0]),2);
    ibz_set(&(mat[1][1]),3);
    ibz_set(&(mat[2][2]),3);
    ibz_set(&(mat[1][3]),1);
    ibz_set(&(mat[3][3]),2);
    ibz_set(&(mat[3][2]),1);
    ibz_set(&(mat[0][4]),1);
    if (ibz_4x5_right_ker_mod_prime(&ker,&mat,&prime)){
        res = res || !quat_test_dim4_4x5_vec_5_mul_is_zero_mod(&mat,&ker,&prime);
    } else {
        res = 1;
    }

    // too large kernel
    ibz_set(&prime,5);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 5; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    res = res || ibz_4x5_right_ker_mod_prime(&ker,&mat,&prime);
    ibz_set(&(mat[0][0]),1);
    ibz_set(&(mat[2][2]),2);
    ibz_set(&(mat[3][3]),3);
    res = res || ibz_4x5_right_ker_mod_prime(&ker,&mat,&prime);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_4x5_right_ker_mod_prime failed\n");
    }
    ibz_vec_5_finalize(&ker);
    ibz_mat_4x5_finalize(&mat);
    ibz_finalize(&prime);
    return(res);
}

// returns 1 if product is 0 mod p, 0 otherwise
int quat_test_dim4_4x4_vec_4_mul_is_zero_mod(const ibz_mat_4x4_t *mat, const ibz_vec_4_t *vec, ibz_t *p){
    int res = 1;
    ibz_t  prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    for(int i = 0; i < 4; i++){
        ibz_set(&sum,0);
        for(int j = 0; j < 4; j++){
            ibz_mul(&prod,&((*mat)[i][j]),&((*vec)[j]));
            ibz_add(&sum,&sum,&prod);
            ibz_mod(&sum,&sum,p);
        }
        res = res && ibz_is_zero(&sum);
    }
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    return(res);
}

//int ibz_4x4_right_ker_mod_prime(ibz_vec_4_t *ker, const ibz_mat_4x4_t *mat, const ibz_t *p);
int quat_test_dim4_ibz_4x4_right_ker_mod_prime(){
    int res = 0;
    ibz_t prime;
    ibz_mat_4x4_t mat;
    ibz_vec_4_t ker;
    ibz_mat_4x4_init(&mat);
    ibz_vec_4_init(&ker);
    ibz_init(&prime);

    // kernel has dimension 1
    ibz_set(&prime,5);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    ibz_set(&(mat[0][0]),2);
    ibz_set(&(mat[1][1]),3);
    ibz_set(&(mat[2][2]),3);
    ibz_set(&(mat[2][3]),1);
    ibz_set(&(mat[3][3]),2);
    ibz_set(&(mat[3][2]),1);
    if (ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime)){
        res = res || !quat_test_dim4_4x4_vec_4_mul_is_zero_mod(&mat,&ker,&prime);
    } else {
        res = 1;
    }

    // too large kernel
    ibz_set(&prime,5);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(mat[i][j]),0);
        }
    }
    res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);
    ibz_set(&(mat[1][2]),2);
    ibz_set(&(mat[3][3]),3);
    res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);

    // too small kernel
    ibz_set(&(mat[0][0]),1);
    ibz_set(&(mat[1][1]),2);
    ibz_set(&(mat[2][2]),3);
    ibz_set(&(mat[3][3]),2);
    res = res || ibz_4x4_right_ker_mod_prime(&ker,&mat,&prime);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_4x4_right_ker_mod_prime failed\n");
    }
    ibz_vec_4_finalize(&ker);
    ibz_mat_4x4_finalize(&mat);
    ibz_finalize(&prime);
    return(res);
}

//int ibz_4x4_right_ker_mod_power_of_2(ibz_vec_4_t *ker, const ibz_mat_4x4_t *mat, unsigned short exp);
int quat_test_dim4_ibz_4x4_right_ker_mod_power_of_2(){
    int res = 0;
    int zero = 1;

    ibz_mat_4x4_t mat;
    ibz_vec_4_t ker;
    ibz_vec_4_t prod;
    ibz_t q, r, two;
    ibz_mat_4x4_init(&mat);
    ibz_vec_4_init(&ker);
    ibz_vec_4_init(&prod);
    ibz_init(&q);
    ibz_init(&r);
    ibz_init(&two);
    ibz_set(&two,2);

    // A diagonal matrix with an obvious kernel over the 2-adics
    ibz_mat_4x4_zero(&mat);
    ibz_set(&mat[0][0], 1);
    ibz_set(&mat[1][1], 1);
    ibz_set(&mat[2][1], 1);
    ibz_set(&mat[3][3], 2);

    res |= ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, 60) != 1;

    res |= ibz_cmp_si(&ker[0], 0);
    res |= ibz_cmp_si(&ker[1], 0);
    res |= ibz_cmp_si(&ker[2], 1);
    res |= ibz_cmp_si(&ker[3], 0);

    // such a vector does exist
    ibz_mat_4x4_zero(&mat);
    ibz_set(&mat[0][0], 4);
    ibz_set(&mat[0][1], 21);
    ibz_set(&mat[0][2], 0);
    ibz_set(&mat[0][3], 2);
    ibz_set(&mat[1][0], 82);
    ibz_set(&mat[1][1], 21);
    ibz_set(&mat[1][2], 0);
    ibz_set(&mat[1][3], 22);
    ibz_set(&mat[2][0], 22);
    ibz_set(&mat[2][1], 128);
    ibz_set(&mat[2][2], 22);
    ibz_set(&mat[2][3], 108);
    ibz_set(&mat[3][0], 32);
    ibz_set(&mat[3][1], 10);
    ibz_set(&mat[3][2], 4);
    ibz_set(&mat[3][3], 4);

    if (ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, 6)){
        zero = 1;
        ibz_mat_4x4_eval(&prod,&mat,&ker);
        for (int i = 0; i < 4; i++){
            ibz_div(&q,&r,&(ker[i]),&two);
            zero = zero && ibz_is_zero(&r);
            ibz_pow(&q,&two,6);
            ibz_mod(&(prod[i]),&(prod[i]),&q);
        }
        res |= !quat_alg_coord_is_zero(&prod);
        res |= zero;
    } else {
        res = 1;
    }

    // such a vector does exist (since 0 row)
    ibz_mat_4x4_zero(&mat);
    ibz_set(&mat[0][0], 1);
    ibz_set(&mat[0][1], 1);
    ibz_set(&mat[0][2], 2);
    ibz_set(&mat[1][0], 3);
    ibz_set(&mat[1][1], 9);
    ibz_set(&mat[1][2], 9);
    ibz_set(&mat[2][0], 7);
    ibz_set(&mat[2][1], 1);
    ibz_set(&mat[2][2], 8);
    ibz_set(&mat[2][3], 5);

    if (ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, 6)){
        zero = 1;
        ibz_mat_4x4_eval(&prod,&mat,&ker);
        for (int i = 0; i < 4; i++){
            ibz_div(&q,&r,&(ker[i]),&two);
            zero = zero && ibz_is_zero(&r);
            ibz_pow(&q,&two,6);
            ibz_mod(&(prod[i]),&(prod[i]),&q);
        }
        res |= !quat_alg_coord_is_zero(&prod);
        res |= zero;
    } else {
        res = 1;
    }


    // Vector exists
    ibz_mat_4x4_zero(&mat);
    ibz_set(&mat[0][0], 1);
    ibz_set(&mat[0][1], 1);
    ibz_set(&mat[0][2], 2);
    ibz_set(&mat[1][0], 3);
    ibz_set(&mat[1][1], 9);
    ibz_set(&mat[1][2], 0);
    ibz_set(&mat[2][0], 7);
    ibz_set(&mat[2][1], 1);
    ibz_set(&mat[2][2], 8);
    ibz_set(&mat[2][3], 5);

    res |= !ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, 6);
    // Can check exact match, thanks to Howell form
    res |= ibz_cmp_si(&ker[0], 43) != 0;
    res |= ibz_cmp_si(&ker[1], 7) != 0;
    res |= ibz_cmp_si(&ker[2], 7) != 0;
    res |= ibz_cmp_si(&ker[3], 4) != 0;

    // Vector exists
    ibz_mat_4x4_zero(&mat);
    ibz_set(&mat[0][0], 1);
    ibz_set(&mat[0][1], 1);
    ibz_set(&mat[0][2], 2);
    ibz_set(&mat[1][0], 3);
    ibz_set(&mat[1][1], 9);
    ibz_set(&mat[1][2], 0);
    ibz_set(&mat[1][3], 1);
    ibz_set(&mat[2][0], 7);
    ibz_set(&mat[2][1], 1);
    ibz_set(&mat[2][2], 1);
    ibz_set(&mat[2][3], 5);

    res |= !ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, 6);
    // Can check exact match, thanks to Howell form
    res |= ibz_cmp_si(&ker[0], 31) != 0;
    res |= ibz_cmp_si(&ker[1], 25) != 0;
    res |= ibz_cmp_si(&ker[2], 4) != 0;
    res |= ibz_cmp_si(&ker[3], 2) != 0;

    // such a vector does not exist, since bad kernel
    ibz_mat_4x4_zero(&mat);
    ibz_set(&mat[0][0], 1);
    ibz_set(&mat[0][1], 1);
    ibz_set(&mat[0][2], 2);
    ibz_set(&mat[1][0], 3);
    ibz_set(&mat[1][1], 9);
    ibz_set(&mat[1][2], 0);
    ibz_set(&mat[2][0], 7);
    ibz_set(&mat[2][1], 1);
    ibz_set(&mat[2][2], 8);
    ibz_set(&mat[2][3], 5);
    ibz_set(&mat[3][0], 1);

    res |= ibz_4x4_right_ker_mod_power_of_2(&ker, &mat, 6);


    
    if (res != 0){
        printf("Quaternion unit dim4_ibz_4x4_right_ker_mod_power_of_2 failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_vec_4_finalize(&ker);
    ibz_vec_4_finalize(&prod);
    ibz_finalize(&q);
    ibz_finalize(&r);
    ibz_finalize(&two);
    return(res);
}

//void ibz_mat_4x4_eval(quat_alg_coord_t  *res, const ibz_mat_4x4_t *mat, const quat_alg_coord_t *vec);
int quat_test_dim4_ibz_mat_4x4_eval(){
    int res = 0;
    ibz_mat_4x4_t mat;
    quat_alg_coord_t vec, cmp, vres;
    quat_alg_coord_init(&cmp);
    quat_alg_coord_init(&vres);
    quat_alg_coord_init(&vec);
    ibz_mat_4x4_init(&mat);

    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            ibz_set(&(mat[i][j]), i*j);
        }
        ibz_set(&(vec[i]), i);
    }
    ibz_set(&(cmp[0]), 0);
    ibz_set(&(cmp[1]), 14);
    ibz_set(&(cmp[2]), 28);
    ibz_set(&(cmp[3]), 42);
    ibz_mat_4x4_eval(&vres,&mat,&vec);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(vres[i]),&(cmp[i]));
    }

    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            ibz_set(&(mat[i][j]), i*(j-1)+1);
        }
        ibz_set(&(vec[i]), i*i-2);
    }
    ibz_set(&(cmp[0]), 6);
    ibz_set(&(cmp[1]), 24);
    ibz_set(&(cmp[2]), 42);
    ibz_set(&(cmp[3]), 60);
    ibz_mat_4x4_eval(&vres,&mat,&vec);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(vres[i]),&(cmp[i]));
    }

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_mat_4x4_eval failed\n");
    }
    quat_alg_coord_finalize(&cmp);
    quat_alg_coord_finalize(&vres);
    quat_alg_coord_finalize(&vec);
    ibz_mat_4x4_finalize(&mat);
    return res;
}

//void quat_qf_eval(ibz_t *res, const ibz_mat_4x4_t *qf, const quat_alg_coord_t *coord);
int quat_test_dim4_qf_eval(){
    int res = 0;
    ibz_t ires, cmp;
    ibz_mat_4x4_t qf;
    quat_alg_coord_t vec;
    ibz_init(&cmp);
    ibz_init(&ires);
    quat_alg_coord_init(&vec);
    ibz_mat_4x4_init(&qf);

    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            ibz_set(&(qf[i][j]), i*j);
        }
        ibz_set(&(vec[i]), i);
    }
    ibz_set(&(cmp), 196);
    quat_qf_eval(&ires,&qf,&vec);
    res = res || ibz_cmp(&ires,&cmp);

    for (int i = 0; i <4; i++){
        for (int j = 0; j <4; j++){
            ibz_set(&(qf[i][j]), (i+1)*(j+1)-4);
        }
        ibz_set(&(vec[i]), (i-1)*2-2);
    }
    ibz_set(&(cmp), -4*16);
    quat_qf_eval(&ires,&qf,&vec);
    res = res || ibz_cmp(&ires,&cmp);

    if (res != 0){
        printf("Quaternion unit test dim4_qf_eval failed\n");
    }
    ibz_finalize(&cmp);
    ibz_finalize(&ires);
    quat_alg_coord_finalize(&vec);
    ibz_mat_4x4_finalize(&qf);
    return res;
}

//void ibz_content(ibz_t *content, const quat_alg_coord_t *v);
int quat_test_dim4_ibz_content(){
    int res = 0;
    ibz_t c,cmp;
    quat_alg_coord_t x;
    ibz_init(&c);
    ibz_init(&cmp);
    quat_alg_coord_init(&x);

    ibz_set(&(x[0]),0);
    ibz_set(&(x[1]),0);
    ibz_set(&(x[2]),0);
    ibz_set(&(x[3]),0);
    ibz_set(&cmp,0);
    ibz_content(&c,&x);
    res = res || ibz_cmp(&c,&cmp);

    ibz_set(&(x[0]),5);
    ibz_set(&(x[1]),25);
    ibz_set(&(x[2]),125);
    ibz_set(&(x[3]),30);
    ibz_set(&cmp,5);
    ibz_content(&c,&x);
    res = res || ibz_cmp(&c,&cmp);

    ibz_set(&(x[0]),5);
    ibz_set(&(x[1]),2);
    ibz_set(&(x[2]),125);
    ibz_set(&(x[3]),30);
    ibz_set(&cmp,1);
    ibz_content(&c,&x);
    res = res || ibz_cmp(&c,&cmp);

    ibz_set(&(x[0]),5);
    ibz_set(&(x[1]),-2);
    ibz_set(&(x[2]),125);
    ibz_set(&(x[3]),0);
    ibz_set(&cmp,1);
    ibz_content(&c,&x);
    res = res || ibz_cmp(&c,&cmp);

    ibz_set(&(x[0]),0);
    ibz_set(&(x[1]),-2);
    ibz_set(&(x[2]),0);
    ibz_set(&(x[3]),0);
    ibz_set(&cmp,2);
    ibz_content(&c,&x);
    res = res || ibz_cmp(&c,&cmp);

    if (res != 0){
        printf("Quaternion unit test dim4_ibz_content failed\n");
    }
    quat_alg_coord_finalize(&x);
    ibz_finalize(&c);
    ibz_finalize(&cmp);
    return res;
}


// run all previous tests
int quat_test_dim4(){
    int res = 0;
    printf("\nRunning quaternion tests of matrices, vectors and quadratic forms in dimension 4\n");
    res = res | quat_test_dim4_ibz_mat_4x4_mul();
    res = res | quat_test_dim4_ibz_vec_4_set();
    res = res | quat_test_dim4_ibz_vec_4_copy();
    res = res | quat_test_dim4_ibz_vec_4_negate();
    res = res | quat_test_dim4_ibz_vec_4_linear_combination();
    res = res | quat_test_dim4_ibz_mat_4x4_copy();
    res = res | quat_test_dim4_ibz_mat_4x4_negate();
    res = res | quat_test_dim4_ibz_mat_4x4_transpose();
    res = res | quat_test_dim4_ibz_mat_4x4_zero();
    res = res | quat_test_dim4_ibz_mat_4x4_identity();
    res = res | quat_test_dim4_ibz_mat_4x4_is_identity();
    res = res | quat_test_dim4_ibz_mat_4x4_equal();
    res = res | quat_test_dim4_ibz_mat_4x4_scalar_mul();
    res = res | quat_test_dim4_ibz_mat_4x4_gcd();
    res = res | quat_test_dim4_ibz_mat_4x4_scalar_div();
    res = res | quat_test_dim4_is_hnf();
    res = res | quat_test_dim4_ibz_mat_4x8_hnf_core();
    res = res | quat_test_dim4_ibz_mat_4x4_hnf_mod();
    res = res | quat_test_dim4_ibq_vec_4_copy_ibz();
    res = res | quat_test_dim4_lll_bilinear();
    res = res | quat_test_dim4_gram_schmidt_transposed_with_ibq();
    res = res | quat_test_dim4_lll_verify();
    res = res | quat_test_dim4_ibz_inv_dim4_make_coeff_pmp();
    res = res | quat_test_dim4_ibz_inv_dim4_make_coeff_mpm();
    res = res | quat_test_dim4_ibz_mat_4x4_inv_with_det_as_denom();
    res = res | quat_test_dim4_ibz_4x5_right_ker_mod_prime();
    res = res | quat_test_dim4_ibz_4x4_right_ker_mod_prime();
    res = res | quat_test_dim4_ibz_4x4_right_ker_mod_power_of_2();
    res = res | quat_test_dim4_ibz_mat_4x4_eval();
    res = res | quat_test_dim4_qf_eval();
    res = res | quat_test_dim4_ibz_content();
    return(res);
}
