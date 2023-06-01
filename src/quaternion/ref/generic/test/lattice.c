#include "quaternion_tests.h"

// helper functions

//int quat_lattice_equal(const quat_lattice_t *lat1, const quat_lattice_t *lat2);
int quat_test_lattice_equal(){
    int res = 0;
    quat_lattice_t lat, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&cmp);

    ibz_mat_4x4_identity(&(lat.basis));
    ibz_mat_4x4_identity(&(cmp.basis));
    res = res || !quat_lattice_equal(&lat,&cmp);
    ibz_set(&(lat.denom),5);
    ibz_set(&(cmp.denom),4);
    res = res || quat_lattice_equal(&lat,&cmp);
    ibz_set(&(lat.denom),1);
    ibz_set(&(cmp.denom),-1);
    res = res || !quat_lattice_equal(&lat,&cmp);
    ibz_set(&(lat.denom),3);
    ibz_set(&(cmp.denom),3);
    res = res || !quat_lattice_equal(&lat,&cmp);
    ibz_set(&(lat.basis[0][0]),1);
    ibz_set(&(lat.basis[0][3]),-1);
    ibz_set(&(lat.basis[1][1]),-2);
    ibz_set(&(lat.basis[2][2]),1);
    ibz_set(&(lat.basis[2][1]),1);
    ibz_set(&(lat.basis[3][3]),-3);
    ibz_set(&(lat.denom),6);
    quat_lattice_hnf(&lat);
    ibz_mat_4x4_copy(&(cmp.basis), &(lat.basis));
    ibz_set(&(cmp.denom),6);
    res = res || !quat_lattice_equal(&lat,&cmp);
    ibz_set(&(cmp.denom),-7);
    res = res || quat_lattice_equal(&lat,&cmp);
    ibz_set(&(cmp.denom),6);
    ibz_set(&(cmp.basis[3][3]),165);
    res = res || quat_lattice_equal(&lat,&cmp);

    if (res != 0){
        printf("Quaternion unit test lattice_equal failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&cmp);
    return(res);
}

//void quat_lattice_reduce_denom(quat_lattice_t *reduced, const quat_lattice_t *lat);
int quat_test_lattice_reduce_denom(){
    int res = 0;
    int s;
    quat_lattice_t red, lat, cmp;
    quat_lattice_init(&red);
    quat_lattice_init(&cmp);
    quat_lattice_init(&lat);

    s = 15;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat.basis[i][j]),(i+j)*s);
            ibz_set(&(cmp.basis[i][j]),(i+j));
        }
    }
    ibz_set(&(lat.denom),4*s);
    ibz_set(&(cmp.denom), 4);

    quat_lattice_reduce_denom(&red,&lat);
    res = res || (!ibz_mat_4x4_equal(&(red.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(red.denom),&(cmp.denom));


    quat_lattice_reduce_denom(&lat,&lat);
    res = res || (!ibz_mat_4x4_equal(&(lat.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(lat.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test lattice_reduce_denom failed\n");
    }
    quat_lattice_finalize(&red);
    quat_lattice_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return(res);
}

//void quat_lattice_dual_without_hnf(quat_lattice_t *dual, const quat_lattice_t *lat);
int quat_test_lattice_dual_without_hnf(){
    int res = 0;
    quat_lattice_t lat, dual,cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&dual);
    quat_lattice_init(&cmp);
    // set lattice
    ibz_mat_4x4_zero(&(lat.basis));
    ibz_set(&(lat.basis[0][0]),1);
    ibz_set(&(lat.basis[0][3]),-1);
    ibz_set(&(lat.basis[1][1]),-2);
    ibz_set(&(lat.basis[2][2]),1);
    ibz_set(&(lat.basis[2][1]),1);
    ibz_set(&(lat.basis[3][3]),-3);
    ibz_set(&(lat.denom),6);
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_set(&(cmp.basis[0][0]),6);
    ibz_set(&(cmp.basis[1][1]),3);
    ibz_set(&(cmp.basis[2][2]),6);
    ibz_set(&(cmp.basis[3][3]),2);
    ibz_set(&(cmp.denom),1);
    quat_lattice_hnf(&lat);
    // test whether dual of dual is original lattice, but dual is not.
    quat_lattice_dual_without_hnf(&dual,&lat);
    quat_lattice_hnf(&dual);
    quat_lattice_hnf(&cmp);
    res = res || !quat_lattice_equal(&dual,&cmp);
    res = res || quat_lattice_equal(&dual,&lat);
    quat_lattice_dual_without_hnf(&dual,&dual);
    quat_lattice_hnf(&dual);
    res = res || !quat_lattice_equal(&dual,&lat);
    

    if (res != 0){
        printf("Quaternion unit test lattice_dual_without_hnf failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&dual);
    quat_lattice_finalize(&cmp);
    return(res);
}

//void quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);
int quat_test_lattice_add(){
    int res = 0;
    quat_lattice_t lat1, lat2, cmp, sum;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&sum);
    quat_lattice_init(&cmp);
    ibz_mat_4x4_zero(&(lat1.basis));
    ibz_mat_4x4_zero(&(lat2.basis));
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_set(&(lat1.basis[0][0]),44);
    ibz_set(&(lat1.basis[0][2]),3);
    ibz_set(&(lat1.basis[0][3]),32);
    ibz_set(&(lat2.basis[0][0]),1);
    ibz_set(&(cmp.basis[0][0]),2);
    ibz_set(&(cmp.basis[0][2]),1);
    ibz_set(&(lat1.basis[1][1]),5);
    ibz_set(&(lat2.basis[1][1]),2);
    ibz_set(&(cmp.basis[1][1]),1);
    ibz_set(&(lat1.basis[2][2]),3);
    ibz_set(&(lat2.basis[2][2]),1);
    ibz_set(&(cmp.basis[2][2]),1);
    ibz_set(&(lat1.basis[3][3]),1);
    ibz_set(&(lat2.basis[3][3]),3);
    ibz_set(&(cmp.basis[3][3]),3);
    ibz_set(&(lat1.denom),4);
    ibz_set(&(lat2.denom),6);
    ibz_set(&(cmp.denom),12);

    quat_lattice_add(&sum,&lat1,&lat2);
    res = res || (!ibz_mat_4x4_equal(&(sum.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(sum.denom),&(cmp.denom));

    // same lattices but not under hnf
    ibz_mat_4x4_zero(&(lat1.basis));
    ibz_mat_4x4_zero(&(lat2.basis));
    ibz_set(&(lat1.basis[0][0]),4);
    ibz_set(&(lat1.basis[0][2]),3);
    ibz_set(&(lat2.basis[0][0]),1);
    ibz_set(&(lat2.basis[0][3]),-1);
    ibz_set(&(lat1.basis[1][1]),5);
    ibz_set(&(lat2.basis[1][1]),-2);
    ibz_set(&(lat1.basis[2][2]),3);
    ibz_set(&(lat2.basis[2][2]),1);
    ibz_set(&(lat2.basis[2][1]),1);
    ibz_set(&(lat1.basis[3][3]),7);
    ibz_set(&(lat2.basis[3][3]),-3);
    ibz_set(&(lat1.denom),4);
    ibz_set(&(lat2.denom),6);

    quat_lattice_add(&sum,&lat1,&lat2);
    res = res || (!ibz_mat_4x4_equal(&(sum.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(sum.denom),&(cmp.denom));

    //double in place gives hnf
    ibz_mat_4x4_copy(&(cmp.basis),&lat2.basis);
    ibz_copy(&(cmp.denom),&(lat2.denom));
    quat_lattice_hnf(&cmp);
    quat_lattice_add(&lat2,&lat2,&lat2);
    res = res || (!ibz_mat_4x4_equal(&(lat2.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(lat2.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test lattice_add failed\n");
    }
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&sum);
    quat_lattice_finalize(&cmp);
    return(res);
}

//void quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2);
int quat_test_lattice_intersect(){
    int res = 0;
    quat_lattice_t lat1,lat2, inter, cmp;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&inter);
    quat_lattice_init(&cmp);
    ibz_mat_4x4_zero(&(cmp.basis));
    ibz_mat_4x4_zero(&(lat1.basis));
    ibz_mat_4x4_zero(&(lat2.basis));
    ibz_set(&(lat1.basis[0][0]),4);
    ibz_set(&(lat1.basis[0][2]),3);
    ibz_set(&(lat2.basis[0][0]),1);
    ibz_set(&(lat2.basis[0][3]),-1);
    ibz_set(&(lat1.basis[1][1]),5);
    ibz_set(&(lat2.basis[1][1]),-2);
    ibz_set(&(lat1.basis[2][2]),3);
    ibz_set(&(lat2.basis[2][2]),1);
    ibz_set(&(lat2.basis[2][1]),1);
    ibz_set(&(lat1.basis[3][3]),7);
    ibz_set(&(lat2.basis[3][3]),-3);
    ibz_set(&(lat1.denom),4);
    ibz_set(&(lat2.denom),6);
    quat_lattice_hnf(&lat1);
    quat_lattice_hnf(&lat2);
    
    ibz_set(&(cmp.basis[0][0]),2);
    ibz_set(&(cmp.basis[0][2]),1);
    ibz_set(&(cmp.basis[1][1]),10);
    ibz_set(&(cmp.basis[2][2]),3);
    ibz_set(&(cmp.basis[3][3]),7);
    ibz_set(&(cmp.denom),2);
    quat_lattice_intersect(&inter,&lat1,&lat2);

    res = res || !quat_lattice_equal(&inter,&cmp);
    quat_lattice_intersect(&lat2,&lat1,&lat2);
    res = res || !quat_lattice_equal(&lat2,&cmp);
    ibz_mat_4x4_copy(&(cmp.basis),&(lat1.basis));
    ibz_copy(&(cmp.denom),&(lat1.denom));
    quat_lattice_intersect(&lat1,&lat1,&lat1);
    res = res || !quat_lattice_equal(&lat1,&cmp);

    if (res != 0){
        printf("Quaternion unit test lattice_intersect failed\n");
    }
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&inter);
    quat_lattice_finalize(&cmp);
    return(res);
}

//void quat_lattice_mul(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2, const quat_alg_t *alg);
int quat_test_lattice_mul(){
    int res = 0;
    quat_lattice_t lat1, lat2, cmp, prod;
    quat_alg_t alg;
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&prod);
    quat_lattice_init(&cmp);
    quat_alg_init_set_ui(&alg, 19);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat1.basis[i][j]),0);
            ibz_set(&(lat2.basis[i][j]),0);
            ibz_set(&(cmp.basis[i][j]),0);
        }
    }

    ibz_set(&(lat1.basis[0][0]),44);
    ibz_set(&(lat1.basis[0][2]),3);
    ibz_set(&(lat1.basis[0][3]),32);
    ibz_set(&(lat2.basis[0][0]),1);
    ibz_set(&(cmp.basis[0][0]),1);
    ibz_set(&(lat1.basis[1][1]),5);
    ibz_set(&(lat2.basis[1][1]),2);
    ibz_set(&(cmp.basis[1][1]),1);
    ibz_set(&(lat1.basis[2][2]),3);
    ibz_set(&(lat2.basis[2][2]),1);
    ibz_set(&(cmp.basis[2][2]),1);
    ibz_set(&(lat1.basis[3][3]),1);
    ibz_set(&(lat2.basis[3][3]),3);
    ibz_set(&(cmp.basis[3][3]),1);
    ibz_set(&(lat1.denom),4);
    ibz_set(&(lat2.denom),6);
    ibz_set(&(cmp.denom),24);

    quat_lattice_mul(&prod,&lat1,&lat2,&alg);
    res = res || (!ibz_mat_4x4_equal(&(prod.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(prod.denom),&(cmp.denom));

    // same lattices but not under hnf
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat1.basis[i][j]),0);
            ibz_set(&(lat2.basis[i][j]),0);
        }
    }
    ibz_set(&(lat1.basis[0][0]),4);
    ibz_set(&(lat1.basis[0][2]),3);
    ibz_set(&(lat2.basis[0][0]),1);
    ibz_set(&(lat2.basis[0][3]),-1);
    ibz_set(&(lat1.basis[1][1]),5);
    ibz_set(&(lat2.basis[1][1]),-2);
    ibz_set(&(lat1.basis[2][2]),3);
    ibz_set(&(lat2.basis[2][2]),1);
    ibz_set(&(lat2.basis[2][1]),1);
    ibz_set(&(lat1.basis[3][3]),7);
    ibz_set(&(lat2.basis[3][3]),-3);
    ibz_set(&(lat1.denom),4);
    ibz_set(&(lat2.denom),6);

    quat_lattice_mul(&prod,&lat1,&lat2,&alg);
    res = res || (!ibz_mat_4x4_equal(&(prod.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(prod.denom),&(cmp.denom));

    //double in place gives hnf
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(cmp.basis[i][j]),0);
        }
    }
    ibz_set(&(cmp.basis[0][0]),1);
    ibz_set(&(cmp.basis[1][1]),1);
    ibz_set(&(cmp.basis[2][2]),1);
    ibz_set(&(cmp.basis[3][3]),1);
    ibz_set(&(cmp.denom),36);
    quat_lattice_mul(&lat2,&lat2,&lat2,&alg);
    res = res || (!ibz_mat_4x4_equal(&(lat2.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(lat2.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test lattice_mul failed\n");
    }
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&prod);
    quat_lattice_finalize(&cmp);
    quat_alg_finalize(&alg);
    return(res);
}

//int quat_lattice_contains_without_alg(quat_alg_coord_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t *x);
int quat_test_lattice_contains_without_alg(){
    int res = 0;
    quat_alg_elem_t x;
    quat_alg_coord_t coord, cmp;
    quat_lattice_t lat;
    quat_alg_elem_init(&x);
    quat_alg_coord_init(&coord);
    quat_alg_coord_init(&cmp);
    quat_lattice_init(&lat);

    // lattice 1
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat.basis[i][j]),0);
        }
    }
    ibz_set(&(lat.basis[0][0]),4);
    ibz_set(&(lat.basis[0][2]),3);
    ibz_set(&(lat.basis[1][1]),5);
    ibz_set(&(lat.basis[2][2]),3);
    ibz_set(&(lat.basis[3][3]),7);
    ibz_set(&(lat.denom),4);

    // x 1, should fail
    ibz_set(&(x.denom),3);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),-2);
    ibz_set(&(x.coord[2]),26);
    ibz_set(&(x.coord[3]),9);

    res = res || quat_lattice_contains_without_alg(&coord,&lat,&x);
    // again, but with NULL
    res = res || quat_lattice_contains_without_alg(NULL,&lat,&x);


    // lattice 2
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat.basis[i][j]),0);
        }
    }
    ibz_set(&(lat.basis[0][0]),1);
    ibz_set(&(lat.basis[0][3]),-1);
    ibz_set(&(lat.basis[1][1]),-2);
    ibz_set(&(lat.basis[2][2]),1);
    ibz_set(&(lat.basis[2][1]),1);
    ibz_set(&(lat.basis[3][3]),-3);
    ibz_set(&(lat.denom),6);
    quat_lattice_hnf(&lat);
    // x 1, should succeed
    ibz_set(&(x.denom),3);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),-2);
    ibz_set(&(x.coord[2]),26);
    ibz_set(&(x.coord[3]),9);
    ibz_set(&(cmp[0]),2);
    ibz_set(&(cmp[1]),-2);
    ibz_set(&(cmp[2]),52);
    ibz_set(&(cmp[3]),6);

    res = res || (0==quat_lattice_contains_without_alg(&coord,&lat,&x));

    res = res || ibz_cmp(&(coord[0]),&(cmp[0]));
    res = res || ibz_cmp(&(coord[1]),&(cmp[1]));
    res = res || ibz_cmp(&(coord[2]),&(cmp[2]));
    res = res || ibz_cmp(&(coord[3]),&(cmp[3]));
    // again, but with NULL
    res = res || (0==quat_lattice_contains_without_alg(NULL,&lat,&x));


    if (res != 0){
        printf("Quaternion unit test lattice_contains_without_alg failed\n");
    }
    quat_alg_elem_finalize(&x);
    quat_alg_coord_finalize(&coord);
    quat_alg_coord_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return(res);
}

//int quat_lattice_contains(quat_alg_coord_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t *x, const quat_alg_t *alg);
int quat_test_lattice_contains(){
    int res = 0;
    quat_alg_t alg;
    quat_alg_elem_t x;
    quat_alg_coord_t coord, cmp;
    quat_lattice_t lat;
    quat_alg_init_set_ui(&alg, 103);
    quat_alg_elem_init(&x);
    quat_alg_coord_init(&coord);
    quat_alg_coord_init(&cmp);
    quat_lattice_init(&lat);

    // lattice 1
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat.basis[i][j]),0);
        }
    }
    ibz_set(&(lat.basis[0][0]),4);
    ibz_set(&(lat.basis[0][2]),3);
    ibz_set(&(lat.basis[1][1]),5);
    ibz_set(&(lat.basis[2][2]),3);
    ibz_set(&(lat.basis[3][3]),7);
    ibz_set(&(lat.denom),4);

    // x 1, should fail
    ibz_set(&(x.denom),3);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),-2);
    ibz_set(&(x.coord[2]),26);
    ibz_set(&(x.coord[3]),9);

    res = res || quat_lattice_contains(&coord,&lat,&x,&alg);
    // again, but with NULL
    res = res || quat_lattice_contains(NULL,&lat,&x,&alg);


    // lattice 2
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat.basis[i][j]),0);
        }
    }
    ibz_set(&(lat.basis[0][0]),1);
    ibz_set(&(lat.basis[0][3]),-1);
    ibz_set(&(lat.basis[1][1]),-2);
    ibz_set(&(lat.basis[2][2]),1);
    ibz_set(&(lat.basis[2][1]),1);
    ibz_set(&(lat.basis[3][3]),-3);
    ibz_set(&(lat.denom),6);
    quat_lattice_hnf(&lat);
    // x 1, should succeed
    ibz_set(&(x.denom),3);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),-2);
    ibz_set(&(x.coord[2]),26);
    ibz_set(&(x.coord[3]),9);
    ibz_set(&(cmp[0]),2);
    ibz_set(&(cmp[1]),-2);
    ibz_set(&(cmp[2]),52);
    ibz_set(&(cmp[3]),6);

    res = res || (0==quat_lattice_contains(&coord,&lat,&x,&alg));

    res = res || ibz_cmp(&(coord[0]),&(cmp[0]));
    res = res || ibz_cmp(&(coord[1]),&(cmp[1]));
    res = res || ibz_cmp(&(coord[2]),&(cmp[2]));
    res = res || ibz_cmp(&(coord[3]),&(cmp[3]));
    // again, but with NULL
    res = res || (0==quat_lattice_contains(NULL,&lat,&x,&alg));


    if (res != 0){
        printf("Quaternion unit test lattice_contains failed\n");
    }
    quat_alg_finalize(&alg);
    quat_alg_elem_finalize(&x);
    quat_alg_coord_finalize(&coord);
    quat_alg_coord_finalize(&cmp);
    quat_lattice_finalize(&lat);
    return(res);
}

//void quat_lattice_index(ibz_t *index, const quat_lattice_t *sublat, const quat_lattice_t *overlat);
int quat_test_lattice_index(){
    int res = 0;
    quat_lattice_t sublat,overlat;
    ibz_t index;
    ibz_init(&index);
    quat_lattice_init(&sublat);
    quat_lattice_init(&overlat);

    ibz_mat_4x4_zero(&(sublat.basis));
    ibz_mat_4x4_identity(&(overlat.basis));
    ibz_set(&(overlat.denom),2);
    ibz_set(&(sublat.basis[0][0]),2);
    ibz_set(&(sublat.basis[0][1]),0);
    ibz_set(&(sublat.basis[0][2]),1);
    ibz_set(&(sublat.basis[0][3]),0);
    ibz_set(&(sublat.basis[1][0]),0);
    ibz_set(&(sublat.basis[1][1]),4);
    ibz_set(&(sublat.basis[1][2]),2);
    ibz_set(&(sublat.basis[1][3]),3);
    ibz_set(&(sublat.basis[2][0]),0);
    ibz_set(&(sublat.basis[2][1]),0);
    ibz_set(&(sublat.basis[2][2]),1);
    ibz_set(&(sublat.basis[2][3]),0);
    ibz_set(&(sublat.basis[3][0]),0);
    ibz_set(&(sublat.basis[3][1]),0);
    ibz_set(&(sublat.basis[3][2]),0);
    ibz_set(&(sublat.basis[3][3]),1);
    ibz_set(&(sublat.denom),2);
    quat_lattice_index(&index,&sublat,&overlat);

    res = res || !(ibz_get(&index)==8);

    if (res != 0){
        printf("Quaternion unit test lattice_index failed\n");
    }
    quat_lattice_finalize(&sublat);
    quat_lattice_finalize(&overlat);
    ibz_finalize(&index);
    return(res);
}

//void quat_lattice_hnf(quat_lattice_t *lat);
int quat_test_lattice_hnf(){
    int res = 0;
    quat_lattice_t lat, cmp;
    quat_lattice_init(&lat);
    quat_lattice_init(&cmp);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(lat.basis[i][j]),0);
            ibz_set(&(cmp.basis[i][j]),0);
        }
    }
    ibz_set(&(lat.basis[0][0]),1);
    ibz_set(&(lat.basis[0][3]),-1);
    ibz_set(&(lat.basis[1][1]),-2);
    ibz_set(&(lat.basis[2][2]),1);
    ibz_set(&(lat.basis[2][1]),1);
    ibz_set(&(lat.basis[3][3]),-3);
    ibz_set(&(cmp.basis[0][0]),1);
    ibz_set(&(cmp.basis[1][1]),2);
    ibz_set(&(cmp.basis[2][2]),1);
    ibz_set(&(cmp.basis[3][3]),3);
    ibz_set(&(cmp.denom),6);
    ibz_set(&(lat.denom),6);

    quat_lattice_hnf(&lat);
    res = res || (!ibz_mat_4x4_equal(&(lat.basis),&(cmp.basis)));
    res = res || ibz_cmp(&(lat.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test lattice_hnf failed\n");
    }
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&cmp);
    return(res);
}

//int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int precision);
int quat_test_lattice_lll(){
    int res = 0;
    quat_lattice_t lat, test;
    ibz_mat_4x4_t red;
    ibz_t num, denom, q;
    ibq_t coeff;
    ibz_init(&num);
    ibz_init(&denom);
    ibz_init(&q);
    ibq_init(&coeff);
    ibz_mat_4x4_init(&red);
    quat_lattice_init(&lat);
    quat_lattice_init(&test);

    // set lattice
    ibz_set(&lat.denom, 60);
    ibz_mat_4x4_zero(&(lat.basis));
    ibz_set(&lat.basis[0][0], 3);
    ibz_set(&lat.basis[1][0], 7);
    ibz_set(&lat.basis[0][1], 1);
    ibz_set(&lat.basis[3][1], -6);
    ibz_set(&lat.basis[1][2], 12);
    ibz_set(&lat.basis[2][2], 5);
    ibz_set(&lat.basis[0][3], -19);
    ibz_set(&lat.basis[3][3], 3);
    
    quat_lattice_hnf(&lat);

    ibz_set(&q,103);
    res = res || quat_lattice_lll(&red,&lat,&q,10);
    // test lll reduced
    ibz_set(&num,3);
    ibz_set(&denom,4);
    ibq_set(&coeff,&num,&denom);
    res = res || !quat_dim4_lll_verify(&red,&coeff,&q);
    // test lattice equality
    ibz_copy(&(test.denom),&(lat.denom));
    ibz_mat_4x4_copy(&(test.basis),&(red));
    quat_lattice_hnf(&test);
    res = res || !quat_lattice_equal(&test,&lat);
    
    if (res != 0){
        printf("Quaternion unit test lattice_lll failed\n");
    }
    ibz_finalize(&num);
    ibz_finalize(&denom);
    ibz_finalize(&q);
    ibq_finalize(&coeff);
    ibz_mat_4x4_finalize(&red);
    quat_lattice_finalize(&lat);
    quat_lattice_finalize(&test);
    return(res);
}

//int quat_lattice_random_elem(quat_alg_elem_t *elem, const quat_lattice_t *lattice, unsigned char n)
int quat_test_lattice_random_elem(){
    int res = 0;
    
    quat_lattice_t lat;
    quat_alg_elem_t elem;
    quat_lattice_init(&lat);
    quat_alg_elem_init(&elem);

    ibz_set(&lat.denom, 60);
    ibz_set(&lat.basis[0][0], 3);
    ibz_set(&lat.basis[1][0], 7);
    ibz_set(&lat.basis[0][1], 1);
    ibz_set(&lat.basis[3][1], -6);
    ibz_set(&lat.basis[1][2], 12);
    ibz_set(&lat.basis[2][2], 5);
    ibz_set(&lat.basis[0][3], -19);
    ibz_set(&lat.basis[3][3], 3);
    quat_lattice_hnf(&lat);

    for (int i = 0; i <= 10; i++) {
        int prng = quat_lattice_random_elem(&elem, &lat, 3);
        assert(prng);
        res |= !quat_lattice_contains_without_alg(&elem.coord, &lat, &elem);
        for (int j = 0; j < 4; j++) {
            res |= !(ibz_get(&elem.coord[j]) < (1 << 2));
            res |= !(-(1 << 2) <= ibz_get(&elem.coord[j]));
        }
    }
    
    quat_lattice_finalize(&lat);
    quat_alg_elem_finalize(&elem);
    
    if (res != 0){
        printf("Quaternion unit test lattice_random_elem failed\n");
    }
    return res;
}

//void quat_lattice_right_transporter(quat_lattice_t *trans, const quat_lattice_t *lat1, const quat_lattice_t *lat1)
int quat_test_lattice_right_transporter_hnf(){
    int res = 0;

    quat_alg_t alg;
    quat_lattice_t lat1, lat2, trans, exp;
    quat_alg_init_set_ui(&alg, 103);
    quat_lattice_init(&lat1);
    quat_lattice_init(&lat2);
    quat_lattice_init(&trans);
    quat_lattice_init(&exp);

    ibz_set(&lat1.basis[0][0], 9);
    ibz_set(&lat1.basis[0][1], 4);
    ibz_set(&lat1.basis[0][2], 7);
    ibz_set(&lat1.basis[1][1], 12);
    ibz_set(&lat1.basis[1][3], 6);
    ibz_set(&lat1.basis[2][2], 8);
    ibz_set(&lat1.basis[2][3], 3);
    ibz_set(&lat1.basis[3][3], 3);
    ibz_set(&lat1.denom, 11);

    quat_lattice_mul(&lat2, &lat1, &lat1, &alg);
    quat_lattice_right_transporter(&trans, &lat1, &lat2, &alg);

    ibz_set(&exp.basis[0][0], 1);
    ibz_set(&exp.basis[1][1], 12);
    ibz_set(&exp.basis[1][3], 6);
    ibz_set(&exp.basis[2][2], 8);
    ibz_set(&exp.basis[2][3], 3);
    ibz_set(&exp.basis[3][3], 3);
    ibz_set(&exp.denom, 11);

    res |= !quat_lattice_equal(&trans, &exp);

    ibz_set(&(alg.p),19);
    ibz_set(&(alg.gram[2][2]),19);
    ibz_set(&(alg.gram[3][3]),19);
    ibz_set(&lat1.basis[0][0], -12);
    ibz_set(&lat1.basis[0][1], 60);
    ibz_set(&lat1.basis[0][2], 42);
    ibz_set(&lat1.basis[0][3], 60);
    ibz_set(&lat1.basis[1][0], 9);
    ibz_set(&lat1.basis[1][1], 24);
    ibz_set(&lat1.basis[1][2], 0);
    ibz_set(&lat1.basis[1][3], -21);
    ibz_set(&lat1.basis[2][0], 0);
    ibz_set(&lat1.basis[2][1], -18);
    ibz_set(&lat1.basis[2][2], 27);
    ibz_set(&lat1.basis[2][3], 3);
    ibz_set(&lat1.basis[3][0], 0);
    ibz_set(&lat1.basis[3][1], 0);
    ibz_set(&lat1.basis[3][2], 15);
    ibz_set(&lat1.basis[3][3], 3);
    ibz_set(&lat1.denom, 15);

    quat_lattice_hnf(&lat1);
    quat_lattice_reduce_denom(&lat1, &lat1);
    quat_lattice_mul(&lat2, &lat1, &lat1, &alg);
    quat_lattice_right_transporter(&trans, &lat1, &lat2, &alg);
    ibz_mat_4x4_zero(&(exp.basis));
    ibz_set(&exp.basis[0][0], 2);
    ibz_set(&exp.basis[1][1], 1);
    ibz_set(&exp.basis[2][2], 2);
    ibz_set(&exp.basis[2][3], 1);
    ibz_set(&exp.basis[3][3], 1);
    ibz_set(&exp.denom, 5);

    res |= !quat_lattice_equal(&trans, &exp);
    
    if (res != 0){
        printf("Quaternion unit test lattice_right_transporter_hnf failed\n");
    }
    quat_alg_finalize(&alg);
    quat_lattice_finalize(&lat1);
    quat_lattice_finalize(&lat2);
    quat_lattice_finalize(&trans);
    quat_lattice_finalize(&exp);
    return(res);
}


// run all previous tests
int quat_test_lattice(){
    int res = 0;
    printf("\nRunning quaternion tests of lattice functions\n");
    res = res | quat_test_lattice_equal();
    res = res | quat_test_lattice_reduce_denom();
    res = res | quat_test_lattice_dual_without_hnf();
    res = res | quat_test_lattice_add();
    res = res | quat_test_lattice_intersect();
    res = res | quat_test_lattice_mul();
    res = res | quat_test_lattice_contains_without_alg();
    res = res | quat_test_lattice_contains();
    res = res | quat_test_lattice_index();
    res = res | quat_test_lattice_hnf();
    res = res | quat_test_lattice_lll();
    res = res | quat_test_lattice_random_elem();
    res = res | quat_test_lattice_right_transporter_hnf();
    return(res);
}
