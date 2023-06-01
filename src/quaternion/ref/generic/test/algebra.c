#include "quaternion_tests.h"

// tests of internal helper functions

//static inline void quat_alg_init_set_ui(quat_alg_t *alg, unsigned int p);
int quat_test_init_set_ui(){
    int res = 0;
    int p = 5;
    quat_alg_t alg;
    ibz_mat_4x4_t cmp;
    ibz_mat_4x4_init(&cmp);
    quat_alg_init_set_ui(&alg, p);
    res = res || (p != ibz_get(&(alg.p)));
    ibz_mat_4x4_identity(&cmp);
    ibz_set(&(cmp[2][2]),p);
    ibz_set(&(cmp[3][3]),p);
    res = res || !ibz_mat_4x4_equal(&cmp,&(alg.gram));
    if (res != 0){
        printf("Quaternion unit test alg_init_set_ui failed\n");
    }
    quat_alg_finalize(&alg);
    ibz_mat_4x4_finalize(&cmp);
    return(res);
}

//void quat_alg_coord_add(quat_alg_coord_t *res, const quat_alg_coord_t *a, const quat_alg_coord_t *b);
int quat_test_alg_coord_add(){
    int res = 0;
    quat_alg_coord_t a, b, c, cmp;
    quat_alg_coord_init(&a);
    quat_alg_coord_init(&b);
    quat_alg_coord_init(&c);
    quat_alg_coord_init(&cmp);


    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),-2);
    ibz_set(&(a[2]),7);
    ibz_set(&(a[3]),199);
    ibz_set(&(b[0]),-6);
    ibz_set(&(b[1]),2);
    ibz_set(&(b[2]),67);
    ibz_set(&(b[3]),-22);
    ibz_set(&(cmp[0]),-5);
    ibz_set(&(cmp[1]),0);
    ibz_set(&(cmp[2]),74);
    ibz_set(&(cmp[3]),177);
    quat_alg_coord_add(&c,&a,&b);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c[i]),&(cmp[i]));
    }

    ibz_set(&(a[0]),-122);
    ibz_set(&(a[1]),0);
    ibz_set(&(a[2]),-7);
    ibz_set(&(a[3]),1889);
    ibz_set(&(b[0]),-6);
    ibz_set(&(b[1]),2);
    ibz_set(&(b[2]),67);
    ibz_set(&(b[3]),-1889);
    ibz_set(&(cmp[0]),-128);
    ibz_set(&(cmp[1]),2);
    ibz_set(&(cmp[2]),60);
    ibz_set(&(cmp[3]),0);
    quat_alg_coord_add(&c,&a,&b);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c[i]),&(cmp[i]));
    }

    ibz_set(&(a[0]),-1);
    ibz_set(&(a[1]),2);
    ibz_set(&(a[2]),-7);
    ibz_set(&(a[3]),19);;
    ibz_set(&(cmp[0]),-2);
    ibz_set(&(cmp[1]),4);
    ibz_set(&(cmp[2]),-14);
    ibz_set(&(cmp[3]),38);
    quat_alg_coord_add(&a,&a,&a);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(a[i]),&(cmp[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_coord_add failed\n");
    }
    quat_alg_coord_finalize(&a);
    quat_alg_coord_finalize(&b);
    quat_alg_coord_finalize(&c);
    quat_alg_coord_finalize(&cmp);
    return(res);
}

//void quat_alg_coord_sub(quat_alg_coord_t *res, const quat_alg_coord_t *a, const quat_alg_coord_t *b);
int quat_test_alg_coord_sub(){
    int res = 0;
    quat_alg_coord_t a, b, c, cmp;
    quat_alg_coord_init(&a);
    quat_alg_coord_init(&b);
    quat_alg_coord_init(&c);
    quat_alg_coord_init(&cmp);

    ibz_set(&(a[0]),1);
    ibz_set(&(a[1]),-2);
    ibz_set(&(a[2]),7);
    ibz_set(&(a[3]),199);
    ibz_set(&(b[0]),-6);
    ibz_set(&(b[1]),2);
    ibz_set(&(b[2]),67);
    ibz_set(&(b[3]),-22);
    ibz_set(&(cmp[0]),7);
    ibz_set(&(cmp[1]),-4);
    ibz_set(&(cmp[2]),-60);
    ibz_set(&(cmp[3]),221);
    quat_alg_coord_sub(&c,&a,&b);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c[i]),&(cmp[i]));
    }

    ibz_set(&(a[0]),-122);
    ibz_set(&(a[1]),0);
    ibz_set(&(a[2]),-7);
    ibz_set(&(a[3]),1889);
    ibz_set(&(b[0]),-6);
    ibz_set(&(b[1]),2);
    ibz_set(&(b[2]),67);
    ibz_set(&(b[3]),-1889);
    ibz_set(&(cmp[0]),-116);
    ibz_set(&(cmp[1]),-2);
    ibz_set(&(cmp[2]),-74);
    ibz_set(&(cmp[3]),3778);
    quat_alg_coord_sub(&c,&a,&b);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c[i]),&(cmp[i]));
    }

    ibz_set(&(a[0]),-1);
    ibz_set(&(a[1]),2);
    ibz_set(&(a[2]),-7);
    ibz_set(&(a[3]),19);;
    ibz_set(&(cmp[0]),0);
    ibz_set(&(cmp[1]),0);
    ibz_set(&(cmp[2]),0);
    ibz_set(&(cmp[3]),0);
    quat_alg_coord_sub(&a,&a,&a);
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(a[i]),&(cmp[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_coord_sub failed\n");
    }
    quat_alg_coord_finalize(&a);
    quat_alg_coord_finalize(&b);
    quat_alg_coord_finalize(&c);
    quat_alg_coord_finalize(&cmp);
    return(res);
}

//void quat_alg_equal_denom(quat_alg_elem_t *res_a, quat_alg_elem_t *res_b, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
int quat_test_alg_equal_denom(){
    int res = 0;
    quat_alg_elem_t a, b, res_a, res_b, cmp_a, cmp_b;
    quat_alg_elem_init(&a);
    quat_alg_elem_init(&b);
    quat_alg_elem_init(&res_a);
    quat_alg_elem_init(&res_b);
    quat_alg_elem_init(&cmp_a);
    quat_alg_elem_init(&cmp_b);

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),67);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),3);
    ibz_set(&(cmp_a.coord[0]),-12);
    ibz_set(&(cmp_a.coord[1]),0);
    ibz_set(&(cmp_a.coord[2]),-7);
    ibz_set(&(cmp_a.coord[3]),19);
    ibz_set(&(cmp_a.denom),9);
    ibz_set(&(cmp_b.coord[0]),-18);
    ibz_set(&(cmp_b.coord[1]),6);
    ibz_set(&(cmp_b.coord[2]),201);
    ibz_set(&(cmp_b.coord[3]),-57);
    ibz_set(&(cmp_b.denom),9);
    quat_alg_equal_denom(&res_a,&res_b,&a,&b);
    res = res || ibz_cmp(&(res_a.denom),&(cmp_a.denom));
    res = res || ibz_cmp(&(res_b.denom),&(cmp_b.denom));
    res = res || ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(res_a.coord[i]),&(cmp_a.coord[i]));
        res = res || ibz_cmp(&(res_b.coord[i]),&(cmp_b.coord[i]));
    }


    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),67);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),6);
    ibz_set(&(cmp_a.coord[0]),-24);
    ibz_set(&(cmp_a.coord[1]),0);
    ibz_set(&(cmp_a.coord[2]),-14);
    ibz_set(&(cmp_a.coord[3]),38);
    ibz_set(&(cmp_a.denom),18);
    ibz_set(&(cmp_b.coord[0]),-18);
    ibz_set(&(cmp_b.coord[1]),6);
    ibz_set(&(cmp_b.coord[2]),201);
    ibz_set(&(cmp_b.coord[3]),-57);
    ibz_set(&(cmp_b.denom),18);
    quat_alg_equal_denom(&res_a,&res_b,&a,&b);
    res = res || ibz_cmp(&(res_a.denom),&(cmp_a.denom));
    res = res || ibz_cmp(&(res_b.denom),&(cmp_b.denom));
    res = res || ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(res_a.coord[i]),&(cmp_a.coord[i]));
        res = res || ibz_cmp(&(res_b.coord[i]),&(cmp_b.coord[i]));
    }

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),6);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),67);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),6);
    ibz_set(&(cmp_a.coord[0]),-12);
    ibz_set(&(cmp_a.coord[1]),0);
    ibz_set(&(cmp_a.coord[2]),-7);
    ibz_set(&(cmp_a.coord[3]),19);
    ibz_set(&(cmp_a.denom),6);
    ibz_set(&(cmp_b.coord[0]),-6);
    ibz_set(&(cmp_b.coord[1]),2);
    ibz_set(&(cmp_b.coord[2]),67);
    ibz_set(&(cmp_b.coord[3]),-19);
    ibz_set(&(cmp_b.denom),6);
    quat_alg_equal_denom(&res_a,&res_b,&a,&b);
    res = res || ibz_cmp(&(res_a.denom),&(cmp_a.denom));
    res = res || ibz_cmp(&(res_b.denom),&(cmp_b.denom));
    res = res || ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(res_a.coord[i]),&(cmp_a.coord[i]));
        res = res || ibz_cmp(&(res_b.coord[i]),&(cmp_b.coord[i]));
    }

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),6);
    ibz_set(&(cmp_a.coord[0]),-12);
    ibz_set(&(cmp_a.coord[1]),0);
    ibz_set(&(cmp_a.coord[2]),-7);
    ibz_set(&(cmp_a.coord[3]),19);
    ibz_set(&(cmp_a.denom),6);
    ibz_set(&(cmp_b.coord[0]),-12);
    ibz_set(&(cmp_b.coord[1]),0);
    ibz_set(&(cmp_b.coord[2]),-7);
    ibz_set(&(cmp_b.coord[3]),19);
    ibz_set(&(cmp_b.denom),6);
    quat_alg_equal_denom(&a,&b,&a,&a);
    res = res || ibz_cmp(&(a.denom),&(a.denom));
    res = res || ibz_cmp(&(b.denom),&(b.denom));
    res = res || ibz_cmp(&(cmp_a.denom),&(cmp_b.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(a.coord[i]),&(a.coord[i]));
        res = res || ibz_cmp(&(b.coord[i]),&(b.coord[i]));
    }
    if (res != 0){
        printf("Quaternion unit test alg_equal_denom failed\n");
    }
    quat_alg_elem_finalize(&a);
    quat_alg_elem_finalize(&b);
    quat_alg_elem_finalize(&res_a);
    quat_alg_elem_finalize(&res_b);
    quat_alg_elem_finalize(&cmp_a);
    quat_alg_elem_finalize(&cmp_b);
    return(res);
}


//Tests of public functions


//void quat_alg_add(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
int quat_test_alg_add(){
    int res = 0;
    quat_alg_elem_t a, b, c, cmp;
    quat_alg_elem_init(&a);
    quat_alg_elem_init(&b);
    quat_alg_elem_init(&c);
    quat_alg_elem_init(&cmp);

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),7);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),3);
    ibz_set(&(cmp.coord[0]),-30);
    ibz_set(&(cmp.coord[1]),6);
    ibz_set(&(cmp.coord[2]),14);
    ibz_set(&(cmp.coord[3]),-38);
    ibz_set(&(cmp.denom),9);
    quat_alg_add(&c, &a, &b);
    res = res || ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),7);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),6);
    ibz_set(&(cmp.coord[0]),-42);
    ibz_set(&(cmp.coord[1]),6);
    ibz_set(&(cmp.coord[2]),7);
    ibz_set(&(cmp.coord[3]),-19);
    ibz_set(&(cmp.denom),18);
    quat_alg_add(&c, &a, &b);
    res = res || ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(cmp.coord[0]),-24);
    ibz_set(&(cmp.coord[1]),0);
    ibz_set(&(cmp.coord[2]),-14);
    ibz_set(&(cmp.coord[3]),38);
    ibz_set(&(cmp.denom),9);
    quat_alg_add(&a, &a, &a);
    res = res || ibz_cmp(&(a.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(a.coord[i]),&(cmp.coord[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_add failed\n");
    }
    quat_alg_elem_finalize(&a);
    quat_alg_elem_finalize(&b);
    quat_alg_elem_finalize(&c);
    quat_alg_elem_finalize(&cmp);
    return(res);
}

//void quat_alg_sub(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b);
int quat_test_alg_sub(){
    int res = 0;
    quat_alg_elem_t a, b, c, cmp;
    quat_alg_elem_init(&a);
    quat_alg_elem_init(&b);
    quat_alg_elem_init(&c);
    quat_alg_elem_init(&cmp);

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),7);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),3);
    ibz_set(&(cmp.coord[0]),-12-3*(-6));
    ibz_set(&(cmp.coord[1]),-3*2);
    ibz_set(&(cmp.coord[2]),-7-3*7);
    ibz_set(&(cmp.coord[3]),19-3*(-19));
    ibz_set(&(cmp.denom),9);
    quat_alg_sub(&c, &a, &b);
    res = res || ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){
        res = res || ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(b.coord[0]),-6);
    ibz_set(&(b.coord[1]),2);
    ibz_set(&(b.coord[2]),7);
    ibz_set(&(b.coord[3]),-19);
    ibz_set(&(b.denom),6);
    ibz_set(&(cmp.coord[0]),-2*12-3*(-6));
    ibz_set(&(cmp.coord[1]),-3*2);
    ibz_set(&(cmp.coord[2]),-2*7-3*7);
    ibz_set(&(cmp.coord[3]),2*19-3*(-19));
    ibz_set(&(cmp.denom),18);
    quat_alg_sub(&a, &a, &b);
    res = res || ibz_cmp(&(a.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(a.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(a.coord[0]),-12);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),-7);
    ibz_set(&(a.coord[3]),19);
    ibz_set(&(a.denom),9);
    ibz_set(&(cmp.coord[0]),0);
    ibz_set(&(cmp.coord[1]),0);
    ibz_set(&(cmp.coord[2]),0);
    ibz_set(&(cmp.coord[3]),0);
    ibz_set(&(cmp.denom),9);
    quat_alg_sub(&a, &a, &a);
    res = res || ibz_cmp(&(a.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(a.coord[i]),&(cmp.coord[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_sub failed\n");
    }
    quat_alg_elem_finalize(&a);
    quat_alg_elem_finalize(&b);
    quat_alg_elem_finalize(&c);
    quat_alg_elem_finalize(&cmp);
    return(res);
}

//void quat_alg_mul(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b, const quat_alg_t *alg);
int quat_test_alg_mul(){
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 7);
    quat_alg_elem_t a, b, c, cmp;
    quat_alg_elem_init(&a);
    quat_alg_elem_init(&b);
    quat_alg_elem_init(&c);
    quat_alg_elem_init(&cmp);

    ibz_set(&(a.coord[0]),152);
    ibz_set(&(a.coord[1]),57);
    ibz_set(&(a.coord[2]),190);
    ibz_set(&(a.coord[3]),28);
    ibz_set(&(a.denom),76);
    ibz_set(&(b.coord[0]),165);
    ibz_set(&(b.coord[1]),35);
    ibz_set(&(b.coord[2]),231);
    ibz_set(&(b.coord[3]),770);
    ibz_set(&(b.denom),385);
    ibz_set(&(cmp.coord[0]),-435065);
    ibz_set(&(cmp.coord[1]),993549);
    ibz_set(&(cmp.coord[2]),23552);
    ibz_set(&(cmp.coord[3]),128177);
    ibz_set(&(cmp.denom),29260);
    quat_alg_mul(&c, &a, &b, &alg);
    res = res || ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }
    ibz_set(&(alg.p),11);
    ibz_set(&(cmp.coord[0]),-696865);
    ibz_set(&(cmp.coord[1]),1552877);
    ibz_set(&(cmp.denom),29260);
    quat_alg_mul(&c, &a, &b, &alg);
    res = res || ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(alg.p),7);
    ibz_set(&(a.coord[0]),1);
    ibz_set(&(a.coord[1]),1);
    ibz_set(&(a.coord[2]),1);
    ibz_set(&(a.coord[3]),1);
    ibz_set(&(a.denom),2);
    ibz_set(&(cmp.coord[0]),-14);
    ibz_set(&(cmp.coord[1]),2);
    ibz_set(&(cmp.coord[2]),2);
    ibz_set(&(cmp.coord[3]),2);
    ibz_set(&(cmp.denom),4);
    quat_alg_mul(&c, &a, &a, &alg);
    res = res || ibz_cmp(&(c.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(c.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(alg.p),7);
    ibz_set(&(a.coord[0]),1);
    ibz_set(&(a.coord[1]),1);
    ibz_set(&(a.coord[2]),1);
    ibz_set(&(a.coord[3]),1);
    ibz_set(&(a.denom),2);
    ibz_set(&(cmp.coord[0]),-14);
    ibz_set(&(cmp.coord[1]),2);
    ibz_set(&(cmp.coord[2]),2);
    ibz_set(&(cmp.coord[3]),2);
    ibz_set(&(cmp.denom),4);
    quat_alg_mul(&a, &a, &a, &alg);
    res = res || ibz_cmp(&(a.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(a.coord[i]),&(cmp.coord[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_mul failed\n");
    }
    quat_alg_elem_finalize(&a);
    quat_alg_elem_finalize(&b);
    quat_alg_elem_finalize(&c);
    quat_alg_elem_finalize(&cmp);
    quat_alg_finalize(&alg);
    return(res);
}

//void quat_alg_norm(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_t *alg);
int quat_test_alg_norm(){
    int res = 0;
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 11);
    quat_alg_elem_t a;
    ibq_t norm, cmp;
    ibz_t num, denom;
    quat_alg_elem_init(&a);
    ibq_init(&norm);
    ibq_init(&cmp);
    ibz_init(&num);
    ibz_init(&denom);

    ibz_set(&(alg.p),11);
    ibz_set(&(a.coord[0]),1);
    ibz_set(&(a.coord[1]),5);
    ibz_set(&(a.coord[2]),7);
    ibz_set(&(a.coord[3]),2);
    ibz_set(&(a.denom),2);
    ibz_set(&num,609);
    ibz_set(&denom,4);
    ibq_set(&cmp,&num,&denom);
    quat_alg_norm(&norm,&a,&alg);
    res = res || (ibq_cmp(&norm,&cmp));

    // same vector, not reduced
    ibz_set(&(alg.p),11);
    ibz_set(&(a.coord[0]),2);
    ibz_set(&(a.coord[1]),10);
    ibz_set(&(a.coord[2]),14);
    ibz_set(&(a.coord[3]),4);
    ibz_set(&(a.denom),4);
    ibz_set(&num,609);
    ibz_set(&denom,4);
    ibq_set(&cmp,&num,&denom);
    quat_alg_norm(&norm,&a,&alg);
    res = res || (ibq_cmp(&norm,&cmp));

    ibz_set(&(alg.p),11);
    ibz_set(&(a.coord[0]),152);
    ibz_set(&(a.coord[1]),57);
    ibz_set(&(a.coord[2]),190);
    ibz_set(&(a.coord[3]),28);
    ibz_set(&(a.denom),76);
    ibz_set(&num,432077);
    ibz_set(&denom,5776);
    ibq_set(&cmp,&num,&denom);
    quat_alg_norm(&norm,&a,&alg);
    res = res || (ibq_cmp(&norm,&cmp));

    ibz_set(&(alg.p),11);
    ibz_set(&(a.coord[0]),0);
    ibz_set(&(a.coord[1]),12);
    ibz_set(&(a.coord[2]),35);
    ibz_set(&(a.coord[3]),49);
    ibz_set(&(a.denom),28);
    ibz_set(&num,20015);
    ibz_set(&denom,392);
    ibq_set(&cmp,&num,&denom);
    quat_alg_norm(&norm,&a,&alg);
    res = res || (ibq_cmp(&norm,&cmp));

    ibz_set(&(alg.p),7);
    ibz_set(&(a.coord[0]),152);
    ibz_set(&(a.coord[1]),57);
    ibz_set(&(a.coord[2]),190);
    ibz_set(&(a.coord[3]),28);
    ibz_set(&(a.denom),76);
    ibz_set(&num,284541);
    ibz_set(&denom,5776);
    ibq_set(&cmp,&num,&denom);
    quat_alg_norm(&norm,&a,&alg);
    res = res || (ibq_cmp(&norm,&cmp));

    if (res != 0){
        printf("Quaternion unit test alg_norm failed\n");
    }
    quat_alg_elem_finalize(&a);
    ibq_finalize(&norm);
    ibq_finalize(&cmp);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    quat_alg_finalize(&alg);
    return(res);
}

//void quat_alg_trace(quat_alg_elem_t *res, const quat_alg_elem_t *a);
int quat_test_alg_trace(){
    int res = 0;
    quat_alg_elem_t a;
    ibq_t trace, cmp;
    ibz_t num, denom;
    quat_alg_elem_init(&a);
    ibq_init(&trace);
    ibq_init(&cmp);
    ibz_init(&num);
    ibz_init(&denom);

    ibz_set(&(a.coord[0]),152);
    ibz_set(&(a.coord[1]),57);
    ibz_set(&(a.coord[2]),190);
    ibz_set(&(a.coord[3]),28);
    ibz_set(&(a.denom),76);
    ibz_set(&num,4);
    ibz_set(&denom,1);
    ibq_set(&cmp,&num,&denom);
    quat_alg_trace(&trace,&a);
    res = res || (ibq_cmp(&trace,&cmp));

    ibz_set(&(a.coord[0]),0);
    ibz_set(&(a.coord[1]),12);
    ibz_set(&(a.coord[2]),35);
    ibz_set(&(a.coord[3]),49);
    ibz_set(&(a.denom),28);
    ibz_set(&num,0);
    ibz_set(&denom,1);
    ibq_set(&cmp,&num,&denom);
    quat_alg_trace(&trace,&a);
    res = res || (ibq_cmp(&trace,&cmp));

    ibz_set(&(a.coord[0]),5);
    ibz_set(&(a.coord[1]),0);
    ibz_set(&(a.coord[2]),0);
    ibz_set(&(a.coord[3]),0);
    ibz_set(&(a.denom),2);
    ibz_set(&num,10);
    ibz_set(&denom,2);
    ibq_set(&cmp,&num,&denom);
    quat_alg_trace(&trace,&a);
    res = res || (ibq_cmp(&trace,&cmp));

    if (res != 0){
        printf("Quaternion unit test alg_trace failed\n");
    }
    quat_alg_elem_finalize(&a);
    ibq_finalize(&trace);
    ibq_finalize(&cmp);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    return(res);
}


//void quat_alg_scalar(quat_alg_elem_t *elem, const ibz_t *numerator, const ibz_t *denominator);
int quat_test_alg_scalar(){
    int res = 0;
    quat_alg_elem_t elem, cmp;
    ibz_t denom, num;
    quat_alg_elem_init(&cmp);
    quat_alg_elem_init(&elem);
    ibz_init(&num);
    ibz_init(&denom);

    ibz_set(&num, 1);
    ibz_set(&denom, 1);
    ibz_set(&(cmp.coord[0]), 1);
    ibz_set(&(cmp.coord[1]), 0);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), 0);
    ibz_set(&(cmp.denom), 1);
    quat_alg_scalar(&elem,&num,&denom);
    res = res || ibz_cmp(&(elem.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(elem.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&num, 5);
    ibz_set(&denom, 9);
    ibz_set(&(cmp.coord[0]), 5);
    ibz_set(&(cmp.coord[1]), 0);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), 0);
    ibz_set(&(cmp.denom), 9);
    quat_alg_scalar(&elem,&num,&denom);
    res = res || ibz_cmp(&(elem.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(elem.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&num, -125);
    ibz_set(&denom, 25);
    ibz_set(&(cmp.coord[0]), -125);
    ibz_set(&(cmp.coord[1]), 0);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), 0);
    ibz_set(&(cmp.denom), 25);
    quat_alg_scalar(&elem,&num,&denom);
    res = res || ibz_cmp(&(elem.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(elem.coord[i]),&(cmp.coord[i]));
    }
    if (res != 0){
        printf("Quaternion unit test alg_scalar failed\n");
    }
    quat_alg_elem_finalize(&elem);
    quat_alg_elem_finalize(&cmp);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    return(res);
}

//void quat_alg_conj(quat_alg_elem_t *conj, const quat_alg_elem_t *x);
int quat_test_alg_conj(){
    int res = 0;
    quat_alg_elem_t a, conj, cmp;
    quat_alg_elem_init(&cmp);
    quat_alg_elem_init(&conj);
    quat_alg_elem_init(&a);

    ibz_set(&(a.coord[0]), 0);
    ibz_set(&(a.coord[1]), 0);
    ibz_set(&(a.coord[2]), 0);
    ibz_set(&(a.coord[3]), 7);
    ibz_set(&(a.denom), 25);
    ibz_set(&(cmp.coord[0]), 0);
    ibz_set(&(cmp.coord[1]), 0);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), -7);
    ibz_set(&(cmp.denom), 25);
    quat_alg_conj(&conj,&a);
    res = res || ibz_cmp(&(conj.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(conj.coord[i]),&(cmp.coord[i]));
    }

    ibz_set(&(a.coord[0]), -125);
    ibz_set(&(a.coord[1]), 2);
    ibz_set(&(a.coord[2]), 0);
    ibz_set(&(a.coord[3]), -30);
    ibz_set(&(a.denom), 25);
    ibz_set(&(cmp.coord[0]), -125);
    ibz_set(&(cmp.coord[1]), -2);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), 30);
    ibz_set(&(cmp.denom), 25);
    quat_alg_conj(&conj,&a);
    res = res || ibz_cmp(&(conj.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(conj.coord[i]),&(cmp.coord[i]));
    }
    if (res != 0){
        printf("Quaternion unit test alg_conj failed\n");
    }

    quat_alg_elem_finalize(&a);
    quat_alg_elem_finalize(&conj);
    quat_alg_elem_finalize(&cmp);
    return(res);
}

//void quat_alg_make_primitive(quat_alg_coord_t *primitive_x, ibz_t *content, const quat_alg_elem_t *x, const quat_order_t *order, const quat_alg_t *alg);
int quat_test_alg_make_primitive(){
    int res = 0;
    quat_alg_elem_t x;
    quat_alg_t alg;
    quat_order_t order;
    quat_alg_coord_t prim, x_coord_in_order;
    ibz_t cmp_cnt, cnt;
    quat_alg_elem_init(&x);
    quat_alg_init_set_ui(&alg, 19);
    quat_alg_coord_init(&x_coord_in_order);
    quat_alg_coord_init(&prim);
    quat_order_init(&order);
    ibz_init(&cmp_cnt);
    ibz_init(&cnt);


    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(order.basis[i][j]),0);
        }
    }
    ibz_set(&(order.basis[0][0]),1);
    ibz_set(&(order.basis[0][3]),-1);
    ibz_set(&(order.basis[1][1]),-2);
    ibz_set(&(order.basis[2][2]),1);
    ibz_set(&(order.basis[2][1]),1);
    ibz_set(&(order.basis[3][3]),-3);
    ibz_set(&(order.denom),6);
    quat_lattice_hnf(&order);
    // x=1, should succeed if order
    ibz_set(&(x.denom),1);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),0);
    ibz_set(&(x.coord[2]),0);
    ibz_set(&(x.coord[3]),0);
    // is it an order?
    res = res || (0==quat_lattice_contains(&prim,&order,&x,&alg));

    // actual test
    ibz_set(&(x.denom),6);
    ibz_set(&(x.coord[0]),2);
    ibz_set(&(x.coord[1]),-4);
    ibz_set(&(x.coord[2]),26);
    ibz_set(&(x.coord[3]),18);

    res = res || (0==quat_lattice_contains(&x_coord_in_order,&order,&x,&alg));
    ibz_content(&cmp_cnt,&x_coord_in_order);
    quat_alg_make_primitive(&prim,&cnt,&x,&order,&alg);
    res = res || ibz_cmp(&cnt,&cmp_cnt);
    ibz_content(&cmp_cnt,&prim);
    res = res || !ibz_is_one(&cmp_cnt);
    //multiply by cnt, and compare to x in order
    //assumes contains is correct
    for(int i = 0; i < 4; i++){
        ibz_mul(&cmp_cnt,&cnt,&(prim[i]));
        res = res || ibz_cmp(&cmp_cnt,&(x_coord_in_order[i]));
    }

    //on a primitive element
    ibz_set(&(x.denom),6);
    ibz_set(&(x.coord[0]),2);
    ibz_set(&(x.coord[1]),-4);
    ibz_set(&(x.coord[2]),5);
    ibz_set(&(x.coord[3]),18);

    res = res || (0==quat_lattice_contains(&x_coord_in_order,&order,&x,&alg));
    ibz_content(&cmp_cnt,&x_coord_in_order);
    quat_alg_make_primitive(&prim,&cnt,&x,&order,&alg);
    res = res || ibz_cmp(&cnt,&cmp_cnt);
    ibz_content(&cmp_cnt,&prim);
    res = res || !ibz_is_one(&cmp_cnt);
    // this x is is primitive, so cnt should be 1
    res = res || !ibz_is_one(&cnt);
    //multiply by cnt, and compare to x in order
    //assumes contains is correct
    for(int i = 0; i < 4; i++){
        ibz_mul(&cmp_cnt,&cnt,&(prim[i]));
        res = res || ibz_cmp(&cmp_cnt,&(x_coord_in_order[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_make_primitive failed\n");
    }
    quat_alg_elem_finalize(&x);
    quat_alg_finalize(&alg);
    quat_alg_coord_finalize(&x_coord_in_order);
    quat_alg_coord_finalize(&prim);
    quat_order_finalize(&order);
    ibz_finalize(&cmp_cnt);
    ibz_finalize(&cnt);
    return(res);
}

//int quat_alg_is_primitive(const quat_alg_elem_t *x, const quat_order_t *order, const quat_alg_t *alg):
int quat_test_alg_is_primitive(){
    int res = 0;
    ibz_t cnt;
    quat_alg_elem_t x;
    quat_alg_t alg;
    quat_order_t order;
    quat_alg_coord_t coord;
    quat_alg_elem_init(&x);
    quat_alg_init_set_ui(&alg, 19);
    quat_order_init(&order);
    quat_alg_coord_init(&coord);
    ibz_init(&cnt);



    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            ibz_set(&(order.basis[i][j]),0);
        }
    }
    ibz_set(&(order.basis[0][0]),1);
    ibz_set(&(order.basis[0][3]),-1);
    ibz_set(&(order.basis[1][1]),-2);
    ibz_set(&(order.basis[2][2]),1);
    ibz_set(&(order.basis[2][1]),1);
    ibz_set(&(order.basis[3][3]),-3);
    ibz_set(&(order.denom),6);
    quat_lattice_hnf(&order);
    // x=1, should succeed if order
    ibz_set(&(x.denom),1);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),0);
    ibz_set(&(x.coord[2]),0);
    ibz_set(&(x.coord[3]),0);
    // is it an order?
    res = res || (0==quat_lattice_contains(&coord,&order,&x,&alg));

    // actual test, this is not primitive
    ibz_set(&(x.denom),6);
    ibz_set(&(x.coord[0]),2);
    ibz_set(&(x.coord[1]),-4);
    ibz_set(&(x.coord[2]),26);
    ibz_set(&(x.coord[3]),18);
    // x not primitive
    res = res || (0==quat_lattice_contains(&coord,&order,&x,&alg));
    ibz_content(&cnt,&coord);
    res = res || ibz_is_one(&cnt);
    res = res || quat_alg_is_primitive(&x,&order,&alg);
    
    // actual test, this is primitive
    ibz_set(&(x.denom),6);
    ibz_set(&(x.coord[0]),2);
    ibz_set(&(x.coord[1]),-4);
    ibz_set(&(x.coord[2]),5);
    ibz_set(&(x.coord[3]),18);
    // x is primitive
    res = res || (0==quat_lattice_contains(&coord,&order,&x,&alg));
    ibz_content(&cnt,&coord);
    res = res || !ibz_is_one(&cnt);
    res = res || !quat_alg_is_primitive(&x,&order,&alg);
    if (res != 0){
        printf("Quaternion unit test alg_is_primitive failed\n");
    }
    quat_alg_elem_finalize(&x);
    quat_alg_finalize(&alg);
    quat_order_finalize(&order);
    quat_alg_coord_finalize(&coord);
    ibz_finalize(&cnt);
    return(res);
}

//void quat_alg_normalize(quat_alg_elem_t *x);
int quat_test_alg_normalize(){
    int res = 0;
    quat_alg_elem_t x, cmp;
    ibz_t gcd;
    quat_alg_elem_init(&x);
    quat_alg_elem_init(&cmp);
    ibz_init(&gcd);

    // sign change
    ibz_set(&(x.coord[0]), -125);
    ibz_set(&(x.coord[1]), 2);
    ibz_set(&(x.coord[2]), 0);
    ibz_set(&(x.coord[3]), -30);
    ibz_set(&(x.denom), -25);
    ibz_set(&(cmp.coord[0]), 125);
    ibz_set(&(cmp.coord[1]), -2);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), 30);
    ibz_set(&(cmp.denom), 25);
    quat_alg_normalize(&x);
    res = res || ibz_cmp(&(x.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(x.coord[i]),&(cmp.coord[i]));
    }
    //divide by gcd
    ibz_set(&(x.coord[0]), -36);
    ibz_set(&(x.coord[1]), 18);
    ibz_set(&(x.coord[2]), 0);
    ibz_set(&(x.coord[3]), -300);
    ibz_set(&(x.denom), 48);
    ibz_set(&(cmp.coord[0]), -6);
    ibz_set(&(cmp.coord[1]), 3);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), -50);
    ibz_set(&(cmp.denom), 8);
    quat_alg_normalize(&x);
    res = res || ibz_cmp(&(x.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(x.coord[i]),&(cmp.coord[i]));
    }
    //divide by gcd
    ibz_set(&(x.coord[0]), -36);
    ibz_set(&(x.coord[1]), 18);
    ibz_set(&(x.coord[2]), 0);
    ibz_set(&(x.coord[3]), -300);
    ibz_set(&(x.denom), -6);
    ibz_set(&(cmp.coord[0]), 6);
    ibz_set(&(cmp.coord[1]), -3);
    ibz_set(&(cmp.coord[2]), 0);
    ibz_set(&(cmp.coord[3]), 50);
    ibz_set(&(cmp.denom), 1);
    quat_alg_normalize(&x);
    res = res || ibz_cmp(&(x.denom),&(cmp.denom));
    for (int i = 0; i<4; i++){;
        res = res || ibz_cmp(&(x.coord[i]),&(cmp.coord[i]));
    }

    if (res != 0){
        printf("Quaternion unit test alg_normalize failed\n");
    }
    quat_alg_elem_finalize(&x);
    quat_alg_elem_finalize(&cmp);
    ibz_finalize(&gcd);
    return(res);
}

//int quat_alg_elem_is_zero(const quat_alg_elem_t *x);
int quat_test_alg_elem_is_zero(){
    int res = 0;
    quat_alg_elem_t x;
    quat_alg_elem_init(&x);
    ibz_set(&(x.denom),1);
    ibz_set(&(x.coord[0]),0);
    ibz_set(&(x.coord[1]),0);
    ibz_set(&(x.coord[2]),0);
    ibz_set(&(x.coord[3]),0);
    res = res | (1-quat_alg_elem_is_zero(&x));
    ibz_set(&(x.denom),56865);
    res = res | (1-quat_alg_elem_is_zero(&x));
    ibz_set(&(x.denom),0);
    // maybe failure should be accepted here, but according to doc, this is still 0
    res = res | (1-quat_alg_elem_is_zero(&x));
    ibz_set(&(x.coord[3]),1);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.denom),56865);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[3]),-1);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[2]),1);
    ibz_set(&(x.coord[3]),0);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[2]),-20);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[1]),1);
    ibz_set(&(x.coord[2]),0);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[1]),-50000);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[0]),1);
    ibz_set(&(x.coord[1]),0);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[0]),-90000);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[0]),0);
    ibz_set(&(x.coord[1]),-500);
    ibz_set(&(x.coord[2]),20);
    ibz_set(&(x.coord[3]),0);
    res = res | quat_alg_elem_is_zero(&x);
    ibz_set(&(x.coord[0]),19);
    ibz_set(&(x.coord[1]),-500);
    ibz_set(&(x.coord[2]),20);
    ibz_set(&(x.coord[3]),-2);
    res = res | quat_alg_elem_is_zero(&x);
    if (res != 0){
        printf("Quaternion unit test alg_elem_is_zero failed\n");
    }
    quat_alg_elem_finalize(&x);
    return(res);
}

//int quat_alg_coord_is_zero(const quat_alg_coord_t *x);
int quat_test_alg_coord_is_zero(){
    int res = 0;
    quat_alg_coord_t x;
    quat_alg_coord_init(&x);
    ibz_set(&(x[0]),0);
    ibz_set(&(x[1]),0);
    ibz_set(&(x[2]),0);
    ibz_set(&(x[3]),0);
    res = res | (1-quat_alg_coord_is_zero(&x));
    ibz_set(&(x[3]),1);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[3]),-1);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[2]),1);
    ibz_set(&(x[3]),0);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[2]),-20);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[1]),1);
    ibz_set(&(x[2]),0);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[1]),-50000);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[0]),1);
    ibz_set(&(x[1]),0);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[0]),-90000);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[0]),0);
    ibz_set(&(x[1]),-500);
    ibz_set(&(x[2]),20);
    ibz_set(&(x[3]),0);
    res = res | quat_alg_coord_is_zero(&x);
    ibz_set(&(x[0]),19);
    ibz_set(&(x[1]),-500);
    ibz_set(&(x[2]),20);
    ibz_set(&(x[3]),-2);
    res = res | quat_alg_coord_is_zero(&x);
    if (res != 0){
        printf("Quaternion unit test alg_coord_is_zero failed\n");
    }
    quat_alg_coord_finalize(&x);
    return(res);
}


//void quat_alg_elem_set(quat_alg_elem_t *elem, const ibz_t *denom, const ibz_t *coord0,const ibz_t *coord1,const ibz_t *coord2,const ibz_t *coord3){
int quat_test_alg_elem_copy_ibz(){
    int res = 0;
    ibz_t a,b,c,d,q;
    quat_alg_elem_t elem;
    quat_alg_elem_init(&elem);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&q);
    ibz_set(&a,1);
    ibz_set(&b,2);
    ibz_set(&c,3);
    ibz_set(&d,4);
    ibz_set(&q,5);
    quat_alg_elem_copy_ibz(&elem,&q,&a,&b,&c,&d);
    res = res || ibz_cmp(&(elem.coord[0]),&a);
    res = res || ibz_cmp(&(elem.coord[1]),&b);
    res = res || ibz_cmp(&(elem.coord[2]),&c);
    res = res || ibz_cmp(&(elem.coord[3]),&d);
    res = res || ibz_cmp(&(elem.denom),&q);

    if (res != 0){
        printf("Quaternion unit test alg_elem_copy_ibz failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&q);
    quat_alg_elem_finalize(&elem);
    return(res);
}


//void quat_alg_elem_set(quat_alg_elem_t *elem, int64_t denom, int64_t coord0, int64_t coord1, int64_t coord2, int64_t coord3)
int quat_test_alg_elem_set(){
    int res = 0;
    quat_alg_elem_t elem;
    quat_alg_elem_init(&elem);
    quat_alg_elem_set(&elem,5,1,2,3,4);
    res = res || (ibz_get(&(elem.coord[0])) != 1);
    res = res || (ibz_get(&(elem.coord[1])) != 2);
    res = res || (ibz_get(&(elem.coord[2])) != 3);
    res = res || (ibz_get(&(elem.coord[3])) != 4);
    res = res || (ibz_get(&(elem.denom)) != 5);

    if (res != 0){
        printf("Quaternion unit test alg_elem_set failed\n");
    }
    quat_alg_elem_finalize(&elem);
    return(res);
}

//void quat_alg_elem_mul_by_scalar(quat_alg_elem_t *res, const ibz_t *scalar, const quat_alg_elem_t *elem);
int quat_test_alg_elem_mul_by_scalar(){
    int res = 0;
    ibz_t scalar;
    quat_alg_elem_t elem, cmp, prod;
    ibz_init(&scalar);
    quat_alg_elem_init(&elem);
    quat_alg_elem_init(&prod);
    quat_alg_elem_init(&cmp);

    ibz_set(&scalar,6);
    ibz_set(&(elem.denom),2);
    ibz_set(&(elem.coord[0]),2);
    ibz_set(&(elem.coord[1]),-4);
    ibz_set(&(elem.coord[2]),5);
    ibz_set(&(elem.coord[3]),25);
    ibz_set(&(cmp.denom),2);
    ibz_set(&(cmp.coord[0]),12);
    ibz_set(&(cmp.coord[1]),-24);
    ibz_set(&(cmp.coord[2]),30);
    ibz_set(&(cmp.coord[3]),150);
    
    quat_alg_elem_mul_by_scalar(&prod,&scalar,&elem);
    res = res || ibz_cmp(&(prod.coord[0]),&(cmp.coord[0]));
    res = res || ibz_cmp(&(prod.coord[1]),&(cmp.coord[1]));
    res = res || ibz_cmp(&(prod.coord[2]),&(cmp.coord[2]));
    res = res || ibz_cmp(&(prod.coord[3]),&(cmp.coord[3]));
    res = res || ibz_cmp(&(prod.denom),&(cmp.denom));
    // denom should not be modified
    res = res || ibz_cmp(&(prod.denom),&(elem.denom));

    ibz_set(&scalar,-3);
    ibz_set(&(cmp.denom),2);
    ibz_set(&(cmp.coord[0]),-6);
    ibz_set(&(cmp.coord[1]),12);
    ibz_set(&(cmp.coord[2]),-15);
    ibz_set(&(cmp.coord[3]),-75);
    
    quat_alg_elem_mul_by_scalar(&prod,&scalar,&elem);
    res = res || ibz_cmp(&(prod.coord[0]),&(cmp.coord[0]));
    res = res || ibz_cmp(&(prod.coord[1]),&(cmp.coord[1]));
    res = res || ibz_cmp(&(prod.coord[2]),&(cmp.coord[2]));
    res = res || ibz_cmp(&(prod.coord[3]),&(cmp.coord[3]));
    res = res || ibz_cmp(&(prod.denom),&(cmp.denom));

    if (res != 0){
        printf("Quaternion unit test alg_elem_mul_by_scalar failed\n");
    }
    ibz_finalize(&scalar);
    quat_alg_elem_finalize(&elem);
    quat_alg_elem_finalize(&prod);
    quat_alg_elem_finalize(&cmp);
    return(res);
}


//void quat_alg_rightmul_mat(ibz_mat_4x4_t *mulmat, const quat_alg_elem_t *a, const quat_alg_t *alg);
int quat_test_alg_rightmul_mat(){
    int res = 0;
    ibz_mat_4x4_t mat;
    quat_alg_t alg;
    quat_alg_elem_t elem, prod, test, prodmat;
    quat_alg_init_set_ui(&alg,19);
    ibz_mat_4x4_init(&mat);
    quat_alg_elem_init(&prod);
    quat_alg_elem_init(&prodmat);
    quat_alg_elem_init(&test);
    quat_alg_elem_init(&elem);

    quat_alg_elem_set(&elem, 2,4,2,6,8);
    quat_alg_elem_set(&elem, 12,9,81,27,-6);
    quat_alg_rightmul_mat(&mat,&elem,&alg);
    ibz_mat_4x4_eval(&(prodmat.coord),&mat,&(test.coord));
    ibz_mul(&(prodmat.denom),&(test.denom),&(elem.denom));
    quat_alg_mul(&prod,&test,&elem,&alg);
    quat_alg_sub(&test,&prod,&prod);
    res = res || !quat_alg_elem_is_zero(&test);

    if (res != 0){
        printf("Quaternion unit test alg_rightmul_mat failed\n");
    }
    ibz_mat_4x4_finalize(&mat);
    quat_alg_elem_finalize(&prod);
    quat_alg_elem_finalize(&prodmat);
    quat_alg_elem_finalize(&test);
    quat_alg_elem_finalize(&elem);
    quat_alg_finalize(&alg);
    return(res);
}

// run all previous tests
int quat_test_algebra(){
    int res = 0;
    printf("\nRunning quaternion tests of algebra operations\n");
    res = res | quat_test_init_set_ui();
    res = res | quat_test_alg_coord_add();
    res = res | quat_test_alg_coord_sub();
    res = res | quat_test_alg_equal_denom();
    res = res | quat_test_alg_add();
    res = res | quat_test_alg_sub();
    res = res | quat_test_alg_mul();
    res = res | quat_test_alg_norm();
    res = res | quat_test_alg_trace();
    res = res | quat_test_alg_scalar();
    res = res | quat_test_alg_conj();
    res = res | quat_test_alg_make_primitive();
    res = res | quat_test_alg_is_primitive();
    res = res | quat_test_alg_normalize();
    res = res | quat_test_alg_elem_is_zero();
    res = res | quat_test_alg_coord_is_zero();
    res = res | quat_test_alg_elem_copy_ibz();
    res = res | quat_test_alg_elem_set();
    res = res | quat_test_alg_elem_mul_by_scalar();
    res = res | quat_test_alg_rightmul_mat();
    return(res);
}
