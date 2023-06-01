#include "quaternion_tests.h"

// void ibz_vec_2_set(ibz_vec_2_t *vec, int a0, int a1);
int quat_test_dim2_ibz_vec_2_set()
{
    int res = 0;
    ibz_vec_2_t vec;
    ibz_vec_2_init(&vec);
    ibz_vec_2_set(&vec, 2, 5);
    res = res || !(ibz_get(&(vec[0])) == 2);
    res = res || !(ibz_get(&(vec[1])) == 5);
    if (res != 0)
    {
        printf("Quaternion unit test dim2_ibz_vec_2_set failed\n");
    }
    ibz_vec_2_finalize(&vec);
    return (res);
}

// void ibz_mat_2x2_set(ibz_mat_2x2_t *mat, int a00, int a01, int a10, int a11);
int quat_test_dim2_ibz_mat_2x2_set()
{
    int res = 0;
    ibz_mat_2x2_t mat;
    ibz_mat_2x2_init(&mat);
    ibz_mat_2x2_set(&mat, 2, 7, -1, 5);
    res = res || !(ibz_get(&(mat[0][0])) == 2);
    res = res || !(ibz_get(&(mat[0][1])) == 7);
    res = res || !(ibz_get(&(mat[1][0])) == -1);
    res = res || !(ibz_get(&(mat[1][1])) == 5);
    if (res != 0)
    {
        printf("Quaternion unit test dim2_ibz_mat_2x2_set failed\n");
    }
    ibz_mat_2x2_finalize(&mat);
    return (res);
}

//void ibz_mat_2x2_det_from_ibz(ibz_t *det, const ibz_t *a11, const ibz_t *a12, const ibz_t *a21, const ibz_t *a22);
int quat_test_dim2_ibz_mat_2x2_det_from_ibz(){
    int res = 0;
    ibz_t det,cmp,a11,a12,a21,a22;
    ibz_init(&a11);
    ibz_init(&a12);
    ibz_init(&a21);
    ibz_init(&a22);
    ibz_init(&det);
    ibz_init(&cmp);
    
    ibz_set(&a11,1);
    ibz_set(&a12,0);
    ibz_set(&a21,0);
    ibz_set(&a22,1);
    ibz_set(&cmp,1);
    ibz_mat_2x2_det_from_ibz(&det,&a11,&a12,&a21,&a22);
    res = res || ibz_cmp(&cmp,&det);

    ibz_set(&a11,2);
    ibz_set(&a12,3);
    ibz_set(&a21,1);
    ibz_set(&a22,-2);
    ibz_set(&cmp,-7);
    ibz_mat_2x2_det_from_ibz(&det,&a11,&a12,&a21,&a22);
    res = res || ibz_cmp(&cmp,&det);

    ibz_set(&a11,0);
    ibz_set(&a12,3);
    ibz_set(&a21,-1);
    ibz_set(&a22,0);
    ibz_set(&cmp,3);
    ibz_mat_2x2_det_from_ibz(&det,&a11,&a12,&a21,&a22);
    res = res || ibz_cmp(&cmp,&det);

    ibz_set(&a11,2);
    ibz_set(&cmp,0);
    ibz_mat_2x2_det_from_ibz(&a11,&a11,&a11,&a11,&a11);
    res = res || ibz_cmp(&cmp,&a11);


    if (res != 0){
        printf("Quaternion unit test dim2_ibz_mat_2x2_det_from_ibz failed\n");
    }
    ibz_finalize(&a11);
    ibz_finalize(&a12);
    ibz_finalize(&a21);
    ibz_finalize(&a22);
    ibz_finalize(&cmp);
    ibz_finalize(&det);
    return res;
}

// void ibz_mat_2x2_eval(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec);
int quat_test_dim2_ibz_mat_2x2_eval()
{
    int res = 0;
    ibz_vec_2_t vec, ret, cmp;
    ibz_mat_2x2_t mat;
    ibz_vec_2_init(&vec);
    ibz_mat_2x2_init(&mat);
    ibz_vec_2_init(&ret);
    ibz_vec_2_init(&cmp);

    ibz_mat_2x2_set(&mat, 1, -1, 2, 4);
    ibz_vec_2_set(&vec, 1, -1);
    ibz_vec_2_set(&cmp, 2, -2);
    ibz_mat_2x2_eval(&ret, &mat, &vec);
    res = res || ibz_cmp(&(ret[0]), &(cmp[0]));
    res = res || ibz_cmp(&(ret[1]), &(cmp[1]));

    ibz_mat_2x2_set(&mat, 2, -2, 1, 3);
    ibz_vec_2_set(&vec, 2, 4);
    ibz_vec_2_set(&cmp, -4, 14);
    ibz_mat_2x2_eval(&vec, &mat, &vec);
    res = res || ibz_cmp(&(vec[0]), &(cmp[0]));
    res = res || ibz_cmp(&(vec[1]), &(cmp[1]));

    if (res != 0)
    {
        printf("Quaternion unit test dim2_ibz_mat_2x2_eval failed\n");
    }
    ibz_vec_2_finalize(&vec);
    ibz_mat_2x2_finalize(&mat);
    ibz_vec_2_finalize(&ret);
    ibz_vec_2_finalize(&cmp);
    return (res);
}

// modular 2x2 operations

// void ibz_2x2_mul_mod(ibz_mat_2x2_t *prod, const ibz_mat_2x2_t *mat_a, const ibz_mat_2x2_t *mat_b, const ibz_t *m);
int quat_test_dim2_ibz_mat_2x2_mul_mod()
{
    int res = 0;
    ibz_t m;
    ibz_mat_2x2_t a, b, cmp, prod;
    ibz_init(&m);
    ibz_mat_2x2_init(&a);
    ibz_mat_2x2_init(&b);
    ibz_mat_2x2_init(&prod);
    ibz_mat_2x2_init(&cmp);

    ibz_set(&m, 7);
    ibz_mat_2x2_set(&a, 2, -2, 1, 3);
    ibz_mat_2x2_set(&b, 5, 3, 4, 1);
    ibz_mat_2x2_set(&cmp, 2, 4, 3, 6);
    ibz_2x2_mul_mod(&prod, &a, &b, &m);
    res = res || ibz_cmp(&(prod[0][0]), &(cmp[0][0]));
    res = res || ibz_cmp(&(prod[0][1]), &(cmp[0][1]));
    res = res || ibz_cmp(&(prod[1][0]), &(cmp[1][0]));
    res = res || ibz_cmp(&(prod[1][1]), &(cmp[1][1]));
    ibz_mat_2x2_set(&cmp, 6, 6, 2, 2);
    ibz_2x2_mul_mod(&prod, &b, &a, &m);
    res = res || ibz_cmp(&(prod[0][0]), &(cmp[0][0]));
    res = res || ibz_cmp(&(prod[0][1]), &(cmp[0][1]));
    res = res || ibz_cmp(&(prod[1][0]), &(cmp[1][0]));
    res = res || ibz_cmp(&(prod[1][1]), &(cmp[1][1]));

    ibz_set(&m, 12);
    ibz_mat_2x2_set(&a, 2, 7, 1, -2);
    ibz_mat_2x2_set(&cmp, 11, 0, 0, 11);
    ibz_2x2_mul_mod(&a, &a, &a, &m);
    res = res || ibz_cmp(&(a[0][0]), &(cmp[0][0]));
    res = res || ibz_cmp(&(a[0][1]), &(cmp[0][1]));
    res = res || ibz_cmp(&(a[1][0]), &(cmp[1][0]));
    res = res || ibz_cmp(&(a[1][1]), &(cmp[1][1]));

    if (res != 0)
    {
        printf("Quaternion unit test dim2_ibz_mat_2x2_mul_mod failed\n");
    }
    ibz_mat_2x2_finalize(&a);
    ibz_mat_2x2_finalize(&b);
    ibz_mat_2x2_finalize(&prod);
    ibz_mat_2x2_finalize(&cmp);
    ibz_finalize(&m);
    return (res);
}

// int ibz_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m);
int quat_test_dim2_ibz_mat_2x2_inv_mod()
{
    int res = 0;
    ibz_t m;
    ibz_mat_2x2_t a, inv, id, prod;
    ibz_init(&m);
    ibz_mat_2x2_init(&a);
    ibz_mat_2x2_init(&inv);
    ibz_mat_2x2_init(&prod);
    ibz_mat_2x2_init(&id);
    ibz_mat_2x2_set(&id, 1, 0, 0, 1);

    // inverse exists
    ibz_set(&m, 7);
    ibz_mat_2x2_set(&a, 2, -3, 1, 3);
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
    ibz_set(&m, 12);
    ibz_mat_2x2_set(&a, 2, 7, 1, -2);
    ibz_mat_2x2_set(&inv, 2, 7, 1, -2);
    if (ibz_2x2_inv_mod(&inv, &inv, &m))
    {
        ibz_2x2_mul_mod(&prod, &a, &inv, &m);
        res = res || ibz_cmp(&(prod[0][0]), &(id[0][0]));
        res = res || ibz_cmp(&(prod[0][1]), &(id[0][1]));
        res = res || ibz_cmp(&(prod[1][0]), &(id[1][0]));
        res = res || ibz_cmp(&(prod[1][1]), &(id[1][1]));
    }
    else
    {
        res = 1;
    }

    // no inverse
    ibz_set(&m, 25);
    ibz_mat_2x2_set(&a, 2, -2, -1, 1);
    res = res || ibz_2x2_inv_mod(&inv, &a, &m);
    ibz_set(&m, 7);
    ibz_mat_2x2_set(&a, 2, 3, 1, -2);
    res = res || ibz_2x2_inv_mod(&inv, &a, &m);
    ibz_set(&m, 25);
    ibz_mat_2x2_set(&a, 2, 1, 1, -2);
    res = res || ibz_2x2_inv_mod(&inv, &a, &m);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_ibz_mat_2x2_inv_mod failed\n");
    }
    ibz_mat_2x2_finalize(&a);
    ibz_mat_2x2_finalize(&inv);
    ibz_mat_2x2_finalize(&prod);
    ibz_mat_2x2_finalize(&id);
    ibz_finalize(&m);
    return (res);
}


// cvp helper functions

// int quat_dim2_lattice_contains(ibz_mat_2x2_t *basis, ibz_t *coord1, ibz_t *coord2);
int quat_test_dim2_lattice_contains()
{
    int res = 0;
    ibz_mat_2x2_t basis;
    ibz_t c1, c2;
    ibz_mat_2x2_init(&basis);
    ibz_init(&c1);
    ibz_init(&c2);
    ibz_mat_2x2_set(&basis, -1, 2, 5, 0);
    ibz_set(&c1, 0);
    ibz_set(&c2, 0);
    res = res || !quat_dim2_lattice_contains(&basis, &c1, &c2);
    ibz_set(&c1, 1);
    ibz_set(&c2, 5);
    res = res || !quat_dim2_lattice_contains(&basis, &c1, &c2);
    ibz_set(&c1, 1);
    ibz_set(&c2, 4);
    res = res || quat_dim2_lattice_contains(&basis, &c1, &c2);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_contains failed\n");
    }
    ibz_mat_2x2_finalize(&basis);
    ibz_finalize(&c1);
    ibz_finalize(&c2);
    return (res);
}

// void quat_dim2_lattice_norm(ibz_t *norm, const ibz_t *coord1, const ibz_t *coord2, const ibz_t *norm_q)
int quat_test_dim2_lattice_norm()
{
    int res = 0;
    ibz_t norm, cmp, q, a, b;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&q);
    ibz_init(&norm);
    ibz_init(&cmp);
    ibz_set(&a, 1);
    ibz_set(&b, 2);
    ibz_set(&q, 3);
    ibz_set(&cmp, 1 * 1 + 3 * 2 * 2);
    quat_dim2_lattice_norm(&norm, &a, &b, &q);
    res = res || ibz_cmp(&norm, &cmp);
    ibz_set(&a, 7);
    ibz_set(&b, -2);
    ibz_set(&q, 5);
    ibz_set(&cmp, 7 * 7 + 5 * 2 * 2);
    quat_dim2_lattice_norm(&norm, &a, &b, &q);
    res = res || ibz_cmp(&norm, &cmp);
    ibz_set(&cmp, 7 * 7 + 7 * 7 * 7);
    quat_dim2_lattice_norm(&a, &a, &a, &a);
    res = res || ibz_cmp(&a, &cmp);
    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_norm failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&q);
    ibz_finalize(&norm);
    ibz_finalize(&cmp);
    return (res);
}

// void quat_dim2_lattice_bilinear(ibz_t *res, const ibz_t *v11, const ibz_t *v12,const ibz_t *v21, const ibz_t *v22, const ibz_t *norm_q);
int quat_test_dim2_lattice_bilinear()
{
    int res = 0;
    ibz_t prod, cmp, q, a, b, c, d;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&d);
    ibz_init(&q);
    ibz_init(&prod);
    ibz_init(&cmp);
    ibz_set(&a, 1);
    ibz_set(&b, 2);
    ibz_set(&c, 3);
    ibz_set(&d, 4);
    ibz_set(&q, 3);
    ibz_set(&cmp, 1 * 3 + 3 * 2 * 4);
    quat_dim2_lattice_bilinear(&prod, &a, &b, &c, &d, &q);
    res = res || ibz_cmp(&prod, &cmp);
    ibz_set(&a, 7);
    ibz_set(&b, -2);
    ibz_set(&c, 4);
    ibz_set(&d, 7);
    ibz_set(&q, 5);
    ibz_set(&cmp, 7 * 4 - 5 * 2 * 7);
    quat_dim2_lattice_bilinear(&prod, &a, &b, &c, &d, &q);
    res = res || ibz_cmp(&prod, &cmp);
    ibz_set(&cmp, 7 * 7 + 7 * 7 * 7);
    quat_dim2_lattice_bilinear(&a, &a, &a, &a, &a, &a);
    res = res || ibz_cmp(&a, &cmp);
    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_bilinear failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&d);
    ibz_finalize(&q);
    ibz_finalize(&prod);
    ibz_finalize(&cmp);
    return (res);
}

// according to conditions from https://cseweb.ucsd.edu/classes/wi12/cse206A-a/lec3.pdf
int quat_dim2_lattice_verify_lll(const ibz_mat_2x2_t *mat, const ibz_t *norm_q){
    ibz_t prod, norm_a, norm_bstar, norm_b, b_ab, b_abstar, b_bbstar,p01_denom,p00_denom, norm_p00,norm_p01;
    ibz_vec_2_t p00,p01,bstar;
    ibz_vec_2_init(&p00);
    ibz_vec_2_init(&p01);
    ibz_vec_2_init(&bstar);
    ibz_init(&prod);
    ibz_init(&norm_a);
    ibz_init(&norm_bstar);
    ibz_init(&norm_b);
    ibz_init(&b_ab);
    ibz_init(&b_bbstar);
    ibz_init(&b_abstar);
    ibz_init(&p00_denom);
    ibz_init(&p01_denom);
    ibz_init(&norm_p00);
    ibz_init(&norm_p01);
    quat_dim2_lattice_norm(&norm_a,&((*mat)[0][0]),&((*mat)[1][0]),norm_q);
    quat_dim2_lattice_norm(&norm_b,&((*mat)[0][1]),&((*mat)[1][1]),norm_q);
    quat_dim2_lattice_bilinear(&b_ab,&((*mat)[0][0]),&((*mat)[1][0]),&((*mat)[0][1]),&((*mat)[1][1]),norm_q);
    ibz_mul(&(bstar[0]),&((*mat)[0][1]),&norm_b);
    ibz_mul(&prod,&((*mat)[0][0]),&b_ab);
    ibz_sub(&(bstar[0]),&(bstar[0]),&prod);
    ibz_mul(&(bstar[1]),&((*mat)[1][1]),&norm_b);
    ibz_mul(&prod,&((*mat)[1][0]),&b_ab);
    ibz_sub(&(bstar[1]),&(bstar[1]),&prod);
    quat_dim2_lattice_norm(&norm_bstar,&(bstar[0]),&(bstar[1]),norm_q);
    quat_dim2_lattice_bilinear(&b_bbstar,&(bstar[0]),&(bstar[1]),&((*mat)[0][1]),&((*mat)[1][1]),norm_q);
    quat_dim2_lattice_bilinear(&b_abstar,&((*mat)[0][0]),&((*mat)[1][0]),&(bstar[0]),&(bstar[1]),norm_q);
    // first projection
    ibz_mul(&(p00[0]),&((*mat)[0][0]),&norm_bstar);
    ibz_mul(&prod,&(bstar[0]),&b_abstar);
    ibz_add(&(p00[0]),&(p00[0]),&prod);
    ibz_mul(&(p00[1]),&((*mat)[1][0]),&norm_bstar);
    ibz_mul(&prod,&(bstar[1]),&b_abstar);
    ibz_add(&(p00[1]),&(p00[1]),&prod);
    ibz_copy(&p00_denom,&norm_bstar);
    // second projection
    ibz_mul(&(p01[0]),&((*mat)[0][0]),&b_ab);
    ibz_mul(&(p01[0]),&(p01[0]),&norm_bstar);
    ibz_mul(&prod,&(bstar[0]),&b_bbstar);
    ibz_mul(&prod,&prod,&norm_a);
    ibz_add(&(p01[0]),&(p01[0]),&prod);
    ibz_mul(&(p01[1]),&((*mat)[1][0]),&b_ab);
    ibz_mul(&(p01[1]),&(p01[1]),&norm_bstar);
    ibz_mul(&prod,&(bstar[1]),&b_bbstar);
    ibz_mul(&prod,&prod,&norm_a);
    ibz_add(&(p01[1]),&(p01[1]),&prod);
    ibz_mul(&p01_denom,&norm_a,&norm_bstar);
    // compute norms
    quat_dim2_lattice_norm(&norm_p00,&(p00[0]),&(p00[1]),norm_q);
    quat_dim2_lattice_norm(&norm_p01,&(p01[0]),&(p01[1]),norm_q);
    // compare on same denom
    ibz_mul(&norm_p00,&norm_p00,&norm_a);
    ibz_mul(&norm_p00,&norm_p00,&norm_a);
    ibz_mul(&p00_denom,&p00_denom,&norm_a);
    int res = (ibz_cmp(&norm_p00,&norm_p01)<=0);
    // also check 2times bilinear against norm
    ibz_set(&prod,2);
    ibz_mul(&prod,&b_ab,&prod);
    res = res && (ibz_cmp(&prod,&norm_a)< 0);

    ibz_finalize(&prod);
    ibz_finalize(&norm_a);
    ibz_finalize(&norm_bstar);
    ibz_finalize(&norm_b);
    ibz_finalize(&norm_p00);
    ibz_finalize(&norm_p01);
    ibz_finalize(&b_ab);
    ibz_finalize(&b_abstar);
    ibz_finalize(&b_bbstar);
    ibz_finalize(&p00_denom);
    ibz_finalize(&p01_denom);
    ibz_vec_2_finalize(&p00);
    ibz_vec_2_finalize(&p01);
    ibz_vec_2_finalize(&bstar);
    return(res);
}

// void quat_dim2_lattice_short_basis(ibz_mat_2x2_t *reduced, const ibz_mat_2x2_t *basis, const ibz_t *norm_q);
int quat_test_dim2_lattice_short_basis()
{
    int res = 0;
    ibz_mat_2x2_t basis, cmp, red;
    ibz_t prod, sum, norm_q, bound;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_init(&norm_q);
    ibz_init(&bound);
    ibz_mat_2x2_init(&basis);
    ibz_mat_2x2_init(&cmp);
    ibz_mat_2x2_init(&red);
    // first test
    ibz_mat_2x2_set(&basis, 3, 5, 7, 9);
    ibz_set(&norm_q, 2);
    quat_dim2_lattice_short_basis(&red, &basis, &norm_q);
    // check second basis vector larger (or at least equal) than 1st
    res = res || (ibz_cmp(&sum, &prod) > 0);
    // check mutual inclusion of lattices
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][0]), &(basis[1][0])));
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][1]), &(basis[1][1])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][0]), &(red[1][0])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][1]), &(red[1][1])));
    // check bilinear form value is small
    ibz_set(&bound, 12);
    quat_dim2_lattice_bilinear(&sum, &(red[0][0]), &(red[0][1]), &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    // check norm smaller than original
    ibz_set(&bound, 12);
    quat_dim2_lattice_norm(&sum, &(red[0][0]), &(red[1][0]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    quat_dim2_lattice_norm(&prod, &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&prod, &bound) > 0);
    // check lll
    res = res || !quat_dim2_lattice_verify_lll(&red,&norm_q);

    // second test
    ibz_mat_2x2_set(&basis, 48, 4, 81, 9);
    ibz_set(&norm_q, 9);
    quat_dim2_lattice_short_basis(&red, &basis, &norm_q);
    // check second basis vector larger (or at least equal) than 1st
    res = res || (ibz_cmp(&sum, &prod) > 0);
    // check mutual inclusion of lattices
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][0]), &(basis[1][0])));
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][1]), &(basis[1][1])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][0]), &(red[1][0])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][1]), &(red[1][1])));
    // check bilinear form value is small
    ibz_set(&bound, 50 * 100);
    quat_dim2_lattice_bilinear(&sum, &(red[0][0]), &(red[0][1]), &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    // check norm smaller than original
    ibz_set(&bound, 50 * 50);
    quat_dim2_lattice_norm(&sum, &(red[0][0]), &(red[1][0]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    quat_dim2_lattice_norm(&prod, &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&prod, &bound) > 0);
    // check lll
    res = res || !quat_dim2_lattice_verify_lll(&red,&norm_q);

    // third test
    ibz_mat_2x2_set(&basis, -1, 2, 3, 2);
    ibz_set(&norm_q, 2);
    quat_dim2_lattice_short_basis(&red, &basis, &norm_q);
    // check second basis vector larger (or at least equal) than 1st
    res = res || (ibz_cmp(&sum, &prod) > 0);
    // check mutual inclusion of lattices
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][0]), &(basis[1][0])));
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][1]), &(basis[1][1])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][0]), &(red[1][0])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][1]), &(red[1][1])));
    // check bilinear form value is small
    ibz_set(&bound, 12);
    quat_dim2_lattice_bilinear(&sum, &(red[0][0]), &(red[0][1]), &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    // check norm smaller than original
    ibz_set(&bound, 12);
    quat_dim2_lattice_norm(&sum, &(red[0][0]), &(red[1][0]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    quat_dim2_lattice_norm(&prod, &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&prod, &bound) > 0);
    // check lll
    res = res || !quat_dim2_lattice_verify_lll(&red,&norm_q);

    //4th test
    ibz_set(&norm_q, 12165);
    ibz_mat_2x2_set(&basis, 364, 0, 1323546, 266606);
    quat_dim2_lattice_short_basis(&red, &basis,&norm_q);
    // check second basis vector larger (or at least equal) than 1st
    res = res || (ibz_cmp(&sum, &prod) > 0);
    // check mutual inclusion of lattices
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][0]), &(basis[1][0])));
    res = res || (!quat_dim2_lattice_contains(&red, &(basis[0][1]), &(basis[1][1])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][0]), &(red[1][0])));
    res = res || (!quat_dim2_lattice_contains(&basis, &(red[0][1]), &(red[1][1])));
    // check bilinear form value is small
    quat_dim2_lattice_norm(&bound, &(basis[0][1]), &(basis[1][1]),&norm_q);
    quat_dim2_lattice_bilinear(&sum, &(red[0][0]), &(red[0][1]), &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    // check norm smaller than original
    quat_dim2_lattice_norm(&bound, &(basis[0][1]), &(basis[1][1]),&norm_q);
    quat_dim2_lattice_norm(&sum, &(red[0][0]), &(red[1][0]), &norm_q);
    res = res || (ibz_cmp(&sum, &bound) > 0);
    quat_dim2_lattice_norm(&prod, &(red[0][1]), &(red[1][1]), &norm_q);
    res = res || (ibz_cmp(&prod, &bound) > 0);
    // check lll
    res = res || !quat_dim2_lattice_verify_lll(&red,&norm_q);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_short_basis failed\n");
    }
    ibz_finalize(&prod);
    ibz_finalize(&sum);
    ibz_finalize(&norm_q);
    ibz_finalize(&bound);
    ibz_mat_2x2_finalize(&basis);
    ibz_mat_2x2_finalize(&cmp);
    ibz_mat_2x2_finalize(&red);
    return res;
}

//void quat_dim2_lattice_get_coefficient_with_orthogonalisation(ibz_t *res, const ibz_t *a0,const ibz_t *a1,const ibz_t *b0,const ibz_t *b1,const ibz_t *t0,const ibz_t *t1,const ibz_t *norm_q); 
int quat_test_dim2_lattice_get_coefficient_with_orthogonalisation(){
    int res = 0;
    ibz_t out, q, cmp;
    ibz_vec_2_t a, b, t;
    ibz_init(&out);
    ibz_init(&q);
    ibz_vec_2_init(&a);
    ibz_vec_2_init(&b);
    ibz_vec_2_init(&t);
    ibz_init(&cmp);

    ibz_vec_2_set(&b,364,203);
    ibz_vec_2_set(&a,136,-606);
    ibz_vec_2_set(&t,4312, 1122);
    ibz_set(&q,1);
    ibz_set(&cmp,2);
    quat_dim2_lattice_get_coefficient_with_orthogonalisation(&out,&(a[0]),&(a[1]),&(b[0]),&(b[1]),&(t[0]),&(t[1]),&q);
    res = res || ibz_cmp(&cmp,&out);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_get_coefficient_with_orthogonalisation failed\n");
    }
    ibz_finalize(&out);
    ibz_finalize(&q);
    ibz_vec_2_finalize(&a);
    ibz_vec_2_finalize(&b);
    ibz_vec_2_finalize(&t);
    ibz_finalize(&cmp);
    return(res);
}

// void quat_dim2_lattice_closest_vector(ibz_vec_2_t *target_minus_closest, ibz_vec_2_t *closest_coords_in_basis, const ibz_mat_2x2_t *reduced_basis, const ibz_vec_2_t *target, const ibz_t *norm_q);
// void quat_dim2_lattice_closest_vector(ibz_vec_2_t *target_minus_closest, ibz_vec_2_t *closest_coords_in_basis, const ibz_mat_2x2_t *reduced_basis, const ibz_vec_2_t *target, const ibz_t *norm_q);
int quat_test_dim2_lattice_closest_vector()
{
    int res = 0;
    ibz_vec_2_t target, coords, diff, cmp_c, cmp_d;
    ibz_t q;
    ibz_mat_2x2_t basis;
    ibz_mat_2x2_init(&basis);
    ibz_init(&q);
    ibz_vec_2_init(&diff);
    ibz_vec_2_init(&target);
    ibz_vec_2_init(&coords);
    ibz_vec_2_init(&cmp_c);
    ibz_vec_2_init(&cmp_d);

    ibz_vec_2_set(&target, 1216765, 879777);
    ibz_set(&q, 12165);
    ibz_mat_2x2_set(&basis, 364, 0, 1323546, 266606);
    quat_dim2_lattice_short_basis(&basis, &basis,&q);
    res = res || !quat_dim2_lattice_verify_lll(&basis,&q);

    quat_dim2_lattice_norm(&(cmp_d[0]),&(basis[0][0]),&(basis[1][0]),&q);
    quat_dim2_lattice_norm(&(cmp_d[1]),&(basis[0][1]),&(basis[1][1]),&q);
    quat_dim2_lattice_closest_vector(&diff, &coords, &basis, &target, &q);
    quat_dim2_lattice_norm(&(cmp_c[0]),&(diff[0]),&(diff[1]),&q);
    // test if diff has short norm
    res = res || (ibz_cmp(&(cmp_d[0]), &(cmp_c[0]))<0);
    res = res || (ibz_cmp(&(cmp_d[1]), &(cmp_c[0]))<0);// test coord in std basis + diff gives target
    ibz_mul(&(cmp_c[0]), &(basis[0][0]), &(coords[0]));
    ibz_mul(&(cmp_d[0]), &(basis[0][1]), &(coords[1]));
    ibz_add(&(cmp_c[0]), &(cmp_c[0]), &(cmp_d[0]));
    ibz_mul(&(cmp_c[1]), &(basis[1][0]), &(coords[0]));
    ibz_mul(&(cmp_d[1]), &(basis[1][1]), &(coords[1]));
    ibz_add(&(cmp_c[1]), &(cmp_c[1]), &(cmp_d[1]));
    ibz_add(&(cmp_d[1]), &(diff[1]), &(cmp_c[1]));
    ibz_add(&(cmp_d[0]), &(diff[0]), &(cmp_c[0]));
    res = res || ibz_cmp(&(target[0]), &(cmp_d[0]));
    res = res || ibz_cmp(&(target[1]), &(cmp_d[1]));

    ibz_set_from_str(&(target[0]),"24331603791763515691729329496192502951056898161849272968170597636916968011043891967931576902713229532952543949167870260",10);
    ibz_set_from_str(&(target[1]),"320864169455876864410246920205472260746024223953596490297261861428291460634916982031632479909031611896628831462307645903",10);
    ibz_set_from_str(&(basis[0][0]),"83823929950128301594915602697228677079260482500022964929273",10);
    ibz_set_from_str(&(basis[0][1]),"0",10);
    ibz_set_from_str(&(basis[1][0]),"26601400243014239200180507899382755112580169588401175904342542265512938583814884922817611265843666255091702318534544478",10);
    ibz_set_from_str(&(basis[1][1]),"702645123228401649030937397460831326065409163240113467483219315690011711473572222698591407715178198721852976513892308529",10);
    ibz_set(&q, 1);
    quat_dim2_lattice_short_basis(&basis, &basis,&q);
    res = res || !quat_dim2_lattice_verify_lll(&basis,&q);

    quat_dim2_lattice_norm(&(cmp_d[0]),&(basis[0][0]),&(basis[1][0]),&q);
    quat_dim2_lattice_norm(&(cmp_d[1]),&(basis[0][1]),&(basis[1][1]),&q);
    quat_dim2_lattice_closest_vector(&diff, &coords, &basis, &target, &q);
    quat_dim2_lattice_norm(&(cmp_c[0]),&(diff[0]),&(diff[1]),&q);
    quat_dim2_lattice_bilinear(&(cmp_c[1]),&(basis[0][0]),&(basis[1][0]),&(basis[0][1]),&(basis[1][1]),&q);
    quat_dim2_lattice_bilinear(&(cmp_c[1]),&(basis[0][0]),&(basis[1][0]),&(diff[0]),&(diff[1]),&q);
    quat_dim2_lattice_bilinear(&(cmp_c[1]),&(basis[0][1]),&(basis[1][1]),&(diff[0]),&(diff[1]),&q);
    // test if diff has short norm
    res = res || (ibz_cmp(&(cmp_d[0]), &(cmp_c[0]))<0);
    res = res || (ibz_cmp(&(cmp_d[1]), &(cmp_c[0]))<0);
    // test coord in std basis + diff gives target
    ibz_mul(&(cmp_c[0]), &(basis[0][0]), &(coords[0]));
    ibz_mul(&(cmp_d[0]), &(basis[0][1]), &(coords[1]));
    ibz_add(&(cmp_c[0]), &(cmp_c[0]), &(cmp_d[0]));
    ibz_mul(&(cmp_c[1]), &(basis[1][0]), &(coords[0]));
    ibz_mul(&(cmp_d[1]), &(basis[1][1]), &(coords[1]));
    ibz_add(&(cmp_c[1]), &(cmp_c[1]), &(cmp_d[1]));
    ibz_add(&(cmp_d[1]), &(diff[1]), &(cmp_c[1]));
    ibz_add(&(cmp_d[0]), &(diff[0]), &(cmp_c[0]));
    res = res || ibz_cmp(&(target[0]), &(cmp_d[0]));
    res = res || ibz_cmp(&(target[1]), &(cmp_d[1]));

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_closest_vector failed\n");
    }
    ibz_mat_2x2_finalize(&basis);
    ibz_finalize(&q);
    ibz_vec_2_finalize(&diff);
    ibz_vec_2_finalize(&target);
    ibz_vec_2_finalize(&coords);
    ibz_vec_2_finalize(&cmp_c);
    ibz_vec_2_finalize(&cmp_d);
    return (res);
}

// void quat_dim2_lattice_get_qf_on_lattice(ibz_t *qf_a, ibz_t *qf_b,ibz_t *qf_c, const ibz_mat_2x2_t *basis ,const ibz_t *norm_q);
int quat_test_dim2_lattice_get_qf_on_lattice()
{
    int res = 0;
    ibz_t a, b, c, q, cmp_a, cmp_b, cmp_c;
    ibz_mat_2x2_t basis;
    ibz_mat_2x2_init(&basis);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&c);
    ibz_init(&q);
    ibz_init(&cmp_a);
    ibz_init(&cmp_b);
    ibz_init(&cmp_c);
    ibz_mat_2x2_set(&basis, -3, 2, 1, 2);
    ibz_set(&q, 2);
    ibz_set(&cmp_a, 11);
    ibz_set(&cmp_b, -4);
    ibz_set(&cmp_c, 12);
    quat_dim2_lattice_get_qf_on_lattice(&a, &b, &c, &basis, &q);
    res = res || ibz_cmp(&a, &cmp_a);
    res = res || ibz_cmp(&b, &cmp_b);
    res = res || ibz_cmp(&c, &cmp_c);
    ibz_mat_2x2_set(&basis, 5, 7, -1, 1);
    ibz_set(&q, 20);
    ibz_set(&cmp_a, 5 * 5 + 20);
    ibz_set(&cmp_b, 2 * (5 * 7 - 20));
    ibz_set(&cmp_c, 7 * 7 + 20);
    quat_dim2_lattice_get_qf_on_lattice(&q, &b, &c, &basis, &q);
    res = res || ibz_cmp(&q, &cmp_a);
    res = res || ibz_cmp(&b, &cmp_b);
    res = res || ibz_cmp(&c, &cmp_c);
    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_get_qf_on_lattice failed\n");
    }
    ibz_mat_2x2_finalize(&basis);
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&c);
    ibz_finalize(&q);
    ibz_finalize(&cmp_a);
    ibz_finalize(&cmp_b);
    ibz_finalize(&cmp_c);
    return (res);
}

// int quat_dim2_lattice_test_cvp_condition(quat_alg_elem_t* elem, const ibz_vec_2_t* vec, const void* params);
int quat_test_dim2_lattice_test_cvp_condition()
{
    int res = 0;
    ibz_vec_2_t vec;
    quat_alg_elem_t elem;
    ibz_t p;
    void *params = (void *)&p;
    quat_alg_elem_init(&elem);
    ibz_vec_2_init(&vec);
    ibz_init(&p);
    ibz_vec_2_set(&vec, 5, 7);
    ibz_set(&p, 5);
    res = res || !quat_dim2_lattice_test_cvp_condition(&elem, &vec, params);
    res = res || ibz_cmp(&(elem.coord[0]), &(vec[0]));
    res = res || ibz_cmp(&(elem.coord[1]), &(vec[1]));
    res = res || ibz_cmp(&(elem.coord[2]), &(vec[0]));
    res = res || ibz_cmp(&(elem.coord[3]), &(vec[1]));
    ibz_vec_2_set(&vec, 5, 8);
    ibz_set(&p, 5);
    res = res || quat_dim2_lattice_test_cvp_condition(&elem, &vec, params);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_test_cvp_condition failed\n");
    }
    quat_alg_elem_finalize(&elem);
    ibz_vec_2_finalize(&vec);
    ibz_finalize(&p);
    return (res);
}

// int quat_dim2_lattice_bound_and_condition(quat_alg_elem_t *res, const ibz_t *x, const ibz_t *y, int (*condition)(quat_alg_elem_t* , const ibz_vec_2_t*, const void*), const void* params, const ibz_vec_2_t* target_minus_closest, const ibz_mat_2x2_t *lat_basis, const ibz_t *norm_q, const ibz_t *norm_bound);
int quat_test_dim2_lattice_bound_and_condition()
{
    int res = 0;
    quat_alg_elem_t out;
    ibz_t x, y, q, bound, p;
    ibz_vec_2_t cmp, diff;
    ibz_mat_2x2_t basis;
    quat_alg_elem_init(&out);
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&q);
    ibz_init(&p);
    ibz_init(&bound);
    ibz_vec_2_init(&cmp);
    ibz_vec_2_init(&diff);
    ibz_mat_2x2_init(&basis);
    void *params = (void *)&p;

    // should pass
    ibz_mat_2x2_set(&basis, -3, 2, 1, 2);
    ibz_set(&q, 2);
    ibz_set(&x, -1);
    ibz_set(&y, 3);
    ibz_set(&p, 3);
    ibz_set(&bound, 1000);
    ibz_vec_2_set(&diff, 1, 0);
    // out contains target - closest - short = diff - short
    ibz_vec_2_set(&cmp, -8, -5);
    if (quat_dim2_lattice_bound_and_condition(&out, &x, &y, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound))
    {
        res = res || ibz_cmp(&(cmp[0]), &(out.coord[0]));
        res = res || ibz_cmp(&(cmp[1]), &(out.coord[1]));
    }
    else
    {
        res = 1;
    }
    // fail due to bound
    ibz_set(&bound, 10);
    res = res || quat_dim2_lattice_bound_and_condition(&out, &x, &y, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound);

    // fail due to condition
    ibz_set(&bound, 1000);
    ibz_set(&p, 2);
    res = res || quat_dim2_lattice_bound_and_condition(&out, &x, &y, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_bound_and_condition failed\n");
    }
    quat_alg_elem_finalize(&out);
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&q);
    ibz_finalize(&p);
    ibz_finalize(&bound);
    ibz_vec_2_finalize(&cmp);
    ibz_vec_2_finalize(&diff);
    ibz_mat_2x2_finalize(&basis);
    return (res);
}

// int quat_dim2_lattice_qf_value_bound_generation(ibz_t *res, const ibz_t *num_a, const ibz_t *denom_a, const ibz_t *num_b, const ibz_t *denom_b){
int quat_test_dim2_lattice_qf_value_bound_generation()
{
    int res = 0;
    ibz_t sqrt_num, sqrt_denom, num_b, denom_b, cmp, bound;
    ibz_init(&sqrt_denom);
    ibz_init(&sqrt_num);
    ibz_init(&denom_b);
    ibz_init(&num_b);
    ibz_init(&cmp);
    ibz_init(&bound);

    ibz_set(&sqrt_denom, 4);
    ibz_set(&sqrt_num, 10);
    ibz_set(&denom_b, 2);
    ibz_set(&num_b, 5);
    ibz_set(&cmp, 5);
    quat_dim2_lattice_qf_value_bound_generation(&bound, &sqrt_num, &sqrt_denom, &num_b, &denom_b);
    res = res || ibz_cmp(&bound, &cmp);

    ibz_set(&sqrt_denom, 6);
    ibz_set(&sqrt_num, 9);
    ibz_set(&denom_b, 2);
    ibz_set(&num_b, 3);
    ibz_set(&cmp, 4);
    quat_dim2_lattice_qf_value_bound_generation(&bound, &sqrt_num, &sqrt_denom, &num_b, &denom_b);
    res = res || ibz_cmp(&bound, &cmp);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_qf_value_bound_generation failed\n");
    }
    ibz_finalize(&sqrt_denom);
    ibz_finalize(&sqrt_num);
    ibz_finalize(&denom_b);
    ibz_finalize(&num_b);
    ibz_finalize(&cmp);
    ibz_finalize(&bound);
    return res;
}

// int quat_dim2_lattice_qf_enumerate_short_vec(quat_alg_elem_t *res,int (*condition)(quat_alg_elem_t*, const ibz_vec_2_t*, const void*), const void *params, const ibz_vec_2_t *target_minus_closest, const ibz_mat_2x2_t *lat_basis, const ibz_t *norm_q,const ibz_t *norm_bound, const int max_tries);
int quat_test_dim2_lattice_qf_enumerate_short_vec()
{
    int res = 0;
    quat_alg_elem_t out;
    ibz_t q, bound, p;
    ibz_vec_2_t diff;
    ibz_mat_2x2_t basis;
    int max_tries = 1000;
    quat_alg_elem_init(&out);
    void *params = (void *)&p;
    ibz_init(&q);
    ibz_init(&p);
    ibz_init(&bound);
    ibz_vec_2_init(&diff);
    ibz_mat_2x2_init(&basis);
    ibz_set(&p, 3);

    // should pass
    ibz_mat_2x2_set(&basis, -3, 2, 1, 2);
    ibz_set(&q, 2);
    ibz_set(&p, 3);
    ibz_set(&bound, 200);
    max_tries = 20000;
    ibz_vec_2_set(&diff, 2, 1);
    if (quat_dim2_lattice_qf_enumerate_short_vec(&out, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound, max_tries))
    {
        res = res || quat_dim2_lattice_bound_and_condition(&out, &(out.coord[0]), &(out.coord[1]), &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound);
    }
    else
    {
        res = 1;
    }

    // should passibz_vec_2_set(&target, 1216765, 879777);
    ibz_mat_2x2_set(&basis, 364, 0, 1323546, 266606);
    ibz_set(&q, 12165);
    ibz_set(&p, 3);
    ibz_set(&bound, 1024*64);
    ibz_mul(&bound,&bound,&bound);
    max_tries = 100;
    ibz_vec_2_set(&diff, -18287, -155);
    if (quat_dim2_lattice_qf_enumerate_short_vec(&out, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound, max_tries))
    {
        res = res || quat_dim2_lattice_bound_and_condition(&out, &(out.coord[0]), &(out.coord[1]), &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound);
    }
    else
    {
        res = 1;
    }

    // should fail
    ibz_mat_2x2_set(&basis, -3, 2, 1, 2);
    ibz_set(&q, 2);
    ibz_set(&p, 2);
    ibz_set(&bound, 200);
    max_tries = 2000;
    ibz_vec_2_set(&diff, 2, 1);
    res = res || quat_dim2_lattice_qf_enumerate_short_vec(&out, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound, max_tries);

    // should fail
    ibz_mat_2x2_set(&basis, -3, 2, 1, 2);
    ibz_set(&q, 2);
    ibz_set(&p, 2);
    ibz_set(&bound, 50);
    max_tries = 2000;
    ibz_vec_2_set(&diff, 2, 1);
    res = res || quat_dim2_lattice_qf_enumerate_short_vec(&out, &quat_dim2_lattice_test_cvp_condition, params, &diff, &basis, &q, &bound, max_tries);

    if (res != 0)
    {
        printf("Quaternion unit test dim2_lattice_qf_enumerate_short_vec failed\n");
    }
    quat_alg_elem_finalize(&out);
    ibz_finalize(&q);
    ibz_finalize(&p);
    ibz_finalize(&bound);
    ibz_vec_2_finalize(&diff);
    ibz_mat_2x2_finalize(&basis);
    return (res);
}

// int quat_2x2_lattice_enumerate_cvp_filter(quat_alg_elem_t *res, const ibz_mat_2x2_t *lat_basis, const ibz_vec_2_t *target,unsigned int qf, unsigned int dist_bound, int (*condition)(quat_alg_elem_t* , const ibz_vec_2_t*, const void*), const void* params, unsigned int max_tries);
int quat_test_dim2_2x2_lattice_enumerate_cvp_filter()
{
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

    // solution should exist
    ibz_mat_2x2_set(&basis, -3, 2, 1, 2);
    ibz_vec_2_set(&target, -3, 3);
    ibz_set(&p, 3);
    q = 2;
    dist_bound = 10;
    max_tries = 20000;
    if (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries))
    {
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
    else
    {
        res = 1;
    }

    // solution should exist
    ibz_mat_2x2_set(&basis, -1, 2, 3, 2);
    ibz_vec_2_set(&target, -3, 3);
    ibz_set(&p, 3);
    q = 2;
    dist_bound = 10;
    max_tries = 20000;
    if (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries))
    {
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
    else
    {
        res = 1;
    }

    // solution should not exist because of lower dist_bound
    ibz_mat_2x2_set(&basis, -1, 2, 1, 2);
    ibz_vec_2_set(&target, -3, 3);
    ibz_set(&p, 5);
    q = 3;
    dist_bound = 1;
    max_tries = 1000;
    res = res || (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries));

    // solution should not exist because of impossible condition
    ibz_mat_2x2_set(&basis, -1, 2, 1, 2);
    ibz_vec_2_set(&target, -3, 3);
    ibz_set(&p, 2);
    q = 3;
    dist_bound = 20;
    max_tries = 1000;
    res = res || (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries));

    // should succeed, but close to bound
    ibz_vec_2_set(&target, 1216765, 879777);
    q = 12165;
    dist_bound = 32;
    ibz_set(&p, 3);
    max_tries = 100;
    ibz_mat_2x2_set(&basis, 364, 0, 1323546, 266606);
    if (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries))
    {
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
    else
    {
        res = 1;
    }
    
    ibz_vec_2_set(&target, 1216765, -223409);
    q = 12;
    dist_bound = 25;
    ibz_set(&p, 3);
    max_tries = 50;
    ibz_mat_2x2_set(&basis, 364, 2, -132346, 26606);
    if (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries))
    {
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
    else
    {
        res = 1;
    }
    
    //large passing test
    ibz_set_from_str(&(target[0]),"24331603791763515691729329496192502951056898161849272968170597636916968011043891967931576902713229532952543949167870260",10);
    ibz_set_from_str(&(target[1]),"320864169455876864410246920205472260746024223953596490297261861428291460634916982031632479909031611896628831462307645903",10);
    ibz_set_from_str(&(basis[0][0]),"83823929950128301594915602697228677079260482500022964929273",10);
    ibz_set_from_str(&(basis[0][1]),"0",10);
    ibz_set_from_str(&(basis[1][0]),"26601400243014239200180507899382755112580169588401175904342542265512938583814884922817611265843666255091702318534544478",10);
    ibz_set_from_str(&(basis[1][1]),"702645123228401649030937397460831326065409163240113467483219315690011711473572222698591407715178198721852976513892308529",10);
    q =  1;
    dist_bound = 600;
    ibz_set(&p, 3);
    max_tries = 1000;
    if (quat_2x2_lattice_enumerate_cvp_filter(&cvp_res, &basis, &target, q, dist_bound, &quat_dim2_lattice_test_cvp_condition, params, max_tries))
    {
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
    else
    {
        res = 1;
    }

    if (res != 0)
    {
        printf("Quaternion unit test dim2_2x2_lattice_enumerate_cvp_filter failed\n");
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

// run all previous tests
int quat_test_dim2()
{
    int res = 0;
    printf("\nRunning quaternion tests of functions for matrices, vectors and lattices in dimension 2\n");
    res = res | quat_test_dim2_ibz_vec_2_set();
    res = res | quat_test_dim2_ibz_mat_2x2_set();
    res = res | quat_test_dim2_ibz_mat_2x2_det_from_ibz();
    res = res | quat_test_dim2_ibz_mat_2x2_eval();
    res = res | quat_test_dim2_ibz_mat_2x2_mul_mod();
    res = res | quat_test_dim2_ibz_mat_2x2_inv_mod();
    res = res | quat_test_dim2_lattice_test_cvp_condition();
    res = res | quat_test_dim2_lattice_contains();
    res = res | quat_test_dim2_lattice_norm();
    res = res | quat_test_dim2_lattice_bilinear();
    res = res | quat_test_dim2_lattice_short_basis();
    res = res | quat_test_dim2_lattice_get_coefficient_with_orthogonalisation();
    res = res | quat_test_dim2_lattice_closest_vector();
    res = res | quat_test_dim2_lattice_get_qf_on_lattice();
    res = res | quat_test_dim2_lattice_bound_and_condition();
    res = res | quat_test_dim2_lattice_qf_value_bound_generation();
    res = res | quat_test_dim2_lattice_qf_enumerate_short_vec();
    res = res | quat_test_dim2_2x2_lattice_enumerate_cvp_filter();
    return (res);
}
