#include <quaternion.h>
#include "internal.h"
#include <stdio.h>

// internal helpers, also for other files
void
ibz_vec_2_set(ibz_vec_2_t *vec, int a0, int a1)
{
    ibz_set(&((*vec)[0]), a0);
    ibz_set(&((*vec)[1]), a1);
}
void
ibz_mat_2x2_set(ibz_mat_2x2_t *mat, int a00, int a01, int a10, int a11)
{
    ibz_set(&((*mat)[0][0]), a00);
    ibz_set(&((*mat)[0][1]), a01);
    ibz_set(&((*mat)[1][0]), a10);
    ibz_set(&((*mat)[1][1]), a11);
}

void
ibz_mat_2x2_copy(ibz_mat_2x2_t *copy, const ibz_mat_2x2_t *copied)
{
    ibz_copy(&((*copy)[0][0]), &((*copied)[0][0]));
    ibz_copy(&((*copy)[0][1]), &((*copied)[0][1]));
    ibz_copy(&((*copy)[1][0]), &((*copied)[1][0]));
    ibz_copy(&((*copy)[1][1]), &((*copied)[1][1]));
}

void
ibz_mat_2x2_add(ibz_mat_2x2_t *sum, const ibz_mat_2x2_t *a, const ibz_mat_2x2_t *b)
{
    ////FILE *fp = fopen("/mnt/c/Users/WonKim/the-sqisign/quaternion_elements_values_lvl1.txt", "a");

    ibz_add(&((*sum)[0][0]), &((*a)[0][0]), &((*b)[0][0]));
    //gmp_fprintf(fp, "%Zd\n", *(sum)[0][0]);
    ibz_add(&((*sum)[0][1]), &((*a)[0][1]), &((*b)[0][1]));
    //gmp_fprintf(fp, "%Zd\n", (*sum)[0][1]);
    ibz_add(&((*sum)[1][0]), &((*a)[1][0]), &((*b)[1][0]));
    //gmp_fprintf(fp, "%Zd\n", (*sum)[1][0]);
    ibz_add(&((*sum)[1][1]), &((*a)[1][1]), &((*b)[1][1]));
    ////gmp_fprintf(fp, "%Zd\n", (*sum[1][1]));

    //fclose(fp);
}

void
ibz_mat_2x2_det_from_ibz(ibz_t *det, const ibz_t *a11, const ibz_t *a12, const ibz_t *a21, const ibz_t *a22)
{
    //FILE *fp = fopen("/mnt/c/Users/WonKim/the-sqisign/quaternion_elements_values_lvl1.txt", "a");

    ibz_t prod;
    ibz_init(&prod);
    ibz_mul(&prod, a12, a21);
    //gmp_fprintf(fp, "%Zd\n", prod);
    ibz_mul(det, a11, a22);
    //gmp_fprintf(fp, "%Zd\n", *det);
    ibz_sub(det, det, &prod);
    //gmp_fprintf(fp, "%Zd\n", *det);
    ibz_finalize(&prod);

    //fclose(fp);
}

void
ibz_mat_2x2_eval(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec)
{
    //FILE *fp = fopen("/mnt/c/Users/WonKim/the-sqisign/quaternion_elements_values_lvl1.txt", "a");

    ibz_t prod;
    ibz_vec_2_t matvec;
    ibz_init(&prod);
    ibz_vec_2_init(&matvec);
    ibz_mul(&prod, &((*mat)[0][0]), &((*vec)[0]));
    //gmp_fprintf(fp, "%Zd\n", prod);
    ibz_copy(&(matvec[0]), &prod);
    ibz_mul(&prod, &((*mat)[0][1]), &((*vec)[1]));
    //gmp_fprintf(fp, "%Zd\n", prod);
    ibz_add(&(matvec[0]), &(matvec[0]), &prod);
    //gmp_fprintf(fp, "%Zd\n", matvec[0]);
    ibz_mul(&prod, &((*mat)[1][0]), &((*vec)[0]));
    //gmp_fprintf(fp, "%Zd\n", prod);
    ibz_copy(&(matvec[1]), &prod);
    ibz_mul(&prod, &((*mat)[1][1]), &((*vec)[1]));
    //gmp_fprintf(fp, "%Zd\n", prod);
    ibz_add(&(matvec[1]), &(matvec[1]), &prod);
    //gmp_fprintf(fp, "%Zd\n", matvec[1]);
    ibz_copy(&((*res)[0]), &(matvec[0]));
    ibz_copy(&((*res)[1]), &(matvec[1]));
    ibz_finalize(&prod);
    ibz_vec_2_finalize(&matvec);

    //fclose(fp);
}

// modular 2x2 operations

void
ibz_2x2_mul_mod(ibz_mat_2x2_t *prod, const ibz_mat_2x2_t *mat_a, const ibz_mat_2x2_t *mat_b, const ibz_t *m)
{
    //FILE *fp = fopen("/mnt/c/Users/WonKim/the-sqisign/quaternion_elements_values_lvl1.txt", "a");

    ibz_t mul;
    ibz_mat_2x2_t sums;
    ibz_init(&mul);
    ibz_mat_2x2_init(&sums);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_set(&(sums[i][j]), 0);
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                ibz_mul(&mul, &((*mat_a)[i][k]), &((*mat_b)[k][j]));
                //gmp_fprintf(fp, "%Zd\n", mul);
                ibz_add(&(sums[i][j]), &(sums[i][j]), &mul);
                //gmp_fprintf(fp, "%Zd\n", sums[i][j]);
                ibz_mod(&(sums[i][j]), &(sums[i][j]), m);
            }
        }
    }
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_copy(&((*prod)[i][j]), &(sums[i][j]));
        }
    }
    ibz_finalize(&mul);
    ibz_mat_2x2_finalize(&sums);

    //fclose(fp);
}

int
ibz_mat_2x2_inv_mod(ibz_mat_2x2_t *inv, const ibz_mat_2x2_t *mat, const ibz_t *m)
{
    //FILE *fp = fopen("/mnt/c/Users/WonKim/the-sqisign/quaternion_elements_values_lvl1.txt", "a");

    ibz_t det, prod;
    ibz_init(&det);
    ibz_init(&prod);
    ibz_mul(&det, &((*mat)[0][0]), &((*mat)[1][1]));
    //gmp_fprintf(fp, "%Zd\n", det);
    ibz_mod(&det, &det, m);
    ibz_mul(&prod, &((*mat)[0][1]), &((*mat)[1][0]));
    //gmp_fprintf(fp, "%Zd\n", prod);
    ibz_sub(&det, &det, &prod);
    //gmp_fprintf(fp, "%Zd\n", det);
    ibz_mod(&det, &det, m);
    int res = ibz_invmod(&det, &det, m);
    // return 0 matrix if non invertible determinant
    ibz_set(&prod, res);
    ibz_mul(&det, &det, &prod);
    //gmp_fprintf(fp, "%Zd\n", det);
    // compute inverse
    ibz_copy(&prod, &((*mat)[0][0]));
    ibz_copy(&((*inv)[0][0]), &((*mat)[1][1]));
    ibz_copy(&((*inv)[1][1]), &prod);
    ibz_neg(&((*inv)[1][0]), &((*mat)[1][0]));
    ibz_neg(&((*inv)[0][1]), &((*mat)[0][1]));
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            ibz_mul(&((*inv)[i][j]), &((*inv)[i][j]), &det);
            //gmp_fprintf(fp, "%Zd\n", (*inv)[i][j]);
            ibz_mod(&((*inv)[i][j]), &((*inv)[i][j]), m);
        }
    }
    ibz_finalize(&det);
    ibz_finalize(&prod);
    return (res);

    //fclose(fp);
}
