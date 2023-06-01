#include <quaternion.h>
#include "internal.h"

// RED(k,l) sub-algorithm
static void RED(ibz_mat_4x4_t *basis, ibq_t (*u)[4][4], ibz_t (*H)[4][4], int k, int l) {
    ibq_t tmp, tmp2;
    ibz_t q, tibz, num, den, r;
    ibq_init(&tmp);
    ibq_init(&tmp2);
    ibz_init(&q);
    ibz_init(&tibz);
    ibz_init(&num);
    ibz_init(&den);
    ibz_init(&r);
    ibz_set(&num,1);
    ibz_set(&den,2);
    ibq_set(&tmp,&num,&den);
    ibz_set(&num,0);
    ibz_set(&den,0);

    // if |u_{k,l}| <= 0.5, terminate
    ibq_abs(&tmp2, &((*u)[k][l]));
    if (ibq_cmp(&tmp2, &tmp) <= 0)
        goto end;

    // q <- floor(0.5 + u_{k,l})
    ibq_add(&tmp, &tmp, &((*u)[k][l]));

    ibq_num(&num, &tmp);
    ibq_denom(&den, &tmp);
    //FDIV was used, needs reeimplementation
    ibz_div_floor(&q,&r, &num, &den);
    //ibq_floor(tmp, tmp);
    //ibz_set_f(q, tmp);

    // b_k = b_k - q*b_l
    for (int i = 0; i < 4; ++i) {
        ibz_mul(&tibz, &q, &((*basis)[l][i]));
        ibz_sub(&((*basis)[k][i]), &((*basis)[k][i]), &tibz);
    }

    // H_k = H_k - q*H_l
    for (int i = 0; i < 4; ++i) {
        ibz_mul(&tibz, &q, &((*H)[l][i]));
        ibz_sub(&((*H)[k][i]), &((*H)[k][i]), &tibz);
    }

    // u_{k,j} = u_{k,l}-q
    ibq_set(&tmp2, &q,&ibz_const_one);
    ibq_sub(&((*u)[k][l]), &((*u)[k][l]), &tmp2);

    // forall_i \in 1..l-1: u_{k,i} = u_{k,i} - q*u_{l,i}
    for (int i = 0; i <= l-1; ++i) {
        ibq_mul(&tmp, &tmp2, &((*u)[l][i]));
        ibq_sub(&((*u)[k][i]), &((*u)[k][i]), &tmp);
    }

end:
    ibq_finalize(&tmp);
    ibq_finalize(&tmp2);
    ibz_finalize(&q);
    ibz_finalize(&tibz);
    ibz_finalize(&num);
    ibz_finalize(&den);
    ibz_finalize(&r);
}

// SWAP(k) sub-algorithm
static void SWAP(ibz_mat_4x4_t *basis, ibq_t (*u)[4][4], ibz_t (*H)[4][4], ibq_t (*B)[4], ibq_t (*bStar)[4][4], int k, int kmax) {
    ibq_t tmp, tmp2, tmp3, u_tmp, B_tmp, b[4];
    ibq_init(&tmp);
    ibq_init(&tmp2);
    ibq_init(&tmp3);
    ibq_init(&u_tmp);
    ibq_init(&B_tmp);

    for (int i = 0; i < 4; ++i) {
        ibq_init(&(b[i]));
    }

    // swap b_k and b_{k-1}
    for (int i = 0; i < 4; ++i) {
        ibz_swap(&((*basis)[k][i]), &((*basis)[k-1][i]));
    }

    // swap H_k and H_{k-1}
    for (int i = 0; i < 4; ++i) {
        ibz_swap(&((*H)[k][i]), &((*H)[k-1][i]));
    }

    if (k > 1) {
        // swap u_{k,j} and u_{k-1,j}
        for (int j = 0; j <= k - 2; ++j) {
            ibq_swap(&((*u)[k][j]), &((*u)[k-1][j]));
        }
    }

    // u = u_{k,k-1}
    ibq_copy(&u_tmp, &((*u)[k][k - 1]));

    // B = B_k + u^2*B_{k-1}
    ibq_mul(&B_tmp, &u_tmp, &u_tmp);
    ibq_mul(&B_tmp, &B_tmp, &((*B)[k-1]));
    ibq_add(&B_tmp, &((*B)[k]), &B_tmp);

    // u_{k,k-1} = u*B_{k-1} / B
    ibq_mul(&tmp, &u_tmp, &((*B)[k-1]));
    ibq_div(&((*u)[k][k-1]), &tmp, &B_tmp);

    // b = bSTAR_{k-1}
    for (int i = 0; i < 4; ++i) {
        ibq_copy(&(b[i]), &((*bStar)[k-1][i]));
    }
    // bSTAR_{k-1}=bSTAR_k+u*b
    for (int i = 0; i < 4; ++i) {
        ibq_mul(&tmp, &u_tmp, &(b[i]));
        ibq_add(&((*bStar)[k-1][i]), &((*bStar)[k][i]), &tmp);
    }
    // bSTAR_k = -u_{k,k-1}*bSTAR_k+(B_k/B)*b
    ibq_div(&tmp2, &((*B)[k]), &B_tmp); // B_k/B
    ibq_neg(&tmp, &((*u)[k][k-1]));
    for (int i = 0; i < 4; ++i) {
        ibq_mul(&((*bStar)[k][i]), &tmp, &((*bStar)[k][i]));
        ibq_mul(&tmp3, &tmp2, &(b[i]));
        ibq_add(&((*bStar)[k][i]), &((*bStar)[k][i]), &tmp3);
    }

    // B_k = B_{k-1}*B_k/B
    ibq_mul(&((*B)[k]), &((*B)[k-1]), &((*B)[k]));
    ibq_div(&((*B)[k]), &((*B)[k]), &B_tmp);

    // B_{k-1} = B
    ibq_copy(&((*B)[k-1]), &B_tmp);

    for (int i = k+1; i <= kmax; ++i) {
        // t = u_{i,k}
        ibq_copy(&tmp, &((*u)[i][k]));

        // u_{i,k} = u_{i,k-1} - u*t
        ibq_mul(&((*u)[i][k]), &u_tmp, &tmp);
        ibq_sub(&((*u)[i][k]), &((*u)[i][k-1]), &((*u)[i][k]));

        // u_{i,k-1} = t + u_{k,k-1}*u_{i,k}
        ibq_mul(&tmp2, &((*u)[k][k-1]), &((*u)[i][k]));
        ibq_add(&((*u)[i][k-1]), &tmp, &tmp2);
    }

    ibq_finalize(&tmp);
    ibq_finalize(&tmp2);
    ibq_finalize(&tmp3);
    ibq_finalize(&u_tmp);
    ibq_finalize(&B_tmp);
    for (int i = 0; i < 4; ++i) {
        ibq_finalize(&(b[i]));
    }

}

// m1[0]*m2[0] + m1[1]*m2[1] + q*(m1[2]*m2[2] + m1[3]*m2[3])
static void dotproduct_row(ibz_t *mul, const ibz_mat_4x4_t *m1, const ibz_mat_4x4_t *m2, const ibz_t *q, int m1j, int m2j) {
    ibz_set(mul, 0);
    ibz_t tmp1, tmp2;
    ibz_init(&tmp1);
    ibz_init(&tmp2);
    for (int i = 0; i < 2; ++i) {
        ibz_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
        ibz_add(mul, mul, &tmp1);
    }
    for (int i = 2; i < 4; ++i) {
        ibz_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
        ibz_add(&tmp2, &tmp2, &tmp1);
    }
    ibz_mul(&tmp2, &tmp2, q);
    ibz_add(mul, mul, &tmp2);

    ibz_finalize(&tmp1);
    ibz_finalize(&tmp2);
}

static void dotproduct_zr_row(ibq_t *mul, const ibz_mat_4x4_t *m1, const ibq_t (*m2)[4][4], const ibz_t *q, int m1j, int m2j) {
    ibq_set(mul, &ibz_const_zero, &ibz_const_one);
    ibq_t tmp1, tmp2;
    ibq_init(&tmp1);
    ibq_init(&tmp2);
    for (int i = 0; i < 2; ++i) {
        ibq_set(&tmp1, &((*m1)[m1j][i]), &ibz_const_one);
        ibq_mul(&tmp1, &tmp1, &((*m2)[m2j][i]));
        ibq_add(mul, mul, &tmp1); 
    }
    for (int i = 2; i < 4; ++i) {
        ibq_set(&tmp1, &((*m1)[m1j][i]), &ibz_const_one);
        ibq_mul(&tmp1, &tmp1, &((*m2)[m2j][i]));
        ibq_add(&tmp2, &tmp2, &tmp1); 
    }
    ibq_set(&tmp1, q, &ibz_const_one);
    ibq_mul(&tmp2, &tmp2, &tmp1);
    ibq_add(mul, mul, &tmp2);

    ibq_finalize(&tmp1);
    ibq_finalize(&tmp2);
}

static void dotproduct_rr_row(ibq_t *mul, const ibq_t (*m1)[4][4], const ibq_t (*m2)[4][4], const ibz_t *q, int m1j, int m2j) {
    //ibq_set(mul, 0);
    ibq_set(mul, &ibz_const_zero, &ibz_const_one);
    ibq_t tmp1, tmp2;
    ibq_init(&tmp1);
    ibq_init(&tmp2);
    for (int i = 0; i < 2; ++i) {
        ibq_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
        ibq_add(mul, mul, &tmp1);
    }
    for (int i = 2; i < 4; ++i) {
        ibq_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
        ibq_add(&tmp2, &tmp2, &tmp1);
    }
    ibq_set(&tmp1, q, &ibz_const_one);
    ibq_mul(&tmp2, &tmp2, &tmp1);
    ibq_add(mul, mul, &tmp2);

    ibq_finalize(&tmp1);
    ibq_finalize(&tmp2);
}

static void mul_row(ibq_t (*mul)[4][4], const ibq_t *a, const ibq_t (*m)[4][4], int j) {
    for (int i = 0; i < 4; ++i) {
        ibq_mul(&((*mul)[j][i]), a, &((*m)[j][i]));
    }
}

static void add_row(ibz_mat_4x4_t *add, const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b, int j, int aj, int bj) {
    for (int i = 0; i < 4; ++i) {
        ibz_add(&((*add)[j][i]), &((*a)[aj][i]), &((*b)[bj][i]));
    }
}

static void sub_row(ibq_t (*add)[4][4], const ibq_t (*a)[4][4], const ibq_t (*b)[4][4], int j, int aj, int bj) {
    for (int i = 0; i < 4; ++i) {
        ibq_sub(&((*add)[j][i]), &((*a)[aj][i]), &((*b)[bj][i]));
    }
}

/// @brief LLL reduction on 4-dimensional lattice
/// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number Theory"
/// @param red 
/// @param lattice 
/// @return 
int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int precision) {
    (void) precision;
    int ret = 0;
    ibz_mat_4x4_t basis;
    ibq_t bStar[4][4];
    ibq_t bStar_tmp[4][4];
    ibq_t tmp;
    ibz_t tmp_z;
    ibz_t den;
    ibz_t num;
    ibq_t cnst;
    ibq_t u[4][4];
    ibz_t H[4][4]; // -> I_4
    ibq_t B[4];
    ibq_init(&tmp);
    ibz_init(&tmp_z);
    ibz_init(&den);
    ibz_init(&num);
    ibq_init(&cnst);
    for (int i = 0; i < 4; ++i)
        ibq_init(&(B[i]));

    ibz_mat_4x4_init(&basis);
    ibz_mat_4x4_transpose(&basis, &(lattice->basis));

    // Step 1: Initialize: ...
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            ibq_init(&(u[i][j]));
            ibq_init(&(bStar[i][j]));
            ibq_init(&(bStar_tmp[i][j]));
            // bSTAR_1 = b_1 (we copy all)
            if (i == j){
                ibz_init(&(H[i][j]));
                ibz_set(&(H[i][j]), 1);
            }
            else {
                ibz_init(&(H[i][j]));
            }
        }
    }
    int k = 1, kmax = 0;
    // bStar_1 = b_1
    for (int i = 0; i < 4; ++i)
        ibq_set(&(bStar[0][i]), &(basis[0][i]), &ibz_const_one);
    // B_1 = b_1 * b_1
    dotproduct_row(&tmp_z, &basis, &basis, q, 0, 0);
    ibq_set(&(B[0]), &tmp_z, &ibz_const_one);
    ibz_set(&num,99);
    ibz_set(&den,100);
    ibq_set(&cnst,&num,&den);

    while (k < 4) {
        // Step 2: Incremental Gram-Schmidt
        // if (k <= kmax) -> we can omit..
        if (k > kmax) {
            kmax = k;
            for (int i = 0; i < 4; ++i) {
                ibq_set(&(bStar[k][i]), &(basis[k][i]), &ibz_const_one);
            }
            for (int j = 0; j <= k-1; ++j) {
                // bStar_k = b_k -> already done initially
                // nop
                // u_{k,j} = b_k*bSTAR_j/B_j
                dotproduct_zr_row(&tmp, &basis, &bStar, q, k, j);
                ibq_div(&(u[k][j]), &tmp, &(B[j]));
                // bStar_k = bStar_k - u_{k,j}*bStar_j
                mul_row(&bStar_tmp, &(u[k][j]), &bStar, j);
                sub_row(&bStar, &bStar, &bStar_tmp, k, k, j);
            }
            // B_k = bStar_k*bStar_k
            dotproduct_rr_row(&(B[k]), &bStar, &bStar, q, k, k);
            if (ibq_is_zero(&(B[k]))) {
                // b_i did not form a basis, terminate with error
                ret = -1;
                goto err;
            }
        }

        while(1) {
            // Step 3: Test LLL condition
            RED(&basis, &u, &H, k, k - 1);
            // If B_k < (0.75 - u_{k,k-1}^2)*B_{k-1}
            ibq_mul(&tmp, &(u[k][k-1]), &(u[k][k-1]));
            ibq_sub(&tmp, &cnst, &tmp);
            ibq_mul(&tmp, &tmp, &(B[k-1]));
            if (ibq_cmp(&(B[k]), &tmp) < 0) {
                SWAP(&basis, &u, &H, &B, &bStar, k, kmax);
                k = (k - 1 > 1 ? k - 1 : 1);
            } else {
                for (int l = k - 2; l >= 0; --l) {
                    RED(&basis, &u, &H, k, l);
                }
                k++;
                break;
            }
        }
    }
    ibz_mat_4x4_transpose(red, &basis);

err:
    ibq_finalize(&tmp);
    ibz_finalize(&tmp_z);
    ibz_finalize(&num);
    ibz_finalize(&den);
    ibq_finalize(&cnst);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            ibq_finalize(&(u[i][j]));
            ibz_finalize(&(H[i][j]));
            ibq_finalize(&(bStar[i][j]));
            ibq_finalize(&(bStar_tmp[i][j]));
        }
    }
    for (int i = 0; i < 4; ++i)
        ibq_finalize(&(B[i]));
    ibz_mat_4x4_finalize(&basis);
    return ret;
}
