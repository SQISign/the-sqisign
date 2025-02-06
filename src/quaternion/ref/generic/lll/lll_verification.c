#include "lll_internals.h"

// functions to verify lll
void
quat_lll_set_ibq_parameters(ibq_t *delta, ibq_t *eta)
{
    ibz_t num, denom;
    ibz_init(&num);
    ibz_init(&denom);
    ibq_set(delta, &ibz_const_one, &ibz_const_two);
    ibz_set(&num, EPSILON_NUM);
    ibz_set(&denom, EPSILON_DENOM);
    ibq_set(eta, &num, &denom);
    ibq_add(eta, eta, delta);
    ibz_set(&num, DELTA_NUM);
    ibz_set(&denom, DELTA_DENOM);
    ibq_set(delta, &num, &denom);
    ibz_finalize(&num);
    ibz_finalize(&denom);
}

void
ibq_vec_4_copy_ibz(ibq_vec_4_t *vec, const ibz_t *coeff0, const ibz_t *coeff1, const ibz_t *coeff2, const ibz_t *coeff3)
{
    ibz_t one;
    ibz_init(&one);
    ibz_set(&one, 1);
    ibq_set(&((*vec)[0]), coeff0, &one);
    ibq_set(&((*vec)[1]), coeff1, &one);
    ibq_set(&((*vec)[2]), coeff2, &one);
    ibq_set(&((*vec)[3]), coeff3, &one);
    ibz_finalize(&one);
}

void
quat_lll_bilinear(ibq_t *b, const ibq_vec_4_t *vec0, const ibq_vec_4_t *vec1, const ibz_t *q)
{
    ibq_t sum, prod, norm_q;
    ibz_t one;
    ibz_init(&one);
    ibz_set(&one, 1);
    ibq_init(&sum);
    ibq_init(&prod);
    ibq_init(&norm_q);
    ibq_set(&norm_q, q, &one);

    ibq_mul(&sum, &((*vec0)[0]), &((*vec1)[0]));
    ibq_mul(&prod, &((*vec0)[1]), &((*vec1)[1]));
    ibq_add(&sum, &sum, &prod);
    ibq_mul(&prod, &((*vec0)[2]), &((*vec1)[2]));
    ibq_mul(&prod, &prod, &norm_q);
    ibq_add(&sum, &sum, &prod);
    ibq_mul(&prod, &((*vec0)[3]), &((*vec1)[3]));
    ibq_mul(&prod, &prod, &norm_q);
    ibq_add(b, &sum, &prod);

    ibz_finalize(&one);
    ibq_finalize(&sum);
    ibq_finalize(&prod);
    ibq_finalize(&norm_q);
}

void
quat_lll_gram_schmidt_transposed_with_ibq(ibq_mat_4x4_t *orthogonalised_transposed,
                                          const ibz_mat_4x4_t *mat,
                                          const ibz_t *q)
{
    ibq_mat_4x4_t work;
    ibq_vec_4_t vec;
    ibq_t norm, b, coeff, prod;
    ibq_init(&norm);
    ibq_init(&coeff);
    ibq_init(&prod);
    ibq_init(&b);
    ibq_mat_4x4_init(&work);
    ibq_vec_4_init(&vec);
    // transpose the input matrix to be able to work on vectors
    for (int i = 0; i < 4; i++) {
        ibq_vec_4_copy_ibz(&(work[i]), &((*mat)[0][i]), &((*mat)[1][i]), &((*mat)[2][i]), &((*mat)[3][i]));
    }

    for (int i = 0; i < 4; i++) {
        quat_lll_bilinear(&norm, &(work[i]), &(work[i]), q);
        ibq_inv(&norm, &norm);
        for (int j = i + 1; j < 4; j++) {
            ibq_vec_4_copy_ibz(&vec, &((*mat)[0][j]), &((*mat)[1][j]), &((*mat)[2][j]), &((*mat)[3][j]));
            quat_lll_bilinear(&b, &(work[i]), &vec, q);
            ibq_mul(&coeff, &norm, &b);
            for (int k = 0; k < 4; k++) {
                ibq_mul(&prod, &coeff, &(work[i][k]));
                ibq_sub(&(work[j][k]), &(work[j][k]), &prod);
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibq_copy(&((*orthogonalised_transposed)[i][j]), &(work[i][j]));
        }
    }

    ibq_finalize(&norm);
    ibq_finalize(&coeff);
    ibq_finalize(&prod);
    ibq_finalize(&b);
    ibq_mat_4x4_finalize(&work);
    ibq_vec_4_finalize(&vec);
}

int
quat_lll_verify(const ibz_mat_4x4_t *mat, const ibq_t *delta, const ibq_t *eta, const quat_alg_t *alg)
{
    int res = 1;
    ibq_mat_4x4_t orthogonalised_transposed;
    ibq_vec_4_t tmp_vec;
    ibq_t div, tmp, mu, two, norm, b;
    ibz_t mu2_floored, num, denom;
    ibq_mat_4x4_init(&orthogonalised_transposed);
    ibq_vec_4_init(&tmp_vec);
    ibq_init(&div);
    ibq_init(&tmp);
    ibq_init(&norm);
    ibq_init(&b);
    ibq_init(&mu);
    ibq_init(&two);
    ibz_init(&mu2_floored);
    ibz_init(&num);
    ibz_init(&denom);
    ibz_set(&num, 2);
    ibz_set(&denom, 1);
    ibq_set(&two, &num, &denom);

    quat_lll_gram_schmidt_transposed_with_ibq(&orthogonalised_transposed, mat, &(alg->p));
    // check small bilinear products/norms
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < i; j++) {
            ibq_vec_4_copy_ibz(&tmp_vec, &((*mat)[0][i]), &((*mat)[1][i]), &((*mat)[2][i]), &((*mat)[3][i]));
            quat_lll_bilinear(&b, &(orthogonalised_transposed[j]), &tmp_vec, &(alg->p));
            quat_lll_bilinear(&norm, &(orthogonalised_transposed[j]), &(orthogonalised_transposed[j]), &(alg->p));
            ibq_inv(&tmp, &norm);
            ibq_mul(&mu, &b, &tmp);
            ibq_abs(&mu, &mu);
            // compare to eta
            res = res && (ibq_cmp(&mu, eta) <= 0);
        }
    }
    for (int i = 1; i < 4; i++) {
        ibq_vec_4_copy_ibz(&tmp_vec, &((*mat)[0][i]), &((*mat)[1][i]), &((*mat)[2][i]), &((*mat)[3][i]));
        quat_lll_bilinear(&b, &(orthogonalised_transposed[i - 1]), &tmp_vec, &(alg->p));
        quat_lll_bilinear(&norm, &(orthogonalised_transposed[i - 1]), &(orthogonalised_transposed[i - 1]), &(alg->p));
        ibq_inv(&tmp, &norm);
        ibq_mul(&mu, &b, &tmp);
        // tmp is mu^2
        ibq_mul(&tmp, &mu, &mu);
        // mu is delta-mu^2
        ibq_sub(&mu, delta, &tmp);
        quat_lll_bilinear(&tmp, &(orthogonalised_transposed[i]), &(orthogonalised_transposed[i]), &(alg->p));
        // get (delta-mu^2)norm(i-1)
        ibq_mul(&div, &norm, &mu);
        res = res && (ibq_cmp(&tmp, &div) >= 0);
    }
    ibq_mat_4x4_finalize(&orthogonalised_transposed);
    ibq_vec_4_finalize(&tmp_vec);
    ibq_finalize(&div);
    ibq_finalize(&norm);
    ibq_finalize(&b);
    ibq_finalize(&tmp);
    ibq_finalize(&mu);
    ibq_finalize(&two);
    ibz_finalize(&mu2_floored);
    ibz_finalize(&num);
    ibz_finalize(&denom);
    return (res);
}
