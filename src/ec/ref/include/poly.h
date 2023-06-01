#ifndef _POLY_H_
#define _POLY_H_

#include <fp2.h>

typedef fp2_t *poly; // Polynomials are arrays of coeffs over Fq, lowest degree first

void poly_mul(poly h, const poly f, const int lenf, const poly g, const int leng);
void poly_mul_low(poly h, const int n, const poly f, const int lenf, const poly g, const int leng);
void poly_mul_middle(poly h, const poly g, const int leng, const poly f, const int lenf);
void poly_mul_selfreciprocal(poly h, const poly g, const int leng, const poly f, const int lenf);

void product_tree(poly H[], int DEG[], const int root, const poly F[], const int LENF, const int n);
void product_tree_LENFeq2(poly H[], int DEG[], const int root, const fp2_t F[][2], const int n);
void product_tree_LENFeq3(poly H[], int DEG[], const int root, const fp2_t F[][3], const int n);
void product_tree_selfreciprocal(poly H[], int DEG[], const int root, const poly F[], const int LENF, const int n);
void product_tree_selfreciprocal_LENFeq3(poly H[], int DEG[], const int root, const fp2_t F[][3], const int n);
void clear_tree(poly H[], const int root, const int n);

void product(fp2_t *c, const fp2_t F[], const int n);

void reciprocal(poly h, fp2_t *c, const poly f, const int lenf, const int n);
void poly_redc(poly h, const poly g, const int leng, const poly f, const int lenf,const poly f_inv, const fp2_t c);
void reciprocal_tree(poly *R, fp2_t *A, const int leng, const poly H[], const int DEG[], const int root, const int n);
void multieval_unscaled(fp2_t REM[], const poly g, const int leng, const poly R[], const fp2_t A[], const poly H[], const int DEG[], const int root, const int n);
void multieval_scaled(fp2_t REM[], const poly G, const poly H[], const int DEG[], const int root, const int n);

#endif /* _POLY_H */
