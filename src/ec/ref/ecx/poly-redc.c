#define _POLY_MUL_REDC_H_
#include "poly.h"
#include <assert.h>

void reciprocal(poly h, fp2_t *c, const poly f, const int lenf, const int n){
  
  // Writes a polynomial to h and a field element to c such that f*h = c mod x^n
  // REQUIRES h to have space for n terms
  // NOT responsible for terms in h beyond h[n-1]

  int i;

  // Case when f needs to be padded with zeroes
  if(n > lenf)
  {
    fp2_t fpad[n];
    for(i = 0; i < lenf; i++)
      fp2_copy(&fpad[i], &f[i]);
    for(i = lenf; i < n; i++)
      fp2_set(&fpad[i], 0);
    reciprocal(h, c, fpad, n, n);
    return;
  }

  // Trivial case
  if(n == 0)
  {
    fp2_set(&*c, 0);
    return;
  }

  // Case n = 1
  if(n == 1)
  {
    fp2_copy(&*c, &f[0]);
    fp_mont_setone(h[0].re);fp_set(h[0].im,0);
    return;
  }

  // Case n = 2
  if(n == 2)
  {
    fp2_sqr(&*c, &f[0]);
    fp2_copy(&h[0], &f[0]);
    fp2_neg(&h[1], &f[1]);
    return;
  }

  // Case n = 3
  if(n == 3)
  {
    fp2_t t0, t1;

    fp2_sqr(&t0, &f[1]);
    fp2_mul(&t1, &f[0], &f[2]);
    fp2_sub(&t1, &t1, &t0);
    fp2_mul(&t1, &t1, &f[0]);

    reciprocal(h, c, f, 2, 2);
    fp2_mul(&h[0], &h[0], &*c);
    fp2_mul(&h[1], &h[1], &*c);
    fp2_neg(&h[2], &t1);
    fp2_sqr(&*c, &*c);
    return;
  }

  // Case n = 4
  if(n == 4)
  {
    fp2_t t0, t1, t2, t3, g[2];

    reciprocal(g, &t3, f, 2, 2);
    fp2_sqr(&t0, &f[1]);
    fp2_mul(&t1, &g[0], &f[2]);
    fp2_mul(&t2, &g[0], &f[3]);
    fp2_mul(&h[1], &g[1], &f[2]);
    fp2_sub(&t0, &t1, &t0);
    fp2_add(&t1, &t2, &h[1]);
    fp2_mul(&t2, &t0, &g[0]);
    fp2_mul(&h[1], &t0, &g[1]);
    fp2_mul(&h[3], &t1, &g[0]);
    fp2_add(&h[3], &h[1], &h[3]);
    
    fp2_mul(&h[0], &g[0], &t3);
    fp2_mul(&h[1], &g[1], &t3);
    fp2_neg(&h[2], &t2);
    fp2_neg(&h[3], &h[3]);
    fp2_sqr(&*c, &t3);
    return;
  }


  // General case
  // Compute the reciprocal g mod x^m for m = ceil(n/2)
  // Then f*g-c is multiple of x^m so we only care about terms from m to n-1
  const int m = n - (n>>1);
  fp2_t g[m], t[m], t0;

  reciprocal(g, &t0, f, lenf, m);
  poly_mul_middle(t, g, m, f, n);
  poly_mul_low(t, n-m, g, m, &(t[2*m-n]), n-m);
  for(i = 0; i < m; i++)
    fp2_mul(&h[i], &g[i], &t0);
  for(i = m; i < n; i++)
    fp2_neg(&h[i], &t[i-m]);
  fp2_sqr(&*c, &t0);
  return;
}


void poly_redc(poly h, const poly g, const int leng, const poly f, const int lenf,//
	       const poly f_rev_inv, const fp2_t c)
{
  // Computes h(x) =  a * g(x) mod f(x) for some scalar a, writting lenf-1 terms to h.
  // REQUIRES an inverse f_rev_inv such that f_rev*f_rev_inv = c mod x^(leng-lenf+1),
  // where f_rev is the polynomial with the coefficients of f listed in reverse order.
  // The scalar a is equal to c, except for special cases:
  //    - If leng<lenf (no reduction needed) then a = 1
  //    - If lenf = leng = 2, then a = f[1] 
  //    - If lenf = leng = 3, then a = f[2] 
  //    - If lenf=2, leng=3 then a = 2*f[1]^2
  //
  // REQUIRES h to have space for lenf-1 terms
  // NOT responsible for terms in h beyond h[lenf-2]

  int i;
  
  // Case without reduction
  if(leng < lenf)
  {
    for(i = 0; i < leng; i++)
      fp2_copy(&h[i], &g[i]);
    for(i = leng; i < lenf-1; i++)
      fp2_set(&h[i], 0);
    return;
  }

  // Small cases for f linear
  if(lenf == 2)
  {
    if(leng == 2)
    {
      fp2_t t0;
      fp2_mul(&t0, &g[0], &f[1]);
      fp2_mul(&h[0], &g[1], &f[0]);
      fp2_sub(&h[0], &t0, &h[0]);
      return;
    }
    
    if(leng == 3)
    {
      fp2_t f0f1, f02, f12;
      fp2_sqr(&f02, &f[0]);
      fp2_sqr(&f12, &f[1]);
      fp2_sub(&f0f1, &f[0], &f[1]);
      fp2_sqr(&f0f1, &f0f1);
      fp2_sub(&f0f1, &f0f1, &f02);
      fp2_sub(&f0f1, &f0f1, &f12);
      fp2_add(&f02, &f02, &f02);
      fp2_add(&f12, &f12, &f12);
      fp2_mul(&f02, &f02, &g[2]);
      fp2_mul(&f12, &f12, &g[0]);
      fp2_mul(&f0f1, &f0f1, &g[1]);
      fp2_add(&h[0], &f02, &f12);
      fp2_add(&h[0], &h[0], &f0f1);
      return;
    }
  }

  // Small case for f cuadratic
  if(lenf == 3 && leng == 3)
  {
    fp2_t f2g1, f2g0, f1g2;
    fp2_mul(&f2g1, &g[1], &f[2]);
    fp2_mul(&f2g0, &g[0], &f[2]);
    fp2_mul(&f1g2, &g[2], &f[1]);
    fp2_mul(&h[0], &g[2], &f[0]);
    fp2_sub(&h[0], &f2g0, &h[0]);
    fp2_sub(&h[1], &f2g1, &f1g2);
    return;
  }

  // General case
  fp2_t g_reversed[leng], Q[leng - lenf + 1], Q_reversed[leng - lenf + 1];
  
  for(i = 0; i < leng; i++)
    fp2_copy(&g_reversed[i], &g[leng-1-i]);

  poly_mul_low(Q, leng-lenf+1, f_rev_inv, leng-lenf+1, g_reversed, leng-lenf+1);

  for(i = 0; i < leng - lenf + 1; i++)
    fp2_copy(&Q_reversed[i], &Q[leng - lenf - i]);

  poly_mul_low(g_reversed, lenf-1, Q_reversed, leng-lenf+1, f, lenf);

  for(i = 0; i < lenf-1; i++)
  {
    fp2_mul(&h[i], &g[i], &c);
    fp2_sub(&h[i], &h[i], &g_reversed[i]);
  }
  return;
}


void reciprocal_tree(poly *R, fp2_t *A, const int leng, const poly H[], const int DEG[],//
		     const int root, const int n)
{
  // Given a product tree H with degrees tree DEG rooted at root and generated 
  // by n polynomials, writes the reverse-reciprocal polynomials to R and field elements 
  // to A such that Rev(H[i])*R[i] = A[i] mod x^(N) for all nodes but the leaves.
  // The mod is N = deg(parent)-deg(self) for inner nodes, or N = leng - deg(root) for the root.
  //
  // REQUIRES that leng >= DEG[0] and that R,A have enough space for the tree (see product_tree)

  if(n == 0)
    return;

  const int parent = (root-1) >> 1;
  const int brother = root - 1 + 2*(root & 1);
  int lenr;

  if(root > 0)
    lenr = DEG[parent] - DEG[root];
  else
    lenr = leng - DEG[root];
  
  R[root] = malloc(sizeof(fp2_t)*lenr);
  
  // ----------------------------------
  // base cases determined by poly_redc
  if(n == 1)
    return;


  // case for computing  g mod f when len(f), len(g) = 3
  if (DEG[root] == 2 && lenr == 1)
  {
    reciprocal_tree(R, A, lenr-1, H, DEG, 2*root+1, n-(n>>1));
    reciprocal_tree(R, A, lenr-1, H, DEG, 2*root+2, n>>1);
    return;
  }
  
  // ----------------------------------

  int i;
  
  // When the parent's inverse was calculated to a smaller modulus, need to invert from scratch
  if(root == 0 || leng < lenr)
  {
    for(i = 0; i < lenr && i < DEG[root]+1; i++)
      fp2_copy(&R[root][i], &H[root][DEG[root]-i]);
    for(i = DEG[root]+1; i < lenr; i++){
      fp2_set(&R[root][i], 0);
    }
    reciprocal(R[root], &(A[root]), R[root], lenr, lenr);
  }
  else
  {
  // When parent's inverse was to a greater/equal modulus, this inverse can be obtained from it
    for(i = 0; i < lenr; i++)
      fp2_copy(&R[root][i], &H[brother][DEG[brother]-i]);
    poly_mul_low(R[root], lenr, R[parent], leng, R[root], lenr);
    fp2_copy(&A[root], &A[parent]);
  }

  // Now move on to the children
  reciprocal_tree(R, A, lenr-1, H, DEG, 2*root+1, n-(n>>1));
  reciprocal_tree(R, A, lenr-1, H, DEG, 2*root+2, n>>1);
  return;
}


void multieval_unscaled(fp2_t REM[], const poly g, const int leng, const poly R[], const fp2_t A[],//
		const poly H[], const int DEG[], const int root, const int n)
{
  // Given the product tree H and reciprocal tree R,A generated by f_0, ... , f_{n-1},
  // with corresponding degrees tree DEG[] and rooted at root,  writes the constant term 
  // of c_i*g mod f_i to REM[i]. The constants c_i are unspecified, but are a function
  // only of leng and f_0,...,f_{n-1} so they cancel out when taking the ratios of
  // remainders of different g's of the same length.
  //
  // REQUIRES REM to have space for n terms

  if(n == 0)
    return;
  
  fp2_t g_mod[DEG[root]];
  poly_redc(g_mod, g, leng, H[root], DEG[root]+1, R[root], A[root]);

  if(n == 1)
  {
    fp2_copy(&REM[0], &g_mod[0]);
    return;
  }
  
  multieval_unscaled(REM, g_mod, DEG[root], R, A, H, DEG, 2*root+1, n-(n>>1));
  multieval_unscaled(&(REM[n-(n>>1)]), g_mod, DEG[root], R, A, H, DEG, 2*root+2, n>>1);
  return;
}


void multieval_scaled(fp2_t REM[], const poly G, const poly H[], //
			   const int DEG[], const int root, const int n)
{
  // Given the product tree H generated by LINEAR f_0,...,f_{n-1} rooted at root and with
  // corresponding degrees tree DEG, writes the constant term of c_i * g mod f_i(x) to REM[i]
  // The constants c_i are unspecified but are only a function of leng and f_0,...,f_{n-1},
  // so they cancel out when taking the ratios of remainders of different g's of the same length.
  //
  // REQUIRES REM to have space for n terms and n > 1
  // Also REQUIRES G = rev((rev(g mod F)) * F_rev_inv mod x^deg(F)-1) where F = H[root]
  // and F_rev_inv is its reverse's reciprocal mod x^deg(F)

  if(root == 0)
  {
    if(n == 1)
    {
      fp2_copy(&REM[0], &G[DEG[root]-1]);
      return;
    }
    else
    {
      multieval_scaled(REM, G, H, DEG, 2*root+1, n-(n>>1));
      multieval_scaled(&(REM[n-(n>>1)]), G, H, DEG, 2*root+2, n>>1);
      return;
    }
  }
    
  const int parent = (root-1) >> 1;
  const int brother = root - 1 + 2*(root & 1);
  const int uncle = parent - 1 + 2*(parent & 1);
  fp2_t fg[DEG[brother]+1];

  if(root > 2)
    poly_mul_middle(fg, H[brother], DEG[brother]+1, G, DEG[uncle]+1);
  else
    poly_mul_middle(fg, H[brother], DEG[brother]+1, G, DEG[0]);
    
  
  if(n == 1)
  {
    fp2_copy(&REM[0], &fg[DEG[brother]]);
    return;
  }

  multieval_scaled(REM, fg, H, DEG, 2*root+1, n-(n>>1));
  multieval_scaled(&(REM[n-(n>>1)]), fg, H, DEG, 2*root+2, n>>1);
  return;
}
