#include "poly.h"
#include <assert.h>
#include <stdio.h>
#define nmax 32

bool fp2_isequal(fp2_t a, fp2_t b){
    return fp_is_equal(a.re, b.re) && fp_is_equal(a.im, b.im);
}

// VERY NOT SECURE (testing only)
void fp2_random(fp2_t *a){
    for(int i = 0; i < NWORDS_FIELD; i++){
        a->re[i] = rand();
        a->im[i] = rand();
    }
    // Normalize
    fp2_t one;
    fp_mont_setone(one.re);fp_set(one.im,0);
    fp2_mul(&*a, &*a, &one);
    // Update seed
    srand((unsigned) a->re[0]);
}

int main(){
  fp2_t fp2_0, fp2_1;
  fp2_set(&fp2_0, 0);
  fp_mont_setone(fp2_1.re);fp_set(fp2_1.im,0);

  int lenf, leng, n, e, iteration, array_size, tree_size, i, root, brother, *DEG, LENF;
  poly f, g, h, f_rev, f_rev_inv, *F, *H, *R, g1, g2, REM1, REM2, G1, G2, G1_rev, G2_rev, R0;
  fp2_t c, *A, *C, ratio, A0;
  
  f_rev_inv = 0;
  
// TEST FOR RECIPROCAL
  for(lenf = 1; lenf < nmax; lenf++)
  {  
    printf("[%3d%%] Testing reciprocals", 100 * lenf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random poly
    f = malloc(sizeof(fp2_t)*lenf);
    for(e = 0; e < lenf; e++)
      fp2_random(&f[e]);

    for(n = 1; n < nmax; n++)
    {
      // Get the reciprocal and multiply them
      h = malloc(sizeof(fp2_t)*n);
      memset(h, 0, sizeof(fp2_t)*n);
      reciprocal(h, &c, f, lenf, n);
      poly_mul_low(h, n, f, lenf, h, n);

      // Compare with expected
      assert(fp2_isequal(h[0],c));
      for(e = 1;  e < n; e++)
	assert(fp2_is_zero(&h[e]));
      free(h);
    }
    free(f); 
  }
  printf("[%3d%%] Tested reciprocals:\t\tNo errors!\n", 100 * lenf / nmax);
  
  

  // TEST FOR REDUCTION
  for(lenf = 2; lenf < nmax; lenf++)
  {
    printf("[%3d%%] Testing polynomial reduction", 100 * lenf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");

    // Get random poly for the mod
    f = malloc(sizeof(fp2_t)*lenf);
    f_rev = malloc(sizeof(fp2_t)*lenf);
    for(e = 0; e < lenf; e++)
    {
      fp2_random(&f[e]);
      fp2_copy(&f_rev[lenf-1-e], &f[e]);
    }

    for(leng = 1; leng < nmax; leng++)
    {
      // Get random poly to reduce
      g = malloc(sizeof(fp2_t)*leng);
      for(e = 0; e < leng; e++){
	fp2_random(&g[e]);
      }

      // Get reverse-inverse mod x^(leng-lenf+1)
      if(leng >= lenf)
      {
	f_rev_inv = malloc(sizeof(fp2_t)*(leng-lenf+1));
	reciprocal(f_rev_inv, &c, f_rev, lenf, leng-lenf+1);
      }
      else{
	fp_mont_setone(c.re);fp_set(c.im,0);
      }
	
      // Compute the reduction
      h = malloc(sizeof(fp2_t)*(lenf-1));
      poly_redc(h, g, leng, f, lenf, f_rev_inv, c);

      // Reduce manually
      int leng_red = leng;
      fp2_t scale, f_e;
      while(leng_red >= lenf)
      {
	fp2_copy(&scale, &f[lenf-1]);
	fp2_inv(&scale);
	fp2_mul(&scale, &scale, &g[leng_red-1]);
	for(e = 0; e < lenf; e++)
	  {
	    fp2_mul(&f_e, &f[e], &scale);
	    fp2_sub(&g[e+leng_red-lenf], &g[e+leng_red-lenf], &f_e);
	  }
	leng_red--;
      }

      // Rescale manual result
      if( leng < lenf){
	      fp_mont_setone(scale.re);fp_set(scale.im,0);
      }
      else
	if(lenf == 2 && leng == 3)
	{
	  fp2_sqr(&scale, &f[1]);
	  fp2_add(&scale, &scale, &scale);
	}
	else
	  fp2_copy(&scale, &c);
      for(e = 0; e < leng_red; e++)
	fp2_mul(&g[e], &g[e], &scale);
     

      // Comapre results
      for(e = leng_red-1; e >= 0; e--)
	      assert(fp2_isequal(h[e], g[e]));
      for(e = leng_red; e < lenf-1; e++)
	      assert(fp2_is_zero(&h[e]));
      
      free(g);
      free(h);
      if(leng >= lenf)
	free(f_rev_inv);
    }
    free(f);
    free(f_rev);
  }
  printf("[%3d%%] Tested polynomial reduction:\tNo errors!\n", 100 * lenf / nmax);

  

// TEST FOR RECIPROCAL TREES
  
  for(tree_size = 3; tree_size < nmax; tree_size++)
  {
    printf("[%3d%%] Testing reciprocal tree:\t\tTree size %d out of %d", 100 * tree_size / nmax, tree_size, nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    // Compute size of arrays
    i = 0;
    while((1<<i) < tree_size){
      i++;
    }
    array_size = (1<<(i+2))-1;
    
    DEG = malloc(sizeof(int)*array_size);
    H = malloc(sizeof(poly)*array_size);
    R = malloc(sizeof(poly)*array_size);
    F = malloc(sizeof(poly)*tree_size);
    A = malloc(sizeof(fp2_t)*array_size);
    
    // Get random polys
    LENF = 2;
    for(i = 0; i < tree_size; i++)
    {
      F[i] = malloc(sizeof(fp2_t)*LENF);
      for(e = 0; e < LENF; e++){
	      fp2_random(&F[i][e]);
      }
    }
    
    // Get product tree then reciprocal tree
    product_tree(H, DEG, 0, F, LENF, tree_size);
    leng = DEG[0]+1+(rand() % nmax);
    reciprocal_tree(R, A, leng, H, DEG, 0, tree_size);
    
    // Check the root
    root = 0;
    lenf = leng-DEG[root];
    f = malloc(sizeof(fp2_t)*lenf);
    for(e = 0; e < DEG[root]+1 && e < lenf; e++){
      fp2_copy(&f[e], &H[root][DEG[root]-e]);
    }
    for(e = DEG[root]+1; e < lenf; e++){
      fp2_set(&f[e], 0);
    }
    poly_mul_low(f, lenf, f, lenf, R[root], lenf);
    assert(fp2_isequal(f[0], A[root]));
    for(e = 1; e < lenf; e++){
      assert(fp2_is_zero(&f[e]));
    }
    free(f);
    
    // Perform random walks
    for(iteration = 0; iteration < nmax - tree_size; iteration++)
    {
      root = 0;
      n = tree_size;
      while(n > 1)
      {
	if(rand() & 1)
	{
	  root = 2*root+1;
	  n = n - (n>>1);
	}
	else
	{
	  root = 2*root+2;
	  n = n>>1;
	}
	brother = root - 1 + 2*(root & 1);
	
	// Check current node
	if(DEG[root] > 2)
	{
	  lenf = DEG[brother];
	  f = malloc(sizeof(fp2_t)*lenf);
	  for(e = 0; e < DEG[root]+1 && e < lenf; e++){
	    fp2_copy(&f[e], &H[root][DEG[root]-e]);
    }
	  for(e = DEG[root]+1; e < lenf; e++){
	    fp2_set(&f[e], 0);
    }
	  poly_mul_low(f, lenf, f, lenf, R[root], lenf);
	  assert(fp2_isequal(f[0], A[root]));
	  for(e = 1; e < lenf; e++){
	    assert(fp2_is_zero(&f[e]));
    }
	  free(f);
	}
      }
    }
    // Clean up
    for(i = 0; i < tree_size; i++)
      free(F[i]);
    clear_tree(H, 0, tree_size);
    clear_tree(R, 0, tree_size);
    free(F);
    free(H);
    free(R);
    free(A);
    free(DEG);
  }
  printf("[%3d%%] Tested reciprocal tree:\t\tNo errors!\n", 100 * tree_size / nmax);
  
  

  // TEST FOR REMAINDERS
  for(tree_size = 2; tree_size < nmax; tree_size++)
  {
    printf("[%3d%%] Testing batched remainders:\t\tTree size %d out of %d", 100 * tree_size / nmax, tree_size, nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    // Compute size of arrays
    i = 0;
    while((1<<i) < tree_size)
      i++;
    array_size = (1<<(i+2))-1;
    
    DEG = malloc(sizeof(int)*array_size);
    H = malloc(sizeof(poly)*array_size);
    R = malloc(sizeof(poly)*array_size);
    F = malloc(sizeof(poly)*tree_size);
    A = malloc(sizeof(fp2_t)*array_size);
    REM1 = malloc(sizeof(fp2_t)*array_size);
    REM2 = malloc(sizeof(fp2_t)*array_size);
    C = malloc(sizeof(fp2_t)*tree_size);
    
    // Get random polys
    LENF = 2;
    for(i = 0; i < tree_size; i++)
    {
      F[i] = malloc(sizeof(fp2_t)*LENF);
      for(e = 0; e < LENF; e++)
	fp2_random(&F[i][e]);
    }
    
    // Get product tree, reciprocal tree, and remainders
    product_tree(H, DEG, 0, F, LENF, tree_size);
    leng = DEG[0]+1+(rand() % nmax);
    g1 = malloc(sizeof(fp2_t)*leng);
    g2 = malloc(sizeof(fp2_t)*leng);
    for(e = 0; e < leng; e++)
    {
      fp2_random(&g1[e]);
      fp2_random(&g2[e]);
    }
    reciprocal_tree(R, A, leng, H, DEG, 0, tree_size);
    multieval_unscaled(REM1, g1, leng, R, (const fp2_t*)A, H, DEG, 0, tree_size);
    multieval_unscaled(REM2, g2, leng, R, (const fp2_t*)A, H, DEG, 0, tree_size);
    
    for(i = 0; i < tree_size; i++)
    {
      // Get ratio of the remainder
      fp2_inv(&REM1[i]);
      fp2_mul(&ratio, &REM1[i], &REM2[i]);
      
      // Compute remainders manually
      f_rev = malloc(sizeof(fp2_t)*LENF);
      f_rev_inv = malloc(sizeof(fp2_t)*(leng-LENF+1));
      h = malloc(sizeof(fp2_t)*(LENF-1));
      for(e = 0; e < LENF; e++)
	fp2_copy(&f_rev[e], &F[i][LENF-1-e]);
      reciprocal(f_rev_inv, &c, f_rev, LENF, leng-LENF+1);
      poly_redc(h, g1, leng, F[i], LENF, f_rev_inv, c);
      fp2_copy(&REM1[i], &h[0]);
      poly_redc(h, g2, leng, F[i], LENF, f_rev_inv, c);
      fp2_copy(&REM2[i], &h[0]);
      free(f_rev);
      free(f_rev_inv);
      free(h);

      // Compare results
      fp2_inv(&REM1[i]);
      fp2_mul(&REM1[i], &REM1[i], &REM2[i]);
      assert(fp2_isequal(REM1[i], ratio));
    }
		 
    // Clean up
    for(i = 0; i < tree_size; i++)
      free(F[i]);
    free(g1);
    free(g2);
    clear_tree(H, 0, tree_size);
    clear_tree(R, 0, tree_size);
    free(F);
    free(H);
    free(R);
    free(A);
    free(DEG);
    free(REM1);
    free(REM2);
    free(C);
  } 
  printf("[%3d%%] Tested batched remainders:\tNo errors!\n", 100 * tree_size / nmax);
  


// TEST FOR SCALED REMAINDER TREE
  for(tree_size = 1; tree_size < nmax; tree_size++)
  {
    printf("[%3d%%] Testing scaled remainder tree:\tTree size %d out of %d", 100 * tree_size / nmax, tree_size, nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    // Compute size of arrays
    i = 0;
    while((1<<i) < tree_size)
      i++;
    array_size = (1<<(i+2))-1;
    
    DEG = malloc(sizeof(int)*array_size);
    H = malloc(sizeof(poly)*array_size);
    F = malloc(sizeof(poly)*tree_size);
    REM1 = malloc(sizeof(fp2_t)*array_size);
    REM2 = malloc(sizeof(fp2_t)*array_size);
    
    // Get random polys
    LENF = 2;
    for(i = 0; i < tree_size; i++)
    {
      F[i] = malloc(sizeof(fp2_t)*LENF);
      for(e = 0; e < LENF; e++)
	fp2_random(&F[i][e]);
    }
    
    // Get random polys to reduce
    product_tree(H, DEG, 0, F, LENF, tree_size);
    leng = DEG[0]+1+(rand() % nmax);
    g1 = malloc(sizeof(fp2_t)*leng);
    g2 = malloc(sizeof(fp2_t)*leng);
    for(e = 0; e < leng; e++)
    {
      fp2_random(&g1[e]);
      fp2_random(&g2[e]);
    }

    // Get the required initial nodes
    G1 = malloc(sizeof(fp2_t)*DEG[0]);
    G2 = malloc(sizeof(fp2_t)*DEG[0]);
    G1_rev = malloc(sizeof(fp2_t)*DEG[0]);
    G2_rev = malloc(sizeof(fp2_t)*DEG[0]);
    R0 = malloc(sizeof(fp2_t)*(leng));
    f_rev = malloc(sizeof(fp2_t)*(DEG[0]+1));
    for(e = 0; e < DEG[0]+1; e++)
      fp2_copy(&f_rev[e], &H[0][DEG[0]-e]);
    if( DEG[0] > leng-DEG[0])
      reciprocal(R0, &A0, f_rev, DEG[0]+1, DEG[0]);
    else
      reciprocal(R0, &A0, f_rev, DEG[0]+1, leng-DEG[0]);
    poly_redc(G1, g1, leng, H[0], DEG[0]+1, R0, A0);
    poly_redc(G2, g2, leng, H[0], DEG[0]+1, R0, A0);
    for(e = 0; e < DEG[0]; e++)
    {
      fp2_copy(&G1_rev[e], &G1[DEG[0]-1-e]);
      fp2_copy(&G2_rev[e], &G2[DEG[0]-1-e]);
    }
    poly_mul_middle(G1_rev, G1_rev, DEG[0], R0, DEG[0]);
    poly_mul_middle(G2_rev, G2_rev, DEG[0], R0, DEG[0]);
    for(e = 0; e < DEG[0]; e++)
    {
      fp2_copy(&G1[e], &G1_rev[DEG[0]-1-e]);
      fp2_copy(&G2[e], &G2_rev[DEG[0]-1-e]);
    }
    free(G1_rev);free(G2_rev);free(R0);free(f_rev);

    // Compute the scaled remainder trees
    multieval_scaled(REM1, G1, H, DEG, 0, tree_size);
    multieval_scaled(REM2, G2, H, DEG, 0, tree_size);
    
    for(i = 0; i < tree_size; i++)
    {
      // Get ratio of the remainder
      fp2_inv(&REM1[i]);
      fp2_mul(&ratio, &REM1[i], &REM2[i]);

      // Compute remainders manually
      f_rev = malloc(sizeof(fp2_t)*LENF);
      f_rev_inv = malloc(sizeof(fp2_t)*(leng-LENF+1));
      h = malloc(sizeof(fp2_t)*(LENF-1));
      for(e = 0; e < LENF; e++)
	fp2_copy(&f_rev[e], &F[i][LENF-1-e]);
      reciprocal(f_rev_inv, &c, f_rev, LENF, leng-LENF+1);
      poly_redc(h, g1, leng, F[i], LENF, f_rev_inv, c);
      fp2_copy(&REM1[i], &h[0]);
      poly_redc(h, g2, leng, F[i], LENF, f_rev_inv, c);
      fp2_copy(&REM2[i], &h[0]);
      free(f_rev);free(f_rev_inv);free(h);

      // Compare results
      fp2_inv(&REM1[i]);
      fp2_mul(&REM1[i], &REM1[i], &REM2[i]);
      assert(fp2_isequal(REM1[i], ratio));
    }
		 
    // Clean up
    for(i = 0; i < tree_size; i++)
      free(F[i]);
    free(F);free(g1);free(g2);free(G1);free(G2);
    clear_tree(H, 0, tree_size);free(H);free(DEG);
    free(REM1);free(REM2);
  } 
  printf("[%3d%%] Tested scaled remainder tree:\tNo errors!\n", 100 * tree_size / nmax);
  
  printf("-- All tests passed.\n");
}
