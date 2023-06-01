#include <poly.h>
#include <assert.h>
#include <stdio.h>

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

void slow_mul(poly h, poly f, int lenf, poly g, int leng){
  // Computes h = f*g by school method

  fp2_t a, b;
  int nf, ng, e;
  int lenh = lenf + leng - 1;
  
  if(lenh <= 0){
    return;
  }
  
  fp2_t fg[lenh];
  
  if (leng > lenf){
    slow_mul(h, g, leng, f, lenf);
    return;
  }
  
  for(e = 0; e < lenh; e++){

    if (lenf - 1 < e){
      nf = lenf - 1;
    }
    else{
      nf = e;
    }

    ng = e - nf;
    fp2_set(&a, 0);
    while( (ng < leng) & (nf >= 0) ){
      fp2_mul(&b, &f[nf], &g[ng]);
      fp2_add(&a, &a, &b);
      nf--;
      ng++;
    }
    fp2_copy(&fg[e], &a);
  }
  for(e = 0; e < lenh; e++){
    fp2_copy(&h[e], &fg[e]);
  }
  return;
}



int main(){
  fp2_t fp2_0, fp2_1;
  #define nmax 16
  int nf, ng, n, e;
        fp2_set(&fp2_0, 0);
        fp_mont_setone(fp2_1.re);fp_set(fp2_1.im,0); 
  
  //TEST MULTIPLICATION BY 0
  
  for(nf = 2; nf < nmax; nf++){
    fp2_t f[nf], h[nf-1];

    printf("[%3d%%] Testing multiplication by 0", 100 * nf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    for(e = 0; e < nf; e++){
      fp2_random(&f[e]);
    }
    poly_mul(h, f, nf, f, 0);
    for(e = 0; e < nf-1; e++){
      assert(fp2_is_zero(&h[e])==1);
    }
    poly_mul(h, f, 0, f, nf);
    for(e = 0; e < nf-1; e++){
      assert(fp2_is_zero(&h[e])==1);
    }
  }
  printf("[%3d%%] Tested multiplication by 0:\t\tNo errors!\n", 100 * nf / nmax);

  
  
  //TEST FOR f, g, h DISJOINT MEMORY SPACES
  
  for(nf = 1; nf < nmax; nf++){
    
    printf("[%3d%%] Testing multiplication", 100 * nf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    for(ng = 1; ng < nmax; ng++){
      
      fp2_t f[nf];   //Random length nf poly
      for(e = 0; e < nf; e++){
	fp2_random(&f[e]);
      }
      
      fp2_t g[ng];  // Random length ng poly
      for(e = 0; e < ng; e++){
	fp2_random(&g[e]);
      }
      
      fp2_t h[nf+ng-1];// Compute product
      poly_mul(h, f, nf, g, ng);

      fp2_t fg[nf+ng-1]; // Compute the product by school method
      slow_mul(fg, f, nf, g, ng);
      
      for(e = 0; e < nf + ng - 1; e++){   // Verify answer term by term
	assert(fp2_isequal(h[e], fg[e])==1);
      }
    }
  }
  printf("[%3d%%] Tested multiplication:\t\t\tNo errors!\n", 100 * nf / nmax);

  

  // TEST FOR f, g CONTIGIOUS AND RESULT SAVED OVER THEM
    
  for(nf = 1; nf < nmax; nf++){
          
    printf("[%3d%%] Testing multiplication in place", 100 * nf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    for(ng = 1; ng < nmax; ng++){
      
      fp2_t h[nf+ng];
      
      //Random length nf poly
      for(e = 0; e < nf; e++){
	fp2_random(&h[e]);
      }
      
      // Random length ng poly
      for(e = 0; e < ng; e++){
	fp2_random(&h[e+nf]);
      }

      // Compute the product
      fp2_t fg[nf+ng-1];
      slow_mul(fg, h, nf, &(h[nf]), ng); // School method
      poly_mul(h, h, nf, &(h[nf]), ng); // Karatsuba method


      for(e = 0; e < nf + ng - 1; e++){   // Verify answer term by term
	assert(fp2_isequal(h[e], fg[e])==1);
      }
    }
  }
    printf("[%3d%%] Tested multiplication in place:\t\tNo errors!\n", 100 * nf / nmax);

    
    
  //TEST FOR MULTIPLICATION MOD X^N BY 0
    
  for(nf = 2; nf < nmax; nf++){
    fp2_t f[nf];
    
    printf("[%3d%%] Testing mul mod x^n by 0", 100 * nf / nmax);
    fflush(stdout);
    printf("\r\x1b[K");
    
    for(e = 0; e < nf; e++){
      fp2_random(&f[e]);
    }
    
    for(n = 1; n < nmax; n++){
      fp2_t h[n];
      poly_mul_low(h, n, f, nf, f, 0);
      for(e = 0; e < n; e++){
	assert(fp2_is_zero(&h[e])==1);
      }
      poly_mul_low(h, n, f, 0, f, nf);
      for(e = 0; e < n; e++){
	assert(fp2_is_zero(&h[e])==1);
      }
    }
  }
  printf("[%3d%%] Tested mul mod x^n by 0:\t\t\tNo errors!\n", 100 * nf / nmax);

  
  
  //TEST FOR MULTIPLICATION MOD X^N
    
    for(nf = 1; nf < nmax; nf++){
    
      printf("[%3d%%] Testing mul mod x^n", 100 * nf / nmax);
      fflush(stdout);
      printf("\r\x1b[K");
      
      for(ng = 1; ng < nmax; ng++){

	fp2_t f[nf], g[ng], fg[nf+ng-1];
	poly h;

	//Get random polynomials
	for(e = 0; e < nf; e++){
	  fp2_random(&f[e]);
	}
	for(e = 0; e < ng; e++){
	  fp2_random(&g[e]);
	}
	
	//Save regular result to fg
	slow_mul(fg, f, nf, g, ng);

	//Compute result mod x^n
	for(n = 1; n < 2*nmax; n++){
	  h = malloc(sizeof(fp2_t)*n);
	  poly_mul_low(h, n, f, nf, g, ng);

	  //Compare with expected
	  e = 0;
	  while(e < nf+ng-1 && e < n){
	    assert(fp2_isequal(h[e], fg[e]) == 1);
	    e++;
	  }
	  while(e < n){
	    assert(fp2_is_zero(&h[e]) == 1);
	    e++;
	  }
	  free(h);
	}
      }
    }
    printf("[%3d%%] Tested mul mod x^n:\t\t\tNo errors!\n", 100 * nf / nmax);

  
     
  //TEST FOR POLY_MUL_MIDDLE
    
    for(nf = 1; nf < 2*nmax; nf+=1){
      fp2_t f[nf];
      
      printf("[%3d%%] Testing poly_mul_middle", 100 * nf / (2*nmax));
      fflush(stdout);
      printf("\r\x1b[K");
      
      for(ng = (nf+1)>>1; ng < (nf+1)-((nf+1)>>1); ng++){
	// This runs from floor((nf+1)/2) to ceil((nf+1)/2)
	fp2_t g[ng];
	for(e = 0; e < nf; e++){
	  fp2_random(&f[e]);
	}
	for(e = 0; e < ng; e++){
	  fp2_random(&g[e]);
	}
	
	fp2_t h[nf+ng-1];
	slow_mul(h, g, ng, f, nf);
	poly_mul_middle(g, g, ng, f, nf);
      
	for(e = 0; e < ng; e++){
	  assert(fp2_isequal(h[e+nf-ng], g[e])==1);
	}
      }
    }
    printf("[%3d%%] Tested poly_mul_middle:\t\t\tNo errors!\n", 100 * nf / (2*nmax));

  
  // TEST FOR SELF RECIPROCAL MULTIPLICATION
    for(nf = 1; nf < nmax; nf++){

      printf("[%3d%%] Testing self reciprocal mul", 100 * nf / nmax);
      fflush(stdout);
      printf("\r\x1b[K");

      for(ng = 1; ng < nmax; ng++){
      
	fp2_t f[nf], g[ng], h[nf+ng-1], fg[nf+ng-1];

	// Get random palyndromes
	for(e = 0; e < (nf>>1); e++){
	  fp2_random(&f[e]);
	  fp2_copy(&f[nf-1-e], &f[e]);
	}
	if(nf & 1){
	  fp2_random(&f[nf>>1]);
	}

	for(e = 0; e < (ng>>1); e++){
	  fp2_random(&g[e]);
	  fp2_copy(&g[ng-1-e], &g[e]);
	}
	if(ng & 1){
	  fp2_random(&g[ng>>1]);
	} 

	// Compute products
	poly_mul_selfreciprocal(h, g, ng, f, nf);
	slow_mul(fg, g, ng, f, nf);

	// Compare
	for(e = 0; e < nf+ng-1; e++){
	  assert(fp2_isequal(fg[e], h[e])==1);
	}
      }
    }		 
    printf("[%3d%%] Tested self reciprocal mul:\t\tNo errors!\n", 100 * nf / nmax);

  // TEST FOR PRODUCT TREES
    int tree_size, iteration, i;
    int  len, *DEG, LENF;
    poly *H, *F, h;
    
    for(tree_size = 1; tree_size < nmax; tree_size++){

      printf("[%3d%%] Testing product tree:\t\t\tSize %d out of %d", 100 * tree_size / nmax, tree_size, nmax-1);
      fflush(stdout);
      printf("\r\x1b[K");

      i = 0;
      while((1<<i) < tree_size){
	i++;
      }
      DEG = malloc(sizeof(int)*((1<<(i+2))-1));
      H = malloc(sizeof(poly)*((1<<(i+2))-1));
      F = malloc(sizeof(poly)*tree_size);
      h = malloc(sizeof(fp2_t)*(nmax+1)*tree_size);

      for(iteration = 0; iteration < nmax + 1 - tree_size ; iteration++){

	// Generate random list of polynomials
	LENF = (rand() % nmax)+1;
	for(i = 0; i < tree_size; i++){
	  F[i] = malloc(sizeof(fp2_t)*LENF);
	  for(e = 0; e < LENF; e++){
	    fp2_random(&F[i][e]);
	  }
	}
	product_tree(H, DEG, 0, F, LENF, tree_size);
	
	// Build product of all polynomials manually
	len = LENF;
	
	//for(e = 0; e < LENF[0]; e++){
	for(e = 0; e < LENF; e++){
	  fp2_copy(&h[e], &F[0][e]);
	}
	for(i = 1; i < tree_size; i++){
	  poly_mul(h, h, len, F[i], LENF);
	  len += LENF-1;
	}

	// Compare to root
	assert (len == DEG[0]+1);
	for(e = 0; e < len; e++){
	  assert(fp2_isequal(H[0][e], h[e])==1);
	}
      clear_tree(H, 0, tree_size);
      for(i = 0; i < tree_size; i++){
	free(F[i]);
      }

      }
      free(DEG);
      free(H);
      free(F); 
      free(h);
    }
    printf("[%3d%%] Tested product tree:\t\t\tNo errors!\n", 100 * tree_size / nmax);
    
  // TEST FOR SELF RECIPROCAL PRODUCT TREES
    
    for(tree_size = 1; tree_size < nmax; tree_size++){

      printf("[%3d%%] Testing selfreciprocal product tree:\tSize %d out of %d", 100 * tree_size / nmax, tree_size, nmax-1);
      fflush(stdout);
      printf("\r\x1b[K");

      i = 0;
      while((1<<i) < tree_size){
	i++;
      }
      DEG = malloc(sizeof(int)*((1<<(i+2))-1));
      H = malloc(sizeof(poly)*((1<<(i+2))-1));
      F = malloc(sizeof(poly)*tree_size);
      h = malloc(sizeof(fp2_t)*(nmax+1)*tree_size);

      for(iteration = 0; iteration < nmax + 1 - tree_size ; iteration++){

	// Generate random list of polynomials
	LENF = (rand() % nmax)+1;;
	for(i = 0; i < tree_size; i++){
	  F[i] = malloc(sizeof(fp2_t)*LENF);
	  for(e = 0; e < (LENF>>1); e++){
	    fp2_random(&F[i][e]);
	    fp2_copy(&F[i][LENF-1-e], &F[i][e]);
	  }
	  if(LENF & 1){
	  	fp2_random(&F[i][(LENF>>1)]);
	  }
	}
	product_tree_selfreciprocal(H, DEG, 0, F, LENF, tree_size);
	
	// Build product of all polynomials manually
	len = LENF;
	for(e = 0; e < LENF; e++){
	  fp2_copy(&h[e], &F[0][e]);
	}
	for(i = 1; i < tree_size; i++){
	  poly_mul(h, h, len, F[i], LENF);
	  len += LENF-1;
	}

	// Compare to root
	assert (len == DEG[0]+1);
	for(e = 0; e < len; e++){
	  assert(fp2_isequal(H[0][e], h[e])==1);
	}
      clear_tree(H, 0, tree_size);
      for(i = 0; i < tree_size; i++){
	free(F[i]);
      }

      }
      free(DEG);
      free(H);
      free(F); 
      free(h);
    }
    printf("[%3d%%] Tested selfreciprocal product tree:\tNo errors!\n", 100 * tree_size / nmax);
    
    printf("-- All tests passed.\n");
    return 0;
}
  
