#define _POLY_MUL_REDC_H_
#include <poly.h>
#include <assert.h>

void poly_mul(poly h, const poly f, const int lenf, const poly g, const int leng)
{
	// Karatsuba multiplication of polynomials h = f*g
	// REQUIRES  that h have enough space for lenf + leng - 1 terms
	// NOT responsible for terms in h beyond h[lenf + leng - 2]

	int i;
	const int lenh = lenf + leng - 1;

	if (lenf < leng)
	{
		poly_mul(h, g, leng, f, lenf);
		return;
	}

	// Resulting degree is negative
	if(lenh <= 0)
		return;

	// Multiplication by 0 polynomial
	if(leng == 0)
	{
		for(i = 0; i < lenf-1; i++)
			fp2_set(&h[i], 0);
		return;
	}

	// Multiplication by constant
	if(leng == 1)
	{
	        fp2_t fc[lenf], gc; // copy of f and g
		for(i = 0; i < lenf; i++)
			fp2_copy(&fc[i], &f[i]);

		fp2_copy(&gc, &g[0]);
		for (i = 0; i < lenh; i++)
			fp2_mul(&h[i], &fc[i], &gc);

		
		return;
	}

	// At this point we ensure lenf >= leng >= 2
  
	// Case when lenf = leng = 2
	if (lenf == 2)
	{
	        fp2_t t0, t1, t2;
		fp2_add(&t1, &f[0], &f[1]);
		fp2_add(&t2, &g[0], &g[1]);
		fp2_mul(&t0, &f[0], &g[0]);
		fp2_mul(&h[2], &f[1], &g[1]);
		fp2_mul(&t1, &t1, &t2);
		fp2_sub(&t1, &t1, &t0);
		fp2_sub(&h[1], &t1, &h[2]);
		fp2_copy(&h[0], &t0);
		return;
	}
  
	// Cases for f cuadratic
	if (lenf == 3)
	{
		// Case when f is cuadratic and g linear
		if (leng == 2)
		{
		  fp2_t t0, t1, t2, t3;
			fp2_mul(&t0, &f[0], &g[0]);
			fp2_mul(&t2, &f[1], &g[1]);
			fp2_add(&t3, &f[0], &f[1]);
			fp2_add(&t1, &g[0], &g[1]);
			fp2_mul(&t1, &t1, &t3);
			fp2_sub(&t1, &t1, &t2);
			fp2_sub(&t1, &t1, &t0);
			fp2_mul(&t3, &f[2], &g[1]);
			fp2_mul(&h[2], &f[2], &g[0]);
			fp2_add(&h[2], &h[2], &t2);
		        fp2_copy(&h[0], &t0);
			fp2_copy(&h[1], &t1);
			fp2_copy(&h[3], &t3);
			return;
		}
   
		// Case when f,g both cuadratic
		if (leng == 3)
		{
		        fp2_t fg_high[3], t0, t1, t2;

			poly_mul(fg_high, &(f[1]), 2, &(g[1]), 2);

			fp2_add(&t0, &f[0], &f[1]);
			fp2_add(&t1, &g[0], &g[1]);
			fp2_mul(&t1, &t0, &t1);
			fp2_add(&t0, &f[0], &f[2]);
			fp2_add(&t2, &g[0], &g[2]);
			fp2_mul(&t2, &t0, &t2);
			
			fp2_mul(&h[0], &f[0], &g[0]);

			fp2_sub(&t1, &t1, &h[0]);
			fp2_sub(&h[1], &t1, &fg_high[0]);
			fp2_sub(&t2, &t2, &h[0]);
			fp2_sub(&t2, &t2, &fg_high[2]);
			fp2_add(&h[2], &t2, &fg_high[0]);
			fp2_copy(&h[3], &fg_high[1]);
			fp2_copy(&h[4], &fg_high[2]);
      
			return;
		}
	}

	// At this point we ensure lenf >= 4 and lenf >= leng
	
	const int nf = lenf >> 1;
	const int mf = lenf - nf;
	const int mg = leng - nf;

	// Case when g is half the size of f  
	if (leng <= nf)
	{
		fp2_t fg_low[nf + leng -1], fg_high[mf + leng -1];
    
		poly_mul(fg_low, f, nf, g, leng);
		poly_mul(fg_high, &(f[nf]), mf, g, leng);
		for(i = 0; i < nf; i++)
			fp2_copy(&h[i], &fg_low[i]);
    
		for(i = nf; i < nf + leng - 1; i++)
			fp2_add(&h[i], &fg_high[i-nf], &fg_low[i]);

		for(i = nf + leng - 1; i < lenh; i++)
			fp2_copy(&h[i], &fg_high[i - nf]);

		return;
	}
  

	// All other cases
	// Strategy is to spli f and g, into
	// f(x) = f0(x) + x^nf*f1(x) and g(x) = g0(x) + x^nf*g1(x) where nf = floor(lenf/2)
	// such that f(x)*g(x) = fg_low(x) + x^nf*fg_mid(x) + x^2nf*fg_high(x) where:
	//   fg_low(x) = f0(x)*g0(x)
	//   fg_mid(x) = f0(x)*g1(x) + f1(x)*g0(x)
	//   fg_high(x) = f1(x)*g1(x)
  
	fp2_t f_mid[mf], g_mid[mf], fg_low[2*nf-1], fg_mid[2*mf-1], fg_high[mf+mg-1];
  
	for(i = 0; i < nf; i++)
		fp2_add(&f_mid[i], &f[i], &f[i+nf]);

	if(lenf & 1)
		fp2_copy(&f_mid[nf], &f[lenf-1]);
  
	i = 0;
	while(i < nf && i < mg)
	{
		fp2_add(&g_mid[i], &g[i], &g[nf+i]);
		i++;
	}
	while(i < nf)
	{
		fp2_copy(&g_mid[i], &g[i]);
		i++;
	}
	
	poly_mul(fg_low, f, nf, g, nf);
	poly_mul(fg_high, &(f[nf]), mf, &(g[nf]), mg);
	
	if((lenf & 1) && (mg == mf))
	{
		fp2_copy(&g_mid[nf], &g[leng - 1]);
		poly_mul(fg_mid, f_mid, mf, g_mid, mf);
	}
	else
	{
	        poly_mul(fg_mid, f_mid, mf, g_mid, nf);
	}


	for(i = 0; i < mf + mg - 1; i++)
		fp2_sub(&fg_mid[i], &fg_mid[i], &fg_high[i]);

	for(i = 0; i < 2*nf-1; i++)
		fp2_sub(&fg_mid[i], &fg_mid[i], &fg_low[i]);

	for(i = 0; i < nf; i++)
		fp2_copy(&h[i], &fg_low[i]);

	for(i = nf; i < 2*nf-1; i++)
		fp2_add(&h[i], &fg_low[i], &fg_mid[i-nf]);

	fp2_copy(&h[2*nf-1], &fg_mid[nf-1]);
  
	for(i = 2*nf; i < 2*nf+mf-1; i++)
		fp2_add(&h[i], &fg_mid[i-nf], &fg_high[i-2*nf]);

	for(i = 2*nf+mf-1; i < lenh; i++)
		fp2_copy(&h[i], &fg_high[i-2*nf]);

	return;
}

void poly_mul_low(poly h, const int n, const poly f, const int lenf, const poly g, const int leng)
{
	// Karatsuba multiplication of polynomials h = f*g mod x^n
	// REQUIRES that h have enough space for n terms
	// NOT responsible for terms in h beyond h[n-1]

	// Ensure f is the one with highest degree
	if(leng > lenf)
	{
		poly_mul_low(h, n, g, leng, f, lenf);
		return;
	}

	// Cleave g and f are too big
	if(leng > n)
	{
		poly_mul_low(h, n, f, n, g, n);
		return;
	}

	// Cleave f is too big
	if(lenf > n)
	{
		poly_mul_low(h, n, f, n, g, leng);
		return;
	}
  
	// Case when resulting degree is too small
	if(n == 0)
		return;

	// Multiplication by zero polynomial
	if( (leng == 0) || (lenf == 0) )
	{
		int i;
		for(i = 0; i < n; i++)
			fp2_set(&h[i], 0);

		return;
	}
  
	// Multiplication mod x
	if(n == 1)
	{
		fp2_mul(&h[0], &f[0], &g[0]);
		return;
	}
  
	// Case when no reduction is necessary
	if(n >= lenf + leng - 1)
	{
		int i;
		poly_mul(h, f, lenf, g, leng);
		for(i = lenf + leng - 1; i < n; i++)
			fp2_set(&h[i], 0);

		return;
	}

	// Multiplication by a constant
	if(leng == 1)
	{
		fp2_t fc[n], gc; // Copies for f,g
		int i;
		for(i = 0; i < n; i++)
			fp2_copy(&fc[i], &f[i]);

		fp2_copy(&gc, &g[0]);
		for(i = 0; i < n; i++)
			fp2_mul(&h[i], &fc[i], &gc);

		return;
	}


	// Multiplication mod x^2 of two linear polynomials
	if(n == 2)
	{
		fp2_t fg[2];

		fp2_mul(&fg[1], &f[1], &g[0]);
		fp2_mul(&fg[0], &f[0], &g[1]);
		fp2_add(&fg[1], &fg[1], &fg[0]);
		fp2_mul(&h[0], &f[0], &g[0]);
		fp2_copy(&h[1], &fg[1]);
		return;
	}
  
	// Cases for multiplication mod x^3
	if(n == 3)
	{

		// Multiplication mod x^3 of linear by cuadratic or linear
		if(leng == 2)
		{
		  fp2_t t0, t1, t2;
			fp2_add(&t0, &f[0], &f[1]);
			fp2_add(&t1, &g[0], &g[1]);
			fp2_mul(&t1, &t0, &t1);
			fp2_mul(&t2, &f[0], &g[0]);
			fp2_sub(&t1, &t1, &t2);
			fp2_mul(&t0, &f[1], &g[1]);
			fp2_sub(&t1, &t1, &t0);
			fp2_mul(&h[2], &f[2], &g[0]);
			fp2_add(&h[2], &t0, &h[2]);
			fp2_copy(&h[1], &t1);
			fp2_copy(&h[0], &t2);
			return;
		}

		// Multiplication mod x^3 of two cuadratic polynomials
		if(leng == 3)
		{
		        fp2_t t0, t1, t2, t3, t4;
			fp2_mul(&t0, &f[0], &g[0]);
			fp2_mul(&t1, &f[1], &g[1]);
			fp2_add(&t2, &f[0], &f[1]);
			fp2_add(&t3, &g[0], &g[1]);
			fp2_mul(&t2, &t2, &t3);
			fp2_sub(&t2, &t2, &t0);
			fp2_sub(&t2, &t2, &t1);

			fp2_add(&t3, &f[0], &f[2]);
			fp2_add(&t4, &g[0], &g[2]);
			fp2_mul(&t3, &t3, &t4);
			fp2_mul(&t4, &f[2], &g[2]);
			fp2_sub(&t3, &t3, &t4);
			fp2_sub(&t3, &t3, &t0);
			fp2_add(&h[2], &t3, &t1);

			fp2_copy(&h[0], &t0);
			fp2_copy(&h[1], &t2);
			return;
		}
	}

	// Special case for small values, currently only used for n=4, leng=4
	if((n==4) && (leng==4))
	{
		int i, j, k;
		int S = n;
		const int nf = n >> 1;
		const int nc = n - nf;
		fp2_t t0, t1;
    
		for(i = 0; i < nf; i++)
		      S = S + n - 1 - 2*i;

		fp2_t c[S];
    
		for(i = 0; i < n; i++)
			fp2_mul(&c[i], &f[i], &g[i]);
		
		k = n;
		for(i = 0; i < nf; i++)
		{
			for(j = 0; j < n-2*i-1; j++)
			{
				fp2_add(&t0, &f[i], &f[i+j+1]);
				fp2_add(&t1, &g[i], &g[i+j+1]);
				fp2_mul(&c[k], &t0, &t1);
				fp2_add(&t0, &c[i], &c[i+j+1]);
				fp2_sub(&c[k], &c[k], &t0);
				k++;
			}
		}

		fp2_copy(&c[n-1], &c[0]);
		for(i = 1; i < nf; i++)
		{
			for(j = 1; j < n-2*i; j++)
				fp2_add(&c[n+2*i-1+j], &c[n+2*i-1+j], &c[(1+i)*n-i*i-1+j]);
		}

		for(i = 1; i < nc; i++)
			fp2_add(&c[n+2*i-1], &c[n+2*i-1], &c[i]);

		for(i = n-1; i < 2*n-1; i++)
			fp2_copy(&h[i-n+1], &c[i]);
		
		return;
	}
  

	// All other cases (n >= 4 and reduction is necessary for f*g but not f nor g):
  
	// Strategy is to split f, g into odd and even powers. E.j. f(x) = f0(x^2) + x*f1(x^2)
	// Such that f*g(x) = fg_0(x^2) + x*fg_01(x^2) + x^2*fg_1(x^2) where:
	//    fg_0(x) = f0(x)*g0(x) mod x^ceil(n/2)
	//    fg_01(x) = (f0(x)g1(x)+f1(x)g0(x)) mod x^floor(n/2)
	//    fg_1(x) = f1(x)g1(x) mod x^ceil(n/2-1)
  
	int i;
	const int l1 = n >> 1;
	const int l0 = n - l1;

	//Split f
	const int lenf1 = lenf >> 1;
	const int lenf0 = lenf - lenf1;
	fp2_t f0[lenf0], f1[lenf1];
	for(i = 0; i < lenf1; i++)
	{
		fp2_copy(&f0[i], &f[2*i]);
		fp2_copy(&f1[i], &f[2*i+1]);
	}

	if(lenf0 > lenf1)
		fp2_copy(&f0[lenf0-1], &f[lenf - 1]);

	//Split g
	const int leng1 = leng >> 1;
	const int leng0 = leng - leng1;
	fp2_t g0[leng0], g1[leng1];
	for(i = 0; i < leng1; i++)
	{
		fp2_copy(&g0[i], &g[2*i]);
		fp2_copy(&g1[i], &g[2*i+1]);
	}
	
	if(leng0 > leng1)
		fp2_copy(&g0[leng0-1], &g[leng - 1]);
  
	//Compute f01 = f0 + f1
	fp2_t f01[lenf0];
	for(i = 0; i < lenf1; i++)
		fp2_add(&f01[i], &f0[i], &f1[i]);

	if(lenf0 > lenf1)
		fp2_copy(&f01[lenf0-1], &f0[lenf0-1]);

	//Compute g01 = g0 + g1
	fp2_t g01[leng0];
	for(i = 0; i < leng1; i++)
		fp2_add(&g01[i], &g0[i], &g1[i]);

	if(leng0 > leng1)
		fp2_copy(&g01[leng0-1], &g0[leng0-1]);
  
	//Compute fg_0(x)
	int n0 = lenf0 + leng0 - 1;
	if(n0 > l0)
		n0 = l0;
	
	fp2_t fg_0[n0];
	poly_mul_low(fg_0, n0, f0, lenf0, g0, leng0);

	//Compute fg_1(x)
	int n1 = lenf1 + leng1 - 1;
	if(n1 > l1)
		n1 = l1;

	fp2_t fg_1[n1];
	poly_mul_low(fg_1, n1, f1, lenf1, g1, leng1);

	//Compute fg_01(x) = f_01(x) * g_01(x) mod x^floor(n/2);
	int n01 = lenf0 + leng0 -1;
	if(n01 > l1)
		n01 = l1;

	fp2_t fg_01[n01];
	poly_mul_low(fg_01, n01, f01, lenf0, g01, leng0);
  
	//Substractic fg_0 and fg_1 from fg_01
	//Note that the sizes satisfy n1 <= n01 <= n0
	//Result can be stored directly to h already
	
	i = 0;
	while(i < n1)
	{
	        fp2_sub(&fg_01[i], &fg_01[i], &fg_0[i]);
		fp2_sub(&h[2*i+1], &fg_01[i], &fg_1[i]);
		i++;
	}
	while(i < n01)
	{
	        fp2_sub(&h[2*i+1], &fg_01[i], &fg_0[i]);
		i++;
	}

	//Combine the result for fg_1 and fg_0
	fp2_copy(&h[0], &fg_0[0]);
	for(i = 1; i < n1; i++)
	        fp2_add(&h[2*i], &fg_0[i], &fg_1[i-1]);
	if(2*n1 < n)
	        fp2_add(&h[2*n1], &fg_0[n1], &fg_1[n1-1]);

	return;
}


// The following two function compute the middle section of a product and are based
// on "The Middle Product Algorithm, I." by G. Hanrot. M. Quercia, and P. Zimmermann
void quasi_poly_mul_middle(poly h, const poly g, const int leng, const poly f, const int lenf)
{
 	// Computes the middle part of the product, sepcifically fg[lenf-leng:lenf],
	// for the special case of lenf = 2*leng-1
	// REQUIRES that h have space for leng terms and that lenf = 2*leng-1
	// NOT responsible for terms in h beyond h[leng-1]
	assert(lenf >= (2*leng - 1));

	if(leng == 0)
		return;

	if(leng == 1)
	{
		fp2_mul(&h[0], &g[0], &f[0]);
		return;
	}
  
	int i;
	const int leng0 = leng >> 1;
	const int leng1 = leng - leng0;
	const int lenF = 2*leng1-1;  
	fp2_t F[lenF], G[leng1];
 
	fp2_t A[leng1];
	for(i = 0; i < lenF; i++)
		fp2_add(&F[i], &f[i], &f[i+leng1]);

	quasi_poly_mul_middle(A, &(g[leng0]), leng1, F, lenF);

	fp2_t B[leng1];
	if(leng & 1)
	{
		fp2_copy(&G[0], &g[leng0]);
		for(i = 0; i < leng0; i++)
			fp2_sub(&G[i+1], &g[leng1+i], &g[i]);

	}
	else
	{
		for(i = 0; i < leng0; i++)
			fp2_sub(&G[i], &g[leng1+i], &g[i]);

	}
	
	quasi_poly_mul_middle(B, G, leng1, &(f[leng1]), lenF);

	fp2_t C[leng0];
	for(i = 0; i < 2*leng0-1; i++)
		fp2_add(&F[i], &f[i+leng1], &f[i+2*leng1]);

	quasi_poly_mul_middle(C, g, leng0, F, 2*leng0-1);

	for(i = 0; i < leng1; i++)
		fp2_sub(&h[i], &A[i], &B[i]);

	for(i = 0; i < leng0; i++)
		fp2_add(&h[i+leng1], &C[i], &B[i]);

	return;
}

void poly_mul_middle(poly h, const poly g, const int leng, const poly f, const int lenf)
{
	// Computes the middle part of the product, sepcifically fg[lenf-leng:lenf], for lenf >= leng
	// Note that this is equivalent to the highest leng terms of fg mod x^(lenf)
	// REQUIRES that h have space for leng terms and that lenf > leng
	// NOT responsible for terms in h beyond h[leng-1]

	int i;
  
	// Case of deg(f) odd and deg(g) = floor(deg(f)/2), ie lengths n and 2n
	if ( (leng == lenf >> 1) && !(lenf & 1) )
	{
		// f1 = f[1:]+[0]
		fp2_t f1[lenf]; 
		for(i = 0; i < lenf - 1; i++)
			fp2_copy(&f1[i], &f[i+1]);

		fp2_set(&f1[lenf-1], 0);
		quasi_poly_mul_middle(h, g, leng, f1, lenf);
		return;
	}

	// Case of deg(f) odd and deg(g) = ceil(deg(f)/2), ie lengths n and 2n-2
	if( (leng == (lenf>>1)+1) && !(lenf & 1) )
	{
		// f1 = [0]+f[:]
		fp2_t f1[lenf+1];
		fp2_set(&f1[0], 0);
		for(i = 0; i < lenf; i++)
			fp2_copy(&f1[i+1], &f[i]);

		quasi_poly_mul_middle(h, g, leng, f1, lenf+1);
		return;
	}

	// Case of deg(f) even and deg(g) = deg(f)/2, ie lengths n and 2n-1
	if( (leng == (lenf>>1)+1) && (lenf & 1) )
	{
		quasi_poly_mul_middle(h, g, leng, f, lenf);
		return;
	}

  
	// Unbalanced case, for deg(g) > ceil(deg(f)/2) or deg(g) < floor(deg(f)/2)

	if(leng == 0)
		return;

	const int lenF0 = lenf-leng;
	fp2_t F0[lenF0], G[leng];
	for(i = 0; i < lenF0; i++)
		fp2_copy(&F0[i], &f[lenF0-1-i]);

	for( i = 0; i < leng; i++)
		fp2_copy(&G[i], &g[leng-1-i]);

	fp2_t fg_low[leng-1], fg_low_reverse[leng-1];
	poly_mul_low(fg_low, leng-1, F0, lenF0, G, leng);
	for(i = 0; i < leng-1; i++)
		fp2_copy(&fg_low_reverse[i], &fg_low[leng-2-i]);

	fp2_t fg_high[leng];
	poly_mul_low(fg_high, leng, &(f[lenF0]), leng, g, leng);

	for(i = 0; i < leng-1; i++)
		fp2_add(&h[i], &fg_low_reverse[i], &fg_high[i]);

	fp2_copy(&h[leng-1], &fg_high[leng-1]);
	return;
}


void poly_mul_selfreciprocal(poly h, const poly g, const int leng, const poly f, const int lenf)
{
	// Computes the product h= f*g of self reciprocal polynomials of same length
	// Note: a polynomial is self reciprocal if f(x) = x^deg(f) * f(1/x), or equivalently
	//       if the list of coefficients is a palyndrome.
	// REQUIRES that h have space for lenf+leng-1 entries
	// NOT responsible for terms in h beyond h[lenf+leng-1]

	// Case for same length
	if(lenf==leng)
	{
		if(lenf == 0)
			return;

		if(lenf == 1)
		{
			fp2_mul(&h[0], &g[0], &f[0]);
			return;
		}

		if(lenf == 2)
		{
			fp2_mul(&h[0], &g[0], &f[0]);
			fp2_add(&h[1], &h[0], &h[0]);
			fp2_copy(&h[2], &h[0]);
			return;
		}

		if(lenf == 3)
		{
			fp2_t t0, t1;
			fp2_add(&t0, &g[0], &g[1]);
			fp2_add(&t1, &f[0], &f[1]);
			fp2_mul(&t1, &t0, &t1);
			fp2_mul(&t0, &g[0], &f[0]);
			fp2_mul(&h[2], &g[1], &f[1]);
			fp2_add(&h[2], &h[2], &t0);
			fp2_copy(&h[0], &t0);
			fp2_sub(&h[1], &t1, &h[2]);
			fp2_add(&h[2], &h[2], &h[0]);
			fp2_copy(&h[3], &h[1]);
			fp2_copy(&h[4], &h[0]);
			return;
		}

		if(lenf == 4)
		{
			fp2_t t0, t1;
			fp2_add(&t0, &g[0], &g[1]);
			fp2_add(&t1, &f[0], &f[1]);
			fp2_mul(&t1, &t0, &t1);
			fp2_mul(&t0, &g[0], &f[0]);
			fp2_mul(&h[3], &g[1], &f[1]);
			fp2_copy(&h[2], &t1);
			fp2_copy(&h[0], &t0);
			fp2_sub(&h[2], &h[2], &h[0]);
			fp2_sub(&h[1], &h[2], &h[3]);
			fp2_add(&h[3], &h[3], &h[0]);
			fp2_add(&h[3], &h[3], &h[3]);
			fp2_copy(&h[4], &h[2]);
			fp2_copy(&h[5], &h[1]);
			fp2_copy(&h[6], &h[0]);
			return;
		}

		if(lenf == 5)
		{
			fp2_t t0, t1, t2, t3, t4;
			fp2_sub(&t1, &g[1], &g[0]);
			fp2_sub(&t0, &f[0], &f[1]);
			fp2_mul(&t1, &t1, &t0);
			fp2_sub(&t2, &g[2], &g[0]);
			fp2_sub(&t0, &f[0], &f[2]);
			fp2_mul(&t2, &t0, &t2);
			fp2_sub(&t3, &g[2], &g[1]);
			fp2_sub(&t0, &f[1], &f[2]);
			fp2_mul(&t3, &t3, &t0);
			fp2_mul(&t0, &g[1], &f[1]);
			fp2_mul(&t4, &g[2], &f[2]);
			fp2_mul(&h[0], &g[0], &f[0]);
			fp2_copy(&h[1], &t1);
			fp2_copy(&h[2], &t2);
			fp2_copy(&h[3], &t3);
			fp2_add(&h[5], &t0, &h[0]);
			fp2_add(&h[1], &h[1], &h[5]);
			fp2_add(&h[3], &h[3], &h[1]);
			fp2_add(&h[4], &h[5], &t4);
			fp2_add(&h[4], &h[4], &h[5]);
			fp2_add(&h[2], &h[2], &h[5]);
			fp2_add(&h[2], &h[2], &t4);
			fp2_add(&h[3], &h[3], &t0);
			fp2_add(&h[3], &h[3], &t4);
			fp2_copy(&h[5], &h[3]);
			fp2_copy(&h[6], &h[2]);
			fp2_copy(&h[7], &h[1]);
			fp2_copy(&h[8], &h[0]);
			return;
		}

		// General case for same lengths
		if(lenf & 1)
		{
			// Odd length
			// Strategy is to split f, g into odd and even powers. E.j. f(x) = f0(x^2) + x*f1(x^2)
			// Such that f*g(x) = h0(x^2) + x*h01(x^2) + x^2*h1(x^2) where:
			//    h0(x) = f0(x)*g0(x) 
			//    h01(x) = (f0(x)g1(x)+f1(x)g0(x))
			//    h1(x) = f1(x)g1(x) 

			int i;
			const int len1 = lenf >> 1;
			const int len0 = len1 + 1;

			fp2_t g0[len0], f0[len0], g1[len1], f1[len1];
			for(i = 0; i < len1; i++)
			{
				fp2_copy(&g0[i], &g[2*i]);
				fp2_copy(&f0[i], &f[2*i]);
				fp2_copy(&g1[i], &g[2*i+1]);
				fp2_copy(&f1[i], &f[2*i+1]);
			}

			fp2_copy(&g0[len1], &g[2*len1]);
			fp2_copy(&f0[len1], &f[2*len1]);

			fp2_t h0[2*len0-1];
			poly_mul_selfreciprocal(h0, g0, len0, f0, len0);

			fp2_t h1[2*len1-1];
			poly_mul_selfreciprocal(h1, g1, len1, f1, len1);

			fp2_t h01[2*len0-1];
			for(i = 0; i < len1; i++)
			{
				fp2_add(&g0[i], &g0[i], &g1[i]);
				fp2_add(&f0[i], &f0[i], &f1[i]);
				fp2_add(&g0[i+1], &g0[i+1], &g1[i]);
				fp2_add(&f0[i+1], &f0[i+1], &f1[i]);
			}
			poly_mul_selfreciprocal(h01, g0, len0, f0, len0);
       
			// Mixing results for the odd degree part
			for(i = 0; i < 2*len0-1; i++)
				fp2_sub(&h01[i], &h01[i], &h0[i]);

			for(i = 0; i < 2*len1-1; i++)
				fp2_sub(&h01[i], &h01[i], &h1[i]);

			for(i = 0; i < 2*len1-1; i++)
			{
				fp2_sub(&h01[i+1], &h01[i+1], &h1[i]);
				fp2_sub(&h01[i+1], &h01[i+1], &h1[i]);
			}
			for(i = 0; i < 2*len1-1; i++)
				fp2_sub(&h01[i+2], &h01[i+2], &h1[i]);

			fp2_copy(&h[1], &h01[0]);
			for(i = 1; i < 2*len0-2; i++)
				fp2_sub(&h[2*i+1], &h01[i], &h[2*i-1]);

			// Mixing results for the even degree parts
			fp2_copy(&h[0], &h0[0]);
			for(i = 1; i < 2*len1; i++)
				fp2_add(&h[2*i], &h0[i], &h1[i-1]);

			fp2_copy(&h[4*len1], &h0[2*len1]);
		  	return;
		}
		else
		{
			// Even length
			// Strategy is to spli f and g, into
			// f(x) = f0(x) + x^nf*f1(x) and g(x) = g0(x) + x^nf*g1(x) where nf = floor(lenf/2)
			// such that f(x)*g(x) = fg_low(x) + x^nf*fg_mid(x) + x^2nf*fg_high(x) where:
			//   fg_low(x) = f0(x)*g0(x)
			//   fg_mid(x) = f0(x)*g1(x) + f1(x)*g0(x) <-[reverse of f0g1]
			//   fg_high(x) = f1(x)*g1(x)  <-[reverse of fg_low]
	
	  		int i;
			const int half = leng >> 1;
			fp2_t h1[leng-1];
      
			poly_mul(h, g, half, f, half);
			poly_mul(h1, g, half, &(f[half]), half);

			fp2_set(&h[leng-1], 0);
			for(i = 0; i < half; i++)
			{
				fp2_add(&h[half+i], &h[half+i], &h1[i]);
				fp2_add(&h[half+i], &h[half+i], &h1[leng-2-i]);
			}

			for(i = 0; i < leng-1; i++)
				fp2_copy(&h[2*leng-2-i], &h[i]);

			return;
		}
	}
	else
	{
		// Case for different lengths
		const int m = (lenf+leng) >> 1;
		int i;

		poly_mul_low(h, m, g, leng, f, lenf);
		for(i = m; i < lenf+leng-1; i++)
			fp2_copy(&h[i], &h[lenf+leng-2-i]);

		return;
	}
}


void product_tree(poly H[], int DEG[], const int root, const poly F[], const int LENF, const int n)
{
	// Given an array F containing n polynomials of the same length LENF,
	// Writes the product tree of the polynomials to the tree rooted at H[root].
	//
	// The tree is represented by an array where, for node H[i], the left child is H[2*i+1] and
	// the right child is H[2*i+2]. The root of the complete tree is intended to be H[0].
	//
	// Also writes to DEG the tree containing the degree of the corresponding polynomials.
	//
	// REQUIRES H and DEG to have enough space for 2^(k+2)-1 terms where k = ceil(log2(n))

	int i;
	if(n == 1)
	{
		H[root] = malloc(sizeof(fp2_t)*LENF);
		for(i = 0; i < LENF; i++)
			fp2_copy(&H[root][i], &F[0][i]);

		DEG[root] = LENF-1;
		return;
	}
  
	const int n1 = n >> 1;
	const int n0 = n - n1;
	const int left = 2*root+1;
	const int right = left + 1;
	product_tree(H, DEG, left, F, LENF, n0);
	product_tree(H, DEG, right, &(F[n0]), LENF, n1);
	DEG[root] = DEG[left] + DEG[right];
	H[root] = malloc(sizeof(fp2_t)*(DEG[left]+DEG[right]+1));
	poly_mul(H[root], H[left], DEG[left]+1, H[right], DEG[right]+1);
	return;
}



void product_tree_LENFeq2(poly H[], int DEG[], const int root, const fp2_t F[][2], const int n)
{
	// Same as product tree but allows F to be given as a double array with LENF fixed to 2
	int i;
	if(n == 1)
	{
		H[root] = malloc(sizeof(fp2_t)*2);
		for(i = 0; i < 2; i++)
			fp2_copy(&H[root][i], &F[0][i]);

		DEG[root] = 1;
		return;
	}
  
	const int n1 = n >> 1;
	const int n0 = n - n1;
	const int left = 2*root+1;
	const int right = left + 1;
	product_tree_LENFeq2(H, DEG, left, F, n0);
	product_tree_LENFeq2(H, DEG, right, &(F[n0]), n1);
	DEG[root] = DEG[left] + DEG[right];
	H[root] = malloc(sizeof(fp2_t)*(DEG[left]+DEG[right]+1));
	poly_mul(H[root], H[left], DEG[left]+1, H[right], DEG[right]+1);
	return;
}



void product_tree_LENFeq3(poly H[], int DEG[], const int root, const fp2_t F[][3], const int n)
{
	// Same as product tree but allows F to be given as a double array with LENF fixed to 3
	int i;
	if(n == 1)
	{
		H[root] = malloc(sizeof(fp2_t)*3);
		for(i = 0; i < 3; i++)
			fp2_copy(&H[root][i], &F[0][i]);

		DEG[root] = 3-1;
		return;
	}
  
	const int n1 = n >> 1;
	const int n0 = n - n1;
	const int left = 2*root+1;
	const int right = left + 1;
	product_tree_LENFeq3(H, DEG, left, F, n0);
	product_tree_LENFeq3(H, DEG, right, &(F[n0]), n1);
	DEG[root] = DEG[left] + DEG[right];
	H[root] = malloc(sizeof(fp2_t)*(DEG[left]+DEG[right]+1));
	poly_mul(H[root], H[left], DEG[left]+1, H[right], DEG[right]+1);
	return;
}


void product_tree_selfreciprocal(poly H[], int DEG[], const int root, const poly F[], const int LENF, const int n)
{
	// Same as product_tree but for selfreciprocal inputs

	int i;
	if(n == 1)
	{
		H[root] = malloc(sizeof(fp2_t)*LENF);
		for(i = 0; i < LENF; i++)
			fp2_copy(&H[root][i], &F[0][i]);

		DEG[root] = LENF-1;
		return;
	}
  
	const int n1 = n >> 1;
	const int n0 = n - n1;
	const int left = 2*root+1;
	const int right = left + 1;
	product_tree_selfreciprocal(H, DEG, left, F, LENF, n0);
	product_tree_selfreciprocal(H, DEG, right, &(F[n0]), LENF, n1);
	DEG[root] = DEG[left] + DEG[right];
	H[root] = malloc(sizeof(fp2_t)*(DEG[left]+DEG[right]+1));
	poly_mul_selfreciprocal(H[root], H[left], DEG[left]+1, H[right], DEG[right]+1);
	return;
}

void product_tree_selfreciprocal_LENFeq3(poly H[], int DEG[], const int root, const fp2_t F[][3], const int n)
{
	// Same as product_tree_selfreciprocal but allows F to be given as a double array with LENF fixed to 3

	int i;
	if(n == 1)
	{
		H[root] = malloc(sizeof(fp2_t)*3);
		for(i = 0; i < 3; i++)
			fp2_copy(&H[root][i], &F[0][i]);

		DEG[root] = 3-1;
		return;
	}
  
	const int n1 = n >> 1;
	const int n0 = n - n1;
	const int left = 2*root+1;
	const int right = left + 1;
	product_tree_selfreciprocal_LENFeq3(H, DEG, left, F, n0);
	product_tree_selfreciprocal_LENFeq3(H, DEG, right, &(F[n0]), n1);
	DEG[root] = DEG[left] + DEG[right];
	H[root] = malloc(sizeof(fp2_t)*(DEG[left]+DEG[right]+1));
	poly_mul_selfreciprocal(H[root], H[left], DEG[left]+1, H[right], DEG[right]+1);
	return;
}

void clear_tree(poly H[], const int root, const int n)
{
	// Frees all dynamic memory in a polynomial tree rooted at H[root] and generated by n leafs

	if(n == 1)
	{
		free(H[root]);
		return;
	}

	clear_tree(H, 2*root+1, n-(n>>1));
	clear_tree(H, 2*root+2, n>>1);
	free(H[root]);
	return;
}



void product(fp2_t *c, const fp2_t F[], const int n)
{
	// Given an array F of n constant polynomials, writes their product to c
	int i;

	fp_mont_setone((*c).re);fp_set((*c).im,0);

	// Empty list must returns 1
	if (n == 0)
		return;

	// At this step, we ensure there is at least one element in the list
	fp2_copy(&*c, &F[0]);
	for(i = 1; i < n; i++)
		fp2_mul(&*c, &*c, &F[i]);

	return;
}
