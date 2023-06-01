#include "isog.h"
#include "curve_extras.h"
#include <assert.h>

int sI, sJ, sK;	// Sizes of each current I, J, and K	

fp2_t I[sI_max][2],		// I plays also as the linear factors of the polynomial h_I(X)
			EJ_0[sJ_max][3], EJ_1[sJ_max][3];	// To be used in xisog y xeval

ec_point_t J[sJ_max], K[sK_max];		// Finite subsets of the kernel
fp2_t XZJ4[sJ_max],		// -4* (Xj * Zj) for each j in J, and x([j]P) = (Xj : Zj)
    rtree_A[(1 << (ceil_log_sI_max+2)) - 1],		// constant multiple of the reciprocal tree computation
    A0;			// constant multiple of the reciprocal R0

poly ptree_hI[(1 << (ceil_log_sI_max+2)) - 1],		// product tree of h_I(X)
     rtree_hI[(1 << (ceil_log_sI_max+2)) - 1],		// reciprocal tree of h_I(X)
     ptree_EJ[(1 << (ceil_log_sJ_max+2)) - 1];		// product tree of E_J(X)
     
fp2_t R0[2*sJ_max + 1];		// Reciprocal of h_I(X) required in the scaled remainder tree approach

int deg_ptree_hI[(1 << (ceil_log_sI_max+2)) - 1],	// degree of each noed in the product tree of h_I(X)
    deg_ptree_EJ[(1 << (ceil_log_sJ_max+2)) - 1];	// degree of each node in the product tree of E_J(X)

fp2_t leaves[sI_max];		// leaves of the remainder tree, which are required in the Resultant computation

// -----------------------------------------------------------
// -----------------------------------------------------------
// Traditional Kernel Point computation (KPs)

// Kernel computation required in tye degree-4 isogeny evaluation
void kps_4(ec_point_t const P)
{
	fp2_sub(&K[1].x, &P.x, &P.z);
	fp2_add(&K[2].x, &P.x, &P.z);
	fp2_sqr(&K[0].x, &P.z);
	fp2_add(&K[0].z, &K[0].x, &K[0].x);
	fp2_add(&K[0].x, &K[0].z, &K[0].z);
}

void eds2mont(ec_point_t* P)
{
	fp2_t t;
	fp2_add(&t, &(P->z), &(P->x));
	fp2_sub(&(P->z), &(P->z), &(P->x));
	fp2_copy(&(P->x), &t);
}


// Differential doubling in Twisted Edwards model
void ydbl(ec_point_t* Q, ec_point_t* const P, ec_point_t const* A)
{
	fp2_t t_0, t_1, X, Z;

	fp2_sqr(&t_0, &(P->x));
	fp2_sqr(&t_1, &(P->z));
	fp2_mul(&Z, &(A->z), &t_0);
	fp2_mul(&X, &Z, &t_1);
	fp2_sub(&t_1, &t_1, &t_0);
	fp2_mul(&t_0, &(A->x), &t_1);
	fp2_add(&Z, &Z, &t_0);
	fp2_mul(&Z, &Z, &t_1);

	fp2_sub(&(Q->x), &X, &Z);
	fp2_add(&(Q->z), &X, &Z);
}

// Differential addition in Twisted Edwards model
void yadd(ec_point_t* R, ec_point_t* const P, ec_point_t* const Q, ec_point_t* const PQ)
{
	fp2_t a, b, c, d, X, Z;

	fp2_mul(&a, &(P->z), &(Q->x));
	fp2_mul(&b, &(P->x), &(Q->z));
	fp2_add(&c, &a, &b);
	fp2_sub(&d, &a, &b);
	fp2_sqr(&c, &c);
	fp2_sqr(&d, &d);

	fp2_add(&a, &(PQ->z), &(PQ->x));
	fp2_sub(&b, &(PQ->z), &(PQ->x));
	fp2_mul(&X, &b, &c);
	fp2_mul(&Z, &a, &d);

	fp2_sub(&(R->x), &X, &Z);
	fp2_add(&(R->z), &X, &Z);
}

// tvelu formulae
void kps_t(uint64_t const i, ec_point_t const P, ec_point_t const A)
{
	int j;
	int d = ((int)TORSION_ODD_PRIMES[i] - 1) / 2;

	// Mapping the input point x(P), which belongs to a 
	// Montogmery curve model, into its Twisted Edwards 
	// representation y(P)
	fp2_sub(&K[0].x, &P.x, &P.z);
	fp2_add(&K[0].z, &P.x, &P.z);
	ydbl(&K[1], &K[0], &A);				// y([2]P)

	for (j = 2; j < d; j++)
		yadd(&K[j], &K[j - 1], &K[0], &K[j - 2]);	// y([j+1]P)
}

// -----------------------------------------------------------
// -----------------------------------------------------------
// Kernel Point computation (KPs) used in velu SQRT
void kps_s(uint64_t const i, ec_point_t const P, ec_point_t const A)
{
	// =================================================================================
	assert(TORSION_ODD_PRIMES[i] > gap);	// Ensuring velusqrt is used for l_i > gap
	// The optimal bounds must corresponds to sI, sJ, and sK

	sI = sizeI[i];	// Size of I
	sJ = sizeJ[i];	// Size of J
	sK = sizeK[i];	// Size of K
	assert(sI >= sJ);	// Ensuring #I >= #J
	assert(sK >= 0);	// Recall, it must be that #K >= 0
	assert(sJ > 1);		// ensuring sI >= sJ > 1
	// =================================================================================
	
	// Now, we can proceed by the general case

	int j;

	// --------------------------------------------------
	// Computing [j]P for each j in {1, 3, ..., 2*sJ - 1}
	ec_point_t P2, P4;
	copy_point(&J[0], &P);				//    x(P)
	// Next computations are required for allowing the use of the function get_A()
	fp2_mul(&XZJ4[0], &J[0].x, &J[0].z);					//   Xj*Zj
	fp2_add(&XZJ4[0], &XZJ4[0], &XZJ4[0]);					//  2Xj*Zj
	fp2_add(&XZJ4[0], &XZJ4[0], &XZJ4[0]);					//  4Xj*Zj
	fp2_neg(&XZJ4[0], &XZJ4[0]);					// -4Xj*Zj
	xDBLv2(&P2, &P, &A);					// x([2]P)
	xADD(&J[1], &P2, &J[0], &J[0]);			// x([3]P)
	// Next computations are required for allowing the use of the function get_A()
	fp2_mul(&XZJ4[1], &J[1].x, &J[1].z);					//   Xj*Zj
	fp2_add(&XZJ4[1], &XZJ4[1], &XZJ4[1]);					//  2Xj*Zj
	fp2_add(&XZJ4[1], &XZJ4[1], &XZJ4[1]);					//  4Xj*Zj
	fp2_neg(&XZJ4[1], &XZJ4[1]);					// -4Xj*Zj
	for (j = 2; j < sJ; j++)
	{
		xADD(&J[j], &J[j - 1], &P2, &J[j - 2]);	// x([2*j + 1]P)
		// Next computations are required for allowing the use of the function get_A()
		fp2_mul(&XZJ4[j], &J[j].x, &J[j].z);					//   Xj*Zj
		fp2_add(&XZJ4[j], &XZJ4[j], &XZJ4[j]);					//  2Xj*Zj
		fp2_add(&XZJ4[j], &XZJ4[j], &XZJ4[j]);					//  4Xj*Zj
		fp2_neg(&XZJ4[j], &XZJ4[j]);					// -4Xj*Zj
	};

	// ----------------------------------------------------------
	// Computing [i]P for i in { (2*sJ) * (2i + 1) : 0 <= i < sI}
	// and the linear factors of h_I(W)
	ec_point_t Q, Q2, tmp1, tmp2;
	int bhalf_floor= sJ >> 1;
	int bhalf_ceil = sJ - bhalf_floor;
	xDBLv2(&P4, &P2, &A);								// x([4]P)
	swap_points(&P2, &P4, -(uint64_t)(sJ % 2));								// x([4]P) <--- coditional swap ---> x([2]P)
	xADD(&Q, &J[bhalf_ceil], &J[bhalf_floor - 1], &P2);	// Q := [2b]P
	swap_points(&P2, &P4, -(uint64_t)(sJ % 2));								// x([4]P) <--- coditional swap ---> x([2]P)

	// .............................................
	xDBLv2(&Q2, &Q, &A);					// x([2]Q)
	xADD(&tmp1, &Q2, &Q, &Q);	// x([3]Q)
	fp2_neg(&I[0][0], &Q.x);
	fp2_copy(&I[0][1], &Q.z);
	fp2_neg(&I[1][0], &tmp1.x);
	fp2_copy(&I[1][1], &tmp1.z);
	copy_point(&tmp2, &Q);
	
	for (j = 2; j < sI; j++){
		xADD(&tmp2, &tmp1, &Q2, &tmp2);	// x([2*j + 1]Q)
		fp2_neg(&I[j][0], &tmp2.x);
		fp2_copy(&I[j][1], &tmp2.z);
		swap_points(&tmp1, &tmp2, -(uint64_t)1);
	}


	// ----------------------------------------------------------------
	// Computing [k]P for k in { 4*sJ*sI + 1, ..., l - 6, l - 4, l - 2}
	// In order to avoid BRANCHES we make allways copy in K[0] and K[1]
	// by assuming that these entries are only used when sK >= 1 and 
	// sK >= 2, respectively.

	//if (sK >= 1)
	copy_point(&K[0], &P2);				//       x([l - 2]P) = x([2]P)
	//if (sK >= 2)
	copy_point(&K[1], &P4);				//       x([l - 4]P) = x([4]P)
	
	for (j = 2; j < sK; j++)
		xADD(&K[j], &K[j - 1], &P2, &K[j - 2]);	// x([l - 2*(j+1)]P) = x([2 * (j+1)]P)

	// ----------------------------------------------------------------
	//                   ~~~~~~~~               ~~~~~~~~
	//                    |    |                 |    |
	// Computing h_I(W) = |    | (W - x([i]P)) = |    | (Zi * W - Xi) / Zi where x([i]P) = Xi/Zi
	//                    i in I                 i in I
	// In order to avoid costly inverse computations in fp, we are gonna work with projective coordinates

	product_tree_LENFeq2(ptree_hI, deg_ptree_hI, 0, I, sI);				// Product tree of hI
	if (!scaled)
	{
		// (unscaled) remainder tree approach
		reciprocal_tree(rtree_hI, rtree_A, 2*sJ + 1, ptree_hI, deg_ptree_hI, 0, sI);	// Reciprocal tree of hI
	}
	else
	{
		// scaled remainder tree approach
		fp2_t f_rev[sI_max + 1];
		for (j = 0; j < (sI + 1); j++)
			fp2_copy(&f_rev[j], &ptree_hI[0][sI - j]);

		if (sI > (2*sJ - sI + 1))
			reciprocal(R0, &A0, f_rev, sI + 1, sI);
		else
			reciprocal(R0, &A0, f_rev, sI + 1, 2*sJ - sI + 1);
	};
}

void kps_clear(int i){
		if (TORSION_ODD_PRIMES[i] > gap)
		{
			if (!scaled)
				clear_tree(rtree_hI, 0, sizeI[i]);
			clear_tree(ptree_hI, 0, sizeI[i]);
		}
}