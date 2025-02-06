
#ifndef TEST_EXTRAS_H
#define TEST_EXTRAS_H

#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <encoded_sizes.h>
#include <ec.h>
#include <fp.h>
#include <fp2.h>

#define PASSED 0
#define FAILED 1

// Generating a pseudo-random field element in [0, p-1]
void fp_random_test(fp_t *a);

// Generating a pseudo-random element in GF(p^2)
void fp2_random_test(fp2_t *a);

// Generating a random projective x-only point
void ec_random_test(ec_point_t *P, const ec_curve_t *curve);

// Generating a random projective x-only point and normalizing it
void ec_random_normalized_test(ec_point_t *P, const ec_curve_t *curve);

// Point difference
void projective_difference_point(ec_point_t *PQ, const ec_point_t *P, const ec_point_t *Q, const ec_curve_t *curve);

// xDBL
void xDBL(ec_point_t *Q, const ec_point_t *P, const ec_point_t *AC);

// Double-and-add
extern void xDBLADD(ec_point_t *R,
                    ec_point_t *S,
                    const ec_point_t *P,
                    const ec_point_t *Q,
                    const ec_point_t *PQ,
                    const ec_point_t *A24,
                    const bool A24_normalized);

#endif
