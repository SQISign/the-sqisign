#ifndef TEDWARDS_H
#define TEDWARDS_H

#include <fp2.h>
#include "ec.h"

// a*x^2+y^2=1+d*x^2*y^2

typedef struct ted_point_t {
    fp2_t x;
    fp2_t y;
    fp2_t z;
    fp2_t t; // t = x*y/z
} ted_point_t;

void ted_init(ted_point_t* P);
bool is_ted_equal(ted_point_t const* P1, ted_point_t const* P2);
void copy_ted_point(ted_point_t* P, ted_point_t const* Q);

void ted_neg(ted_point_t* Q, ted_point_t const* P);
void ted_dbl(ted_point_t* Q, ted_point_t const* P, ec_curve_t const* E);
void ted_add(ted_point_t* S, ted_point_t const* P, ted_point_t const* Q, ec_curve_t const* E);

void mont_to_ted(ec_curve_t* E, ec_curve_t const* A);
void mont_to_ted_point(ted_point_t* Q, ec_point_t const* P, ec_curve_t const* A);
void ted_to_mont_point(ec_point_t* Q, ted_point_t const* P);

#endif
