#ifndef CURVE_EXTRAS_H
#define CURVE_EXTRAS_H

#include "ec.h"
#include "torsion_constants.h"

typedef struct jac_point_t {
    fp2_t x;
    fp2_t y;
    fp2_t z;
} jac_point_t;

bool ec_is_zero(ec_point_t const* P);
void copy_point(ec_point_t* P, ec_point_t const* Q);
void swap_points(ec_point_t* P, ec_point_t* Q, const digit_t option);
void ec_init(ec_point_t* P);
void xDBLv2(ec_point_t* Q, ec_point_t const* P, ec_point_t const* A24);
void xDBLADD(ec_point_t* R, ec_point_t* S, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24);
void xDBLMUL(ec_point_t* S, ec_point_t const* P, digit_t const* k, ec_point_t const* Q, digit_t const* l, ec_point_t const* PQ, ec_curve_t const* curve);
void xDBL(ec_point_t* Q, ec_point_t const* P, ec_point_t const* AC);
void xMUL(ec_point_t* Q, ec_point_t const* P, digit_t const* k, ec_curve_t const* curve);
void xDBLMUL(ec_point_t* S, ec_point_t const* P, digit_t const* k, ec_point_t const* Q, digit_t const* l, ec_point_t const* PQ, ec_curve_t const* curve);

#define is_point_equal ec_is_equal
#define xADD ec_add

#endif

