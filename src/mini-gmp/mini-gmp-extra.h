#ifndef __MINI_GMP_EXTRA_H__
#define __MINI_GMP_EXTRA_H__

#include "mini-gmp.h"

// Random number generation functions
typedef void* gmp_randstate_t;
typedef gmp_randstate_t *gmp_randstate_ptr;
typedef const gmp_randstate_t *gmp_randstate_srcptr;

void gmp_randinit_mt(gmp_randstate_ptr);
void gmp_randclear(gmp_randstate_ptr);
void gmp_randseed_ui(gmp_randstate_ptr, unsigned long int);

void mpz_urandomb(mpz_ptr, gmp_randstate_ptr, mp_bitcnt_t);
void mpz_urandomm(mpz_ptr, gmp_randstate_ptr, mpz_srcptr);

// Number-theoretic functions
int mpz_legendre(const mpz_t a, const mpz_t p);
void mpz_nextprime(mpz_ptr, mpz_srcptr);

#endif
