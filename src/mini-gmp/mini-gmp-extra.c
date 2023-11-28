#include <stddef.h>
#include "mini-gmp.h"
#include "mini-gmp-extra.h"
#include <rng.h>

void gmp_randinit_mt(gmp_randstate_ptr state) {
    (void)state;
}

void gmp_randclear(gmp_randstate_ptr state) {
    (void)state;
}

void gmp_randseed_ui(gmp_randstate_ptr state, unsigned long int seed) {
    (void)state;
    (void)seed;
}

void mpz_urandomb(mpz_ptr rop, gmp_randstate_ptr state, mp_bitcnt_t n) {
    (void)state;

    int num_limbs = (n + 8*sizeof(mp_limb_t) - 1)/(8*sizeof(mp_limb_t));
    n %= 8*sizeof(mp_limb_t);

    mp_limb_t* limbs = mpz_limbs_write(rop, num_limbs);

    randombytes((unsigned char*)limbs, num_limbs*sizeof(mp_limb_t));
    limbs[num_limbs - 1] &= ((mp_limb_t)1 << n) - 1;

    mpz_limbs_finish(rop, num_limbs);
}

void mpz_urandomm(mpz_ptr rop, gmp_randstate_ptr state, mpz_srcptr n) {
    (void)state;

    // Use rejection sampling on integers in [0, 2^log2(n) - 1]
    size_t n_bits = mpz_sizeinbase(n, 2);

    do {
        mpz_urandomb(rop, state, n_bits);
    } while (mpz_cmp(rop, n) >= 0);
}

int mpz_legendre(const mpz_t a, const mpz_t p) {
    int ret = -1;
    mpz_t tmp;

    mpz_init(tmp);

    mpz_sub_ui(tmp, p, 1);
    mpz_fdiv_q_2exp(tmp, tmp, 1);

    mpz_powm(tmp, a, tmp, p);

    if (mpz_cmp_ui(tmp, 1) == 0)
    {
        ret = 1;
    }
    else if (mpz_cmp_ui(tmp, 0) == 0)
    {
        ret = 0;
    }

    mpz_clear(tmp);

    return ret;
}

void mpz_nextprime(mpz_ptr rop, mpz_srcptr op) {
    // Note: this is only used for test code, so it's not performance critical. Should this function be used in the main
    // library code, some improvements are admissible, such as sieving by small primes.
    if (mpz_odd_p(op)) {
        mpz_add_ui(rop, rop, 2);
    } else {
        mpz_add_ui(rop, rop, 1);
    }

    while (!mpz_probab_prime_p(rop, 30)) {
        mpz_add_ui(rop, rop, 2);
    }
}
