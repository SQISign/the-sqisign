#include <stddef.h>
#include "mini-gmp.h"
#include "mini-gmp-extra.h"
#include "tutil.h"
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
    n = ((n - 1) % (8*sizeof(mp_limb_t))) + 1;

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
    // Based on https://eprint.iacr.org/2023/1261.pdf

    mpz_t f, g, tmp, tmp2;
    sdigit_t delta, fp, gp, mask;
    int i, j, t = 0, d, n, k;
    int a_bits = mpz_sizeinbase(a, 2);
    int p_bits = mpz_sizeinbase(p, 2);

    // d = max(a_bits, p_bits)
    d = a_bits - p_bits;
    d = a_bits - ((d >> (8*sizeof(int) - 1)) & d);

    n = (45907 * d + 26313)/19929;

    mpz_init_set(f, p);
    mpz_init_set(g, a);
    mpz_init(tmp);
    mpz_init(tmp2);

    // Unlike the paper, we do not check any of the conditions for the Kronecker symbol; this is just a Legendre symbol
    // calculation, and as such, p is assumed to be positive and an odd prime

    delta = 0;

    for (i = 0; i < (n + (RADIX - 2) - 1)/(RADIX - 2); i++) {
        sdigit_t u = 0, ai = 1, bi = 0, ci = 0, di = 1;

        fp = mpz_get_ui(f);
        gp = mpz_get_ui(g);

        fp *= mpz_sgn(f);
        gp *= mpz_sgn(g);

        for (j = 0; j < RADIX - 2; j++) {
            digit_t y = fp, m0, m1, d0, t0, t1, t2;

            d0 = ~(delta >> (RADIX - 1));
            m1 = -(gp & 1);
            m0 = d0 & m1;

            t0 = (fp ^ d0) - d0;
            t1 = (ci ^ d0) - d0;
            t2 = (di ^ d0) - d0;

            gp += t0 & m1;
            ai += t1 & m1;
            bi += t2 & m1;

            fp += gp & m0;
            ci += ai & m0;
            di += bi & m0;

            gp >>= 1;
            ci <<= 1;
            di <<= 1;

            delta = (delta ^ m0) + 1;

            u += ((y & fp) ^ (fp >> 1)) & 2;
            u += (u & 1) ^ -(ci >> (RADIX - 1));
        }

        mpz_mul_si(tmp, g, ai);
        mpz_mul_si(tmp2, f, bi);
        mpz_add(tmp2, tmp2, tmp);
        mpz_fdiv_q_2exp(tmp2, tmp2, RADIX - 2);

        mpz_mul_si(tmp, g, ci);
        mpz_mul_si(f, f, di);
        mpz_add(f, f, tmp);
        mpz_fdiv_q_2exp(f, f, RADIX - 2);

        mpz_set(g, tmp2);

        t = (t + u) & 3;
        t = (t + ((t & 1) ^ ((mpz_sgn(f) < 0) ? 1 : 0))) % 4;
    }

    t = (t + (t & 1)) & 3;

    mpz_clear(f);
    mpz_clear(g);
    mpz_clear(tmp);
    mpz_clear(tmp2);

    // fp' = abs(fp)
    mask = fp >> (RADIX - 1);
    fp = (fp + mask) ^ mask;

    // mask = (fp' == 1) ? 0 : -1;
    fp--;
    mask = (fp | -fp) >> (RADIX - 1);

    // return (fp' == 1) ? (1 - t) : 0;
    return (1 - t) & ~mask;
}

void mpz_nextprime(mpz_ptr rop, mpz_srcptr op) {
    // Note: this is only used for test code, so it's not performance critical. Should this function be used in the main
    // library code, some improvements are admissible, such as sieving by small primes.
    if (mpz_odd_p(op)) {
        mpz_add_ui(rop, op, 2);
    } else {
        mpz_add_ui(rop, op, 1);
    }

    while (!mpz_probab_prime_p(rop, 30)) {
        mpz_add_ui(rop, rop, 2);
    }
}
