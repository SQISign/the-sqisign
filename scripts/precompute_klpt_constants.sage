#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(9, 8, print_message=True):
    exit('')

################################################################

from parameters import f, p, Tpls, Tmin
negl = 2**-64   #TODO optimize

################################################################

logp = ceil(log(p, 2))
logT = ceil(log(Tpls*Tmin, 2))
tors2val = (p+1).valuation(2)

defs = dict()

# lideal_equiv

defs['KLPT_equiv_bound_coeff'] = ceil((log(negl, 1-2/logp) ** (1/4) - 1) / 2) + 2
assert (1 - 2/logp) ** ((2 * defs['KLPT_equiv_bound_coeff'] + 1) ** 4) <= negl

defs['KLPT_equiv_num_iter'] = (2 * defs['KLPT_equiv_bound_coeff'] + 1) ** 4

defs['KLPT_primality_num_iter'] = ceil(-log(negl, 4))

# signing KLPT

defs['KLPT_signing_klpt_length'] = f * ceil (ceil((log(negl, 2) / -2) + 15/4*logp + 25)/f)
assert 2**(-2 * (defs['KLPT_signing_klpt_length'] - 15/4*logp - 25)) <= negl

defs['KLPT_signing_num_gamma_trial'] = ceil(log(negl, 2) / -1)
assert 2 ** ( - defs['KLPT_signing_num_gamma_trial']) <= negl

defs['KLPT_gamma_exponent_interval_size'] = 0

defs['KLPT_gamma_exponent_center_shift'] = ceil(log(log(negl, 1-1/logp) + defs['KLPT_signing_num_gamma_trial'], 2) + defs['KLPT_gamma_exponent_interval_size'])
assert (1 - 1/logp) ** (2**(defs['KLPT_gamma_exponent_center_shift'] - defs['KLPT_gamma_exponent_interval_size']) - defs['KLPT_signing_num_gamma_trial']) <= negl

defs['KLPT_repres_num_gamma_trial'] = 2**(defs['KLPT_gamma_exponent_center_shift'] + defs['KLPT_gamma_exponent_interval_size'])

defs['KLPT_signing_number_strong_approx'] = ceil(log(1/64, 1-4/13/logp))
assert (1 - 4/13/logp) ** defs['KLPT_signing_number_strong_approx'] <= 1/64

# keygen KLPT

defs['KLPT_random_prime_attempts'] = 64

defs['KLPT_secret_key_prime_size'] = ceil(logp / 4)

defs['KLPT_keygen_length'] =   f* ceil ( ceil(log(negl, 2) / -2 + 5/2*logp -25 ) / f)
assert 2 ** (-2 * (defs['KLPT_keygen_length'] - 5/2*logp +25)) <= negl

defs['KLPT_keygen_num_gamma_trial'] = ceil(log(negl, 2) / -1)

defs['KLPT_eichler_smallnorm_bitsize'] = ceil(1/2*logp - 4/3*( logT - 5/4*logp))

defs['KLPT_keygen_number_strong_approx'] = ceil(log(1/64, 1-2/5/logp))
assert (1 - 2/5/logp) ** defs['KLPT_keygen_number_strong_approx'] <= 1/64

# Eichler

defs['KLPT_eichler_number_mu_norm'] = ceil((logT - 5/4*logp) / log(3,2))

defs['KLPT_eichler_strong_approx_log_margin'] = 2

defs['KLPT_eichler_num_equiv_ideal'] = ceil(logp / 10)

defs['KLPT_eichler_number_strong_approx'] = ceil(10 * logp)

# signature response

defs['SQISIGN_response_attempts'] = 64

# signature isogeny degrees

defs['SQISIGN_random_length'] = 0
defs['SQISIGN_signing_total_length'] = defs['KLPT_signing_klpt_length']
defs['SQISIGN_signing_length'] = ZZ(defs['SQISIGN_signing_total_length'] / tors2val)
defs['SQISIGN_keygen_length'] = ZZ(defs['KLPT_keygen_length'] / tors2val)

# prime data for Cornacchia

primes_1mod4 = [p for p in primes(100) if p%4==1]
prod_primes_3mod4 = prod(p for p in primes(100) if p%4==3)

################################################################

from cformat import Ibz, Object, ObjectFormatter

objs = ObjectFormatter([
        Object('short[]', 'SMALL_PRIMES_1MOD4', [int(v) for v in primes_1mod4]),
        Object('ibz_t', 'PROD_SMALL_PRIMES_3MOD4', Ibz(prod_primes_3mod4)),
    ])

################################################################

with open('include/klpt_constants.h','w') as hfile:
    with open('klpt_constants.c','w') as cfile:
        print(f'#include <intbig.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <klpt_constants.h>', file=cfile)

        for k,v in defs.items():
            v = ZZ(v)
            print(f'#define {k} {v}', file=hfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

