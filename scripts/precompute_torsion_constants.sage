#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(9, 8, print_message=True):
    exit('')

################################################################

from parameters import p, B, f, Tpls, Tmin, Dcom, Dchall

################################################################

Lpls = sorted(set(Tpls.prime_factors()) - {2})
Epls = [Tpls.valuation(l) for l in Lpls]

Lmin = sorted(set(Tmin.prime_factors()) - {2})
Emin = [Tmin.valuation(l) for l in Lmin]

tors2part = (p+1).p_primary_part(2)
tors3part = (p+1).p_primary_part(3)
tors23part = tors2part * tors3part

defs = {
        'TORSION_2POWER_BYTES': (int(tors2part).bit_length() + 7) // 8,
        'TORSION_3POWER_BYTES': (int(tors3part).bit_length() + 7) // 8,
        'TORSION_23POWER_BYTES': (int(tors23part).bit_length() + 7) // 8,
    }

from cformat import Ibz, Object, ObjectFormatter

objs = ObjectFormatter([
        Object('digit_t', 'TORSION_PLUS_EVEN_POWER', int(f)),
        Object('digit_t[]', 'TORSION_ODD_PRIMES', Lpls + Lmin),
        Object('digit_t[]', 'TORSION_ODD_POWERS', Epls + Emin),
        Object('digit_t[]', 'TORSION_PLUS_ODD_PRIMES', Lpls),      # TODO deduplicate?
        Object('size_t[]', 'TORSION_PLUS_ODD_POWERS', Epls),        # TODO deduplicate?
        Object('digit_t[]', 'TORSION_MINUS_ODD_PRIMES', Lmin),     # TODO deduplicate?
        Object('size_t[]', 'TORSION_MINUS_ODD_POWERS', Emin),       # TODO deduplicate?
        Object('size_t[]', 'DEGREE_COMMITMENT_POWERS', [Dcom.valuation(l) for l in Lpls+Lmin]), #FIXME should be ec_degree_odd_t
        Object('ibz_t', 'CHARACTERISTIC', Ibz(p)),
        Object('ibz_t', 'TORSION_ODD', Ibz(Tpls * Tmin)),
        Object('ibz_t[]', 'TORSION_ODD_PRIMEPOWERS', [Ibz(l^e) for Tpm in (Tpls,Tmin) for l,e in Tpm.factor()]),
        Object('ibz_t', 'TORSION_ODD_PLUS', Ibz(Tpls)),
        Object('ibz_t', 'TORSION_ODD_MINUS', Ibz(Tmin)),
        Object('ibz_t', 'TORSION_PLUS_2POWER', Ibz(tors2part)),
        Object('ibz_t', 'TORSION_PLUS_3POWER', Ibz(tors3part)),
        Object('ibz_t', 'TORSION_PLUS_23POWER', Ibz(tors23part)),
        Object('ibz_t', 'DEGREE_COMMITMENT', Ibz(Dcom)),
        Object('ibz_t', 'DEGREE_COMMITMENT_PLUS', Ibz(gcd(Dcom, Tpls))),
        Object('ibz_t', 'DEGREE_COMMITMENT_MINUS', Ibz(gcd(Dcom, Tmin))),
        Object('ibz_t', 'DEGREE_CHALLENGE', Ibz(Dchall)),
    ])

with open('include/torsion_constants.h','w') as hfile:
    with open('torsion_constants.c','w') as cfile:
        print(f'#include <intbig.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <torsion_constants.h>', file=cfile)

        for k,v in defs.items():
            print(f'#define {k} {v}', file=hfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

