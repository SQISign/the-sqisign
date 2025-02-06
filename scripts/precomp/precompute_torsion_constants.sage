#!/usr/bin/env sage
proof.all(False)  # faster

################################################################

from parameters import p

################################################################

tors2part = (p+1).p_primary_part(2)
lambda_security = round(p.bit_length() / 128) * 64
N_sec = next_prime(1 << 4*lambda_security)
N_com = N_sec

defs = {
        'TORSION_2POWER_BYTES': (tors2part.bit_length() + 7) // 8,
    }

from cformat import Ibz, Object, ObjectFormatter

objs = ObjectFormatter([
        Object('ibz_t', 'TWO_TO_SECURITY_BITS', Ibz(1 << lambda_security)),  # lambda_security = SECURITY_BITS (128, 192, 256)
        Object('ibz_t', 'TORSION_PLUS_2POWER', Ibz(tors2part)),
        Object('ibz_t', 'SEC_DEGREE', Ibz(N_sec)),
        Object('ibz_t', 'COM_DEGREE', Ibz(N_com)),
    ])

with open('include/torsion_constants.h','w') as hfile:
    with open('torsion_constants.c','w') as cfile:
        print(f'#include <quaternion.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <torsion_constants.h>', file=cfile)

        for k,v in defs.items():
            print(f'#define {k} {v}', file=hfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

