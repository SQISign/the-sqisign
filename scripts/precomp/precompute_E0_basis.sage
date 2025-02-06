#!/usr/bin/env sage
proof.all(False)  # faster

################################################################

from parameters import p, f

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')

assert (1 << f).divides(p + 1)
Fp2.<i> = GF((p,2), modulus=[1,0,1])
E0 = EllipticCurve(Fp2, [1, 0])

from torsion_basis import even_torsion_basis_E0
P, Q = even_torsion_basis_E0(E0, f)

################################################################

from cformat import FpEl, Object, ObjectFormatter

def Fp2_to_list(el):
    return [FpEl(int(c), p, True) for c in Fp2(el)]

objs = ObjectFormatter([
        Object('fp2_t', 'BASIS_E0_PX', Fp2_to_list(P.x())),
        Object('fp2_t', 'BASIS_E0_QX', Fp2_to_list(Q.x())),
    ])

################################################################

with open('include/e0_basis.h','w') as hfile:
    with open('e0_basis.c','w') as cfile:
        print(f'#include <fp2.h>', file=hfile)
        print(f'#include <e0_basis.h>', file=cfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)
