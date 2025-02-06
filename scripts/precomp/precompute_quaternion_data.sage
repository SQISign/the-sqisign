#!/usr/bin/env sage
proof.all(False)  # faster



from maxorders import p, orders

from cformat import Ibz, Object, ObjectFormatter

# Prime of same size than p for random ideal of fixed norm
bitlength_p = int(p).bit_length()
prime_cofactor = next_prime((2^(bitlength_p)))

algobj = [Ibz(p)]

objs = \
    [
        [
            # basis (columns)
            [
                Ibz(mat.denominator()),
                [[Ibz(v) for v in vs]
                    for vs in mat.transpose()*mat.denominator()],
            ],
            # sqrt(-q)
            [
                Ibz(mat.denominator()),
                [Ibz(c) for c in ii*mat.denominator()],
            ],
            # sqrt(-p)
            [
                Ibz(1),
                [Ibz(c) for c in (0,0,1,0)]
            ],
            q
        ]
        for q,_,mat,ii,_,_ in orders
    ]

idlobjs = \
    [
        [
            # basis (columns)
            [
                Ibz(idl.denominator()),
                [[Ibz(v) for v in vs]
                    for vs in idl.transpose()*idl.denominator()],
            ],
            # norm
            Ibz(abs(idl.row_space(ZZ).intersection((ZZ^4).submodule([[1,0,0,0]])).basis()[0][0])),
            # left order
            '&MAXORD_O0',
        ]
        for _,_,mat,_,idl,_ in orders
    ]

gammaobjs = \
    [
        [
            Ibz(gamma.denominator()),
            list(map(Ibz, gamma * gamma.denominator())),
        ]
        for _,_,_,_,_,gamma in orders
    ]

objs = ObjectFormatter([
        Object('ibz_t', 'QUAT_prime_cofactor', Ibz(prime_cofactor)),
        Object('quat_alg_t', 'QUATALG_PINFTY', algobj),
        Object('quat_p_extremal_maximal_order_t[]', 'EXTREMAL_ORDERS', objs),
        Object('quat_left_ideal_t[]', 'CONNECTING_IDEALS', idlobjs),  # ideal corresponding to an isogeny from E0 which acts as identity w.r.t. the basis_even
        Object('quat_alg_elem_t[]', 'CONJUGATING_ELEMENTS', gammaobjs), # elements γ such that each I has right order γ O₁ γ^-1
    ])

with open('include/quaternion_data.h','w') as hfile:
    with open('quaternion_data.c','w') as cfile:
        print(f'#include <quaternion.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <quaternion_data.h>', file=cfile)

        #FIXME this should eventually go away?
        print(f'#define MAXORD_O0 (EXTREMAL_ORDERS->order)', file=hfile)
        print(f'#define STANDARD_EXTREMAL_ORDER (EXTREMAL_ORDERS[0])', file=hfile)
        print(f'#define NUM_ALTERNATE_EXTREMAL_ORDERS {len(orders)-1}', file=hfile)
        print(f'#define ALTERNATE_EXTREMAL_ORDERS (EXTREMAL_ORDERS+1)', file=hfile)
        print(f'#define ALTERNATE_CONNECTING_IDEALS (CONNECTING_IDEALS+1)', file=hfile)
        print(f'#define ALTERNATE_CONJUGATING_ELEMENTS (CONJUGATING_ELEMENTS+1)', file=hfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

