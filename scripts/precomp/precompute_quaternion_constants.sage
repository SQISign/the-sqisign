#!/usr/bin/env sage
proof.all(False)  # faster

################################################################

from parameters import p
negl = 2**-64

################################################################

logp = ceil(log(p, 2))
loglogp = ceil(log(logp,2))
tors2val = (p+1).valuation(2)

defs = dict()

# RepresentInteger data
small = ceil(log(negl, 2) / -1)
assert 2**-small <= negl

add_shift = ceil(log(log(negl, 1-1/(64*logp)), 2))
assert (1 - 1/(64*logp)) ** (2**(add_shift)) <= negl

defs['QUAT_primality_num_iter'] = ceil(-log(negl, 4))
defs['QUAT_repres_bound_input'] = add_shift

# Equivalent ideal data
defs['QUAT_equiv_bound_coeff'] = 2**(1 + add_shift//4)

# Find_uv constants
m = 2 + floor((logp - tors2val) / 4)
defs['FINDUV_box_size'] = m
defs['FINDUV_cube_size'] = (2 * m + 1)**4 - 1

################################################################

with open('include/quaternion_constants.h','w') as hfile:
    print(f'#include <quaternion.h>', file=hfile)

    for k,v in defs.items():
        v = ZZ(v)
        print(f'#define {k} {v}', file=hfile)


