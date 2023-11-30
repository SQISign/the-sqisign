#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')
    
################################################################

# The following benchmarks are required for computing optimal strategies.
# The costs can be specified in arbitrary units and can be obtained by
# runing src/ec/ref/lvl1/test/mont.test or left as "None" to use a
# generic estimate.
# ++++++++++ 
# p2 =  None   # cost of xdbl
# q4 = None    # cost of xeval4

##SqiSign lvl 1 costs
p2 = 1766
q4 = 2452

##SqiSign lvl 3 costs
# p2 = 4071
# q4 = 5224

##SqiSign lvl 5 costs
# p2 = 7755
# q4 = 9847
# ++++++++++

################################################################

from parameters import p, f, Tpls, Tmin
g = valuation(Tpls, 3)
T = Tpls * Tmin
Lpls = [x[0] for x in factor(Tpls)]
Lmin = [x[0] for x in factor(Tmin)]
Epls = [x[1] for x in factor(Tpls)]
Emin = [x[1] for x in factor(Tmin)]
L = Lpls + Lmin
E = Epls + Emin
################################################################

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')
if 3 not in Lpls or 3 in Lmin:
    raise NotImplementedError('power of 3 must be in the positive torsion')

################################################################

from math import log, ceil, floor, sqrt

def strategy(n, p, q):
    S = { 1: [] }
    C = { 1: 0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)),
                      key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n]

def ijk(l):
    i,j,k = ceil(sqrt((l-1)/4)), floor((l-1)/4/ceil(sqrt((l-1)/4))), int((l-1)/2 - 2*ceil(sqrt((l-1)/4))*floor((l-1)/4/ceil(sqrt((l-1)/4))))
    assert i < 2**16 and j < 2**16 and k < 2**16
    return i,j,k

bL = '{'+', '.join([str(ceil(log(l,2))) for l in L])+'}'
strategy4 = strategy(f//2-1, 2*p2, q4)
strategy4 = '{'+', '.join([str(s) for s in strategy4])+'}'
sizeI = []
sizeJ = []
sizeK = []
for l in L:
    I,J,K = ijk(l)
    sizeI.append(I)
    sizeJ.append(J)
    sizeK.append(K)
sImax = max(sizeI)
sJmax = max(sizeJ)
sKmax = max(sizeK)
sizeI = '{'+', '.join([str(s) for s in sizeI])+'}'
sizeJ = '{'+', '.join([str(s) for s in sizeJ])+'}'
sizeK = '{'+', '.join([str(s) for s in sizeK])+'}'

################################################################

FpEls = {
    'TWOpF':2**f,
    'TWOpFm1': 2**(f-1),
    'THREEpE': 3**(g//2),
    'THREEpF': 3**g,
    'THREEpFdiv2': (3**g)//2,
    'p_cofactor_for_2f': (p+1)//(2**f),
    'p_cofactor_for_3g': (p+1)//(3**g),
    'p_cofactor_for_6fg': (p+1)//(2**f*3**g)
}

from cformat import FpEl, ObjectStatic, ObjectFormatter
objs = ObjectFormatter(
    [ObjectStatic('digit_t', f'{k}[NWORDS_ORDER]', FpEl(v, p))
    for k,v in FpEls.items()]
    )

################################################################


with open('include/ec_params.h', 'w') as hfile:

    print('#ifndef EC_PARAMS_H', file=hfile)
    print('#define EC_PARAMS_H', file=hfile)
    print('',file=hfile)
    print('#include <fp_constants.h>', file=hfile)
    print('',file=hfile)
    print(f'#define POWER_OF_2 {f}', file=hfile)
    print(f'#define POWER_OF_3 {g}', file=hfile)
    print('',file=hfile)
    print('#define scaled 1', file=hfile)
    print('#define gap 83', file=hfile)
    print('',file=hfile)
    print(f'#define P_LEN {len(Lpls)}', file=hfile)
    print(f'#define M_LEN {len(Lmin)}', file=hfile)
    print('',file=hfile)
    print('static digit_t p_plus_minus_bitlength[P_LEN + M_LEN] =\n\t'+bL+';', file=hfile)
    print('',file=hfile)
    print(f'static digit_t STRATEGY4[] =\n\t'+strategy4+';', file=hfile)
    print('',file=hfile)
    print(f'static digit_t sizeI[] =\n\t'+sizeI+';', file=hfile)
    print(f'static digit_t sizeJ[] =\n\t'+sizeJ+';', file=hfile)
    print(f'static digit_t sizeK[] =\n\t'+sizeK+';', file=hfile)
    print('',file=hfile)
    print(f'#define sI_max {sImax}', file=hfile)
    print(f'#define sJ_max {sJmax}', file=hfile)
    print(f'#define sK_max {max(sKmax, 41)}', file=hfile)
    print('',file=hfile)
    print(f'#define ceil_log_sI_max {ceil(log(sImax,2))}', file=hfile)
    print(f'#define ceil_log_sJ_max {ceil(log(sJmax,2))}', file=hfile)
    print('',file=hfile)
    objs.implementation(file=hfile)
    print('',file=hfile)
    print(f'#define P_COFACTOR_FOR_2F_BITLENGTH {ceil(log((p+1)//(2**f),2))}', file=hfile)
    print(f'#define P_COFACTOR_FOR_3G_BITLENGTH {ceil(log((p+1)//(3**g),2))}', file=hfile)
    print(f'#define P_COFACTOR_FOR_6FG_BITLENGTH {ceil(log((p+1)//(2**f*3**g),2))}', file=hfile)
    print('',file=hfile)
    print('#endif', file=hfile)

