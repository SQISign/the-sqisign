#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')

################################################################

from parameters import p

################################################################

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')

################################################################

from math import ceil, log

with open('include/fp_constants.h', 'w') as hfile:

    print('#ifndef FP_CONSTANTS_H', file=hfile)
    print('#define FP_CONSTANTS_H', file=hfile)
    print('',file=hfile)
    print('#if 0',file=hfile)
    for RADIX in (16,32,64):
        print('',file=hfile)
        print(f'#elif 8*DIGIT_LEN == {RADIX}',file=hfile)
        print(f'#define NWORDS_FIELD {ceil(log(p,2**RADIX))}',file=hfile)
        print(f'#define NWORDS_ORDER {ceil(log(p+1,2**RADIX))}',file=hfile)
        print(f'#define BITS {RADIX*ceil(log(p,2**RADIX))}',file=hfile)
        print(f'#define LOG2P {ceil(log(RADIX*ceil(log(p,2**RADIX)),2))}',file=hfile)
    print('',file=hfile)
    print('#endif', file=hfile)
    print('',file=hfile)
    print('#endif', file=hfile)

