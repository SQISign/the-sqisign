#!/usr/bin/env python3

from sage.all import *
from parameters import p, f

if __name__ == '__main__':

    cof = (p+1)//(2**f)

    from cformat import Object, ObjectFormatter

    obj_cof = ObjectFormatter(
        [
            Object('digit_t[]', 'p_cofactor_for_2f', [cof]),
        ]
    )

    with open("include/ec_params.h", "w") as hfile:
        with open("ec_params.c", "w") as cfile:

            hfile.write('#ifndef EC_PARAMS_H\n')
            hfile.write('#define EC_PARAMS_H\n')
            hfile.write('\n')

            hfile.write('#include <fp.h>\n')
            cfile.write('#include <ec_params.h>\n')
            hfile.write('\n')

            hfile.write(f'#define TORSION_EVEN_POWER {f}\n')
            hfile.write('\n')

            hfile.write('// p+1 divided by the power of 2\n')
            cfile.write('// p+1 divided by the power of 2\n')
            obj_cof.header(file=hfile)
            obj_cof.implementation(file=cfile)
            hfile.write(f'#define P_COFACTOR_FOR_2F_BITLENGTH {((p+1)//(2**f)).bit_length()}\n')
            hfile.write('\n')
            cfile.write('\n')


            hfile.write('#endif\n')
