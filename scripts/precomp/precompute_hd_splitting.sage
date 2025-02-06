#!/usr/bin/env sage
proof.all(False)  # faster

################################################################

from parameters import p

# Field
Fp2.<i> = GF((p,2), modulus=[1,0,1])

Fp2_constants = [
    [Fp2(0), Fp2(1), Fp2(i), Fp2(-1), Fp2(-i)],
    ["FP2_ZERO", "FP2_ONE", "FP2_I", "FP2_MINUS_ONE", "FP2_MINUS_I"]
]

################################################################

from cformat import FpEl

def Fp2_to_list(el):
    return [FpEl(int(c), p, True) for c in Fp2(el)]

def Fp2_to_name(el):
    return Fp2_constants[1][Fp2_constants[0].index(el)]

################################################################

# Splitting Data

chi_eval = [
    [1,1,1,1],
    [1,-1,1,-1],
    [1,1,-1,-1],
    [1,-1,-1,1]
]

even_indices = [
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [1, 0],
    [1, 2],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 3],
]

splitting_map = {
    (0, 2): [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]],
    (3, 3): [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    (0, 3): [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]],
    (2, 1): [[1, 1, 1, 1], [1, -1, 1, -1], [1, -1, -1, 1], [1, 1, -1, -1]],
    (0, 1): [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, -1, 0, 0]],
    (1, 2): [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
    (2, 0): [[1, 1, 1, 1], [1, -1, 1, -1], [1, -1, -1, 1], [-1, -1, 1, 1]],
    (3, 0): [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [-1, 1, 1, -1]],
    (1, 0): [[1, 1, 1, 1], [1, -1, -1, 1], [1, 1, -1, -1], [-1, 1, -1, 1]],
    (0, 0): [[1, i, 1, i], [1, -i, -1, i], [1, i, -1, -i], [-1, i, -1, i]],
}

# 24 2x2 maps used for normalization. Applying a uniform matrix to a level 2
# theta null point ensure that we get a random theta null point in the
# equivalence class of Γ/Γ(2,4)
full_normalization_maps_2d = [
    matrix(2, 2, [1, 0, 0, 1]),
    matrix(2, 2, [0, 1, 1, 0]),
    matrix(2, 2, [1, 0, 0, -1]),
    matrix(2, 2, [0, 1, -1, 0]),

    matrix(2, 2, [1, 1, 1, -1]),
    matrix(2, 2, [1, -1, 1, 1]),
    matrix(2, 2, [1, 1, -1, 1]),
    matrix(2, 2, [-1, 1, 1, 1]),

    matrix(2, 2, [1, 0, 0, i]),
    matrix(2, 2, [0, i, 1, 0]),
    matrix(2, 2, [1, 0, 0, -i]),
    matrix(2, 2, [0, i, -1, 0]),

    matrix(2, 2, [i, 1, 1, i]),
    matrix(2, 2, [1, i, i, 1]),
    matrix(2, 2, [i, 1, -1, -i]),
    matrix(2, 2, [1, i, -i, -1]),

    matrix(2, 2, [1, i, 1, -i]),
    matrix(2, 2, [1, -i, 1, i]),
    matrix(2, 2, [1, i, -1, i]),
    matrix(2, 2, [1, -i, -1, -i]),

    matrix(2, 2, [1, 1, i, -i]),
    matrix(2, 2, [i, -i, 1, 1]),
    matrix(2, 2, [1, 1, -i, i]),
    matrix(2, 2, [i, -i, -1, -1]),
]

# 6 2x2 maps used for normalisation. A subset of the preceding 24 matrices,
# that are sufficient to ensure uniform normalisation (under Γ/Γ^0(4))
# when using the Montgomery model
normalization_maps_2d = [
    matrix(2, 2, [1, 0, 0, 1]),
    matrix(2, 2, [0, 1, 1, 0]),
    matrix(2, 2, [1, 1, 1, -1]),
    matrix(2, 2, [-1, 1, 1, 1]),
    matrix(2, 2, [i, 1, 1, i]),
    matrix(2, 2, [1, i, i, 1]),
]


# Format from dictionary to list of lists
splitting_matrices = []
for ind in even_indices:
    # Create a list of all ten matrices (represented as 4x4 matrices)
    splitting_matrices.append([[Fp2_to_name(Fp2(x)) for x in row] for row in splitting_map[tuple(ind)]])

# 6 4x4 maps constructed from the above
normalization_maps_4d = []
for m in normalization_maps_2d:
    M = m.tensor_product(m, subdivide=False).list()
    matrix_elements_list = list(map(Fp2, M))
    # Reshape into matrix
    matrix_elements = [ matrix_elements_list[i:i+4] for i in range(0, len(matrix_elements_list), 4) ]
    normalization_maps_4d.append([[Fp2_to_name(x) for x in row] for row in matrix_elements])

################################################################

from cformat import Object, ObjectFormatter

objs = ObjectFormatter(
    [
        Object('int[][]', 'EVEN_INDEX', even_indices),
        Object('int[][]', 'CHI_EVAL', chi_eval),
        Object('fp2_t[]', 'FP2_CONSTANTS', list(map(Fp2_to_list, Fp2_constants[0]))),
        Object('precomp_basis_change_matrix_t[]', 'SPLITTING_TRANSFORMS', [[x] for x in splitting_matrices]),
        Object('precomp_basis_change_matrix_t[]', 'NORMALIZATION_TRANSFORMS', [[x] for x in normalization_maps_4d]),
    ]
)

with open("include/hd_splitting_transforms.h", "w") as hfile:
    with open("hd_splitting_transforms.c", "w") as cfile:
        print("#ifndef HD_SPLITTING_H", file=hfile)
        print("#define HD_SPLITTING_H", file=hfile)
        print(f"\n#include <hd.h>", file=hfile)
        print(f"#include <stdint.h>\n", file=hfile)
        print("typedef struct precomp_basis_change_matrix {", file=hfile)
        print("    uint8_t m[4][4];", file=hfile)
        print("} precomp_basis_change_matrix_t;\n", file=hfile)

        print(f"#include <hd_splitting_transforms.h>\n", file=cfile)
        for i in range(len(Fp2_constants[1])):
            print(f"#define {Fp2_constants[1][i]} {i}", file=cfile)
        print("", file=cfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

        print("\n#endif\n", file=hfile)
