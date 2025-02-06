#include "quaternion_tests.h"
#include <stdlib.h>

int
quat_test_input_random_ideal_generation(quat_left_ideal_t *ideals,
                                        int norm_bitsize,
                                        int iterations,
                                        const quat_represent_integer_params_t *params)
{
    int randret = 0;
    ibz_t norm, cofactor;
    ibz_init(&norm);
    ibz_init(&cofactor);
    int size_p = ibz_bitsize(&(params->algebra->p));
    randret = !ibz_generate_random_prime(&cofactor,
                                         0,
                                         ((norm_bitsize + 5) > size_p) * (norm_bitsize + 5) +
                                             ((norm_bitsize + 5) <= size_p) * size_p,
                                         params->primality_test_iterations);
    if ((randret != 0) || ibz_is_zero(&cofactor)) {
        printf("Random generation failed in quat_test_input_random_ideal_generation\n");
        goto fin;
    }

    for (int iter = 0; iter < iterations; iter++) {
        // generate random odd norm
        randret = !ibz_rand_interval_bits(&norm, norm_bitsize - 1);
        ibz_add(&norm, &norm, &norm);
        ibz_add(&norm, &norm, &ibz_const_one);
        if ((randret != 0) || ibz_is_zero(&norm)) {
            printf("Random generation failed in quat_test_input_random_ideal_generation\n");
            goto fin;
        }
        ibz_abs(&norm, &norm);
        // compute ideal
        randret = !quat_sampling_random_ideal_O0_given_norm(&(ideals[iter]), &norm, 0, params, &cofactor);
        if (randret != 0) {
            printf("Random generation failed in quat_test_input_random_ideal_generation\n");
            goto fin;
        }
    }
fin:;
    ibz_finalize(&norm);
    ibz_finalize(&cofactor);
    return (randret);
}

// norms is either of length iterations or NULL
int
quat_test_input_random_ideal_lattice_generation(quat_lattice_t *lattices,
                                                ibz_t *norms,
                                                int norm_bitsize,
                                                int iterations,
                                                const quat_represent_integer_params_t *params)
{
    quat_left_ideal_t *ideals;
    ideals = malloc(iterations * sizeof(quat_left_ideal_t));
    for (int i = 0; i < iterations; i++)
        quat_left_ideal_init(&(ideals[i]));
    int randret = quat_test_input_random_ideal_generation(ideals, norm_bitsize, iterations, params);

    for (int iter = 0; iter < iterations; iter++) {
        ibz_mat_4x4_copy(&(lattices[iter].basis), &(ideals[iter].lattice.basis));
        ibz_copy(&(lattices[iter].denom), &(ideals[iter].lattice.denom));
        if (norms != NULL) {
            ibz_copy(&(norms[iter]), &(ideals[iter].norm));
        }
    }
    for (int i = 0; i < iterations; i++)
        quat_left_ideal_finalize(&(ideals[i]));
    free(ideals);
    return randret;
}

int
quat_test_input_random_lattice_generation(quat_lattice_t *lattices, int bitsize, int iterations, int in_hnf)
{
    ibz_t det;
    ibz_init(&det);
    int randret = 0;
    for (int iter = 0; iter < iterations; iter++) {
        // generate random invertible matrix
        ibz_set(&det, 0);
        while (ibz_is_zero(&det)) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    randret = !ibz_rand_interval_bits(&((lattices[iter]).basis[i][j]), bitsize);
                    if (randret != 0) {
                        printf("Random generation failed in "
                               "quat_test_input_random_lattice_generation\n");
                        goto fin;
                    }
                }
            }
            randret = !ibz_rand_interval_bits(&((lattices[iter]).denom), bitsize);
            if (randret != 0) {
                printf("Random generation failed in quat_test_input_random_lattice_generation\n");
                goto fin;
            }
            ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &((lattices[iter]).basis));
            ibz_mul(&det, &det, &(lattices[iter].denom));
            if (in_hnf && !ibz_is_zero(&det))
                quat_lattice_hnf(&(lattices[iter]));
        }
    }
fin:;
    ibz_finalize(&det);
    return (randret);
}
