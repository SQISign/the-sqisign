#include "quaternion_tests.h"
#include <stdlib.h>

// int quat_lattice_bound_parallelogram(ibz_vec_4_t *box, ibz_mat_4x4_t *U, const ibz_mat_4x4_t *G,
// const ibz_t *radius);
int
quat_test_lat_ball_paralellogram_randomized(int iterations, int bitsize)
{
    int res = 0;

    quat_lattice_t lattice;
    ibz_t radius, length, tmp;
    ibz_mat_4x4_t U, G;
    ibz_vec_4_t box, dbox, x;
    quat_lattice_init(&lattice);
    ibz_vec_4_init(&box);
    ibz_vec_4_init(&dbox);
    ibz_vec_4_init(&x);
    ibz_mat_4x4_init(&U);
    ibz_mat_4x4_init(&G);
    ibz_init(&radius);
    ibz_init(&length);
    ibz_init(&tmp);

    for (int it = 0; it < iterations; it++) {
// Create a random positive definite quadatic form
#ifndef NDEBUG
        int randret = quat_test_input_random_lattice_generation(&lattice, bitsize, 1, 0);
        assert(randret == 0);
#else
        (void)quat_test_input_random_lattice_generation(&lattice, bitsize, 1, 0);
#endif
        ibz_mat_4x4_transpose(&G, &lattice.basis);
        ibz_mat_4x4_mul(&G, &G, &lattice.basis);

// Set radius to 2 × sqrt(lattice volume)
#ifndef NDEBUG
        int ok = ibz_mat_4x4_inv_with_det_as_denom(NULL, &radius, &(lattice.basis));
        assert(ok);
#else
        (void)ibz_mat_4x4_inv_with_det_as_denom(NULL, &radius, &(lattice.basis));
#endif
        ibz_abs(&radius, &radius);
        ibz_sqrt_floor(&radius, &radius);
        ibz_mul(&radius, &radius, &ibz_const_two);

        quat_lattice_bound_parallelogram(&box, &U, &G, &radius);
        for (int i = 0; i < 4; i++) {
            // dbox is a box with sides dbox[i] =  2*box[i] + 1
            ibz_add(&dbox[i], &box[i], &box[i]);
            ibz_add(&dbox[i], &dbox[i], &ibz_const_one);
            // initialize x[i] to the bottom of dbox[i]
            ibz_neg(&x[i], &dbox[i]);
        }

        // Integrate U into the Gram matrix
        ibz_mat_4x4_mul(&G, &U, &G);
        ibz_mat_4x4_transpose(&U, &U);
        ibz_mat_4x4_mul(&G, &G, &U);

        // We treat x[0]...x[4] as a counter, incrementing one by one
        // but skipping values that are inside the parallelogram defined
        // by box.
        while (1) {
            // x is out of the parallelogram, so its length must be
            // greater than the radius
            quat_qf_eval(&length, &G, &x);
            if (ibz_cmp(&length, &radius) <= 0) {
                res = 1;
                break;
            }

            // Increment counter
            ibz_add(&x[0], &x[0], &ibz_const_one);
            // if x[0] just entered the interval
            ibz_add(&tmp, &x[0], &box[0]);
            if (ibz_is_zero(&tmp)) {
                int inbox = 1;
                for (int i = 1; i < 4; i++) {
                    ibz_abs(&tmp, &x[i]);
                    inbox &= ibz_cmp(&tmp, &box[i]) <= 0;
                }
                // if x[1]...x[3] are all in the respective intervals
                // jump straight to the end of x[0]'s interval
                if (inbox)
                    ibz_set(&x[0], 1);
            }

            // if x[0] became positive, loop the counter
            if (ibz_is_one(&x[0])) {
                ibz_neg(&x[0], &dbox[0]);
                int carry = 1;
                for (int i = 1; carry && i < 4; i++) {
                    ibz_add(&x[i], &x[i], &ibz_const_one);
                    if (ibz_cmp(&x[i], &dbox[i]) > 0) {
                        ibz_neg(&x[i], &dbox[i]);
                    } else {
                        carry = 0;
                    }
                }

                // If there still is a carry, we are at the end
                if (carry)
                    break;
            }
        }
    }

    quat_lattice_finalize(&lattice);
    ibz_vec_4_finalize(&box);
    ibz_vec_4_finalize(&dbox);
    ibz_vec_4_finalize(&x);
    ibz_mat_4x4_finalize(&U);
    ibz_mat_4x4_finalize(&G);
    ibz_finalize(&radius);
    ibz_finalize(&length);
    ibz_finalize(&tmp);

    if (res != 0) {
        printf("Quaternion unit test lat_ball_paralellogram_randomized failed\n");
    }
    return (res);
}

// helper which tests quat_lattice_sample_from_ball on given input
int
quat_test_lat_ball_sample_helper(const quat_lattice_t *lat, const ibz_t *radius, const quat_alg_t *alg)
{
    int res = 0;
    quat_alg_elem_t vec;
    ibz_t norm_d, norm_n;
    ibz_init(&norm_d);
    ibz_init(&norm_n);
    quat_alg_elem_init(&vec);
    // check return value
    res |= !quat_lattice_sample_from_ball(&vec, lat, alg, radius);
    // check result is in lattice
    res |= !quat_lattice_contains(NULL, lat, &vec);
    quat_alg_norm(&norm_n, &norm_d, &vec, alg);
    // test that n/d <= r so that n <= rd
    ibz_mul(&norm_d, &norm_d, radius);
    res |= !(ibz_cmp(&norm_n, &norm_d) <= 0);

    ibz_finalize(&norm_d);
    ibz_finalize(&norm_n);
    quat_alg_elem_finalize(&vec);
    return res;
}

// int quat_lattice_sample_from_ball(ibz_vec_4_t *x, const quat_lattice_t *lattice, const quat_alg_t
// *alg, const ibz_t *radius);
int
quat_test_lat_ball_sample_from_ball()
{
    int res = 0;

    ibz_t norm_n, norm_d;
    quat_alg_t alg;
    quat_alg_elem_t vec;
    quat_lattice_t lattice;
    ibz_t radius;

    quat_alg_init_set_ui(&alg, 11);
    ibz_init(&norm_n);
    ibz_init(&norm_d);
    quat_lattice_init(&lattice);
    ibz_init(&radius);
    quat_alg_elem_init(&vec);

    for (int it = 0; it < 3; it++) {
        if (it == 0) {
            ibz_mat_4x4_identity(&(lattice.basis)); // Test inner product
        } else if (it == 1) {
            ibz_set(&lattice.denom, 13); // Test with denominator
        } else {                         // if (it == 2)
            // Test a very much non-orthogonal qf
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    ibz_add(&(lattice.basis[i][j]), &(lattice.basis[i][j]), &ibz_const_one);
        }

        for (int i = 0; i < 100; i++) {
            ibz_set(&radius, i + 1);
            res |= quat_test_lat_ball_sample_helper(&lattice, &radius, &alg);
            if (res != 0)
                break;
        }
    }

    if (res != 0) {
        printf("Quaternion unit test lat_ball_sample_from_ball failed\n");
    }

    quat_alg_finalize(&alg);
    quat_lattice_finalize(&lattice);
    ibz_finalize(&radius);
    quat_alg_elem_finalize(&vec);
    ibz_finalize(&norm_n);
    ibz_finalize(&norm_d);
    return (res);
}

int
quat_test_lat_ball_sample_from_ball_randomized(int iterations, int bitsize)
{
    int res = 0;

    ibz_t norm_n, norm_d;
    quat_alg_t alg;
    quat_alg_elem_t vec;
    quat_lattice_t *lattices;
    ibz_t radius;

    quat_alg_init_set_ui(&alg, 11);
    ibz_init(&norm_n);
    ibz_init(&norm_d);
    lattices = malloc(iterations * sizeof(quat_lattice_t));
    for (int i = 0; i < iterations; i++)
        quat_lattice_init(&(lattices[i]));
    ibz_init(&radius);
    quat_alg_elem_init(&vec);

    int randret = quat_test_input_random_lattice_generation(lattices, bitsize, iterations, 0);

    if (!randret) {
        for (int i = 0; i < iterations; i++) {
#ifndef NDEBUG

            int ok = ibz_mat_4x4_inv_with_det_as_denom(NULL, &radius, &(lattices[i].basis));
            assert(ok);
#else
            (void)ibz_mat_4x4_inv_with_det_as_denom(NULL, &radius, &(lattices[i].basis));
#endif
            ibz_abs(&radius, &radius);
            res |= quat_test_lat_ball_sample_helper(&(lattices[i]), &radius, &alg);
            if (res != 0)
                break;
        }
    }

    if (res != 0) {
        printf("Quaternion unit test lat_ball_sample_from_ball_randomized failed\n");
    }

    quat_alg_finalize(&alg);
    for (int i = 0; i < iterations; i++)
        quat_lattice_finalize(&(lattices[i]));
    free(lattices);
    ibz_finalize(&radius);
    quat_alg_elem_finalize(&vec);
    ibz_finalize(&norm_n);
    ibz_finalize(&norm_d);
    return (res);
}

// run all previous tests
int
quat_test_lat_ball(void)
{
    int res = 0;
    printf("\nRunning quaternion tests for sampling lattice points of bounded norm\n");
    res = res | quat_test_lat_ball_sample_from_ball();
    res = res | quat_test_lat_ball_sample_from_ball_randomized(100, 10);
    res = res | quat_test_lat_ball_paralellogram_randomized(100, 100);
    return (res);
}
