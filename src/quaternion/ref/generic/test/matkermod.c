#include "quaternion_tests.h"

//void ibz_mat_howell(int rows, int cols, ibz_t howell[rows][rows+1], const ibz_t mat[rows][cols], ibz_t *mod)
int quat_test_ibz_mat_howell(){
    int res = 0;

    // Examples from Mulders & Storjohan
    
    ibz_t mod;
    ibz_t stormin[3][3], stormout[3][4], stormtrans[4][4], expected[3][4];
    ibz_init(&mod);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            ibz_init(&stormin[i][j]);
        for (int j = 0; j < 4; j++) {
            ibz_init(&stormout[i][j]);
            ibz_init(&expected[i][j]);
        }
    }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_init(&stormtrans[i][j]);

    ibz_set(&mod, 12);
    ibz_set(&expected[0][2], 4);
    ibz_set(&expected[1][2], 1);
    ibz_set(&expected[2][3], 1);
    
    ibz_set(&stormin[0][0], 4);
    ibz_set(&stormin[1][0], 1);
    ibz_set(&stormin[2][1], 5);
    int zeros = ibz_mat_howell(3, 3, stormout, stormtrans, stormin, &mod);
    res |= zeros != 2;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            res |= ibz_cmp(&stormout[i][j], &expected[i][j]);

    ibz_set(&stormin[0][0], 8);
    ibz_set(&stormin[1][0], 5);
    ibz_set(&stormin[2][0], 5);
    ibz_set(&stormin[1][1], 9);
    ibz_set(&stormin[2][1], 8);
    ibz_set(&stormin[2][2], 10);
    zeros = ibz_mat_howell(3, 3, stormout, stormtrans, stormin, &mod);
    res |= zeros != 2;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            res |= ibz_cmp(&stormout[i][j], &expected[i][j]);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            ibz_finalize(&stormin[i][j]);
        for (int j = 0; j < 4; j++) {
            ibz_finalize(&stormout[i][j]);
            ibz_finalize(&expected[i][j]);
        }
    }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_finalize(&stormtrans[i][j]);

    // An example computed in Pari
    ibz_t in[5][3], out[5][6], exp[5][6], trans[6][6];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 3; j++)
            ibz_init(&in[i][j]);
        for (int j = 0; j < 6; j++) {
            ibz_init(&out[i][j]);
            ibz_init(&exp[i][j]);
        }
    }
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            ibz_init(&trans[i][j]);

    ibz_set(&mod, 60);

    ibz_set(&in[0][0], 1);
    ibz_set(&in[0][1], 2);
    ibz_set(&in[0][2], -3);
    ibz_set(&in[1][0], 4);
    ibz_set(&in[1][1], 5);
    ibz_set(&in[1][2], 6);
    ibz_set(&in[2][0], 7);
    ibz_set(&in[2][1], 8);
    ibz_set(&in[2][2], 19);
    ibz_set(&in[3][0], 10);
    ibz_set(&in[3][1], 11);
    ibz_set(&in[3][2], 12);
    ibz_set(&in[4][0], 13);
    ibz_set(&in[4][1], 14);
    ibz_set(&in[4][2], 15);

    ibz_set(&exp[0][2], 12);
    ibz_set(&exp[0][3], 6);
    ibz_set(&exp[2][3], 10);
    ibz_set(&exp[1][4], 9);
    ibz_set(&exp[2][4], 6);
    ibz_set(&exp[3][4], 3);
    ibz_set(&exp[0][5], 1);
    ibz_set(&exp[1][5], 1);
    ibz_set(&exp[2][5], 1);
    ibz_set(&exp[3][5], 1);
    ibz_set(&exp[4][5], 1);
    
    zeros = ibz_mat_howell(5, 3, out, trans, in, &mod);
    res |= zeros != 2;
    
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 6; j++)
            res |= ibz_cmp(&out[i][j], &exp[i][j]);
    
    if (res != 0){
        printf("Quaternion unit test ibz_mat_howell failed\n");
    }
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 3; j++)
            ibz_finalize(&in[i][j]);
        for (int j = 0; j < 6; j++) {
            ibz_finalize(&out[i][j]);
            ibz_finalize(&exp[i][j]);
        }
    }
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            ibz_finalize(&trans[i][j]);
    ibz_finalize(&mod);
    return res;
}


//void ibz_mat_right_ker_mod(int rows, int cols, ibz_t ker[cols][rows], const ibz_t mat[rows][cols], ibz_t *mod)
int quat_test_ibz_mat_right_ker_mod() {
    int res = 0;

    // An example computed in Pari
    ibz_t mod, a, b, tmp;
    ibz_init(&mod);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&tmp);
    ibz_t mat[5][3], ker[3][5];
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 3; j++) {
            ibz_init(&mat[i][j]);
            ibz_init(&ker[j][i]);
        }
    }

    ibz_set(&mod, 60);

    ibz_set(&mat[0][0], 1);
    ibz_set(&mat[0][1], 2);
    ibz_set(&mat[0][2], -3);
    ibz_set(&mat[1][0], 4);
    ibz_set(&mat[1][1], 5);
    ibz_set(&mat[1][2], 6);
    ibz_set(&mat[2][0], 7);
    ibz_set(&mat[2][1], 8);
    ibz_set(&mat[2][2], 19);
    ibz_set(&mat[3][0], 10);
    ibz_set(&mat[3][1], 11);
    ibz_set(&mat[3][2], 12);
    ibz_set(&mat[4][0], 13);
    ibz_set(&mat[4][1], 14);
    ibz_set(&mat[4][2], 15);

    // self-testing thanks to assertions
    ibz_mat_right_ker_mod(5, 3, ker, mat, &mod);

    // Randomized test
    ibz_set(&mod, 6402373705728000l);
    for (int r = 0; r < 10; r++) {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 3; j++)
                ibz_rand_interval(&mat[i][j], &ibz_const_zero, &mod);
        for (int i = 2; i < 5; i++) {
            ibz_rand_interval(&a, &ibz_const_zero, &mod);
            ibz_rand_interval(&b, &ibz_const_zero, &mod);
            for (int j = 0; j < 3; j++) {
                ibz_mul(&tmp, &mat[0][j], &a);
                ibz_mul(&mat[i][j], &mat[1][1], &b);
                ibz_add(&mat[i][j], &mat[i][j], &tmp);
                ibz_mod(&mat[i][j], &mat[i][j], &mod);
            }
        }

        // self-testing thanks to assertions
        ibz_mat_right_ker_mod(5, 3, ker, mat, &mod);
    }
             
    
    if (res != 0){
        printf("Quaternion unit test ibz_mat_howell failed\n");
    }
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 3; j++) {
            ibz_finalize(&mat[i][j]);
            ibz_finalize(&ker[j][i]);
        }
    }
    ibz_finalize(&mod);
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&tmp);
    return res;
}

// run all previous tests
int quat_test_matkermod(){
    int res = 0;
    printf("\nRunning quaternion tests of matkermod\n");
    res = res | quat_test_ibz_mat_howell();
    res = res | quat_test_ibz_mat_right_ker_mod();
    return(res);
}
