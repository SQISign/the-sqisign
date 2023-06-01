#include "id2iso_tests.h"

static void random_scalar(ibz_t *k)
{
    ibz_pow(k, &ibz_const_two, 123);
    ibz_rand_interval(k, &ibz_const_zero, k);
}

int _id2iso_test_ker2id() {
    int res = 1;

    ibz_vec_2_t vec2, vec3;
    ibz_vec_2_init(&vec2);
    ibz_vec_2_init(&vec3);

    {
        ibz_t gcd;
        ibz_init(&gcd);
        do {
            random_scalar(&vec2[0]);
            random_scalar(&vec2[1]);
            ibz_gcd(&gcd, &vec2[0], &vec2[1]);
        } while (ibz_divides(&gcd, &ibz_const_two));
        ibz_t three;
        ibz_init(&three);
        ibz_set(&three, 3);
        do {
            random_scalar(&vec3[0]);
            random_scalar(&vec3[1]);
            ibz_gcd(&gcd, &vec3[0], &vec3[1]);
        } while (ibz_divides(&gcd, &three));
        ibz_finalize(&three);
        ibz_finalize(&gcd);
    }
//    ibz_printf("vec2 = (%Zd,%Zd)\n", &vec2[0], &vec2[1]);
//    ibz_printf("vec3 = (%Zd,%Zd)\n", &vec3[0], &vec3[1]);

    quat_left_ideal_t I;
    quat_left_ideal_init(&I);
    id2iso_kernel_dlogs_to_ideal(&I, &vec2, &vec3);
//    quat_left_ideal_print(&I);
    quat_left_ideal_finalize(&I);

    //TODO FIXME this really only tests that the function doesn't crash

    return res;
}


int id2iso_test_ker2id() {
    int res = 1;
    printf("\n \nRunning id2iso tests for kernel_dlogs_to_ideal() \n \n");

    for (int i = 0; i < 10; i++) {
        res &= _id2iso_test_ker2id();
    }

    if (!res) {
        printf("ID2ISO unit test kernel_dlogs_to_ideal() failed\n");
    }

    return res;
}
