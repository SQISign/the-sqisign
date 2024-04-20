
#include <ec.h>
#include <quaternion.h>
#include <protocols.h>
#include <inttypes.h>
#include "test_protocols.h"


static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.im[i]);
    printf("\n");
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}

bool curve_is_canonical(ec_curve_t const *E);   // test_protocols.c

int test_commitment()
{
    int res = 1;

    quat_left_ideal_t ideal;
    quat_left_ideal_init(&ideal);

    ec_curve_t E1;
    ec_basis_t E1basis = BASIS_CHALLENGE;

    protocols_commit(&ideal, &E1, &E1basis);

    res &= !ibz_cmp(&ideal.norm, &DEGREE_COMMITMENT);

    res &= curve_is_canonical(&E1);

// quat_left_ideal_print(&ideal);
    res &= !fp2_is_zero(&E1.C);
// curve_print("E1", E1);

    quat_left_ideal_finalize(&ideal);

    return res;
}

