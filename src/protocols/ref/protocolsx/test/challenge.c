
#include <inttypes.h>

#include <ec.h>
#include <quaternion.h>
#include <protocols.h>
#include <endomorphism_action.h>
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


int test_challenge()
{
    int res = 1;

    quat_left_ideal_t ideal_comm, ideal_chall;
    quat_left_ideal_init(&ideal_comm);
    quat_left_ideal_init(&ideal_chall);

    // commitment
    ec_curve_t E1;
    ec_basis_t E1basis = BASIS_CHALLENGE;
    protocols_commit(&ideal_comm, &E1, &E1basis);
// curve_print("E1", E1);

#ifdef RADIX_32
    signature_t sig = {{0}, {0xcafecafe, 0xcafecafe}, {0xcc, {0xdeadbeef}, {0xdeadbeef}}};
#else
    signature_t sig = {{0}, {0xcafecafecafecafe}, {0xcc, {0xdeadbeef}, {0xdeadbeef}}};
#endif
    {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);

        char msg[] = "Hello, world!";
        hash_to_challenge(&vec, &E1, (unsigned char *) msg, strlen(msg));
// ibz_vec_2_print(&vec);

        protocols_challenge(&ideal_chall, &sig, &E1, &E1basis, &vec, NULL);
// quat_left_ideal_print(&ideal_chall);
{
// printf("r = 0x");
if (0) {
for (ssize_t i = sizeof(sig.r)/sizeof(*sig.r)-1; i >= 0; --i)
    printf(HEX_FS, sig.r[i]);
printf("\n");
printf("s.bit2 = %u\n", !!(sig.s.select23 & 1));
printf("s.bit3 = %u\n", !!(sig.s.select23 & 2));
printf("s.scalar2 = 0x");
for (ssize_t i = sizeof(sig.s.scalar2)/sizeof(*sig.s.scalar2)-1; i >= 0; --i)
    printf(HEX_FS, sig.s.scalar2[i]);
printf("\n");
printf("s.scalar3 = 0x");
for (ssize_t i = sizeof(sig.s.scalar3)/sizeof(*sig.s.scalar3)-1; i >= 0; --i)
    printf(HEX_FS, sig.s.scalar3[i]);
printf("\n");
}
}
        ibz_vec_2_finalize(&vec);
    }

    quat_left_ideal_finalize(&ideal_comm);
    quat_left_ideal_finalize(&ideal_chall);

    return res;
}

