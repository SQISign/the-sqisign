
#include <inttypes.h>
#include <locale.h>

#include <ec.h>
#include <quaternion.h>
#include <protocols.h>
#include <endomorphism_action.h>
#include "test_protocols.h"
#include <bench.h>

static __inline__ uint64_t rdtsc(void)
{
    return (uint64_t) cpucycles();
}

static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.im[i]);
    printf("\n");
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}

int test_full_verif() {
    int res = 1;
    public_key_t pk;
    secret_key_t sk;
    secret_key_init(&sk); 

setlocale(LC_NUMERIC, "");
uint64_t t0, t1;


t0 = rdtsc();
    res= res & protocols_keygen(&pk,&sk);
    if (!res) {
        printf("keygen failed \n");
    }
t1 = rdtsc();
// printf("\x1b[34mkeygen took %'" PRIu64 " cycles\x1b[0m\n", t1-t0);

    // printf(" \n\n keygen done \n\n");

    // init the sig
    signature_t sig;
    signature_init(&sig);
t0 = rdtsc();
    res = res & protocols_sign(&sig,&pk,&sk,(unsigned char*)"test",4);
    if (!res) {
        printf("sign failed \n");
    }

    
t1 = rdtsc();
// printf("\x1b[34msigning took %'" PRIu64 " cycles\x1b[0m\n", t1-t0);

    assert(res);
    // printf(" \n signing done \n");


    // printf(" \n starting verification \n ");

// t0 = rdtsc();
    res = res & protocols_verif(&sig,&pk,(unsigned char*)"test",4);
// t1 = rdtsc();
// printf("\x1b[34mverifying took %'" PRIu64 " cycles\x1b[0m\n", t1-t0);

    printf(" \x1b[35mfull\x1b[0m signature was: %s\n\n", res ? "\x1b[32mvalid\x1b[0m" : "\x1b[31minvalid\x1b[0m");


    signature_finalize(&sig);
    secret_key_finalize(&sk);

    return res;
}

bool curve_is_canonical(ec_curve_t const *E);   // test_protocols.c

int test_verif_from_chall()
{
    int res = 1;

    quat_left_ideal_t ideal_comm, ideal_chall;
    quat_left_ideal_init(&ideal_comm);
    quat_left_ideal_init(&ideal_chall);

    // commitment
    ec_curve_t E1;
    ec_basis_t E1basis = BASIS_CHALLENGE;
    protocols_commit(&ideal_comm, &E1, &E1basis);

    // challenge
    char const msg[] = "Hello, world!";
    signature_t sig = {{0}, {0xcafecafecafecafe}, {0xcc, {0xdeadbeef}, {0xdeadbeef}}};
    ec_curve_t E2;
    {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);

        hash_to_challenge(&vec, &E1, (unsigned char *) msg, strlen(msg));
//ibz_vec_2_print(&vec);

        protocols_challenge(&ideal_chall, &sig, &E1, &E1basis, &vec, &E2);

        ibz_vec_2_finalize(&vec);
    }

    ec_point_t dual;
    // now let's check it
    int valid = protocols_verif_from_chall(&sig, &E2, &dual, (unsigned char *) msg, strlen(msg));
    printf("\n\x1b[35mchallenge\x1b[0m signature was: %s\n\n", valid ? "\x1b[32mvalid\x1b[0m" : "\x1b[31minvalid\x1b[0m");
    res &= valid;

    res &= curve_is_canonical(&E1);
    res &= curve_is_canonical(&E2);

    quat_left_ideal_finalize(&ideal_comm);
    quat_left_ideal_finalize(&ideal_chall);

    return res;
}

int test_sign_verif()
{
    int res = 1;

    for (int i =0;i<2;i++) {
        printf("#%d ",i);
        res &= test_full_verif();
    }

    if (!res) {
        printf("\n Keygen + Sign + Verif failed ! \n");
    }
    // res &= test_verif_from_chall(); 
    

    return res;
}
