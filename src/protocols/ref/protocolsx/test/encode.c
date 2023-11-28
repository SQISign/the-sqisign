#include <inttypes.h>

#include <ec.h>
#include <quaternion.h>
#include <protocols.h>
#include "test_protocols.h"

#include "../encode.c"


static int test_ibz_encdec() {
    int res = 1;
    ibz_t nm, nm_rec;
    ibz_init(&nm_rec);

    {
        ibz_set_from_str(&nm, "8892374823748723894723847892374897234897238478923748237", 16);
        int nm_bits = ibz_bitsize(&nm);
        int nm_bytes = (nm_bits + 7) / 8;
        unsigned char *enc = malloc(nm_bytes);

        ibz_encode(enc, &nm, nm_bytes);

        ibz_decode(&nm_rec, enc, nm_bytes);

        res &= (ibz_cmp(&nm, &nm_rec) == 0);

        free(enc);
        ibz_finalize(&nm);
    }

    {
        ibz_set_from_str(&nm, "-123", 16);
        int nm_bits = 23;
        int nm_bytes = (nm_bits + 7) / 8;
        unsigned char *enc = malloc(nm_bytes);

        ibz_encode(enc, &nm, nm_bytes);

        ibz_decode(&nm_rec, enc, nm_bytes);

        res &= (ibz_cmp(&nm, &nm_rec) == 0);

        free(enc);
        ibz_finalize(&nm);
    }

    ibz_finalize(&nm_rec);
    return res;
}


static void print_fp2(const fp2_t d, int len) {
    printf("re: ");
    for (int i = 0; i < len; ++i) {
        printf(HEX_FS, d.re[i]);
    } printf("\n");
    printf("im: ");
    for (int i = 0; i < len; ++i) {
        printf(HEX_FS, d.im[i]);
    } printf("\n");
}
static void print_pk(const public_key_t *pk) {
    print_fp2(pk->E.A, NWORDS_FIELD);
    print_fp2(pk->E.C, NWORDS_FIELD);
}

static int test_with_key(public_key_t *pk, secret_key_t *sk) {
    printf("Encoded bytes: %d\n", PUBLICKEY_BYTES);
    printf("pkenc\n");
    print_pk(pk);
    unsigned char *pkEnc = malloc(PUBLICKEY_BYTES);
    unsigned char *pkEnc2 = malloc(PUBLICKEY_BYTES);

    public_key_t pkdec = { 0 };

    public_key_encode(pkEnc, pk);

    public_key_decode(&pkdec, pkEnc);
    printf("pkdec\n");
    print_pk(&pkdec);

    public_key_encode(pkEnc2, &pkdec);

    int res = !memcmp(pkEnc, pkEnc2, PUBLICKEY_BYTES);

    free(pkEnc);
    free(pkEnc2);
    return res;
}

int test_encode() {
    int res = 1; 

    public_key_t pk;
    secret_key_t sk;
    secret_key_init(&sk); 

    res= res & protocols_keygen(&pk,&sk);
    res = res & test_with_key(&pk, &sk);

    secret_key_finalize(&sk);
    return res;
}
