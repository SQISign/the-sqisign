// SPDX-License-Identifier: Apache-2.0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rng.h>
#include <sig.h>
#include <api.h>

#ifdef ENABLE_CT_TESTING
#include <valgrind/memcheck.h>
#endif

#ifdef ENABLE_CT_TESTING
static void print_hex(const unsigned char *hex, int len) {
    unsigned char *copy  = calloc(len, 1);
    memcpy(copy, hex, len); // make a copy that we can tell valgrind is okay to leak
    VALGRIND_MAKE_MEM_DEFINED(copy, len);

    for (int i = 0; i < len;  ++i) {
        printf("%02x", copy[i]);
    }
    printf("\n");
    free(copy);
}
#else
static void print_hex(const unsigned char *hex, int len) {
    for (int i = 0; i < len;  ++i) {
        printf("%02x", hex[i]);
    }
    printf("\n");
}
#endif


static int test_sqisign() {
    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);
    unsigned char *sig = calloc(CRYPTO_BYTES + 32, 1);

    unsigned char seed[48] = { 0 };
    unsigned char msg[32] = { 0 };
    unsigned long long msglen = 32;

    randombytes_init(seed, NULL, 256);

    printf("Testing Keygen, Sign, Open: %s\n", CRYPTO_ALGNAME);

    int res = sqisign_keypair(pk, sk);
    if (res != 0) {
        res = -1;
        goto err;
    }

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(pk, CRYPTO_PUBLICKEYBYTES);
#endif

    unsigned long long smlen = CRYPTO_BYTES + 32;

    res = sqisign_sign(sig, &smlen, msg, 32, sk);
    if (res != 0) {
        res = -1;
        goto err;
    }

    printf("pk: ");
    print_hex(pk, CRYPTO_PUBLICKEYBYTES);
    printf("sk: ");
    print_hex(sk, CRYPTO_SECRETKEYBYTES);
    printf("sm: ");
    print_hex(sig, smlen);

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(sig, smlen);
#endif

    res = sqisign_open(msg, &msglen, sig, smlen, pk);
    if (res != 0) {
        res = -1;
        goto err;
    }

    sig[0] = ~sig[0];
    res = sqisign_open(msg, &msglen, sig, smlen, pk);
    if (res != 1) {
        res = -1;
        goto err;
    } else {
        res = 0;
    }

err:
    free(pk);
    free(sk);
    free(sig);
    return res;
}

int main(int argc, char *argv[]) {
    int rc = 0;

    rc = test_sqisign();

    if (rc != 0) {
        printf("test failed for %s\n", argv[1]);
    }
    return rc;
}
