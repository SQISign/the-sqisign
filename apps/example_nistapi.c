// SPDX-License-Identifier: Apache-2.0

/**
 * An example to demonstrate how to use SQIsign with the NIST API.
 */

#include <api.h>
#include <mem.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/**
 * Example for SQIsign variant:
 * - crypto_sign_keypair
 * - crypto_sign
 * - crypto_sign_open
 * 
 * @return int return code
 */
static int example_sqisign(void) {

    unsigned long long msglen = 32;
    unsigned long long smlen = CRYPTO_BYTES + msglen;

    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);

    unsigned char *sig = calloc(smlen, 1);

    unsigned char msg[32] = { 0xe };
    unsigned char msg2[32] = { 0 };

    printf("Example with %s\n", CRYPTO_ALGNAME);

    printf("crypto_sign_keypair -> ");
    int res = crypto_sign_keypair(pk, sk);
    if (res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("crypto_sign -> ");
    res = crypto_sign(sig, &smlen, msg, msglen, sk);
    if (res) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        printf("OK\n");
    }

    printf("crypto_sign_open (with correct signature) -> ");
    res = crypto_sign_open(msg2, &msglen, sig, smlen, pk);
    if (res || memcmp(msg, msg2, msglen)) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = 0;
        printf("OK\n");
    }

    printf("crypto_sign_open (with altered signature) -> ");
    sig[0] = ~sig[0];
    memset(msg2, 0, msglen);
    res = crypto_sign_open(msg2, &msglen, sig, smlen, pk);
    if (!res || !memcmp(msg, msg2, msglen)) {
        printf("FAIL\n");
        res = -1;
        goto err;
    } else {
        res = 0;
        printf("OK\n");
    }

err:
    free(pk);
    sqisign_secure_free(sk, CRYPTO_SECRETKEYBYTES);
    free(sig);
    return res;
}

int main(void) {
    return example_sqisign();
}
