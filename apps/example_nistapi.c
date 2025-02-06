// SPDX-License-Identifier: Apache-2.0

/**
 * An example to demonstrate how to use SQIsign with the NIST API.
 */

#include <inttypes.h>
#include <mem.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include <api.h>
#include <rng.h>
#include <bench_test_arguments.h>
#if defined(TARGET_BIG_ENDIAN)
#include <tutil.h>
#endif

static uint32_t rand_u32()
{
    unsigned char buf[4];
    if (randombytes(buf, sizeof(buf)))
        abort();
    return ((uint32_t) buf[3] << 24)
         | ((uint32_t) buf[2] << 16)
         | ((uint32_t) buf[1] <<  8)
         | ((uint32_t) buf[0] <<  0);
}

/**
 * Example for SQIsign variant:
 * - crypto_sign_keypair
 * - crypto_sign
 * - crypto_sign_open
 *
 * @return int return code
 */
static int
example_sqisign(void)
{

    unsigned long long msglen = rand_u32() % 100;
    unsigned long long smlen = CRYPTO_BYTES + msglen;

    unsigned char *sk = calloc(CRYPTO_SECRETKEYBYTES, 1);
    unsigned char *pk = calloc(CRYPTO_PUBLICKEYBYTES, 1);

    unsigned char *sm = calloc(smlen, 1);

    unsigned char msg[msglen], msg2[msglen];

    printf("Example with %s\n", CRYPTO_ALGNAME);

    printf("crypto_sign_keypair -> ");
    int res = crypto_sign_keypair(pk, sk);
    if (res) {
        printf("FAIL\n");
        goto err;
    } else {
        printf("OK\n");
    }

    // choose a random message
    for (size_t i = 0; i < msglen; ++i)
        msg[i] = rand_u32();

    printf("crypto_sign -> ");
    res = crypto_sign(sm, &smlen, msg, msglen, sk);
    if (res) {
        printf("FAIL\n");
        goto err;
    } else {
        printf("OK\n");
    }

    printf("crypto_sign_open (with correct signature) -> ");
    res = crypto_sign_open(msg2, &msglen, sm, smlen, pk);
    if (res || msglen != sizeof(msg) || memcmp(msg, msg2, msglen)) {
        printf("FAIL\n"); // signature was not accepted!?
        goto err;
    } else {
        printf("OK\n");
    }


    // fill with random bytes
    for (size_t i = 0; i < msglen; ++i)
        msg2[i] = rand_u32();

    // let's try a single bit flip
    size_t pos = rand_u32() % smlen;
    sm[pos / 8] ^= 1 << pos % 8;

    res = crypto_sign_open(msg2, &msglen, sm, smlen, pk);

    printf("crypto_sign_open (with altered signature) -> ");
    if (!res) {
        printf("FAIL\n"); // signature was accepted anyway!?
        res = -1;
        goto err;
    }
    else {
        printf("OK\n");
        res = 0;

        if (msglen)
            printf("WARNING: verification failed but the message length was returned nonzero; misuse-prone API\n");

        unsigned char any = 0;
        for (size_t i = 0; i < msglen; ++i)
            any |= msg2[i];
        if (any)
            printf("WARNING: verification failed but the message buffer was not zeroed out; misuse-prone API\n");
    }

err:
    sqisign_secure_free(sk, CRYPTO_SECRETKEYBYTES);
    free(pk);
    free(sm);

    return res;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }
    }

    if (help) {
        printf("Usage: %s [--seed=<seed>]\n", argv[0]);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    return example_sqisign();
}
