// SPDX-License-Identifier: Apache-2.0

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <rng.h>
#include <sig.h>
#include <api.h>
#include <bench_test_arguments.h>
#ifdef TARGET_BIG_ENDIAN
#include <tutil.h>
#endif

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

static int test_sqisign(unsigned long long in_msglen) {
    unsigned char *pk  = calloc(CRYPTO_PUBLICKEYBYTES, 1);
    unsigned char *sk  = calloc(CRYPTO_SECRETKEYBYTES, 1);
    unsigned char *sig = calloc(CRYPTO_BYTES + in_msglen, 1);

    unsigned char *msg = malloc(in_msglen);
    unsigned char *msg_open = malloc(in_msglen);

    unsigned long long msglen = in_msglen;

    randombytes(msg, in_msglen);
    randombytes(msg_open, in_msglen);

    printf("Testing Keygen, Sign, Open: %s\n", CRYPTO_ALGNAME);

    int res = sqisign_keypair(pk, sk);
    if (res != 0) {
        res = -1;
        goto err;
    }

#ifdef ENABLE_CT_TESTING
    VALGRIND_MAKE_MEM_DEFINED(pk, CRYPTO_PUBLICKEYBYTES);
#endif

    unsigned long long smlen = CRYPTO_BYTES + in_msglen;

    res = sqisign_sign(sig, &smlen, msg, in_msglen, sk);
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

    res = sqisign_open(msg_open, &msglen, sig, smlen, pk);
    if (res != 0 || msglen != in_msglen || memcmp(msg_open, msg, msglen)) {
        res = -1;
        goto err;
    }

    randombytes(msg_open, in_msglen);

    sig[0] = ~sig[0];
    res = sqisign_open(msg_open, &msglen, sig, smlen, pk);
    if (res != 1) {
        res = -1;
        goto err;
    } else {
        res = 0;
    }

    // Test with `sm` and `m` as the same buffer
    msglen = in_msglen;
    randombytes(msg, in_msglen);
    memcpy(sig, msg, msglen);

    res = sqisign_sign(sig, &smlen, sig, in_msglen, sk);
    if (res != 0) {
        res = -1;
        goto err;
    }

    randombytes(msg_open, in_msglen);

    res = sqisign_open(msg_open, &msglen, sig, smlen, pk);
    if (res != 0 || msglen != in_msglen || memcmp(msg_open, msg, msglen)) {
        res = -1;
        goto err;
    }


err:
    free(pk);
    free(sk);
    free(sig);
    free(msg);
    free(msg_open);
    return res;
}

int main(int argc, char *argv[]) {
    uint32_t seed[12] = { 0 };
    int help = 0;
    int seed_set = 0;
    int msglen_set = 0;
    int res = 0;
    unsigned long long msglen = 32;

    for (int i = 1; i < argc; i++) {
        unsigned int _msglen;

        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (!msglen_set && sscanf(argv[i], "--msglen=%u", &_msglen) == 1) {
            msglen = (unsigned long long) _msglen;
            msglen_set = 1;
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

    res = test_sqisign(msglen);

    if (res != 0) {
        printf("test failed for %s\n", argv[1]);
    }
    return res;
}
