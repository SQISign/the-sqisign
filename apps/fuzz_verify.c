// SPDX-License-Identifier: Apache-2.0

#include <mem.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <api.h>
#include <rng.h>

#include "encoded_sizes.h"

// Largest signed message the harness will hand to the verifier: a signature plus a 64-byte message.
#define SM_MAX (CRYPTO_BYTES + 64)

typedef struct
{
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sm[SM_MAX];
    size_t smlen;
} testcase_t;

// Crash on purpose, so that the fuzzer records the input as a failing test case. abort() raises SIGABRT, which every
// fuzzer records, and being _Noreturn it cannot be optimised away.
static void
crash()
{
    abort();
}

static int
load_signature(testcase_t *sig, int iter)
{
    // Add some extra bytes because in theory iter might be larger than 1 million, and gcc complains about a possible
    // buffer overflow. In practice the number of testcases should be much smaller than even 1 million, and anyway there
    // are no consequences to such an overflow, as this is just a fuzzing app.
    char filename[sizeof("testcases/" CRYPTO_ALGNAME "/signature000000.bin") + 16];
    snprintf(filename, sizeof(filename), "testcases/%s/signature%06d.bin", CRYPTO_ALGNAME, iter);
    FILE *f = fopen(filename, "rb");

    if (!f) {
        fprintf(stderr, "Can't open file: %s\n", filename);
        return 1;
    }

    if (fread(sig->pk, CRYPTO_PUBLICKEYBYTES, 1, f) != 1) {
        fprintf(stderr, "Can't read public key from file: %s\n", filename);
        fclose(f);
        return 1;
    }

    // The signed message is whatever remains in the file, so a corpus entry can record any length up to SM_MAX.
    sig->smlen = fread(sig->sm, 1, SM_MAX, f);
    if (sig->smlen == 0) {
        fprintf(stderr, "Can't read signature from file: %s\n", filename);
        fclose(f);
        return 1;
    }

    fclose(f);

    return 0;
}

static void
verify_signature(const testcase_t corpus[], int testcases)
{
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char staging[SM_MAX];
    unsigned long long msglen = 0;

    // AFL rewinds the fd between persistent-mode iterations, but not our stdio buffer.
    (void)fseek(stdin, 0, SEEK_SET);
    clearerr(stdin);

    if (fread(pk, CRYPTO_PUBLICKEYBYTES, 1, stdin) != 1) {
        fprintf(stderr, "Error reading public key from stdin\n");
        return;
    }

    // smlen is however many bytes the input supplies, so inputs of any length reach the verifier.
    size_t smlen = fread(staging, 1, sizeof(staging), stdin);

    // Both buffers are sized to exactly the lengths passed to crypto_sign_open(), so a sanitizer catches any read or
    // write past them.
    size_t msgmax = smlen > (size_t)SIGNATURE_BYTES ? smlen - (size_t)SIGNATURE_BYTES : 0;
    unsigned char *sm = malloc(smlen ? smlen : 1);
    unsigned char *msg = malloc(msgmax ? msgmax : 1);

    // A crash on a failed allocation would be recorded by the fuzzer as a verification failure, so return instead.
    if (!sm || !msg) {
        fprintf(stderr, "Memory allocation failed\n");
        free(sm);
        free(msg);
        return;
    }

    memcpy(sm, staging, smlen);

    int accepted = crypto_sign_open(msg, &msglen, sm, smlen, pk) == 0;

    // A signed message shorter than the fixed signature size carries no message and cannot be accepted, whatever its
    // contents.
    if (accepted && smlen < (size_t)SIGNATURE_BYTES)
        crash();

    // An accepted signature must recover exactly the bytes appended to the signature.
    if (accepted && (msglen != (unsigned long long)msgmax || memcmp(msg, sm + SIGNATURE_BYTES, msgmax) != 0))
        crash();

    // The corpus holds known-good inputs, so acceptance must coincide with corpus membership. The length, the key and
    // the signed message all have to match.
    int in_corpus = 0;
    for (int i = 0; i < testcases; ++i)
        if (smlen == corpus[i].smlen && !memcmp(pk, corpus[i].pk, CRYPTO_PUBLICKEYBYTES) &&
            !memcmp(sm, corpus[i].sm, smlen)) {
            in_corpus = 1;
            break;
        }

    if (accepted != in_corpus)
        crash();

    free(sm);
    free(msg);
}

int
main(int argc, char *argv[])
{
    int testcases = 10;

    if (argc == 2) {
        sscanf(argv[1], "--testcases=%d", &testcases);
    }

    testcase_t corpus[testcases];
    for (int i = 0; i < testcases; ++i)
        if (load_signature(&corpus[i], i))
            return 1;

#ifdef __AFL_LOOP
    while (__AFL_LOOP(1000))
        verify_signature(corpus, testcases);
#else
    verify_signature(corpus, testcases);
#endif

    return 0;
}
