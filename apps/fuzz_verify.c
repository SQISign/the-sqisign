// SPDX-License-Identifier: Apache-2.0

#include <mem.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <api.h>
#include <rng.h>

#include "encoded_sizes.h"

typedef struct {
  unsigned char pk[CRYPTO_PUBLICKEYBYTES];
  unsigned char sm[CRYPTO_BYTES + 64];
} signature_t;

static void crash() {
  int *p = 0;
  *p = 0;
}

static int load_signature(signature_t *sig, int iter) {
  char filename[sizeof("testcases/SQIsign_lvl1/signature000000.bin")];
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

  if (fread(sig->sm, CRYPTO_BYTES + 64, 1, f) != 1) {
    fprintf(stderr, "Can't read signature from file: %s\n", filename);
    fclose(f);
    return 1;
  }

  fclose(f);

  return 0;
}

static void verify_signature(signature_t corpus[], int testcases) {
  unsigned long long msglen = 64;
  unsigned long long smlen = CRYPTO_BYTES + msglen;

  unsigned char *pk = calloc(CRYPTO_PUBLICKEYBYTES, 1);
  unsigned char *sm = calloc(smlen, 1);

  unsigned char msg[msglen];

  if (fread(pk, CRYPTO_PUBLICKEYBYTES, 1, stdin) != 1) {
    fprintf(stderr, "Error reading public key from stdin\n");
    free(pk);
    free(sm);
    return;
  }
  if (fread(sm, smlen, 1, stdin) != 1) {
    fprintf(stderr, "Error reading signature from stdin\n");
    free(pk);
    free(sm);
    return;
  }

  int res = crypto_sign_open(msg, &msglen, sm, smlen, pk);
  if (res || msglen != sizeof(msg) || memcmp(msg, sm + SIGNATURE_BYTES, msglen)) {
    // Signature was not accepted -- check if it was in the corpus and, in that case, crash
    for (int i = 0; i < testcases; ++i)
      if (!memcmp(pk, corpus[i].pk, CRYPTO_PUBLICKEYBYTES) || !memcmp(sm, corpus[i].sm, smlen))
        crash();
  } else {
    // Signature was accepted -- check if it was not in the corpus and, in that case, crash
    int in_corpus = 0;
    for (int i = 0; i < testcases; ++i)
      if (!memcmp(pk, corpus[i].pk, CRYPTO_PUBLICKEYBYTES) || !memcmp(sm, corpus[i].sm, smlen)) {
        in_corpus = 1;
        break;
      }

    if (!in_corpus)
      crash();
  }

  free(pk);
  free(sm);
}

int
main(int argc, char *argv[]) {
  int testcases = 10;

  if (argc == 2) {
    sscanf(argv[1], "--testcases=%d", &testcases);
  }

  signature_t corpus[testcases];
  for (int i = 0; i < testcases; ++i)
    if (!load_signature(&corpus[i], i))
      return 1;

#ifdef __AFL_LOOP
  while (__AFL_LOOP(1000))
    verify_signature(corpus, testcases);
#else
  verify_signature(corpus, testcases);
#endif

  return 0;
}
