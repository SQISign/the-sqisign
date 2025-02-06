// pqm4 KAT generator

// SPDX-License-Identifier: Apache-2.0 and Unknown

/*
NIST-developed software is provided by NIST as a public service. You may use,
copy, and distribute copies of the software in any medium, provided that you
keep intact this entire notice. You may improve, modify, and create derivative
works of the software or any portion of the software, and you may copy and
distribute such modifications or works. Modified works should carry a notice
stating that you changed the software and should note the date and nature of any
such change. Please explicitly acknowledge the National Institute of Standards
and Technology as the source of the software.

NIST-developed software is expressly provided "AS IS." NIST MAKES NO WARRANTY OF
ANY KIND, EXPRESS, IMPLIED, IN FACT, OR ARISING BY OPERATION OF LAW, INCLUDING,
WITHOUT LIMITATION, THE IMPLIED WARRANTY OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE, NON-INFRINGEMENT, AND DATA ACCURACY. NIST NEITHER REPRESENTS
NOR WARRANTS THAT THE OPERATION OF THE SOFTWARE WILL BE UNINTERRUPTED OR
ERROR-FREE, OR THAT ANY DEFECTS WILL BE CORRECTED. NIST DOES NOT WARRANT OR MAKE
ANY REPRESENTATIONS REGARDING THE USE OF THE SOFTWARE OR THE RESULTS THEREOF,
INCLUDING BUT NOT LIMITED TO THE CORRECTNESS, ACCURACY, RELIABILITY, OR
USEFULNESS OF THE SOFTWARE.

You are solely responsible for determining the appropriateness of using and
distributing the software and you assume all risks associated with its use,
including but not limited to the risks and costs of program errors, compliance
with applicable laws, damage to or loss of data, programs or equipment, and the
unavailability or interruption of operation. This software is not intended to be
used in any situation where a failure could cause risk of injury or damage to
property. The software developed by NIST employees is not subject to copyright
protection within the United States.
*/

#include "api.h"
#include "rng.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

#define MAX_MARKER_LEN 50

#define KAT_SUCCESS 0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR -3
#define KAT_CRYPTO_FAILURE -4

#define NUM_KATS 2
#define MAX_MSG_LEN 59

void output_header(FILE *fp) {
  const char header[] =
    "// SPDX-License-Identifier: Apache-2.0\n"
    "\n"
    "#ifndef api_h\n"
    "#define api_h\n"
    "\n"
    "#include <stddef.h>\n"
    "#include <sqisign_namespace.h>\n"
    "\n"
    "#define CRYPTO_SECRETKEYBYTES " STRINGIFY(CRYPTO_SECRETKEYBYTES) "\n"
    "#define CRYPTO_PUBLICKEYBYTES " STRINGIFY(CRYPTO_PUBLICKEYBYTES) "\n"
    "#define CRYPTO_BYTES " STRINGIFY(CRYPTO_BYTES) "\n"
    "\n"
    "#define CRYPTO_ALGNAME \"SQIsign_" STRINGIFY(SQISIGN_VARIANT) "\"\n"
    "\n"
    "SQISIGN_API\n"
    "int\n"
    "crypto_sign_keypair(unsigned char *pk, unsigned char *sk);\n"
    "\n"
    "SQISIGN_API\n"
    "int\n"
    "crypto_sign(unsigned char *sm, size_t *smlen,\n"
    "            const unsigned char *m, size_t mlen,\n"
    "            const unsigned char *sk);\n"
    "\n"
    "SQISIGN_API\n"
    "int\n"
    "crypto_sign_open(unsigned char *m, size_t *mlen,\n"
    "                 const unsigned char *sm, size_t smlen,\n"
    "                 const unsigned char *pk);\n"
    "\n"
    "#endif /* api_h */\n";

  fputs(header, fp);
}

void output_rng(FILE *fp) {
  const char rng[] =
    "// SPDX-License-Identifier: Apache-2.0\n"
    "\n"
    "#ifndef rng_h\n"
    "#define rng_h\n"
    "\n"
    "#include \"randombytes.h\"\n"
    "\n"
    "#endif /* rng_h */\n";
  
  fputs(rng, fp);
}

void output_preamble(FILE *fp) {
  const char preamble[] =
    "// SPDX-License-Identifier: Apache-2.0\n"
    "\n"
    "#include <api.h>\n"
    "#include <sig.h>\n"
    "#include <string.h>\n"
    "\n"
    "typedef struct {\n"
    "  size_t mlen;\n"
    "  char msg[" STRINGIFY(MAX_MSG_LEN) "];\n"
    "  size_t smlen;\n"
    "  char sm[" STRINGIFY(MAX_MSG_LEN) " + CRYPTO_BYTES];\n"
    "} SQISign_KAT_t;\n"
    "\n";

  fputs(preamble, fp);
}

void output_pk(FILE *fp, const unsigned char *pk) {
  fprintf(fp, "const char kat_" STRINGIFY(SQISIGN_VARIANT) "_pk[CRYPTO_PUBLICKEYBYTES] = {\n  ");
  for (int i = 0; i < CRYPTO_PUBLICKEYBYTES; i++) {
    fprintf(fp, "0x%02X, ", pk[i]);
  }
  fprintf(fp, "\n};\n\n");
}

void output_message_signature(FILE *fp, const unsigned char *m, unsigned long long mlen, const unsigned char *sm, unsigned long long smlen) {
  fprintf(fp, "  {\n"
              "    .mlen = %llu,\n"
              "    .msg = { ", mlen);
  for (unsigned long long i = 0; i < mlen; i++) {
    fprintf(fp, "0x%02X, ", m[i]);
  }
  fprintf(fp, "},\n"
              "    .smlen = %llu + CRYPTO_BYTES,\n"
              "    .sm = { ", mlen);
  for (unsigned long long i = 0; i < smlen; i++) {
    fprintf(fp, "0x%02X, ", sm[i]);
  }
  fprintf(fp, "},\n"
              "  },\n");
}

void output_implementation(FILE *fp) {
  const char api[] =
    "int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {\n"
    "  memcpy(pk, kat_" STRINGIFY(SQISIGN_VARIANT) "_pk, CRYPTO_PUBLICKEYBYTES);\n"
    "  // We don't need the secret key\n"
    "  memset(sk, 0, CRYPTO_SECRETKEYBYTES);\n"
    "}\n"
    "\n"
    "int crypto_sign(unsigned char *sm, size_t *smlen, const unsigned char *m,\n"
    "                size_t mlen, const unsigned char *sk) {\n"
    "  for (size_t i = 0; i < sizeof(kat_" STRINGIFY(SQISIGN_VARIANT) ") / sizeof(kat_" STRINGIFY(SQISIGN_VARIANT) "[0]); i++) {\n"
    "    if (mlen == kat_" STRINGIFY(SQISIGN_VARIANT) "[i].mlen) {\n"
    "      memcpy(sm, kat_" STRINGIFY(SQISIGN_VARIANT) "[i].sm, kat_" STRINGIFY(SQISIGN_VARIANT) "[i].smlen);\n"
    "      *smlen = kat_" STRINGIFY(SQISIGN_VARIANT) "[i].smlen;\n"
    "      return 0;\n"
    "    }\n"
    "  }\n"
    "\n"
    "  return 1;\n"
    "}\n"
    "\n"
    "int crypto_sign_open(unsigned char *m, size_t *mlen, const unsigned char *sm,\n"
    "                     size_t smlen, const unsigned char *pk) {\n"
    "  unsigned long long mlen_ull = *mlen;\n"
    "  int ret = sqisign_open(m, &mlen_ull, sm, smlen, pk);\n"
    "  if (mlen) {\n"
    "    *mlen = mlen_ull;\n"
    "  }\n"
    "  return ret;\n"
    "}\n";

  fputs(api, fp);
}

int main(void) {
  // pqm4 KATs use only all-zeros messages, one of length 32 and another of length 59;
  // arrays whose dimension are the length of the message are declared the size of the
  // largest of the two, but for the former, only 32 out of 59 bytes are actually used
  unsigned char seed[NUM_KATS][48];
  unsigned char entropy_input[48];
  const unsigned char m[NUM_KATS][MAX_MSG_LEN] = { { 0 }, { 0 } };
  unsigned char sm[NUM_KATS][MAX_MSG_LEN + CRYPTO_BYTES] = { { 0 }, { 0 } };
  unsigned char m1[MAX_MSG_LEN];
  const unsigned long long mlen[2] = { 32, 59 };
  unsigned long long smlen[2], mlen1;
  unsigned char pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
  int ret_val;

  for (int i = 0; i < 48; i++)
    entropy_input[i] = i;

  randombytes_init(entropy_input, NULL, 256);

  // Generate the keypair, shared between both KATs, as required for pqm4
  if ((ret_val = crypto_sign_keypair(pk, sk)) != 0) {
    printf("crypto_sign_keypair returned <%d>\n", ret_val);
    return KAT_CRYPTO_FAILURE;
  }

  // Choose two seeds (for the 32-byte and 59-byte KATs)
  for (int i = 0; i < NUM_KATS; i++)
    randombytes(seed[i], 48);

  // Fill m1 with random bytes. Note that the memcmp check below for a valid signature
  // compares m1 with m[i], but since the KATs use all-zero messages, the comparison
  // may suceed even if m was untouched from a previous iteration. This ensures that
  // memcmp will fail in that case.
  randombytes(m1, MAX_MSG_LEN);

  for (int i = 0; i < NUM_KATS; i++) {
    randombytes_init(seed[i], NULL, 256);

    if ((ret_val = crypto_sign(sm[i], &smlen[i], m[i], mlen[i], sk)) != 0) {
      printf("crypto_sign returned <%d>\n", ret_val);
      return KAT_CRYPTO_FAILURE;
    }

    if ((ret_val = crypto_sign_open(m1, &mlen1, sm[i], smlen[i], pk)) != 0) {
      printf("crypto_sign_open returned <%d>\n", ret_val);
      return KAT_CRYPTO_FAILURE;
    }

    if (mlen[i] != mlen1) {
      printf(
          "crypto_sign_open returned bad 'mlen': Got <%llu>, expected <%llu>\n",
          mlen1, mlen[i]);
      return KAT_CRYPTO_FAILURE;
    }

    if (memcmp(m, m1, mlen[i])) {
      printf("crypto_sign_open returned bad 'm' value\n");
      return KAT_CRYPTO_FAILURE;
    }

    // Fill m1 with random bytes for the next iteration
    randombytes(m1, MAX_MSG_LEN);
  }  

  // Output rng.h
  FILE *fp = fopen("src/pqm4/sqisign_" STRINGIFY(SQISIGN_VARIANT) "/ref/rng.h", "w");

  if (!fp) {
    printf("Couldn't open rng.h file for writing. Are you in the correct folder?\n");
    return KAT_FILE_OPEN_ERROR;
  }

  output_rng(fp);

  fclose(fp);

  // Output the header file
  fp = fopen("src/pqm4/sqisign_" STRINGIFY(SQISIGN_VARIANT) "/ref/api.h", "w");

  if (!fp) {
    printf("Couldn't open api.h file for writing. Are you in the correct folder?\n");
    return KAT_FILE_OPEN_ERROR;
  }

  output_header(fp);

  fclose(fp);

  // Output the implementation
  fp = fopen("src/pqm4/sqisign_" STRINGIFY(SQISIGN_VARIANT) "/ref/pqm4_api.c", "w");

  if (!fp) {
    printf("Couldn't open pqm4_api.c file for writing. Are you in the correct folder?\n");
    return KAT_FILE_OPEN_ERROR;
  }

  output_preamble(fp);

  output_pk(fp, pk);

  fprintf(fp, "const SQISign_KAT_t kat_" STRINGIFY(SQISIGN_VARIANT) "[%d] = {\n", NUM_KATS);

  for (int i = 0; i < NUM_KATS; i++) {
    output_message_signature(fp, m[i], mlen[i], sm[i], smlen[i]);
  }

  fprintf(fp, "};\n\n");

  output_implementation(fp);

  fclose(fp);

  return KAT_SUCCESS;
}
