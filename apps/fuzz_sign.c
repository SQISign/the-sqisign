// SPDX-License-Identifier: Apache-2.0

#include <mem.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <api.h>
#include <rng.h>

/**
 * Example for SQIsign variant:
 * - crypto_sign_keypair
 * - crypto_sign
 * - crypto_sign_open
 *
 * @return int return code
 */
static int example_sqisign(int iter) {
  int ret = 0;

  unsigned long long msglen = 64;
  unsigned long long smlen = CRYPTO_BYTES + msglen;

  unsigned char *sk = calloc(CRYPTO_SECRETKEYBYTES, 1);
  unsigned char *pk = calloc(CRYPTO_PUBLICKEYBYTES, 1);

  unsigned char *sm = calloc(smlen, 1);

  unsigned char msg[msglen];

  FILE *f = NULL;

  int res = crypto_sign_keypair(pk, sk);
  if (res) {
    fprintf(stderr, "crypto_sign_keypair -> FAIL\n");
    ret = 1;
    goto end;
  }

  // choose a random message
  randombytes(msg, msglen);

  res = crypto_sign(sm, &smlen, msg, msglen, sk);
  if (res) {
    fprintf(stderr, "crypto_sign -> FAIL\n");
    ret = 1;
    goto end;
  }

  // This string is larger than necessary, but gcc is not smart enough
  // to detect that iter < 1000000 in the snprintf call below
  char filename[sizeof("testcases/SQIsign_lvl1/signature4294967296.bin") + 1];
  if (iter > 999999) {
    fprintf(stderr, "Too many iterations: %d\n", iter);
    ret = 1;
    goto end;
  }
  snprintf(filename, sizeof(filename), "testcases/%s/signature%06d.bin",
           CRYPTO_ALGNAME, iter);
  f = fopen(filename, "wb");
  if (!f) {
    fprintf(stderr,
            "Can't open file: %s (have you created the testcases/%s folder?)\n",
            filename, CRYPTO_ALGNAME);
    ret = 1;
    goto end;
  }

  if (fwrite(pk, CRYPTO_PUBLICKEYBYTES, 1, f) != 1) {
    fprintf(stderr, "Error writing public key to file\n");
    ret = 1;
    goto end;
  }
  if (fwrite(sm, smlen, 1, f) != 1) {
    fprintf(stderr, "Error writing signature to file\n");
    ret = 1;
    goto end;
  }

end:
  if (f)
    fclose(f);

  free(sk);
  free(pk);
  free(sm);

  return ret;
}

// Brief fuzzing tutorial (assumes level 1, but works for other levels)
// Assumes an Intel Linux system
//
// 0. Some OS configurations required for AFL to work:
//    echo core | sudo tee /proc/sys/kernel/core_pattern
//    echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
//
// 1. Install AFL++ (https://aflplus.plus) -- on Ubuntu, install packages
//    afl++, clang, lld, tmux
//
// 2. Configure CMake using AFL's compilers:
//    AFL_HARDEN=1 cmake -DCMAKE_C_COMPILER=afl-clang-lto <other params>
//
// 3. Build
//
// 4. cd to the build/apps folder
//
// 5. Create required folders:
//    mkdir -p testcases/SQIsign_lvl{1,3,5}
//
// 6. Run ./fuzz_sign_lvl1 to create some initial testcases
//
// 7. Run:
//    tmux new-session -s afl1 afl-fuzz -i testcases/SQIsign_lvl1/ -o syncdir/ -D -M fuzz1 -- ./fuzz_verify_lvl1
//
// 8. Optionally run using other cores in the machine (e.g. for 24 cores),
//    for i in $(seq 2 24)
//    do
//      tmux new-session -s afl$i -d afl-fuzz -d -i testcases/SQIsign_lvl1/ -o syncdir/ -S fuzz$i -- ./fuzz_verify_lvl1
//    done
//
// 9. Attach to a specific instance by running:
//    tmux attach -t afl$i
//
// 10. To get summary statistics for all runs, run:
//     afl-whatsup syncdir/
//
// 11. "Interesting" signatures, in a binary format understood by
//     fuzz_verify_lvl1, will be found in syncdir/fuzz$i/crashes; to
//     reproduce the crash, pipe one of these files to fuzz_verify_lvl1

int
main(int argc, char *argv[]) {
  int testcases = 10;

  unsigned char seed[48];
  randombytes_select(seed, sizeof(seed));
  randombytes_init(seed, NULL, 256);

  if (argc == 2) {
    sscanf(argv[1], "--testcases=%d", &testcases);
  }

  for (int i = 0; i < testcases; ++i) {
    example_sqisign(i);
  }

  return 0;
}
