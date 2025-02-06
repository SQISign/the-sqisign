// SPDX-License-Identifier: Apache-2.0

#include "randombytes_arm64crypto.h"

#include <arm_neon.h>
#include <string.h>

static AES256_CTR_DRBG_struct DRBG_ctx;

static inline uint32_t AES_sbox_x4(uint32_t in) {
  uint8x16_t sbox_val = vreinterpretq_u8_u32(vdupq_n_u32(in));
  sbox_val = vaeseq_u8(sbox_val, vdupq_n_u8(0));

  return vgetq_lane_u32(vreinterpretq_u32_u8(sbox_val), 0);
}

#define ROTR32(x, n) ((x << (32 - n)) | (x >> n))

typedef union {
  uint8_t u8[15][16];
  uint32_t u32[15][4];
} subkeys_t;

static void AES256_key_schedule(uint8_t subkeys[15][16], const uint8_t *key) {
  subkeys_t *sk = (subkeys_t *)subkeys;
  uint8_t rcon = 1;
  uint32_t s;
  int i, j;

  memcpy(&subkeys[0][0], key, 32 * sizeof(uint8_t));

  for (i = 2; i < 14; i += 2) {
    s = AES_sbox_x4(sk->u32[i - 1][3]);
    sk->u32[i][0] = ROTR32(s, 8) ^ rcon ^ sk->u32[i - 2][0];

    for (j = 1; j < 4; j++) {
      sk->u32[i][j] = sk->u32[i][j - 1] ^ sk->u32[i - 2][j];
    }

    s = AES_sbox_x4(sk->u32[i][3]);
    sk->u32[i + 1][0] = s ^ sk->u32[i - 1][0];

    for (j = 1; j < 4; j++) {
      sk->u32[i + 1][j] = sk->u32[i + 1][j - 1] ^ sk->u32[i - 1][j];
    }

    rcon = (rcon << 1) ^ ((rcon >> 7) * 0x11b);
  }

  s = AES_sbox_x4(sk->u32[13][3]);
  sk->u32[14][0] = ROTR32(s, 8) ^ rcon ^ sk->u32[12][0];

  for (j = 1; j < 4; j++) {
    sk->u32[14][j] = sk->u32[14][j - 1] ^ sk->u32[12][j];
  }
}

#define AES256_ECB_XWAYS(ways, vsubkeys, ctr, out)                             \
  do {                                                                         \
    uint8x16_t state[ways];                                                    \
                                                                               \
    for (int j = 0; j < ways; j++) {                                           \
      state[j] = vaeseq_u8(ctr[j], vsubkeys[0]);                               \
      state[j] = vaesmcq_u8(state[j]);                                         \
    }                                                                          \
                                                                               \
    for (int i = 1; i < 13; i++) {                                             \
      for (int j = 0; j < ways; j++) {                                         \
        state[j] = vaeseq_u8(state[j], vsubkeys[i]);                           \
        state[j] = vaesmcq_u8(state[j]);                                       \
      }                                                                        \
    }                                                                          \
                                                                               \
    for (int j = 0; j < ways; j++) {                                           \
      state[j] = vaeseq_u8(state[j], vsubkeys[13]);                            \
      state[j] = veorq_u8(state[j], vsubkeys[14]);                             \
      vst1q_u8(out + j * 16, state[j]);                                        \
    }                                                                          \
  } while (0);

//    subkeys - subkeys for AES-256
//    ctr - a 128-bit plaintext value
//    buffer - a 128-bit ciphertext value
static void AES256_ECB(uint8x16_t vsubkeys[15], uint8x16_t ctr,
                       unsigned char *buffer) {
  AES256_ECB_XWAYS(1, vsubkeys, (&ctr), buffer);
}

// vsubkeys - subkeys for AES-256
// ctr - an array of 3 x 128-bit plaintext value
// buffer - an array of 3 x 128-bit ciphertext value
static void AES256_ECB_x3(uint8x16_t vsubkeys[15], uint8x16_t ctr[3],
                          unsigned char *buffer) {
  AES256_ECB_XWAYS(3, vsubkeys, ctr, buffer);
}

static void bswap128(__uint128_t *x) {
  uint64_t *x64 = (uint64_t *)x;

  uint64_t t = x64[0];
  x64[0] = x64[1];
  x64[1] = t;

  x64[0] = __builtin_bswap64(x64[0]);
  x64[1] = __builtin_bswap64(x64[1]);
}

static void add_to_V(unsigned char V[], int incr) {
  __uint128_t *V128 = (__uint128_t *)V;
  bswap128(V128);
  (*V128) += incr;
  bswap128(V128);
}

static void AES256_CTR_DRBG_Update(unsigned char *provided_data,
                                   uint8x16_t vsubkeys[15], unsigned char *Key,
                                   unsigned char *V) {
  unsigned char temp[48];
  __uint128_t V128, t;
  uint64x2_t vV[3];

  memcpy(&V128, DRBG_ctx.V, sizeof(V128));

  bswap128(&V128);

  for (int j = 0; j < 3; j++) {
    V128++;
    t = V128;
    bswap128(&t);
    vV[j] = vld1q_u64((uint64_t *)&t);
  }

  AES256_ECB_x3(vsubkeys, (uint8x16_t *)vV, temp);

  if (provided_data != NULL)
    for (int i = 0; i < 48; i++)
      temp[i] ^= provided_data[i];
  memcpy(Key, temp, 32);
  memcpy(V, temp + 32, 16);

  add_to_V(DRBG_ctx.V, 1);
}

void randombytes_init_arm64crypto(unsigned char *entropy_input,
                                  unsigned char *personalization_string,
                                  int security_strength) {
  (void)security_strength;

  unsigned char seed_material[48];
  uint8_t subkeys[15][16];
  uint8x16_t vsubkeys[15];

  memcpy(seed_material, entropy_input, 48);
  if (personalization_string)
    for (int i = 0; i < 48; i++)
      seed_material[i] ^= personalization_string[i];
  memset(DRBG_ctx.Key, 0x00, 32);
  memset(DRBG_ctx.V, 0x00, 16);

  AES256_key_schedule(subkeys, DRBG_ctx.Key);
  for (int i = 0; i < 15; i++) {
    vsubkeys[i] = vld1q_u8(subkeys[i]);
  }

  AES256_CTR_DRBG_Update(seed_material, vsubkeys, DRBG_ctx.Key, DRBG_ctx.V);
  DRBG_ctx.reseed_counter = 1;
}

#define WAYS 4

int randombytes_arm64crypto(unsigned char *x, unsigned long long xlen) {
  uint8_t subkeys[15][16];
  unsigned char block[16];
  __uint128_t V[WAYS], Vle[WAYS];
  uint8x16x4_t vV;
  uint8x16_t vsubkeys[15];

  AES256_key_schedule(subkeys, DRBG_ctx.Key);

  for (int j = 0; j < 15; j++) {
    vsubkeys[j] = vld1q_u8(subkeys[j]);
  }

  memcpy(&Vle[0], DRBG_ctx.V, sizeof(Vle[0]));
  V[0] = Vle[0];
  vV.val[0] = vld1q_u8((uint8_t *)&V[0]);
  bswap128(&Vle[0]);
  for (int j = 1; j < WAYS; j++) {
    Vle[j] = Vle[j - 1] + 1;
    V[j] = Vle[j];
    bswap128(&V[j]);
    vV.val[j] = vld1q_u8((uint8_t *)&V[j]);
  }

  int entered_fast_path = (xlen >= WAYS * 16) ? 1 : 0;

  while (xlen >= WAYS * 16) {
    for (int j = 0; j < WAYS; j++) {
      Vle[j] += 4;
    }

    for (int j = 0; j < WAYS; j++) {
      vV.val[j] = vaeseq_u8(vV.val[j], vsubkeys[0]);
      vV.val[j] = vaesmcq_u8(vV.val[j]);
    }

    for (int i = 1; i < 13; i++) {
      for (int j = 0; j < WAYS; j++) {
        vV.val[j] = vaeseq_u8(vV.val[j], vsubkeys[i]);
        vV.val[j] = vaesmcq_u8(vV.val[j]);
      }
    }

    for (int j = 0; j < WAYS; j++) {
      vV.val[j] = vaeseq_u8(vV.val[j], vsubkeys[13]);
      vV.val[j] = veorq_u8(vV.val[j], vsubkeys[14]);
      vst1q_u8(x + j * 16, vV.val[j]);
    }

    for (int j = 0; j < WAYS; j++) {
      V[j] = Vle[j];
      bswap128(&V[j]);
    }

    vV = vld1q_u8_x4((uint8_t *)V);

    x += WAYS * 16;
    xlen -= WAYS * 16;
  }

  if (entered_fast_path && xlen == 0) {
    asm volatile("" : "+r,m"(Vle[3]) : : "memory");
    V[0] = Vle[3] - 4;
    bswap128(&V[0]);
  }

  while (xlen > 0) {
    if (xlen > 16) {
      AES256_ECB(vsubkeys, vld1q_u8((uint8_t *)&V[0]), x);
      x += 16;
      xlen -= 16;

      Vle[0]++;
      V[0] = Vle[0];
      bswap128(&V[0]);
    } else {
      AES256_ECB(vsubkeys, vld1q_u8((uint8_t *)&V[0]), block);
      memcpy(x, block, xlen);
      xlen = 0;
    }
  }

  memcpy(DRBG_ctx.V, &V[0], sizeof(V[0]));

  AES256_CTR_DRBG_Update(NULL, vsubkeys, DRBG_ctx.Key, DRBG_ctx.V);
  DRBG_ctx.reseed_counter++;

  return RNG_SUCCESS;
}

#ifdef RANDOMBYTES_ARM64CRYPTO
int randombytes(unsigned char *random_array, unsigned long long nbytes) {
  int ret = randombytes_arm64crypto(random_array, nbytes);
#ifdef ENABLE_CT_TESTING
  VALGRIND_MAKE_MEM_UNDEFINED(random_array, ret);
#endif
  return ret;
}

void randombytes_init(unsigned char *entropy_input,
                      unsigned char *personalization_string,
                      int security_strength) {
  randombytes_init_arm64crypto(entropy_input, personalization_string,
                               security_strength);
}
#endif
