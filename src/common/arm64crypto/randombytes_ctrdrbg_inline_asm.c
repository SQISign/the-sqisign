// SPDX-License-Identifier: Apache-2.0

#include <arm_neon.h>
#include <string.h>

#include "randombytes_arm64crypto.h"

typedef union {
  uint8_t u8[16];
  uint64_t u64[2];
  __uint128_t u128;
} u128_t;

static AES256_CTR_DRBG_struct DRBG_ctx;

static inline uint32_t AES_sbox_x4(uint32_t in) {
  uint8x16_t sbox_val = vreinterpretq_u8_u32(vdupq_n_u32(in));
  sbox_val = vaeseq_u8(sbox_val, vdupq_n_u8(0));

  return vgetq_lane_u32(vreinterpretq_u32_u8(sbox_val), 0);
}

#define ROTR32(x, n) ((x << (32 - n)) | (x >> n))

typedef union {
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

static void bswap128(u128_t *x) {
  uint64_t t = x->u64[0];
  x->u64[0] = x->u64[1];
  x->u64[1] = t;

  x->u64[0] = __builtin_bswap64(x->u64[0]);
  x->u64[1] = __builtin_bswap64(x->u64[1]);
}

static void incr_V(u128_t *V) {
  bswap128(V);
  V->u128++;
  bswap128(V);
}

static void AES256_CTR_DRBG_Update(const unsigned char *provided_data,
                                   uint8x16_t vsubkeys[15], unsigned char *Key,
                                   unsigned char *V) {
  (void)V;

  unsigned char temp[48];
  u128_t V128, t;
  uint64x2_t vV[3];

  memcpy(&V128, DRBG_ctx.V, sizeof(V128));

  bswap128(&V128);

  for (int j = 0; j < 3; j++) {
    V128.u128++;
    t = V128;
    bswap128(&t);
    vV[j] = vld1q_u64((uint64_t *)&t);
  }

  AES256_ECB_x3(vsubkeys, (uint8x16_t *)vV, temp);

  if (provided_data != NULL)
    for (int i = 0; i < 48; i++)
      temp[i] ^= provided_data[i];
  memcpy(Key, temp, 32);
  memcpy(V128.u8, temp + 32, 16);

  incr_V(&V128);

  memcpy(DRBG_ctx.V, V128.u8, 16);
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
  u128_t V[WAYS], Vle[WAYS];
  uint8x16x4_t vV;
  uint8x16_t vsubkeys[15];

  AES256_key_schedule(subkeys, DRBG_ctx.Key);

  for (int j = 0; j < 15; j++) {
    vsubkeys[j] = vld1q_u8(subkeys[j]);
  }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
  asm("ldp         %[V0l],     %[V0h],  %[DRBG_ctx_V]     \n\t"
      "stp         %[V0l],     %[V0h],    [%[V]     ]     \n\t"
      "rev       %[Vle0h],     %[V0l]                     \n\t"
      "rev       %[Vle0l],     %[V0h]                     \n\t"
      "adds      %[Vle1l],   %[Vle0l],             #1     \n\t"
      "adc       %[Vle1h],   %[Vle0h],            xzr     \n\t"
      "rev         %[V1h],   %[Vle1l]                     \n\t"
      "rev         %[V1l],   %[Vle1h]                     \n\t"
      "stp         %[V1l],     %[V1h],    [%[V], #16]     \n\t"
      "adds      %[Vle2l],   %[Vle0l],             #2     \n\t"
      "adc       %[Vle2h],   %[Vle0h],            xzr     \n\t"
      "rev         %[V2h],   %[Vle2l]                     \n\t"
      "rev         %[V2l],   %[Vle2h]                     \n\t"
      "stp         %[V2l],     %[V2h],    [%[V], #32]     \n\t"
      "adds      %[Vle3l],   %[Vle0l],             #3     \n\t"
      "adc       %[Vle3h],   %[Vle0h],            xzr     \n\t"
      "rev         %[V3h],   %[Vle3l]                     \n\t"
      "rev         %[V3l],   %[Vle3h]                     \n\t"
      "stp         %[V3l],     %[V3h],    [%[V], #48]     \n\t"
      "ld1       { %[vV0].16b, %[vV1].16b, %[vV2].16b, %[vV3].16b }, [%[V]]\n\t"
      "cmp        %[xlen],          #64                   \n\t"
      "b.lo            2f                                 \n\t"
      ".p2align         6                                 \n\t"
      "1:                                                 \n\t"
      "aese    %[vV0].16b,  %[vsk0].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk0].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk0].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk0].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk1].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk1].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk1].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk1].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "adds      %[Vle0l],     %[Vle0l],           #4     \n\t"
      "adc       %[Vle0h],     %[Vle0h],          xzr     \n\t"
      "adds      %[Vle1l],     %[Vle1l],           #4     \n\t"
      "adc       %[Vle1h],     %[Vle1h],          xzr     \n\t"
      "adds      %[Vle2l],     %[Vle2l],           #4     \n\t"
      "adc       %[Vle2h],     %[Vle2h],          xzr     \n\t"
      "adds      %[Vle3l],     %[Vle3l],           #4     \n\t"
      "adc       %[Vle3h],     %[Vle3h],          xzr     \n\t"
      "aese    %[vV0].16b,  %[vsk2].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk2].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk2].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk2].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk3].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk3].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk3].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk3].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "rev         %[V0h],     %[Vle0l]                   \n\t"
      "rev         %[V0l],     %[Vle0h]                   \n\t"
      "rev         %[V1h],     %[Vle1l]                   \n\t"
      "rev         %[V1l],     %[Vle1h]                   \n\t"
      "rev         %[V2h],     %[Vle2l]                   \n\t"
      "rev         %[V2l],     %[Vle2h]                   \n\t"
      "rev         %[V3h],     %[Vle3l]                   \n\t"
      "rev         %[V3l],     %[Vle3h]                   \n\t"
      "aese    %[vV0].16b,  %[vsk4].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk4].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk4].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk4].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk5].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk5].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk5].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk5].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk6].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk6].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk6].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk6].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk7].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk7].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk7].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk7].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk8].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk8].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk8].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk8].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b,  %[vsk9].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b,  %[vsk9].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b,  %[vsk9].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b,  %[vsk9].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "stp         %[V0l],       %[V0h],  [%[V]]          \n\t"
      "stp         %[V1l],       %[V1h],  [%[V], #16]     \n\t"
      "stp         %[V2l],       %[V2h],  [%[V], #32]     \n\t"
      "stp         %[V3l],       %[V3h],  [%[V], #48]     \n\t"
      "aese    %[vV0].16b, %[vsk10].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b, %[vsk10].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b, %[vsk10].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b, %[vsk10].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b, %[vsk11].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b, %[vsk11].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b, %[vsk11].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b, %[vsk11].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b, %[vsk12].16b                   \n\t"
      "aesmc   %[vV0].16b,   %[vV0].16b                   \n\t"
      "aese    %[vV1].16b, %[vsk12].16b                   \n\t"
      "aesmc   %[vV1].16b,   %[vV1].16b                   \n\t"
      "aese    %[vV2].16b, %[vsk12].16b                   \n\t"
      "aesmc   %[vV2].16b,   %[vV2].16b                   \n\t"
      "aese    %[vV3].16b, %[vsk12].16b                   \n\t"
      "aesmc   %[vV3].16b,   %[vV3].16b                   \n\t"
      "aese    %[vV0].16b, %[vsk13].16b                   \n\t"
      "eor     %[vV0].16b,   %[vV0].16b, %[vsk14].16b     \n\t"
      "aese    %[vV1].16b, %[vsk13].16b                   \n\t"
      "eor     %[vV1].16b,   %[vV1].16b, %[vsk14].16b     \n\t"
      "stp        %q[vV0],      %q[vV1],       [%[x]], #32\n\t"
      "aese    %[vV2].16b, %[vsk13].16b                   \n\t"
      "eor     %[vV2].16b,   %[vV2].16b, %[vsk14].16b     \n\t"
      "aese    %[vV3].16b, %[vsk13].16b                   \n\t"
      "eor     %[vV3].16b,   %[vV3].16b, %[vsk14].16b     \n\t"
      "stp        %q[vV2],      %q[vV3],       [%[x]], #32\n\t"
      "sub        %[xlen],      %[xlen],          #64     \n\t"
      "ld1       { %[vV0].16b, %[vV1].16b, %[vV2].16b, %[vV3].16b }, [%[V]]\n\t"
      "cmp        %[xlen],          #64                   \n\t"
      "b.hs            1b                                 \n\t"
      "cbnz       %[xlen],           2f                   \n\t"
      "subs        %[V0h],     %[Vle3l],           #4     \n\t"
      "sbc         %[V0l],     %[Vle3h],          xzr     \n\t"
      "rev         %[V0h],       %[V0h]                   \n\t"
      "rev         %[V0l],       %[V0l]                   \n\t"
      "stp         %[V0l],       %[V0h],       [%[V]]     \n\t"
      "2:                                                 \n\t"
      : [vV0] "=&w"(vV.val[0]), [vV1] "=&w"(vV.val[1]), [vV2] "=&w"(vV.val[2]),
        [vV3] "=&w"(vV.val[3]), [Vle0l] "=&r"(Vle[0].u64[0]),
        [Vle0h] "=&r"(Vle[0].u64[1]), [Vle1l] "=&r"(Vle[1].u64[0]),
        [Vle1h] "=&r"(Vle[1].u64[1]), [Vle2l] "=&r"(Vle[2].u64[0]),
        [Vle2h] "=&r"(Vle[2].u64[1]), [Vle3l] "=&r"(Vle[3].u64[0]),
        [Vle3h] "=&r"(Vle[3].u64[1]), [x] "+r"(x), [xlen] "+r"(xlen),
        [V0l] "=&r"(V[0].u64[0]), [V0h] "=&r"(V[0].u64[1]),
        [V1l] "=&r"(V[1].u64[0]), [V1h] "=&r"(V[1].u64[1]),
        [V2l] "=&r"(V[2].u64[0]), [V2h] "=&r"(V[2].u64[1]),
        [V3l] "=&r"(V[3].u64[0]), [V3h] "=&r"(V[3].u64[1]),
        "=m"(*(unsigned char(*)[64])x), "=m"(*(unsigned char(*)[64])V)
      :
      [vsk0] "w"(vsubkeys[0]), [vsk1] "w"(vsubkeys[1]), [vsk2] "w"(vsubkeys[2]),
      [vsk3] "w"(vsubkeys[3]), [vsk4] "w"(vsubkeys[4]), [vsk5] "w"(vsubkeys[5]),
      [vsk6] "w"(vsubkeys[6]), [vsk7] "w"(vsubkeys[7]), [vsk8] "w"(vsubkeys[8]),
      [vsk9] "w"(vsubkeys[9]), [vsk10] "w"(vsubkeys[10]),
      [vsk11] "w"(vsubkeys[11]), [vsk12] "w"(vsubkeys[12]),
      [vsk13] "w"(vsubkeys[13]), [vsk14] "w"(vsubkeys[14]), [V] "r"(V),
      [DRBG_ctx_V] "m"(DRBG_ctx.V)
      : "cc");
#pragma GCC diagnostic pop

  while (xlen > 0) {
    if (xlen > 16) {
      AES256_ECB(vsubkeys, vld1q_u8((uint8_t *)&V[0]), x);
      x += 16;
      xlen -= 16;

      Vle[0].u128++;
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
