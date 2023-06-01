#ifndef EC_PARAMS_H
#define EC_PARAMS_H

#define POWER_OF_2 145
#define POWER_OF_3 72

static digit_t TWOpF[NWORDS_ORDER] = {0x0, 0x0, 0x20000, 0x0, 0x0, 0x0, 0x0, 0x0}; // 2^g
static digit_t TWOpFm1[NWORDS_ORDER] = {0x0, 0x0, 0x10000, 0x0, 0x0, 0x0, 0x0, 0x0}; // 2^(g-1)
static digit_t THREEpE[NWORDS_ORDER] = {0x2153E468B91C6D1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};    // 3^e
static digit_t THREEpF[NWORDS_ORDER] = {0x1E679735C929F6A1, 0x000456BC60E76C11, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};    // 3^f
static digit_t THREEpFdiv2[NWORDS_ORDER] = {0x8F33CB9AE494FB50, 0x00022B5E3073B608, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};    // Floor(3^f/2)

#define scaled 1	// unscaled (0) or scaled (1) remainder tree approach
#define gap 83

#define P_LEN 7
#define M_LEN 27

static long p_plus_minus_bitlength[P_LEN + M_LEN] ={
        2, 4, 6, 7, 7, 9, 10, 3, 3, 5, 6, 6, 7, 7, 8, 10,
        10, 10, 10, 12, 13, 13, 13, 14, 14, 14, 14, 14,
        15, 15, 16, 16, 16, 19
};

static digit_t p_cofactor_for_2f[6] = {
    0xBA674CA63456D371,
    0xF532FD78514D3C0C,
    0x3E80FF1D62C82506,
    0x73C9015EDF319328,
    0x4DE3460AD801B49B,
    0x00000012ACA35443
};
#define P_COFACTOR_FOR_2F_BITLENGTH 357

static digit_t p_cofactor_for_3g[7] = {
        0x0000000000000000,
        0x0000000000000000,
        0x0F469D6875A20000,
        0x10CAFA623AD43E00,
        0x27F3E11532C58F78,
        0x9BA98880DA1DC006,
        0x0000000000000008

};

#define P_COFACTOR_FOR_3G_BITLENGTH 388


static digit_t p_cofactor_for_6fg[4] = {
        0x1F0007A34EB43AD1,
        0xC7BC08657D311D6A,
        0xE00313F9F08A9962,
        0x00044DD4C4406D0E
};
#define P_COFACTOR_FOR_6FG_BITLENGTH 243

static int STRATEGY4[71] = { 28, 16, 11, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 4, 3, 2, 1, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1};

static int sizeI[P_LEN+M_LEN] = {
        0, 3, 3, 4, 6, 11, 14, 1, 1, 2, 3, 4, 4, 6, 6, 13, 14, 14, 18, 26, 32, 38, 48, 54, 60, 62, 62, 64, 72, 92, 118, 126, 130, 310
};
static int sizeJ[P_LEN+M_LEN] = {
        0, 1, 3, 4, 4, 10, 14, 1, 1, 2, 3, 3, 4, 5, 6, 12, 13, 13, 14, 24, 32, 33, 37, 49, 56, 60, 62, 62, 64, 63, 102, 125, 126, 256
};
static int sizeK[P_LEN+M_LEN] = {
        1, 0, 2, 1, 3, 10, 21, 0, 1, 0, 0, 2, 4, 3, 3, 9, 2, 5, 0, 21, 28, 21, 11, 6, 75, 21, 82, 59, 75, 21, 21, 123, 0, 396
};

#define sI_max 310
#define sJ_max 256
#define sK_max 396

#define ceil_log_sI_max 9
#define ceil_log_sJ_max 9

#endif