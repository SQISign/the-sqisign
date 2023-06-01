#ifndef EC_PARAMS_H
#define EC_PARAMS_H

#include <fp_constants.h>

#define POWER_OF_2 75
#define POWER_OF_3 36

static digit_t TWOpF[NWORDS_ORDER] = {0x0, 0x0800, 0x0, 0x0}; // 2^g
static digit_t TWOpFm1[NWORDS_ORDER] = {0x0, 0x0400, 0x0, 0x0}; // 2^(g-1)
static digit_t THREEpE[NWORDS_ORDER] = {0x0000000017179149, 0x0, 0x0, 0x0};    // 3^e
static digit_t THREEpF[NWORDS_ORDER] = {0x02153E468B91C6D1, 0x0, 0x0, 0x0};    // 3^f
static digit_t THREEpFdiv2[NWORDS_ORDER] = {0x010A9F2345C8E368, 0x0, 0x0, 0x0};    // Floor(3^f/2)

#define scaled 1	// unscaled (0) or scaled (1) remainder tree approach
#define gap 83

#define P_LEN 9
#define M_LEN 19

static digit_t p_plus_minus_bitlength[P_LEN + M_LEN] =
        { 2,5,6,7,7,8,9,10,11,3,4,4,6,7,7,7,8,8,8,8,9,9,9,10,11,11,11,11 };

static digit_t p_cofactor_for_2f[3] = { 0x86e4a593c926aa29,0x318674d50cb0e80e,0x00069c53c50d72bb };
#define P_COFACTOR_FOR_2F_BITLENGTH 179

static digit_t p_cofactor_for_3g[4] = { 0x0000000000000000,0x74f9dace0d9ec800,0x63a25b437f655001,0x0000000000000019 };
#define P_COFACTOR_FOR_3G_BITLENGTH 197

static digit_t p_cofactor_for_6fg[4] = { 0x002E9F3B59C1B3D9,0x032C744B686FECAA };
#define P_COFACTOR_FOR_6FG_BITLENGTH 122

static int STRATEGY4[36] = { 15, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1};

static int sizeI[] = {
        0, 2, 4, 6, 6, 7, 12, 14, 28, 1, 2, 3, 3, 5, 6, 6, 6, 6, 9, 8, 10, 12, 12, 12, 16, 16, 18, 22
};
static int sizeJ[] = {
        0, 2, 3, 4, 4, 7, 10, 13, 17, 1, 1, 1, 3, 4, 4, 4, 5, 5, 6, 7, 9, 8, 10, 12, 16, 16, 16, 22
};
static int sizeK[] = {
        1, 3, 5, 2, 6, 0, 5, 7, 4, 1, 1, 0, 0, 4, 0, 5, 5, 8, 3, 7, 11, 2, 9, 15, 4, 12, 20, 18
};

#define sI_max 28
#define sJ_max 22
#define sK_max 59

#define ceil_log_sI_max 5
#define ceil_log_sJ_max 5

#endif