#ifndef EC_PARAMS_H
#define EC_PARAMS_H

#define POWER_OF_2 97
#define POWER_OF_3 68

static digit_t TWOpF[NWORDS_ORDER] = {0x0, 0x200000000, 0x0, 0x0, 0x0, 0x0}; // 2^g
static digit_t TWOpFm1[NWORDS_ORDER] = {0x0, 0x100000000, 0x0, 0x0, 0x0, 0x0}; // 2^(g-1)
static digit_t THREEpE[NWORDS_ORDER] = {0x003B3FCEF3103289, 0x0, 0x0, 0x0, 0x0, 0x0};    // 3^e
static digit_t THREEpF[NWORDS_ORDER] = {0x58DE83727119CD51, 0x00000DB6794B8C65, 0x0, 0x0, 0x0, 0x0};    // 3^f
static digit_t THREEpFdiv2[NWORDS_ORDER] = {0xAC6F41B9388CE6A8, 0x000006DB3CA5C632, 0x0, 0x0, 0x0, 0x0};    // Floor(3^f/2)

#define scaled 1	// unscaled (0) or scaled (1) remainder tree approach
#define gap 83

#define P_LEN 7
#define M_LEN 21

static digit_t p_plus_minus_bitlength[P_LEN + M_LEN] =
{ 2, 3, 4, 6, 8, 9, 10, 3, 4, 7, 7, 8, 8, 8, 8, 9, 10, 11, 12,
  12, 13, 13, 13, 13, 14, 15, 16, 16};

static digit_t p_cofactor_for_2f[5] = {
0x9AB752342630BA61,
0x48575BA8E3917B34,
0x22E8856B32DE1705,
0x558438D163573025,
0x0000000001EFB777
};
#define P_COFACTOR_FOR_2F_BITLENGTH 281

static digit_t p_cofactor_for_3g[5] = {
        0x0000000000000000,
        0x97E5302200000000,
        0xEE54EDDA13CFE081,
        0xC6BFE309D9209F3F,
        0x000000000000484C
};

#define P_COFACTOR_FOR_3G_BITLENGTH 271


static digit_t p_cofactor_for_6fg[3] = { 
        0x09E7F040CBF29811,
        0xEC904F9FF72A76ED,
        0x00002426635FF184
};
#define P_COFACTOR_FOR_6FG_BITLENGTH 174

static int STRATEGY4[47] = { 19, 12, 7, 4, 2, 1, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 7, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1};

static int sizeI[P_LEN+M_LEN] = {
        0, 1, 2, 3, 6, 14, 14, 1, 3, 5, 5, 7, 8, 8, 10, 11, 14, 19, 26, 30, 35, 40, 44, 44, 52, 88, 114, 114
};
static int sizeJ[P_LEN+M_LEN] = {
        0, 1, 1, 3, 6, 9, 13, 1, 1, 4, 5, 6, 7, 7, 6, 10, 10, 16, 23, 28, 32, 32, 32, 33, 45, 76, 87, 104
};
static int sizeK[P_LEN+M_LEN] = {
        1, 1, 1, 5, 6, 2, 16, 0, 0, 4, 6, 2, 4, 7, 0, 1, 4, 6, 0, 5, 18, 13, 30, 2, 18, 12, 3, 8
};

#define sI_max 114
#define sJ_max 104
#define sK_max 47

#define ceil_log_sI_max 7
#define ceil_log_sJ_max 7

#endif