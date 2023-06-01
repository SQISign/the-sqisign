#ifndef _SDACS_H_
#define _SDACS_H_

static char SDAC_P_0[] = "0";
static char SDAC_P_1[] = "10";
static char SDAC_P_2[] = "100";
static char SDAC_P_3[] = "0100";
static char SDAC_P_4[] = "10000";
static char SDAC_P_5[] = "110000";
static char SDAC_P_6[] = "100000";
static char SDAC_P_7[] = "1100010001";
static char SDAC_P_8[] = "1001010000";
static char SDAC_P_9[] = "0101001000";
static char SDAC_P_10[] = "110110010000";
static char SDAC_P_11[] = "10000000000";
static char SDAC_P_12[] = "1010100001001000";

static char SDAC_M_0[] = "";
static char SDAC_M_1[] = "000";
static char SDAC_M_2[] = "1010";
static char SDAC_M_3[] = "100010";
static char SDAC_M_4[] = "0010000";
static char SDAC_M_5[] = "110000000";
static char SDAC_M_6[] = "1010101010";
static char SDAC_M_7[] = "1010001000";
static char SDAC_M_8[] = "1001000000";
static char SDAC_M_9[] = "0100001000";
static char SDAC_M_10[] ="101101010000"; 
static char SDAC_M_11[] = "100100010010";
static char SDAC_M_12[] = "010100011000";
static char SDAC_M_13[] = "101010000001";
static char SDAC_M_14[] = "010100001000";
static char SDAC_M_15[] = "1101010010000";
static char SDAC_M_16[] = "1001010001010";
static char SDAC_M_17[] = "101001000000101";

static char *SDACs[31] = {
	SDAC_P_0, SDAC_P_1, SDAC_P_2, SDAC_P_3, SDAC_P_4, 
	SDAC_P_5, SDAC_P_6, SDAC_P_7, SDAC_P_8, SDAC_P_9, 
	SDAC_P_10, SDAC_P_11, SDAC_P_12, 
	SDAC_M_0, SDAC_M_1, SDAC_M_2, SDAC_M_3, SDAC_M_4, 
	SDAC_M_5, SDAC_M_6, SDAC_M_7, SDAC_M_8, SDAC_M_9, 
	SDAC_M_10, SDAC_M_11, SDAC_M_12, SDAC_M_13, SDAC_M_14, 
	SDAC_M_15, SDAC_M_16, SDAC_M_17
	};

static int LENGTHS[] =	{
1, 2, 3, 4, 5, 6, 6, 10, 10, 10, 12, 11, 16, 0, 3, 4, 6, 7, 9, 10, 10, 10, 10, 12, 12, 12, 12, 12, 13, 13, 15
	};
#endif