#include <assert.h>
#include "include/fp.h"

#ifdef RADIX_32
const digit_t p[NWORDS_FIELD] =  { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xA6E1FFFF, 0x994C68AD, 0x781974CE, 0xFAF0A29A, 0x4A0DEA65, 0xFE3AC590, 0x26507D01, 0x02BDBE63, 0x6936E792, 0x8C15B003, 0xA8869BC6, 0x00255946 };
const digit_t R2[NWORDS_FIELD] = { 0xC7549CBD, 0x46E4E8A0, 0x43E89EA5, 0xCB993B59, 0x2F1B55C8, 0x545AC09F, 0xACAA06EC, 0x1ADB99DD, 0x55D8B8D4, 0x87994B89, 0x2F9E57C8, 0x2CC2EA62, 0xDAF1003C, 0x2780B5F2, 0x6B8674B8, 0x00169167 };
const digit_t pp[NWORDS_FIELD] = { 0x00000001, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000 };
#elif defined(RADIX_64)
const digit_t p[NWORDS_FIELD] =  { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x994C68ADA6E1FFFF, 0xFAF0A29A781974CE, 0xFE3AC5904A0DEA65, 0x02BDBE6326507D01, 0x8C15B0036936E792, 0x00255946A8869BC6 };
const digit_t R2[NWORDS_FIELD] = { 0x46E4E8A0C7549CBD, 0xCB993B5943E89EA5, 0x545AC09F2F1B55C8, 0x1ADB99DDACAA06EC, 0x87994B8955D8B8D4, 0x2CC2EA622F9E57C8, 0x2780B5F2DAF1003C, 0x001691676B8674B8 };
const digit_t pp[NWORDS_FIELD] = { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 };
#endif

void fp_set(digit_t* x, const digit_t val)
{ // Set field element x = val, where val has wordsize

    x[0] = val;
    for (unsigned int i = 1; i < NWORDS_FIELD; i++) {
        x[i] = 0;
    }
}

bool fp_is_equal(const digit_t* a, const digit_t* b)
{ // Compare two field elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ b[i];

    return (bool)is_digit_zero_ct(r);
}

bool fp_is_zero(const digit_t* a)
{ // Is a field element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void fp_copy(digit_t* out, const digit_t* a)
{
    memcpy(out, a, NWORDS_FIELD*RADIX/8);
}

void fp_neg(digit_t* out, const digit_t* a)
{ // Modular negation, out = -a mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(out[i], borrow, ((digit_t*)p)[i], a[i], borrow);
    }
    fp_sub(out, out, (digit_t*)p);
}

void MUL(digit_t* out, const digit_t a, const digit_t b)
{ // Digit multiplication, digit*digit -> 2-digit result 
  // Inputs: a, b in [0, 2^w-1], where w is the computer wordsize 
  // Output: 0 < out < 2^(2w)-1    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t)*4);            // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t)*4);

    albl = al * bl;
    albh = al * bh;
    ahbl = ah * bl;
    ahbh = ah * bh;
    out[0] = albl & mask_low;                 // out00

    res1 = albl >> (sizeof(digit_t)*4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t)*4);
    out[0] ^= temp << (sizeof(digit_t)*4);    // out01   

    res1 = ahbl >> (sizeof(digit_t)*4);
    res2 = albh >> (sizeof(digit_t)*4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    out[1] = temp & mask_low;                 // out10 
    carry = temp & mask_high;
    out[1] ^= (ahbh & mask_high) + carry;     // out11
}

digit_t mp_shiftr(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision right shift
    digit_t bit_out = x[0] & 1;

    for (unsigned int i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], shift, x[i], RADIX);
    }
    x[nwords-1] >>= shift;
    return bit_out;
}

void mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision left shift

    assert(shift < RADIX*nwords);

    int shift_words = shift / RADIX;

    for (int i = nwords-1; i > shift_words; i--) {
        SHIFTL(x[i-shift_words], x[i-shift_words-1], shift % RADIX, x[i], RADIX);
    }

    x[shift_words] = x[0] << (shift % RADIX);

    for (int i = shift_words - 1; i >= 0; i--) {
        x[i] = 0;
    }
}

static void fp_exp3div4(digit_t* out, const digit_t* a)
{ // Fixed exponentiation out = a^((p-3)/4) mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
  // Requirement: p = 3(mod 4)
    fp_t p_t, acc;
    digit_t bit;

    memcpy((digit_t*)p_t, (digit_t*)p, NWORDS_FIELD*RADIX/8);
    memcpy((digit_t*)acc, (digit_t*)a, NWORDS_FIELD*RADIX/8);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    fp_set(out, 1);
    fp_tomont(out, out);

    for (int i = 0; i < NWORDS_FIELD*RADIX-2; i++) {
        bit = p_t[0] & 1;
        mp_shiftr(p_t, 1, NWORDS_FIELD);
        if (bit == 1) {
            fp_mul(out, out, acc);
        }
        fp_sqr(acc, acc);
    }
}

void fp_inv(digit_t* a)
{ // Modular inversion, out = x^-1*R mod p, where R = 2^(w*nwords), w is the computer wordsize and nwords is the number of words to represent p
  // Input: a=xR in [0, p-1] 
  // Output: out in [0, p-1]. It outputs 0 if the input does not have an inverse  
  // Requirement: Ceiling(Log(p)) < w*nwords
    fp_t t;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_sqr(t, t);
    fp_mul(a, t, a);    // a^(p-2)
}

bool fp_is_square(const digit_t* a)
{ // Is field element a square?
  // Output: out = 0 (false), 1 (true)
    fp_t t, one;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_mul(t, t, a);    // a^((p-1)/2)
    fp_frommont(t, t);
    fp_set(one, 1);

    return fp_is_equal(t, one);
}

void fp_sqrt(digit_t* a)
{ // Square root computation, out = a^((p+1)/4) mod p
    fp_t t;

    fp_exp3div4(t, a);
    fp_mul(a, t, a);    // a^((p+1)/4)
}
