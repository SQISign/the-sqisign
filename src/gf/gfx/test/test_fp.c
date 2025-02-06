#include <stdio.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include "test_utils.h"
#include <rng.h>
#include <bench_test_arguments.h>

bool
fp_test(int iterations)
{ // Tests for the field arithmetic
    bool OK = true;
    int n, passed;
    fp_t a, b, c, d, e, f;

    printf("\n-------------------------------------------------------------------------------------"
           "-------------------\n\n");
    printf("Testing field arithmetic over GF(p): \n\n");

    // Test equality
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);
        fp_add(&b, &a, (fp_t *)&ONE);
        fp_set_zero(&c);

        if (fp_is_equal(&a, &a) == 0) {
            passed = 0;
            break;
        }
        if (fp_is_equal(&a, &b) != 0) {
            passed = 0;
            break;
        }
        if (fp_is_equal(&c, (fp_t *)&ZERO) == 0) {
            passed = 0;
            break;
        }

        if (fp_is_zero((fp_t *)&ZERO) == 0) {
            passed = 0;
            break;
        }
        if (fp_is_zero((fp_t *)&ONE) != 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) equality tests ............................................ PASSED");
    else {
        printf("  GF(p) equality tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Test mul small
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);
        uint32_t val = rand();

        // Multiply by small value
        fp_mul_small(&b, &a, val);

        // Convert and use usual multiplication
        fp_set_small(&c, val);
        fp_mul(&d, &a, &c);

        if (fp_is_equal(&b, &d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) mul small test ............................................ PASSED");
    else {
        printf("  GF(p) mul small tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Test half
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);

        fp_add(&b, &a, &a);
        fp_half(&c, &b);
        if (fp_is_equal(&a, &c) == 0) {
            passed = 0;
            break;
        }

        fp_random_test(&a);

        fp_half(&b, &a);
        fp_add(&c, &b, &b);
        if (fp_is_equal(&a, &c) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) half test ................................................. PASSED");
    else {
        printf("  GF(p) half tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Field division by 3
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);

        fp_add(&b, &a, &a);
        fp_add(&b, &b, &a);
        fp_div3(&c, &b);

        if (fp_is_equal(&a, &c) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&a);
        fp_div3(&d, &a); // d = 0/3
        if (fp_is_zero(&d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) div by 3 tests............................................. PASSED");
    else {
        printf("  GF(p) div by 3 tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Test set small
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_set_one(&a);
        fp_add(&b, &a, &a);
        fp_set_small(&c, (uint32_t)2);
        if (fp_is_equal(&b, &c) == 0) {
            passed = 0;
            break;
        }

        fp_set_one(&a);
        fp_add(&b, &a, &a);
        fp_add(&b, &b, &b);
        fp_set_small(&c, 4);
        if (fp_is_equal(&b, &c) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) set small test ............................................ PASSED");
    else {
        printf("  GF(p) set small... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Field addition
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);
        fp_random_test(&b);
        fp_random_test(&c);
        fp_random_test(&d);

        fp_add(&d, &a, &b);
        fp_add(&e, &d, &c); // e = (a+b)+c
        fp_add(&d, &b, &c);
        fp_add(&f, &d, &a); // f = a+(b+c)
        if (fp_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp_add(&d, &a, &b); // d = a+b
        fp_add(&e, &b, &a); // e = b+a
        if (fp_is_equal(&d, &e) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&b);
        fp_add(&d, &a, &b); // d = a+0
        if (fp_is_equal(&a, &d) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&b);
        fp_neg(&d, &a);
        fp_add(&e, &a, &d); // e = a+(-a)
        if (fp_is_equal(&e, &b) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) addition tests ............................................ PASSED");
    else {
        printf("  GF(p) addition tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Field subtraction
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);
        fp_random_test(&b);
        fp_random_test(&c);
        fp_random_test(&d);

        fp_sub(&d, &a, &b);
        fp_sub(&e, &d, &c); // e = (a-b)-c
        fp_add(&d, &b, &c);
        fp_sub(&f, &a, &d); // f = a-(b+c)
        if (fp_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp_sub(&d, &a, &b); // d = a-b
        fp_sub(&e, &b, &a);
        fp_neg(&e, &e); // e = -(b-a)
        if (fp_is_equal(&d, &e) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&b);
        fp_sub(&d, &a, &b); // d = a-0
        if (fp_is_equal(&a, &d) == 0) {
            passed = 0;
            break;
        }

        fp_sub(&e, &a, &a); // e = a+(-a)
        if (fp_is_zero(&e) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) subtraction tests ......................................... PASSED");
    else {
        printf("  GF(p) subtraction tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Field multiplication
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);
        fp_random_test(&b);
        fp_random_test(&c);
        fp_mul(&d, &a, &b);
        fp_mul(&e, &d, &c); // e = (a*b)*c
        fp_mul(&d, &b, &c);
        fp_mul(&f, &d, &a); // f = a*(b*c)
        if (fp_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp_add(&d, &b, &c);
        fp_mul(&e, &a, &d); // e = a*(b+c)
        fp_mul(&d, &a, &b);
        fp_mul(&f, &a, &c);
        fp_add(&f, &d, &f); // f = a*b+a*c
        if (fp_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp_mul(&d, &a, &b); // d = a*b
        fp_mul(&e, &b, &a); // e = b*a
        if (fp_is_equal(&d, &e) == 0) {
            passed = 0;
            break;
        }

        fp_set_one(&b);
        fp_mul(&d, &a, &b); // d = a*1
        if (fp_is_equal(&a, &d) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&b);
        fp_mul(&d, &a, &b); // d = a*0
        if (fp_is_equal(&b, &d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) multiplication tests ...................................... PASSED");
    else {
        printf("  GF(p) multiplication tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Field squaring
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);

        fp_sqr(&b, &a);     // b = a^2
        fp_mul(&c, &a, &a); // c = a*a
        if (fp_is_equal(&b, &c) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&a);
        fp_sqr(&d, &a); // d = 0^2
        if (fp_is_zero(&d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) squaring tests............................................. PASSED");
    else {
        printf("  GF(p) squaring tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Field inversion
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);

        fp_copy(&b, &a);
        fp_inv(&b);
        fp_mul(&c, &a, &b); // c = a*a^-1
        if (fp_is_equal(&c, (fp_t *)&ONE) == 0) {
            passed = 0;
            break;
        }

        fp_set_zero(&a);
        fp_inv(&a); // c = 0^-1
        if (fp_is_equal(&a, (fp_t *)&ZERO) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p) inversion tests............................................ PASSED");
    else {
        printf("  GF(p) inversion tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Square root and square detection
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp_random_test(&a);

        fp_sqr(&c, &a); // c = a^2
        if (fp_is_square(&c) == 0) {
            passed = 0;
            break;
        }

        fp_sqrt(&c); // c, d = ±sqrt(c)
        fp_neg(&d, &c);
        if ((fp_is_equal(&a, &c) == 0) && (fp_is_equal(&a, &d) == 0)) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  Square root, square tests........................................ PASSED");
    else {
        printf("  Square root, square tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    return OK;
}

int
main(int argc, char *argv[])
{
    uint32_t seed[12] = { 0 };
    int iterations = 1000 * SQISIGN_TEST_REPS;
    int help = 0;
    int seed_set = 0;

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (!seed_set && !parse_seed(argv[i], seed)) {
            seed_set = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>] [--seed=<seed>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for testing; if not "
               "present, uses the default: %d)\n",
               iterations);
        printf("Where <seed> is the random seed to be used; if not present, a random seed is "
               "generated\n");
        return 1;
    }

    if (!seed_set) {
        randombytes_select((unsigned char *)seed, sizeof(seed));
    }

    print_seed(seed);

#if defined(TARGET_BIG_ENDIAN)
    for (int i = 0; i < 12; i++) {
        seed[i] = BSWAP32(seed[i]);
    }
#endif

    randombytes_init((unsigned char *)seed, NULL, 256);

    return !fp_test(iterations);
}
