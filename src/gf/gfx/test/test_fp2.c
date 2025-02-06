#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include "test_utils.h"
#include <rng.h>
#include <bench_test_arguments.h>

bool
fp2_test(int iterations)
{ // Tests for the GF(p^2) arithmetic
    bool OK = true;
    int n, passed;
    fp2_t a, b, c, d, e, f;

    printf("\n-------------------------------------------------------------------------------------"
           "-------------------\n\n");
    printf("Testing arithmetic over GF(p^2): \n\n");

    // Addition in GF(p^2)
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);
        fp2_random_test(&b);
        fp2_random_test(&c);
        fp2_random_test(&d);

        fp2_add(&d, &a, &b);
        fp2_add(&e, &d, &c); // e = (a+b)+c
        fp2_add(&d, &b, &c);
        fp2_add(&f, &d, &a); // f = a+(b+c)
        if (fp2_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp2_add(&d, &a, &b); // d = a+b
        fp2_add(&e, &b, &a); // e = b+a
        if (fp2_is_equal(&d, &e) == 0) {
            passed = 0;
            break;
        }

        fp2_set_zero(&b);
        fp2_add(&d, &a, &b); // d = a+0
        if (fp2_is_equal(&d, &a) == 0) {
            passed = 0;
            break;
        }

        fp2_neg(&d, &a);
        fp2_add(&e, &a, &d); // e = a+(-a)
        if (fp2_is_zero(&e) == 0) {
            passed = 0;
            break;
        }

        // Test add_one special method
        fp2_set_one(&b);
        fp2_add(&e, &a, &b);
        fp2_add_one(&f, &a);
        if (fp2_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p^2) addition tests ............................................ PASSED");
    else {
        printf("  GF(p^2) addition tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Subtraction in GF(p^2)
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);
        fp2_random_test(&b);
        fp2_random_test(&c);
        fp2_random_test(&d);

        fp2_sub(&d, &a, &b);
        fp2_sub(&e, &d, &c); // e = (a-b)-c
        fp2_add(&d, &b, &c);
        fp2_sub(&f, &a, &d); // f = a-(b+c)
        if (fp2_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp2_sub(&d, &a, &b); // d = a-b
        fp2_sub(&e, &b, &a);
        fp2_neg(&e, &e); // e = -(b-a)
        if (fp2_is_equal(&d, &e) == 0) {
            passed = 0;
            break;
        }

        fp2_set_zero(&b);
        fp2_sub(&d, &a, &b); // d = a-0
        if (fp2_is_equal(&d, &a) == 0) {
            passed = 0;
            break;
        }

        fp2_sub(&e, &a, &a); // e = a+(-a)
        if (fp2_is_zero(&e) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p^2) subtraction tests ......................................... PASSED");
    else {
        printf("  GF(p^2) subtraction tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Multiplication in GF(p^2)
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);
        fp2_random_test(&b);
        fp2_random_test(&c);

        fp2_mul(&d, &a, &b);
        fp2_mul(&e, &d, &c); // e = (a*b)*c
        fp2_mul(&d, &b, &c);
        fp2_mul(&f, &d, &a); // f = a*(b*c)
        if (fp2_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp2_add(&d, &b, &c);
        fp2_mul(&e, &a, &d); // e = a*(b+c)
        fp2_mul(&d, &a, &b);
        fp2_mul(&f, &a, &c);
        fp2_add(&f, &d, &f); // f = a*b+a*c
        if (fp2_is_equal(&e, &f) == 0) {
            passed = 0;
            break;
        }

        fp2_mul(&d, &a, &b); // d = a*b
        fp2_mul(&e, &b, &a); // e = b*a
        if (fp2_is_equal(&d, &e) == 0) {
            passed = 0;
            break;
        }

        fp2_set_one(&b);
        fp2_mul(&d, &a, &b); // d = a*1
        if (fp2_is_equal(&a, &d) == 0) {
            passed = 0;
            break;
        }

        fp2_set_zero(&b);
        fp2_mul(&d, &a, &b); // d = a*0
        if (fp2_is_zero(&d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p^2) multiplication tests ...................................... PASSED");
    else {
        printf("  GF(p^2) multiplication tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Test mul small
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);
        uint32_t val = rand();

        // Multiply by a small value
        fp2_mul_small(&b, &a, val);

        // Convert and multiply as a fp val
        fp2_set_small(&c, val);
        fp2_mul(&d, &a, &c);

        // Values should be the same
        if (fp2_is_equal(&b, &d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p^2) mul small test ............................................ PASSED");
    else {
        printf("  GF(p^2) mul small tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Squaring in GF(p^2)
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);
        fp2_sqr(&b, &a);     // b = a^2
        fp2_mul(&c, &a, &a); // c = a*a
        if (fp2_is_equal(&b, &c) == 0) {
            passed = 0;
            break;
        }

        fp2_set_zero(&a);
        fp2_sqr(&d, &a); // d = 0^2
        if (fp2_is_zero(&d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p^2) squaring tests............................................. PASSED");
    else {
        printf("  GF(p^2) squaring tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Inversion in GF(p^2)
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);

        fp2_set_one(&d);
        fp2_copy(&b, &a);
        fp2_inv(&a);
        fp2_mul(&c, &a, &b); // c = a*a^-1
        if (fp2_is_equal(&c, &d) == 0) {
            passed = 0;
            break;
        }

        fp2_set_zero(&a);
        fp2_set_zero(&d);
        fp2_inv(&a); // c = 0^-1
        if (fp2_is_equal(&a, &d) == 0) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  GF(p^2) inversion tests............................................ PASSED");
    else {
        printf("  GF(p^2) inversion tests... FAILED");
        printf("\n");
        return false;
    }
    printf("\n");

    // Square root and square detection in GF(p^2)
    passed = 1;
    for (n = 0; n < iterations; n++) {
        fp2_random_test(&a);

        fp2_sqr(&c, &a); // c = a^2
        if (fp2_is_square(&c) == 0) {
            passed = 0;
            break;
        }

        fp2_copy(&b, &c);
        if (!fp2_sqrt_verify(&b)) {
            passed = 0;
            break;
        }

        fp2_sqrt(&c); // c = a = sqrt(c)
        fp2_neg(&d, &c);
        if ((fp2_is_equal(&a, &c) == 0) && (fp2_is_equal(&a, &d) == 0)) {
            passed = 0;
            break;
        }
    }
    if (passed == 1)
        printf("  Square root, square tests.......................................... PASSED");
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

    return !fp2_test(iterations);
}
