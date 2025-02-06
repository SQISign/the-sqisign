#include "quaternion_tests.h"

// int ibz_generate_random_prime(ibz_t *p, int is3mod4, int bitsize);
int
quat_test_ibz_generate_random_prime()
{
    int res = 0;
    int bitsize, is3mod4, primality_test_attempts;
    ibz_t p;
    ibz_init(&p);
    bitsize = 20;
    primality_test_attempts = 30;
    is3mod4 = 1;
    res = res || !ibz_generate_random_prime(&p, is3mod4, bitsize, primality_test_attempts);
    res = res || (ibz_probab_prime(&p, 20) == 0);
    res = res || (ibz_bitsize(&p) < bitsize);
    res = res || (is3mod4 && (ibz_get(&p) % 4 != 3));
    bitsize = 30;
    is3mod4 = 0;
    res = res || !ibz_generate_random_prime(&p, is3mod4, bitsize, primality_test_attempts);
    res = res || (ibz_probab_prime(&p, 20) == 0);
    res = res || (ibz_bitsize(&p) < bitsize);
    res = res || (is3mod4 && (ibz_get(&p) % 4 != 3));
    is3mod4 = 1;
    res = res || !ibz_generate_random_prime(&p, is3mod4, bitsize, primality_test_attempts);
    res = res || (ibz_probab_prime(&p, 20) == 0);
    res = res || (ibz_bitsize(&p) < bitsize);
    res = res || (is3mod4 && (ibz_get(&p) % 4 != 3));
    if (res) {
        printf("Quaternion unit test ibz_generate_random_prime failed\n");
    }
    ibz_finalize(&p);
    return (res);
}

// int ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p);
int
quat_test_integer_ibz_cornacchia_prime(void)
{
    int res = 0;
    ibz_t x, y, n, prod, c_res, p;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&n);
    ibz_init(&p);
    ibz_init(&prod);
    ibz_init(&c_res);

    ibz_set(&n, 1);
    // there is a solution in these cases
    ibz_set(&p, 5);
    if (ibz_cornacchia_prime(&x, &y, &n, &p)) {
        ibz_mul(&c_res, &x, &x);
        ibz_mul(&prod, &y, &y);
        ibz_mul(&prod, &prod, &n);
        ibz_add(&c_res, &c_res, &prod);
        res = res || ibz_cmp(&p, &c_res);
    } else {
        res = 1;
    }
    ibz_set(&p, 2);
    if (ibz_cornacchia_prime(&x, &y, &n, &p)) {
        ibz_mul(&c_res, &x, &x);
        ibz_mul(&prod, &y, &y);
        ibz_mul(&prod, &prod, &n);
        ibz_add(&c_res, &c_res, &prod);
        res = res || ibz_cmp(&p, &c_res);
    } else {
        res = 1;
    }
    ibz_set(&p, 41);
    if (ibz_cornacchia_prime(&x, &y, &n, &p)) {
        ibz_mul(&c_res, &x, &x);
        ibz_mul(&prod, &y, &y);
        ibz_mul(&prod, &prod, &n);
        ibz_add(&c_res, &c_res, &prod);
        res = res || ibz_cmp(&p, &c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 2);
    ibz_set(&p, 3);
    if (ibz_cornacchia_prime(&x, &y, &n, &p)) {
        ibz_mul(&c_res, &x, &x);
        ibz_mul(&prod, &y, &y);
        ibz_mul(&prod, &prod, &n);
        ibz_add(&c_res, &c_res, &prod);
        res = res || ibz_cmp(&p, &c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 3);
    ibz_set(&p, 7);
    if (ibz_cornacchia_prime(&x, &y, &n, &p)) {
        ibz_mul(&c_res, &x, &x);
        ibz_mul(&prod, &y, &y);
        ibz_mul(&prod, &prod, &n);
        ibz_add(&c_res, &c_res, &prod);
        res = res || ibz_cmp(&p, &c_res);
    } else {
        res = 1;
    }

    ibz_set(&n, 1);
    // there is no solution in these cases
    ibz_set(&p, 7);
    res = res || ibz_cornacchia_prime(&x, &y, &n, &p);
    ibz_set(&p, 3);
    res = res || ibz_cornacchia_prime(&x, &y, &n, &p);
    ibz_set(&n, 3);
    ibz_set(&p, 5);
    res = res || ibz_cornacchia_prime(&x, &y, &n, &p);
    // This should be solved
    ibz_set(&n, 3);
    ibz_set(&p, 3);
    res = res || !ibz_cornacchia_prime(&x, &y, &n, &p);
    if (res != 0) {
        printf("Quaternion unit test integer_ibz_cornacchia_prime failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&n);
    ibz_finalize(&p);
    ibz_finalize(&prod);
    ibz_finalize(&c_res);
    return res;
}

// run all previous tests
int
quat_test_integers(void)
{
    int res = 0;
    printf("\nRunning quaternion tests of integer functions\n");
    res = res | quat_test_ibz_generate_random_prime();
    res = res | quat_test_integer_ibz_cornacchia_prime();
    return (res);
}
