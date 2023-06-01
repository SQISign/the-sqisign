#include "quaternion_tests.h"


//integer helpers
//void ibz_rounded_div(ibz_t *q, const ibz_t *a, const ibz_t *b);
int quat_test_integer_ibz_rounded_div(){
    int res = 0;
    ibz_t q, a, b;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&q);

    // basic tests
    ibz_set(&a,15);
    ibz_set(&b,3);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==5);
    ibz_set(&a,16);
    ibz_set(&b,3);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==5);
    ibz_set(&a,17);
    ibz_set(&b,3);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==6);
    ibz_set(&a,37);
    ibz_set(&b,5);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==7);
    // test sign combination
    ibz_set(&a,149);
    ibz_set(&b,12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==12);
    ibz_set(&a,149);
    ibz_set(&b,-12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-12);
    ibz_set(&a,-149);
    ibz_set(&b,-12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==12);
    ibz_set(&a,-149);
    ibz_set(&b,12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-12);
    ibz_set(&a,151);
    ibz_set(&b,12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==13);
    ibz_set(&a,-151);
    ibz_set(&b,-12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==13);
    ibz_set(&a,151);
    ibz_set(&b,-12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-13);
    ibz_set(&a,-151);
    ibz_set(&b,12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-13);
    //divisibles with sign
    ibz_set(&a,144);
    ibz_set(&b,12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==12);
    ibz_set(&a,-144);
    ibz_set(&b,-12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==12);
    ibz_set(&a,144);
    ibz_set(&b,-12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-12);
    ibz_set(&a,-144);
    ibz_set(&b,12);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-12);
    // tests close to 0
    ibz_set(&a,-12);
    ibz_set(&b,-25);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==0);
    ibz_set(&a,12);
    ibz_set(&b,25);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==0);
    ibz_set(&a,-12);
    ibz_set(&b,25);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==0);
    ibz_set(&a,12);
    ibz_set(&b,-25);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==0);
    ibz_set(&a,-12);
    ibz_set(&b,-23);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==1);
    ibz_set(&a,12);
    ibz_set(&b,23);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==1);
    ibz_set(&a,-12);
    ibz_set(&b,23);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-1);
    ibz_set(&a,12);
    ibz_set(&b,-23);
    ibz_rounded_div(&q,&a,&b);
    res = res || !(ibz_get(&q)==-1);
    // test output equal input
    ibz_set(&a,-151);
    ibz_set(&b,12);
    ibz_rounded_div(&a,&a,&b);
    res = res || !(ibz_get(&a)==-13);
    ibz_set(&a,-151);
    ibz_set(&b,12);
    ibz_rounded_div(&b,&a,&b);
    res = res || !(ibz_get(&b)==-13);
    // test for cmp not returning 1 or -1 or 0
    ibz_set(&a,4292606433816540);
    ibz_set(&b,864673106105940);
    ibz_rounded_div(&b,&a,&b);
    res = res || !(ibz_get(&b)==5);

    if (res != 0){
        printf("Quaternion unit test integer_ibz_rounded_div failed\n");
    }
    ibz_finalize(&a);
    ibz_finalize(&b);
    ibz_finalize(&q);
    return(res);
}

// tests for cornacchia helper functions
//void ibz_complex_mul(ibz_t *re_res, ibz_t *im_res, const ibz_t *re_a, const ibz_t *im_a, const ibz_t *re_b, const ibz_t *im_b);
int quat_test_integer_ibz_complex_mul(){
    int res = 0;
    ibz_t re_res, re_a, re_b, re_cmp, im_res, im_a, im_b, im_cmp;
    ibz_init(&re_res);
    ibz_init(&re_a);
    ibz_init(&re_b);
    ibz_init(&re_cmp);
    ibz_init(&im_res);
    ibz_init(&im_a);
    ibz_init(&im_b);
    ibz_init(&im_cmp);

    ibz_set(&re_a, 1);
    ibz_set(&im_a, 2);
    ibz_set(&re_b, 3);
    ibz_set(&im_b, -4);
    ibz_set(&re_cmp, 11);
    ibz_set(&im_cmp, 2);
    ibz_complex_mul(&re_res,&im_res,&re_a,&im_a,&re_b,&im_b);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    ibz_set(&re_a, -3);
    ibz_set(&im_a, -5);
    ibz_set(&re_b, 2);
    ibz_set(&im_b, 4);
    ibz_set(&re_cmp, 14);
    ibz_set(&im_cmp, -22);
    ibz_complex_mul(&re_res,&im_res,&re_a,&im_a,&re_b,&im_b);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    if (res != 0){
        printf("Quaternion unit test integer_ibz_complex_mul failed\n");
    }
    ibz_finalize(&re_res);
    ibz_finalize(&re_cmp);
    ibz_finalize(&re_a);
    ibz_finalize(&re_b);
    ibz_finalize(&im_res);
    ibz_finalize(&im_cmp);
    ibz_finalize(&im_a);
    ibz_finalize(&im_b);
    return(res);
}

//void ibz_complex_mul_by_complex_power(ibz_t *re_res, ibz_t *im_res, const ibz_t *re_a, const ibz_t *im_a, int64_t exp);
int quat_test_integer_ibz_complex_mul_by_complex_power(){
    int res = 0;
    int64_t exp;
    ibz_t re_res, re_a, re_cmp, im_res, im_a, im_cmp;
    ibz_init(&re_res);
    ibz_init(&re_a);
    ibz_init(&re_cmp);
    ibz_init(&im_res);
    ibz_init(&im_a);
    ibz_init(&im_cmp);

    exp = 0;
    ibz_set(&re_a, 1);
    ibz_set(&im_a, 2);
    ibz_set(&re_res, 3);
    ibz_set(&im_res, -4);
    ibz_set(&re_cmp, 3);
    ibz_set(&im_cmp, -4);
    ibz_complex_mul_by_complex_power(&re_res,&im_res,&re_a,&im_a,exp);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    exp = 1;
    ibz_set(&re_a, 1);
    ibz_set(&im_a, 2);
    ibz_set(&re_res, 3);
    ibz_set(&im_res, -4);
    ibz_set(&re_cmp, 11);
    ibz_set(&im_cmp, 2);
    ibz_complex_mul_by_complex_power(&re_res,&im_res,&re_a,&im_a,exp);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    exp = 2;
    ibz_set(&re_a, 1);
    ibz_set(&im_a, 2);
    ibz_set(&re_res, 3);
    ibz_set(&im_res, -4);
    ibz_set(&re_cmp, 7);
    ibz_set(&im_cmp, 24);
    ibz_complex_mul_by_complex_power(&re_res,&im_res,&re_a,&im_a,exp);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    if (res != 0){
        printf("Quaternion unit test integer_ibz_complex_mul_by_complex_power failed\n");
    }
    ibz_finalize(&re_res);
    ibz_finalize(&re_cmp);
    ibz_finalize(&re_a);
    ibz_finalize(&im_res);
    ibz_finalize(&im_cmp);
    ibz_finalize(&im_a);
    return(res);
}

//int ibz_cornacchia_extended_prime_loop(ibz_t *re_res, ibz_t *im_res, int64_t prime, int64_t val);
int quat_test_integer_ibz_cornacchia_extended_prime_loop(){
    int res = 0;
    int64_t p, re_a, im_a;
    ibz_t re_res, re_cmp, im_res, im_cmp, prod;
    ibz_init(&re_res);
    ibz_init(&re_cmp);
    ibz_init(&im_res);
    ibz_init(&im_cmp);
    ibz_init(&prod);

    p = 5;
    ibz_set(&re_res, 1);
    ibz_set(&im_res, 1);
    ibz_set(&re_cmp, -1);
    ibz_set(&im_cmp, 7);
    ibz_cornacchia_extended_prime_loop(&re_res, &im_res, p, 2);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    p = 7;
    ibz_set(&re_res, -1);
    ibz_set(&im_res, 7);
    ibz_set(&re_cmp, -1);
    ibz_set(&im_cmp, 7);
    ibz_cornacchia_extended_prime_loop(&re_res, &im_res, p, 0);
    res = res || ibz_cmp(&re_res,&re_cmp);
    res = res || ibz_cmp(&im_res,&im_cmp);

    if (res != 0){
        printf("Quaternion unit test integer_ibz_cornacchia_extended_prime_loop failed\n");
    }
    ibz_finalize(&re_res);
    ibz_finalize(&re_cmp);
    ibz_finalize(&im_res);
    ibz_finalize(&im_cmp);
    ibz_finalize(&prod);
    return(res);
}

//int ibz_cornacchia_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p);
int quat_test_integer_ibz_cornacchia_prime(){
    int res = 0;
    ibz_t x,y,n, prod,c_res,p;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&n);
    ibz_init(&p);
    ibz_init(&prod);
    ibz_init(&c_res);

    ibz_set(&n, 1);
    // there is a solution in these cases
    ibz_set(&p, 5);
    if(ibz_cornacchia_prime(&x,&y,&n,&p)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&p,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&p, 2);
    if(ibz_cornacchia_prime(&x,&y,&n,&p)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&p,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&p, 41);
    if(ibz_cornacchia_prime(&x,&y,&n,&p)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&p,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 2);
    ibz_set(&p, 3);
    if(ibz_cornacchia_prime(&x,&y,&n,&p)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&p,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 3);
    ibz_set(&p, 7);
    if(ibz_cornacchia_prime(&x,&y,&n,&p)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&p,&c_res);
    } else {
        res = 1;
    }

    ibz_set(&n, 1);
    // there is no solution in these cases
    ibz_set(&p, 7);
    res = res || ibz_cornacchia_prime(&x,&y,&n,&p);
    ibz_set(&p, 3);
    res = res || ibz_cornacchia_prime(&x,&y,&n,&p);
    ibz_set(&n, 3);
    ibz_set(&p, 5);
    res = res || ibz_cornacchia_prime(&x,&y,&n,&p);
    if (res != 0){
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


//int ibz_cornacchia_special_prime(ibz_t *x, ibz_t *y, const ibz_t *n, const ibz_t *p, const int exp_adjust);
int quat_test_integer_ibz_cornacchia_special_prime(){
    int res = 0;
    int exp_adjust;
    ibz_t x,y,n, prod,c_res,p, cmp, two;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&n);
    ibz_init(&p);
    ibz_init(&prod);
    ibz_init(&c_res);
    ibz_init(&cmp);
    ibz_init(&two);
    ibz_set(&two,2);

    ibz_set(&n, 3);
    // there is a solution in these cases
    ibz_set(&p, 43);
    exp_adjust = 2;
    //x^2 + 3y^2 = 4*43, (5,7) is solution
    if(ibz_cornacchia_special_prime(&x,&y,&n,&p, exp_adjust)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        ibz_pow(&cmp,&two,exp_adjust);
        ibz_mul(&cmp,&cmp,&p);
        res = res || ibz_cmp(&cmp,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 7);
    exp_adjust = 3;
    ibz_set(&p, 11);
    //x^2 + 7y^2 = 8*11, (5,3) is solution
    if(ibz_cornacchia_special_prime(&x,&y,&n,&p, exp_adjust)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        ibz_pow(&cmp,&two,exp_adjust);
        ibz_mul(&cmp,&cmp,&p);
        res = res || ibz_cmp(&cmp,&c_res);
    } else {
        res = 1;
    }
    /* Failing testcase, needs investigation
    ibz_set(&n, 15);
    exp_adjust = 4;
    ibz_set(&p, 31);
    //x^2 + 15y^2 = 16*31, (19,3) is solution
    if(ibz_cornacchia_special_prime(&x,&y,&n,&p, exp_adjust)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        ibz_pow(&cmp,&two,exp_adjust);
        ibz_mul(&cmp,&cmp,&p);
        res = res || ibz_cmp(&cmp,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 3);
    exp_adjust = 2;
    ibz_set(&p, 7);
    //x^2 + 3y^2 = 4*7, (1,3) is solution
    if(ibz_cornacchia_special_prime(&x,&y,&n,&p, exp_adjust)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_mul(&prod,&prod,&n);
        ibz_add(&c_res,&c_res,&prod);
        ibz_pow(&cmp,&two,exp_adjust);
        ibz_mul(&cmp,&cmp,&p);
        res = res || ibz_cmp(&cmp,&c_res);
    } else {
        res = 1;
    }
    */

    // there is no solution in these cases
    ibz_set(&n, 3);
    ibz_set(&p, 11);
    exp_adjust = 3;
    res = res || ibz_cornacchia_special_prime(&x,&y,&n,&p, exp_adjust);

    if (res != 0){
        printf("Quaternion unit test integer_ibz_cornacchia_special_prime failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&n);
    ibz_finalize(&p);
    ibz_finalize(&prod);
    ibz_finalize(&c_res);
    ibz_finalize(&two);
    ibz_finalize(&cmp);
    return res;
}

//int ibz_cornacchia_extended(ibz_t *x, ibz_t *y, const ibz_t *n, const short *prime_list, const int prime_list_length, short primality_test_iterations, const ibz_t *bad_primes_prod); 
int quat_test_integer_ibz_cornacchia_extended(){
    int res = 0;
    ibz_t x,y,n, prod,c_res,bad;
    short primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101};	 	 	 	 	 
    int primes_length = 26;
    short iterations = 20;
    ibz_init(&x);
    ibz_init(&y);
    ibz_init(&n);
    ibz_init(&prod);
    ibz_init(&c_res);
    ibz_init(&bad);
    ibz_set(&bad,3*7*11*19);

    // there is a solution in these cases
    ibz_set(&n, 5);
    if(ibz_cornacchia_extended(&x,&y,&n, primes, 4,iterations, NULL)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&n,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 50);
    if(ibz_cornacchia_extended(&x,&y,&n, primes, 7,iterations,NULL)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&n,&c_res);
    } else {
        res = 1;
    }
    ibz_set(&n, 4100);
    if(ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,&bad)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&n,&c_res);
    } else {
        res = 1;
    }
    // test product of all small primes
    ibz_set(&n, 5*13*17*29*37*41);
    ibz_set(&x, 53*61*73*89*97);
    ibz_mul(&n, &n, &x);
    ibz_set(&x, 404);
    ibz_mul(&n, &n, &x);
    if(ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,NULL)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&n,&c_res);
    } else {
        res = 1;
    }
    // test with large prime
    ibz_set(&n, 1381); // prime and 1 mod 4
    if(ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,&bad)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&n,&c_res);
    } else {
        res = 1;
    }
    // test with large prime part
    ibz_set(&n, 5*13*17*29*37*97);
    ibz_set(&x, 1381); // prime and 1 mod 4
    ibz_mul(&n, &n, &x);
    if(ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,NULL)){
        ibz_mul(&c_res,&x,&x);
        ibz_mul(&prod,&y,&y);
        ibz_add(&c_res,&c_res,&prod);
        res = res || ibz_cmp(&n,&c_res);
    } else {
        res = 1;
    }

    // there is no solution in these cases
    ibz_set(&n, 7);
    res = res || ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,&bad);
    ibz_set(&n, 3);
    res = res || ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,&bad);
    ibz_set(&n, 6);
    res = res || ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,NULL);
    ibz_set(&n, 30);
    res = res || ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,&bad);
    ibz_set(&n, 30*1381);
    res = res || ibz_cornacchia_extended(&x,&y,&n, primes, primes_length,iterations,&bad);

    if (res != 0){
        printf("Quaternion unit test integer_ibz_cornacchia_extended failed\n");
    }
    ibz_finalize(&x);
    ibz_finalize(&y);
    ibz_finalize(&n);
    ibz_finalize(&prod);
    ibz_finalize(&c_res);
    ibz_finalize(&bad);
    return res;
}


// run all previous tests
int quat_test_integers(){
    int res = 0;
    printf("\nRunning quaternion tests of integer functions\n");
    res = res | quat_test_integer_ibz_rounded_div();
    res = res | quat_test_integer_ibz_complex_mul();
    res = res | quat_test_integer_ibz_complex_mul_by_complex_power();
    res = res | quat_test_integer_ibz_cornacchia_extended_prime_loop();
    res = res | quat_test_integer_ibz_cornacchia_prime();
    res = res | quat_test_integer_ibz_cornacchia_special_prime();
    res = res | quat_test_integer_ibz_cornacchia_extended();
    return(res);
}

