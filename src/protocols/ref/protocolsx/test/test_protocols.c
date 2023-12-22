#include <rng.h>
#include <stdio.h>
#include <ec.h>
#include <inttypes.h>

#include "test_protocols.h"




//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf(HEX_FS, b.im[i]);
    printf("\n");
}

void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}

void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}
//XXX


bool curve_is_canonical(ec_curve_t const *E)
{
    ec_curve_t EE;
    ec_isom_t isom;
    ec_curve_normalize(&EE, &isom, E);

    fp2_t lhs, rhs;
    fp2_mul(&lhs, &E->A, &EE.C);
    fp2_mul(&rhs, &E->C, &EE.A);
    return fp2_is_equal(&lhs, &rhs);
}


// run all tests in module
int main(){
    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    // printf("\nRunning encoding tests\n");
    // res &= test_encode();

    printf("\nRunning protocols module unit tests\n \n");

    printf("\nRunning keygen tests\n \n");
    res &= test_keygen();

    printf("\nRunning signature tests\n \n");
    res &= test_commitment();
    res &= test_challenge();
    res &= test_sign_verif();


    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("All tests passed!\n");
    }
    return(!res);
}
