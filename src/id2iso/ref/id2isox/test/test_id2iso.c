
#include "id2iso_tests.h"

// run all tests in module
int main(){
    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    printf("Running id2iso module unit tests\n");

    res = res & id2iso_test_ker2id();
    res = res & id2iso_test_id2ker_even();
    res = res & id2iso_test_id2ker_odd();
    res = res & id2iso_test_id2iso();

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("All tests passed!\n");
    }
    return(!res);
}
