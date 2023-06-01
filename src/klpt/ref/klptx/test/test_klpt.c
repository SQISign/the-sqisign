
#include "klpt_tests.h"

// run all tests in module
int main(){
    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    printf("\nRunning klpt module unit tests\n \n");

    res = res & klpt_test_tools();
    res = res & klpt_test_equiv();
    res = res & klpt_test_eichler();
    res = res & klpt_test_klpt(); 
    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("All tests passed!\n");
    }
    return(!res);
}
