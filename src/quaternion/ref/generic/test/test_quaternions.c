#include "quaternion_tests.h"

// run all tests in module
int main(){
    int res = 0;
    printf("Running quaternion module unit tests\n");
    res = res | quat_test_finit();
    res = res | quat_test_dim4();
    res = res | quat_test_dim2();
    res = res | quat_test_matkermod();
    res = res | quat_test_integers();
    res = res | quat_test_algebra();
    res = res | quat_test_lattice();
    res = res | quat_test_lideal();
    res = res | quat_test_with_randomization();
    if(res != 0){
        printf("\nSome tests failed!\n");
    }
    return(res);
}
