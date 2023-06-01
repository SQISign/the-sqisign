#include <ec.h>
#include <quaternion.h>
#include <protocols.h>
#include "test_protocols.h"

bool curve_is_canonical(ec_curve_t const *E);   // test_protocols.c

int test_inner_keygen() {

    int res = 1;
    public_key_t pk;
    secret_key_t sk;
    secret_key_init(&sk);

    res &= protocols_keygen(&pk, &sk);
    res &= curve_is_canonical(&sk.curve);
    res &= curve_is_canonical(&pk.E);

    secret_key_finalize(&sk);
    return res;
}

int test_keygen() {

    int res = 1;

    for (int i = 0; i<1; i++) {
        res &= test_inner_keygen();
    }

    return 1;
}
