#include <bench.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <ec.h>

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)

/******************************
Util functions
******************************/

int
cmp_u64(const void *v1, const void *v2)
{
    uint64_t x1 = *(const uint64_t *)v1;
    uint64_t x2 = *(const uint64_t *)v2;
    if (x1 < x2) {
        return -1;
    } else if (x1 == x2) {
        return 0;
    } else {
        return 1;
    }
}

void
bench_basis_generation(unsigned int n, int iterations)
{
    int i, j;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[20];

    ec_basis_t basis;
    ec_curve_t curve;
    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Full even torsion generation without hints
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (j = 0; j < iterations; j++) {
            (void)ec_curve_to_basis_2f_to_hint(&basis, &curve, n);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  2^%d torsion generation takes .................................... %" PRIu64 " cycles\n",
           n,
           cycle_runs[4] / (iterations));
}

void
bench_basis_generation_from_hint(unsigned int n, int iterations)
{
    int i, j;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[20];

    ec_basis_t basis;
    ec_curve_t curve;
    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    uint8_t hint = ec_curve_to_basis_2f_to_hint(&basis, &curve, n);

    // Full even torsion generation without hints
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (j = 0; j < iterations; j++) {
            ec_curve_to_basis_2f_from_hint(&basis, &curve, n, hint);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  2^%d torsion generation takes .................................... %" PRIu64 " cycles\n",
           n,
           cycle_runs[4] / (iterations));
}

void
bench_basis(int iterations)
{
    printf("\n-------------------------------------------------------------------------------------"
           "-------------------\n\n");
    printf("Benchmarking E[2^n] basis generation for " STRINGIFY(SQISIGN_VARIANT) ": \n\n");
    bench_basis_generation(TORSION_EVEN_POWER, iterations);
    bench_basis_generation(128, iterations);

    printf("\nBenchmarking E[2^n] basis generation with hint for " STRINGIFY(SQISIGN_VARIANT) ": \n\n");
    bench_basis_generation_from_hint(TORSION_EVEN_POWER, iterations);
    bench_basis_generation_from_hint(128, iterations);
}

int
main(int argc, char *argv[])
{
    int iterations = 100 * SQISIGN_TEST_REPS;
    int help = 0;

#ifndef NDEBUG
    fprintf(stderr,
            "\x1b[31mIt looks like SQIsign was compiled with assertions enabled.\n"
            "This will severely impact performance measurements.\x1b[0m\n");
#endif

    for (int i = 1; i < argc; i++) {
        if (!help && strcmp(argv[i], "--help") == 0) {
            help = 1;
            continue;
        }

        if (sscanf(argv[i], "--iterations=%d", &iterations) == 1) {
            continue;
        }
    }

    if (help || iterations <= 0) {
        printf("Usage: %s [--iterations=<iterations>]\n", argv[0]);
        printf("Where <iterations> is the number of iterations used for benchmarking; if not "
               "present, uses the default: %d)\n",
               iterations);
        return 1;
    }

    cpucycles_init();

    bench_basis(iterations);
    return 0;
}
