// SPDX-License-Identifier: Apache-2.0
#ifndef BENCH_TEST_ARGUMENTS_H__
#define BENCH_TEST_ARGUMENTS_H__

#include <inttypes.h>
#include <stdio.h>
#include <stdint.h>

static int parse_seed(const char *arg, uint32_t *seed)
{
    if (sscanf(arg, "--seed=%u", &seed[0]) == 1)
        return 0;

    if (sscanf(arg, "--seed={ "
        "0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", "
        "0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 ", 0x%" PRIx32 " }",
        &seed[0], &seed[1], &seed[2], &seed[3], &seed[4], &seed[5],
        &seed[6], &seed[7], &seed[8], &seed[9], &seed[10], &seed[11]) == 12)
        return 0;

    return 1;
}

static void print_seed(const uint32_t *seed)
{
    printf("Random seed: \"--seed={ ");
    for (int i = 0; i < 12; i++) {
        printf("0x%08x%s", seed[i], (i < 11) ? ", " : " }\"\n");
    }
}

#endif