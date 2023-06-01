#ifndef ID2ISO_TESTS_H
#define ID2ISO_TESTS_H

#include <klpt.h>
#include <quaternion.h>
#include <intbig.h>
#include <stdio.h>
#include <ec.h>
#include <id2iso.h>
#include <rng.h>
#include <quaternion_data.h>
#include <endomorphism_action.h>

/** @internal
 * @ingroup id2iso_id2iso
 * @defgroup id2iso_tests id2iso module test functions
 * @{
 */

int id2iso_test_ker2id();

int id2iso_test_id2ker_even();

int id2iso_test_id2ker_odd();

int id2iso_test_id2iso();

/** @}
 */


#endif
