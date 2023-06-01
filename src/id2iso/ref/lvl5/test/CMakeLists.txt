set(SOURCE_FILES_ID2ISO_GENERIC_REF_TESTS
    ${ID2ISOX_DIR}/test/ker2id.c
    ${ID2ISOX_DIR}/test/id2ker_even.c
    ${ID2ISOX_DIR}/test/id2ker_odd.c
    ${ID2ISOX_DIR}/test/id2iso.c
    ${ID2ISOX_DIR}/test/test_id2iso.c
)
add_executable(sqisign_test_id2iso_${SVARIANT_LOWER} ${SOURCE_FILES_ID2ISO_GENERIC_REF_TESTS})
target_link_libraries(sqisign_test_id2iso_${SVARIANT_LOWER} ${LIB_ID2ISO_${SVARIANT_UPPER}} ${LIB_KLPT_${SVARIANT_UPPER}} ${LIB_QUATERNION} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_PUBLIC} ${GMP} sqisign_common_sys )
target_include_directories(sqisign_test_id2iso_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_INTBIG} ${INC_EC} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_ID2ISO} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ./include/ )

add_test(sqisign_test_id2iso_${SVARIANT_LOWER} sqisign_test_id2iso_${SVARIANT_LOWER})
