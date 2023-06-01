set(SOURCE_FILES_PROTOCOLS_GENERIC_REF_TESTS
    ${PROTOCOLSX_DIR}/test/commitment.c
    ${PROTOCOLSX_DIR}/test/challenge.c
    ${PROTOCOLSX_DIR}/test/verif.c
    ${PROTOCOLSX_DIR}/test/keygen.c
    ${PROTOCOLSX_DIR}/test/encode.c
    ${PROTOCOLSX_DIR}/test/test_protocols.c
)
add_executable(sqisign_test_protocols_${SVARIANT_LOWER} ${SOURCE_FILES_PROTOCOLS_GENERIC_REF_TESTS})
target_link_libraries(sqisign_test_protocols_${SVARIANT_LOWER} ${LIB_PROTOCOLS_${SVARIANT_UPPER}} ${LIB_ID2ISO_${SVARIANT_UPPER}} ${LIB_EC} ${LIB_KLPT_${SVARIANT_UPPER}} ${LIB_QUATERNION} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_PUBLIC} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}} ${GMP} sqisign_common_sys )
target_include_directories(sqisign_test_protocols_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS} ./include/)

add_test(sqisign_test_protocols_${SVARIANT_LOWER} sqisign_test_protocols_${SVARIANT_LOWER})
