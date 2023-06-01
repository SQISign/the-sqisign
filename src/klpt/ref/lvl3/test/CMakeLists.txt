set(SOURCE_FILES_KLPT_GENERIC_REF_TESTS
    ${KLPTX_DIR}/test/tools.c
    ${KLPTX_DIR}/test/equiv.c
    ${KLPTX_DIR}/test/klpt.c
    ${KLPTX_DIR}/test/eichler.c
    ${KLPTX_DIR}/test/test_klpt.c
)
add_executable(sqisign_test_klpt_${SVARIANT_UPPER} ${SOURCE_FILES_KLPT_GENERIC_REF_TESTS})
target_link_libraries(sqisign_test_klpt_${SVARIANT_UPPER} ${LIB_KLPT_${SVARIANT_UPPER}} ${LIB_QUATERNION} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_PUBLIC} ${GMP} sqisign_common_sys )
target_include_directories(sqisign_test_klpt_${SVARIANT_UPPER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ./include/ )

add_test(sqisign_test_klpt_${SVARIANT_UPPER} sqisign_test_klpt_${SVARIANT_UPPER})