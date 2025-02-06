set(SOURCE_FILES_ID2ISO_GENERIC_REF_TESTS
    ${LVLX_DIR}/test/ker2id.c
    ${LVLX_DIR}/test/test_id2iso.c
    ${LVLX_DIR}/test/test_dim2id2iso.c
)
add_executable(sqisign_test_id2iso_${SVARIANT_LOWER} ${SOURCE_FILES_ID2ISO_GENERIC_REF_TESTS})
target_link_libraries(sqisign_test_id2iso_${SVARIANT_LOWER} ${LIB_ID2ISO_${SVARIANT_UPPER}} sqisign_common_test)
target_include_directories(sqisign_test_id2iso_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_EC} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_HD} ${INC_ID2ISO} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ./include/ )

add_test(sqisign_test_id2iso_${SVARIANT_LOWER} sqisign_test_id2iso_${SVARIANT_LOWER})

set(SOURCE_FILES_ID2ISO_RI_BENCH
    ${LVLX_DIR}/test/represent_integer_benchmarks.c
)
add_executable(sqisign_id2iso_benchmark_represent_integer_${SVARIANT_LOWER} ${SOURCE_FILES_ID2ISO_RI_BENCH})
target_link_libraries(sqisign_id2iso_benchmark_represent_integer_${SVARIANT_LOWER} ${LIB_ID2ISO_${SVARIANT_UPPER}} sqisign_common_test)
target_include_directories(sqisign_id2iso_benchmark_represent_integer_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_EC} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_HD} ${INC_ID2ISO} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ./include/ )
set(BM_BINS ${BM_BINS} sqisign_id2iso_benchmark_represent_integer_${SVARIANT_LOWER} CACHE INTERNAL "List of benchmark executables")
