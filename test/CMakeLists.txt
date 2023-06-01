# SPDX-License-Identifier: Apache-2.0

FOREACH(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
    add_executable(sqisign_test_kat_${SVARIANT} test_kat.c)
    target_link_libraries(sqisign_test_kat_${SVARIANT} sqisign_${SVARIANT_LOWER}_test)
    target_include_directories(sqisign_test_kat_${SVARIANT} PRIVATE ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER} ${INC_PUBLIC} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON} ${INC_ID2ISO} ${INC_PROTOCOLS})

    add_executable(sqisign_bench_${SVARIANT} bench.c)
    target_link_libraries(sqisign_bench_${SVARIANT} sqisign_${SVARIANT_LOWER})
    target_include_directories(sqisign_bench_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})

    add_executable(sqisign_test_scheme_${SVARIANT} test_sqisign.c)
    target_link_libraries(sqisign_test_scheme_${SVARIANT} sqisign_${SVARIANT_LOWER})
    target_include_directories(sqisign_test_scheme_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})

    add_executable(sqisign_test_prof_${SVARIANT} test_sqisign_prof.c)
    target_link_libraries(sqisign_test_prof_${SVARIANT} sqisign_${SVARIANT_LOWER})
    target_include_directories(sqisign_test_prof_${SVARIANT} PUBLIC ${PROJECT_SOURCE_DIR}/src/common ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/nistapi/${SVARIANT_LOWER})

    add_test(sqisign_${SVARIANT}_KAT sqisign_test_kat_${SVARIANT})
    add_test(sqisign_${SVARIANT}_SELFTEST sqisign_test_scheme_${SVARIANT})
ENDFOREACH()

