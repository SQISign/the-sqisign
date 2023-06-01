add_executable(sqisign_test_gf_${SVARIANT_LOWER}_fp test_fp.c test_extras.c)
target_link_libraries(sqisign_test_gf_${SVARIANT_LOWER}_fp ${LIB_GF_${SVARIANT_UPPER}})
target_include_directories(sqisign_test_gf_${SVARIANT_LOWER}_fp PRIVATE ../include ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC})
add_test(sqisign_test_gf_${SVARIANT_LOWER}_fp sqisign_test_gf_${SVARIANT_LOWER}_fp test ${SQISIGN_TEST_REPS})

add_executable(sqisign_test_gf_${SVARIANT_LOWER}_fp2 test_fp2.c test_extras.c)
target_link_libraries(sqisign_test_gf_${SVARIANT_LOWER}_fp2 ${LIB_GF_${SVARIANT_UPPER}})
target_include_directories(sqisign_test_gf_${SVARIANT_LOWER}_fp2 PRIVATE ../include ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC})
add_test(sqisign_test_gf_${SVARIANT_LOWER}_fp2 sqisign_test_gf_${SVARIANT_LOWER}_fp2 test ${SQISIGN_TEST_REPS})