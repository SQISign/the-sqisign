add_executable(sqisign_test_signature_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_signature.c)
target_link_libraries(sqisign_test_signature_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_test)
target_include_directories(sqisign_test_signature_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_executable(sqisign_test_threadsafety_${SVARIANT_LOWER} ${LVLX_DIR}/test/test_threadsafety.c)
target_link_libraries(sqisign_test_threadsafety_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_test pthread)
target_include_directories(sqisign_test_threadsafety_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

add_test(sqisign_test_signature_${SVARIANT_LOWER} sqisign_test_signature_${SVARIANT_LOWER} 3)
add_test(sqisign_test_threadsafety_${SVARIANT_LOWER} sqisign_test_threadsafety_${SVARIANT_LOWER} 3)

add_custom_command(
  TARGET sqisign_test_signature_${SVARIANT_LOWER}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy $<TARGET_FILE:sqisign_test_signature_${SVARIANT_LOWER}> ${CMAKE_BINARY_DIR}/test/sqisign_${SVARIANT_LOWER}
)

add_executable(sqisign_bench_signature_${SVARIANT_LOWER} ${LVLX_DIR}/test/bench_signature.c)
target_link_libraries(sqisign_bench_signature_${SVARIANT_LOWER} ${LIB_SIGNATURE_${SVARIANT_UPPER}} ${LIB_VERIFICATION_${SVARIANT_UPPER}} sqisign_common_sys)
target_include_directories(sqisign_bench_signature_${SVARIANT_LOWER} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_VERIFICATION} ${INC_SIGNATURE})

set(BM_BINS ${BM_BINS} sqisign_bench_signature_${SVARIANT_LOWER} CACHE INTERNAL "List of benchmark executables")
