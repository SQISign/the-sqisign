add_executable(curve-arith.test_${SVARIANT_LOWER} ${LVLX_DIR}/test/curve-arith-test.c ${LVLX_DIR}/test/test_extras.c)
target_include_directories(curve-arith.test_${SVARIANT_LOWER} PUBLIC ${INC_COMMON} ${INC_MP} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ../include ${INC_EC} .)
target_link_libraries(curve-arith.test_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}} sqisign_common_test)

add_executable(biextension.test_${SVARIANT_LOWER} ${LVLX_DIR}/test/biextension-test.c)
target_include_directories(biextension.test_${SVARIANT_LOWER} PUBLIC ${INC_COMMON} ${INC_MP} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ../include ${INC_EC} .)
target_link_libraries(biextension.test_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}} sqisign_common_test)

add_executable(basis-gen.test_${SVARIANT_LOWER} ${LVLX_DIR}/test/basis-gen-test.c)
target_include_directories(basis-gen.test_${SVARIANT_LOWER} PUBLIC ${INC_COMMON} ${INC_MP} ${LVLX_DIR}/test ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ../include ${INC_EC} .)
target_link_libraries(basis-gen.test_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}})

add_test(curve_arith.test_${SVARIANT_LOWER} curve-arith.test_${SVARIANT_LOWER})
add_test(ec_biextension.test_${SVARIANT_LOWER} biextension.test_${SVARIANT_LOWER})
add_test(ec_basis_gen.test_${SVARIANT_LOWER} basis-gen.test_${SVARIANT_LOWER})

add_executable(curve-arith.bench_${SVARIANT_LOWER} ${LVLX_DIR}/test/curve-arith-bench.c ${LVLX_DIR}/test/test_extras.c)
target_include_directories(curve-arith.bench_${SVARIANT_LOWER} PUBLIC ${INC_COMMON} ${INC_MP} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ../include ${INC_EC} .)
target_link_libraries(curve-arith.bench_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}} sqisign_common_sys)

add_executable(biextension.bench_${SVARIANT_LOWER} ${LVLX_DIR}/test/biextension-bench.c)
target_include_directories(biextension.bench_${SVARIANT_LOWER} PUBLIC ${INC_COMMON} ${INC_MP} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ../include ${INC_EC} .)
target_link_libraries(biextension.bench_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}} sqisign_common_sys)

add_executable(basis-gen.bench_${SVARIANT_LOWER} ${LVLX_DIR}/test/basis-gen-bench.c)
target_include_directories(basis-gen.bench_${SVARIANT_LOWER} PUBLIC ${INC_COMMON} ${INC_MP} ${LVLX_DIR}/test ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ../include ${INC_EC} .)
target_link_libraries(basis-gen.bench_${SVARIANT_LOWER} ${LIB_EC_${SVARIANT_UPPER}})

set(BM_BINS ${BM_BINS}
    curve-arith.bench_${SVARIANT_LOWER} basis-gen.bench_${SVARIANT_LOWER} biextension.bench_${SVARIANT_LOWER}
    CACHE INTERNAL "List of benchmark executables")

