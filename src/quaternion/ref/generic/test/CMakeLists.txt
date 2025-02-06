set(SOURCE_FILES_QUATERNION_GENERIC_REF
    intbig.c
    algebra.c
    ideal.c
    dim4.c
    dim2.c
    integers.c
    lattice.c
    lat_ball.c
    finit.c
    mini-gmp.c
    normeq.c
    randomized.c
)

if (NOT GMP_LIBRARY STREQUAL "MINI")
    list(APPEND SOURCE_FILES_QUATERNION_GENERIC_REF ${PROJECT_SOURCE_DIR}/src/mini-gmp/mini-gmp-extra.c)
endif()

set(SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS
    ${SOURCE_FILES_QUATERNION_GENERIC_REF}
    ../hnf/hnf_tests.c
    ../lll/lll_verification.c
    ../lll/lll_tests.c
)
add_executable(sqisign_test_quaternion ${SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS} test_quaternions.c)
target_link_libraries(sqisign_test_quaternion ${LIB_QUATERNION} ${GMP} sqisign_common_sys)
target_include_directories(sqisign_test_quaternion PRIVATE ../internal_quaternion_headers  ${INC_COMMON} ${INC_QUATERNION} ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/mini-gmp)

add_executable(sqisign_bm_quaternion ${SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS} ../lll/lll_benchmarks.c)
target_link_libraries(sqisign_bm_quaternion ${LIB_QUATERNION} ${GMP} sqisign_common_test)
target_include_directories(sqisign_bm_quaternion PRIVATE ../internal_quaternion_headers ${INC_COMMON} ${INC_QUATERNION} ${INC_PUBLIC} ${PROJECT_SOURCE_DIR}/src/mini-gmp)

# MSAN and GMP lead to false positives, see
# https://gmplib.org/list-archives/gmp-bugs/2019-March/004529.html
if(NOT CMAKE_BUILD_TYPE STREQUAL "MSAN")
    add_test(sqisign_test_quaternion sqisign_test_quaternion)
endif()

set(BM_BINS ${BM_BINS} sqisign_bm_quaternion CACHE INTERNAL "List of benchmark executables")
