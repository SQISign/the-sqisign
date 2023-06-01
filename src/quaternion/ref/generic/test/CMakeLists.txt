set(SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS
    algebra.c
    ideal.c
    dim4.c
    dim2.c
    integers.c
    lattice.c
    finit.c
    matkermod.c
    randomized.c
    test_quaternions.c
)
add_executable(sqisign_test_quaternion ${SOURCE_FILES_QUATERNION_GENERIC_REF_TESTS})
target_link_libraries(sqisign_test_quaternion ${LIB_INTBIG} ${LIB_QUATERNION} ${GMP} sqisign_common_sys )
target_include_directories(sqisign_test_quaternion PRIVATE ${INC_INTBIG} ${INC_COMMON} ${INC_QUATERNION} ${INC_PUBLIC} )

# MSAN and GMP lead to false positives, see
# https://gmplib.org/list-archives/gmp-bugs/2019-March/004529.html
if(NOT CMAKE_BUILD_TYPE STREQUAL "MSAN")
    add_test(sqisign_test_quaternion sqisign_test_quaternion)
endif()
