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
    printer.c
    lll/rationals.c
    lll/l2.c
    lll/lll_verification.c
    lll/lll_applications.c
    lll/rationals.c
    normeq.c
    hnf/ibz_division.c
    hnf/hnf_internal.c
    hnf/hnf.c
    test/random_input_generation.c
)

add_library(${LIB_QUATERNION} STATIC ${SOURCE_FILES_QUATERNION_GENERIC_REF})
target_link_libraries(${LIB_QUATERNION} GMP m)
target_include_directories(${LIB_QUATERNION} PRIVATE common ${INC_PUBLIC} ${INC_COMMON} ${INC_QUATERNION} internal_quaternion_headers)
target_compile_options(${LIB_QUATERNION} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
