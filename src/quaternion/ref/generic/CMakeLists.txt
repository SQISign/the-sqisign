set(SOURCE_FILES_QUATERNION_GENERIC_REF
    algebra.c
    ideal.c
    dim4.c
    dim2.c
    integers.c
    lattice.c
    finit.c
    printer.c
    lll.c
    matkermod.c
)

add_library(${LIB_QUATERNION} ${SOURCE_FILES_QUATERNION_GENERIC_REF})
target_include_directories(${LIB_QUATERNION} PRIVATE common ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_QUATERNION})
target_compile_options(${LIB_QUATERNION} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
