
set(SOURCE_FILES_GF_${SVARIANT_UPPER}_BROADWELL
    fp_asm.S fp.c fp2.c
)

add_library(${LIB_GF_${SVARIANT_UPPER}} ${SOURCE_FILES_GF_${SVARIANT_UPPER}_BROADWELL})
target_include_directories(${LIB_GF_${SVARIANT_UPPER}} PRIVATE common ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} include ${PROJECT_SOURCE_DIR}/include ${INC_COMMON})
target_compile_options(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
