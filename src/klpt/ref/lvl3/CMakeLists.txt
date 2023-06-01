set(SOURCE_FILES_KLPT_GENERIC_REF
    ${KLPTX_DIR}/eichler.c
    ${KLPTX_DIR}/equiv.c 
    ${KLPTX_DIR}/klpt.c 
    ${KLPTX_DIR}/tools.c
)

add_library(${LIB_KLPT_${SVARIANT_UPPER}} ${SOURCE_FILES_KLPT_GENERIC_REF})
target_include_directories(${LIB_KLPT_${SVARIANT_UPPER}} PRIVATE common ${INC_PUBLIC} ${INC_COMMON} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_QUATERNION} ${INC_KLPT} )
target_compile_options(${LIB_KLPT_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_subdirectory(test)
