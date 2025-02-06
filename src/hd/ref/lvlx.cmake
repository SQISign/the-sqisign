set(SOURCE_FILES_HD_GENERIC_REF
    ${LVLX_DIR}/hd.c
    ${LVLX_DIR}/theta_structure.c 
    ${LVLX_DIR}/theta_isogenies.c 
)

add_library(${LIB_HD_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_HD_GENERIC_REF})
target_link_libraries(${LIB_HD_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})
target_include_directories(${LIB_HD_${SVARIANT_UPPER}} PRIVATE ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_HD})
target_compile_options(${LIB_HD_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_compile_definitions(${LIB_HD_${SVARIANT_UPPER}} PUBLIC SQISIGN_VARIANT=${SVARIANT_LOWER})
