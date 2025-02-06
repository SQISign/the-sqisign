set(SOURCE_FILES_EC_${SVARIANT_UPPER}_REF
    ${LVLX_DIR}/ec.c
    ${LVLX_DIR}/ec_jac.c
    ${LVLX_DIR}/xisog.c
    ${LVLX_DIR}/xeval.c
    ${LVLX_DIR}/isog_chains.c
    ${LVLX_DIR}/basis.c
    ${LVLX_DIR}/biextension.c
)

add_library(${LIB_EC_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_EC_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_EC_${SVARIANT_UPPER}} PRIVATE ${INC_COMMON} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_PUBLIC} ${INC_MP} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_EC})
target_compile_options(${LIB_EC_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_link_libraries(${LIB_EC_${SVARIANT_UPPER}} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_MP} ${LIB_GF_${SVARIANT_UPPER}})
target_compile_definitions(${LIB_EC_${SVARIANT_UPPER}} PUBLIC SQISIGN_VARIANT=${SVARIANT_LOWER})

add_subdirectory(test)
