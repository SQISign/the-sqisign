set(SOURCE_FILES_GF_${SVARIANT_UPPER}_REF
    ${SOURCE_FILES_GF_SPECIFIC}
    ${LVLX_DIR}/fp.c
    ${LVLX_DIR}/fp2.c
)

add_library(${LIB_GF_${SVARIANT_UPPER}} STATIC ${SOURCE_FILES_GF_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${INC_COMMON} ${PROJECT_SOURCE_DIR}/src/precomp/ref/${SVARIANT_LOWER}/include ${INC_GF} ${INC_MP} ${INC_PUBLIC})
target_compile_options(${LIB_GF_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})
target_link_libraries(${LIB_GF_${SVARIANT_UPPER}} ${LIB_MP})
target_compile_definitions(${LIB_GF_${SVARIANT_UPPER}} PUBLIC SQISIGN_VARIANT=${SVARIANT_LOWER})

add_subdirectory(test)
