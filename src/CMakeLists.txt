add_subdirectory(common)

if(ENABLE_SIGN)
    add_subdirectory(quaternion)
endif()

add_subdirectory(mp)
add_subdirectory(gf)
add_subdirectory(precomp)
add_subdirectory(ec)
add_subdirectory(hd)
add_subdirectory(verification)

if(ENABLE_SIGN)
    add_subdirectory(id2iso)
    add_subdirectory(signature)
endif()

FOREACH(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
    set(SOURCE_FILES_VARIANT sqisign.c)

    # Library for SQIsign variant
    add_library(sqisign_${SVARIANT_LOWER} ${SOURCE_FILES_VARIANT})
    target_link_libraries(sqisign_${SVARIANT_LOWER} PUBLIC
        $<$<BOOL:${ENABLE_SIGN}>:${LIB_SIGNATURE_${SVARIANT_UPPER}}>
        ${LIB_VERIFICATION_${SVARIANT_UPPER}}
        $<$<BOOL:${ENABLE_SIGN}>:${LIB_ID2ISO_${SVARIANT_UPPER}}>
        $<$<BOOL:${ENABLE_SIGN}>:${LIB_QUATERNION}>
        ${LIB_MP}
        ${LIB_GF_${SVARIANT_UPPER}}
        ${LIB_EC_${SVARIANT_UPPER}}
        ${LIB_HD_${SVARIANT_UPPER}}
        ${LIB_PRECOMP_${SVARIANT_UPPER}}
        $<$<BOOL:${ENABLE_SIGN}>:GMP>
        sqisign_common_sys
    )

    target_include_directories(sqisign_${SVARIANT_LOWER} PUBLIC $<$<BOOL:${ENABLE_SIGN}>:${INC_SIGNATURE}> ${INC_VERIFICATION} $<$<BOOL:${ENABLE_SIGN}>:${INC_QUATERNION}> ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_MP} ${INC_EC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_HD} $<$<BOOL:${ENABLE_SIGN}>:${INC_ID2ISO}> ../include PRIVATE common/generic internal)
    target_compile_definitions(sqisign_${SVARIANT_LOWER} PUBLIC SQISIGN_VARIANT=${SVARIANT})

    # Library for SQIsign variant (test)
    add_library(sqisign_${SVARIANT_LOWER}_test ${SOURCE_FILES_VARIANT})
    target_link_libraries(sqisign_${SVARIANT_LOWER}_test PUBLIC
        $<$<BOOL:${ENABLE_SIGN}>:${LIB_SIGNATURE_${SVARIANT_UPPER}}>
        ${LIB_VERIFICATION_${SVARIANT_UPPER}}
        $<$<BOOL:${ENABLE_SIGN}>:${LIB_ID2ISO_${SVARIANT_UPPER}}>
        $<$<BOOL:${ENABLE_SIGN}>:${LIB_QUATERNION}>
        ${LIB_MP}
        ${LIB_GF_${SVARIANT_UPPER}}
        ${LIB_EC_${SVARIANT_UPPER}}
        ${LIB_HD_${SVARIANT_UPPER}}
        ${LIB_PRECOMP_${SVARIANT_UPPER}}
        $<$<BOOL:${ENABLE_SIGN}>:GMP>
        sqisign_common_test
    )

    target_include_directories(sqisign_${SVARIANT_LOWER}_test PUBLIC $<$<BOOL:${ENABLE_SIGN}>:${INC_SIGNATURE}> ${INC_VERIFICATION} $<$<BOOL:${ENABLE_SIGN}>:${INC_QUATERNION}> ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_MP} ${INC_EC} ${INC_GF} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_HD} $<$<BOOL:${ENABLE_SIGN}>:${INC_ID2ISO}> ../include PRIVATE common/generic internal)
    target_compile_definitions(sqisign_${SVARIANT_LOWER}_test PUBLIC SQISIGN_VARIANT=${SVARIANT})

    # Library with NIST API
    set(SOURCE_FILE_NISTAPI nistapi/${SVARIANT_LOWER}/api.c)
    add_library(sqisign_${SVARIANT_LOWER}_nistapi ${SOURCE_FILE_NISTAPI})
    target_link_libraries(sqisign_${SVARIANT_LOWER}_nistapi PUBLIC sqisign_${SVARIANT_LOWER})
    target_include_directories(sqisign_${SVARIANT_LOWER}_nistapi PUBLIC nistapi/${SVARIANT_LOWER} PUBLIC ../include)
    target_compile_definitions(sqisign_${SVARIANT_LOWER}_nistapi PUBLIC SQISIGN_VARIANT=${SVARIANT})

    # Library with NIST API (test)
    add_library(sqisign_${SVARIANT_LOWER}_test_nistapi ${SOURCE_FILE_NISTAPI})
    target_link_libraries(sqisign_${SVARIANT_LOWER}_test_nistapi PUBLIC sqisign_${SVARIANT_LOWER}_test)
    target_include_directories(sqisign_${SVARIANT_LOWER}_test_nistapi PUBLIC nistapi/${SVARIANT_LOWER})
    target_compile_definitions(sqisign_${SVARIANT_LOWER}_test_nistapi PUBLIC SQISIGN_VARIANT=${SVARIANT})
ENDFOREACH()
