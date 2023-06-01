# There are the following dependencies
#     ┌─┬──────┬─┐           ┌─┬────┬─┐            ┌─┬──────┬─┐
#     │ ├──────┤ │           │ ├────┤ │            │ ├──────┤ │
#     │ │Keygen│ │           │ │Sign│ │            │ │Verify│ │
#     │ ├──────┤ │           │ ├────┤ │            │ ├──────┤ │
#     └─┴───┬──┴─┘           └─┴─┬──┴─┘            └─┴───┬──┴─┘
#           │                    │                       │
#           │                    │                       │
#           ├────────────────────┼─────────────────┐     │
#           │                    │                 │     │
#           │                    │                 │     │
#       ┌───▼──┐          ┌──────▼────────┐   ┌────▼─────▼───────────┐
#       │ PRNG ◄────┬─────┤ Iso <-> Ideal ├───►   Elliptic Curves,   │
#       └───▲──┘    │     └──────┬────────┘   │ Pairings & Isogenies │
#           │       │            │            └───▲──────┬───────────┘
#           │       │            │                │      │
#       ┌───┴──┐    │            │                │      │
#       │ KLPT ◄────┘            │     ┌──────────┘      │
#       └───┬──┘                 │     │                 │
#           │                    │     │                 │
# ┌─────────▼─────────┐          │     │                 │
# │ Quaternion orders │          │     │            ┌────▼───┐
# │     and ideals    │          │     │            │ GF(p²) │
# └─────────┬─────────┘          │     │            └────┬───┘
#           │           ┌─┬──────▼─────┴──┬─┐            │
#     ┌─────▼─────┐     │ ├───────────────┤ │      ┌─────▼─────┐
#     │ MP BigInt │     │ │Precomputations│ │      │ FP BigInt │
#     └───────────┘     │ ├───────────────┤ │      └───────────┘
#                       └─┴───────────────┴─┘                    

add_subdirectory(common)
add_subdirectory(intbig)
add_subdirectory(quaternion)
add_subdirectory(precomp)
add_subdirectory(klpt)
add_subdirectory(gf)
add_subdirectory(ec)

add_subdirectory(id2iso)
add_subdirectory(protocols)

FOREACH(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    string(TOUPPER ${SVARIANT} SVARIANT_UPPER)
    set(SOURCE_FILES_VARIANT sqisign.c)
    # Library for SQIsign variant
    add_library(sqisign_${SVARIANT_LOWER} ${SOURCE_FILES_VARIANT})
    target_link_libraries(sqisign_${SVARIANT_LOWER} PUBLIC 
        ${LIB_PROTOCOLS_${SVARIANT_UPPER}} 
        ${LIB_ID2ISO_${SVARIANT_UPPER}} 
        ${LIB_KLPT_${SVARIANT_UPPER}} 
        ${LIB_QUATERNION} 
        ${LIB_PRECOMP_${SVARIANT_UPPER}} 
        ${LIB_INTBIG} 
        ${LIB_GF_${SVARIANT_UPPER}} 
        ${LIB_EC_${SVARIANT_UPPER}} 
        ${GMP} 
        sqisign_common_sys
    )

    target_include_directories(sqisign_${SVARIANT_LOWER} PUBLIC ${INC_PROTOCOLS} ${INC_INTBIG} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_EC} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_KLPT} ${INC_ID2ISO} ../include PRIVATE common/generic internal)
    target_compile_definitions(sqisign_${SVARIANT_LOWER} PUBLIC SQISIGN_VARIANT=${SVARIANT})

    # Library for SQIsign variant (test)
    add_library(sqisign_${SVARIANT_LOWER}_test ${SOURCE_FILES_VARIANT})
    target_link_libraries(sqisign_${SVARIANT_LOWER}_test PUBLIC 
        ${LIB_PROTOCOLS_${SVARIANT_UPPER}} 
        ${LIB_ID2ISO_${SVARIANT_UPPER}} 
        ${LIB_KLPT_${SVARIANT_UPPER}} 
        ${LIB_QUATERNION} 
        ${LIB_PRECOMP_${SVARIANT_UPPER}} 
        ${LIB_INTBIG} 
        ${LIB_GF_${SVARIANT_UPPER}} 
        ${LIB_EC_${SVARIANT_UPPER}} 
        ${GMP} 
        sqisign_common_test
    )

    target_include_directories(sqisign_${SVARIANT_LOWER}_test PUBLIC ${INC_PROTOCOLS} ${INC_INTBIG} ${INC_QUATERNION} ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_EC} ${INC_GF_${SVARIANT_UPPER}} ${INC_COMMON} ${INC_KLPT} ${INC_ID2ISO} ../include PRIVATE common/generic internal)
    target_compile_definitions(sqisign_${SVARIANT_LOWER}_test PUBLIC SQISIGN_VARIANT=${SVARIANT})

    # Library with NIST API
    set(SOURCE_FILE_NISTAPI nistapi/${SVARIANT_LOWER}/api.c)
    add_library(sqisign_${SVARIANT_LOWER}_nistapi ${SOURCE_FILE_NISTAPI})
    target_link_libraries(sqisign_${SVARIANT_LOWER}_nistapi PRIVATE sqisign_${SVARIANT_LOWER})
    target_include_directories(sqisign_${SVARIANT_LOWER}_nistapi PUBLIC nistapi/${SVARIANT_LOWER} PUBLIC ../include)

    # Library with NIST API (test)
    add_library(sqisign_${SVARIANT_LOWER}_test_nistapi ${SOURCE_FILE_NISTAPI})
    target_link_libraries(sqisign_${SVARIANT_LOWER}_test_nistapi PRIVATE sqisign_${SVARIANT_LOWER}_test)
    target_include_directories(sqisign_${SVARIANT_LOWER}_test_nistapi PUBLIC nistapi/${SVARIANT_LOWER})
ENDFOREACH()
