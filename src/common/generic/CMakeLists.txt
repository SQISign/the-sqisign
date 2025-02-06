set(SOURCE_FILES_COMMON_GENERIC
    randombytes_system.c 
    fips202.c 
    mem.c
    tools.c
)

foreach (SQISIGN_COMMON_TARGET sqisign_common_test sqisign_common_sys)
    add_library(${SQISIGN_COMMON_TARGET} STATIC ${SOURCE_FILES_COMMON_GENERIC})
    target_include_directories(${SQISIGN_COMMON_TARGET} PRIVATE include ${INC_PUBLIC})
    target_compile_options(${SQISIGN_COMMON_TARGET} PUBLIC ${C_OPT_FLAGS})
    if (ENABLE_CT_TESTING)
        target_compile_definitions(${SQISIGN_COMMON_TARGET} PUBLIC ENABLE_CT_TESTING)
    endif()
endforeach()
