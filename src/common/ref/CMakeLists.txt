set(SOURCE_FILES_COMMON_TEST_REF
    randombytes_ctrdrbg.c 
    aes_c.c 
)

target_sources(sqisign_common_test PRIVATE ${SOURCE_FILES_COMMON_TEST_REF})
target_include_directories(sqisign_common_test PRIVATE include)
target_compile_definitions(sqisign_common_test PRIVATE RANDOMBYTES_C)

target_compile_definitions(sqisign_common_sys PRIVATE RANDOMBYTES_SYSTEM)
