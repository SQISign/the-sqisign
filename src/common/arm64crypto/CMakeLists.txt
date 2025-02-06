if(CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(SOURCE_FILES_COMMON_ARM64CRYPTO randombytes_ctrdrbg_inline_asm.c)
else()
    set(SOURCE_FILES_COMMON_ARM64CRYPTO randombytes_ctrdrbg.c)
    set_source_files_properties(randombytes_ctrdrbg.c PROPERTIES COMPILE_FLAGS -fno-strict-aliasing)
endif()

foreach(SQISIGN_COMMON_TARGET sqisign_common_test sqisign_common_sys)
    target_sources(${SQISIGN_COMMON_TARGET} PRIVATE ${SOURCE_FILES_COMMON_ARM64CRYPTO})
    target_include_directories(${SQISIGN_COMMON_TARGET} PRIVATE include)
    target_compile_definitions(${SQISIGN_COMMON_TARGET} PRIVATE RANDOMBYTES_ARM64CRYPTO)
    target_compile_options(${SQISIGN_COMMON_TARGET} PRIVATE -march=armv8-a+crypto)
endforeach()

set(SOURCE_FILES_CTRDRBG_TEST_BENCHMARK
    ${SOURCE_FILES_COMMON_ARM64CRYPTO}
    ../ref/aes_c.c
    ../ref/randombytes_ctrdrbg.c
    ../generic/randombytes_system.c
)

add_executable(sqisign_test_ctrdrbg_arm64crypto ${SOURCE_FILES_CTRDRBG_TEST_BENCHMARK} ../generic/test/test_ctrdrbg.c)
target_include_directories(sqisign_test_ctrdrbg_arm64crypto PRIVATE ${INC_PUBLIC} ${INC_COMMON} include ../ref/include)
target_compile_definitions(sqisign_test_ctrdrbg_arm64crypto PRIVATE
    CTRDRBG_TEST_BENCH
    RANDOMBYTES_INIT_PLATFORM=randombytes_init_arm64crypto
    RANDOMBYTES_PLATFORM=randombytes_arm64crypto)
target_compile_options(sqisign_test_ctrdrbg_arm64crypto PRIVATE -march=armv8-a+crypto)

add_test(sqisign_test_ctrdrbg_arm64crypto sqisign_test_ctrdrbg_arm64crypto)

add_executable(sqisign_bench_ctrdrbg_arm64crypto ${SOURCE_FILES_CTRDRBG_TEST_BENCHMARK} ../generic/test/bench_ctrdrbg.c)
target_include_directories(sqisign_bench_ctrdrbg_arm64crypto PRIVATE ${INC_PUBLIC} ${INC_COMMON} include ../ref/include)
target_compile_definitions(sqisign_bench_ctrdrbg_arm64crypto PRIVATE
    CTRDRBG_TEST_BENCH
    RANDOMBYTES_INIT_PLATFORM=randombytes_init_arm64crypto
    RANDOMBYTES_PLATFORM=randombytes_arm64crypto)
target_compile_options(sqisign_bench_ctrdrbg_arm64crypto PRIVATE -march=armv8-a+crypto)

set(BM_BINS ${BM_BINS} sqisign_bench_ctrdrbg_arm64crypto CACHE INTERNAL "List of benchmark executables")
