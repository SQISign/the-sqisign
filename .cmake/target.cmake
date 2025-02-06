# SPDX-License-Identifier: Apache-2.0

include(CheckTypeSize)

function(check_target_feature CODE RUN_RESULT)
    set(TEMP_FILE "${CMAKE_BINARY_DIR}/check_target_feature.c")
    file(WRITE
        ${TEMP_FILE}
        "int main(void) {
            ${CODE}
            return 0;
        }")

    try_run(TEMP_RUN_RESULT TEMP_COMPILE_RESULT ${CMAKE_BINARY_DIR} ${TEMP_FILE})

    set(${RUN_RESULT} ${TEMP_RUN_RESULT} PARENT_SCOPE)
    if (ARGC EQUAL 3)
        set(${ARGV2} ${TEMP_COMPILE_RESULT} PARENT_SCOPE)
    endif()

    file(REMOVE ${TEMP_FILE})
endfunction()

if (${CMAKE_SYSTEM_PROCESSOR} MATCHES "aarch64" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm64")
    add_compile_definitions(TARGET_ARM64)
    set(RADIX 64)

    if (NOT APPLE)
        check_target_feature("asm volatile(\"mrs x0, PMCCNTR_EL0\" : : : \"x0\");" CYCCNT)

        if (CYCCNT STREQUAL "FAILED_TO_RUN")
            message(STATUS "Cycle counter not supported, reverting to fallback measurement")
            add_compile_definitions(NO_CYCLE_COUNTER)
        endif()
    endif()
elseif (${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
    add_compile_definitions(TARGET_ARM)
    set(RADIX 32)
elseif (${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    add_compile_definitions(TARGET_AMD64)
    set(RADIX 64)
elseif (${CMAKE_SYSTEM_PROCESSOR} MATCHES "i386" OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES "i686")
    add_compile_definitions(TARGET_X86)
    set(RADIX 32)
elseif (${CMAKE_SYSTEM_PROCESSOR} MATCHES "^(s390x.*|S390X.*)")
    add_compile_definitions(TARGET_S390X TARGET_BIG_ENDIAN)
    set(RADIX 64)
else()
    add_compile_definitions(TARGET_OTHER)
    set(RADIX 64)
    message("Warning: system architecture not detected, defaulting to 64 bit")
endif()

if (NOT GF_RADIX STREQUAL "AUTO")
    if (NOT((GF_RADIX EQUAL 64) OR (GF_RADIX EQUAL 32)))
        message(FATAL_ERROR "Currently supported options for GF_RADIX: 32 or 64. Aborting")
    endif()
    set(RADIX ${GF_RADIX})
endif()

if (NOT DEFINED SQISIGN_BUILD_TYPE)
    set(SQISIGN_BUILD_TYPE "ref")
endif()

if (RADIX EQUAL 32)
    if (${SQISIGN_BUILD_TYPE} MATCHES "broadwell")
        message(FATAL_ERROR "Broadwell implementation not supported in 32-bit build")
    endif()
else()
    # Testing for unsigned 128-bit integer support
    check_type_size("__uint128_t" uint128_t)
    if (${HAVE_uint128_t} AND (uint128_t EQUAL 16))
        add_compile_definitions(HAVE_UINT128)
    elseif(${SQISIGN_BUILD_TYPE} MATCHES "ref")
        message(WARNING "Compiler/platform does not support unsigned 128-bit integers, falling back to 32-bit build")
        set(RADIX 32)
    endif()
endif()

message(STATUS "Using ${RADIX}-bit radix for gf module")

if (RADIX EQUAL 32)
    add_compile_definitions(RADIX_32)
elseif (RADIX EQUAL 64)
    add_compile_definitions(RADIX_64)
endif()

if (UNIX)
    add_compile_definitions(TARGET_OS_UNIX)
else()
    add_compile_definitions(TARGET_OS_OTHER)
endif()

set(C_OPT_FLAGS "")

if (NOT DEFINED SQISIGN_TEST_REPS)
    set(SQISIGN_TEST_REPS 10)
endif()

add_compile_definitions(SQISIGN_TEST_REPS=${SQISIGN_TEST_REPS})
