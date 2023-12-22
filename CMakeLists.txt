# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.5)
project(SQIsign VERSION 1.0 LANGUAGES C ASM)

set(SQISIGN_SO_VERSION "0")
set(CMAKE_C_STANDARD 99)

include(CTest)

option(ENABLE_STRICT "Build with strict compile options." OFF)
option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
option(ENABLE_GMP_BUILD "Download and build external version of GMP" OFF)
option(ENABLE_MINI_GMP "Build with mini-gmp library" OFF)
option(ENABLE_DOC_TARGET "Enable building API documentation using doxygen" OFF)

if (NOT DEFINED SQISIGN_BUILD_TYPE)
	SET(SQISIGN_BUILD_TYPE "ref")
endif()

if(SQISIGN_BUILD_TYPE STREQUAL "broadwell")
	SET(SVARIANT_S "lvl1")
else()
	SET(SVARIANT_S "lvl1;lvl3;lvl5")
endif()

include(.cmake/flags.cmake)
include(.cmake/sanitizers.cmake)
include(.cmake/target.cmake)
if(ENABLE_DOC_TARGET)
	include(.cmake/target_docs.cmake)
endif()
include(.cmake/gmpconfig.cmake)


set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)
set(SELECT_SQISIGN_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/sqisign_variant.cmake)

set(INC_PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(apps)

if(ENABLE_TESTS)
	enable_testing()
	add_subdirectory(test)
endif()
