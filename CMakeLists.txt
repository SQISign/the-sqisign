# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)
project(SQIsign VERSION 2.0 LANGUAGES C ASM)

set(SQISIGN_SO_VERSION "0")
set(CMAKE_C_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)

option(ENABLE_STRICT "Build with strict compile options." ON)
option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
option(ENABLE_SIGN "Build with sign functionality" ON)
set(GMP_LIBRARY "SYSTEM" CACHE STRING "Which version of GMP to use: SYSTEM, BUILD or MINI")
set(GF_RADIX "AUTO" CACHE STRING "Set the radix for the gf module (currently supported values: 32 or 64), or AUTO.")

set(GLOBAL_PUBLIC_DEFS "")

if (NOT DEFINED SQISIGN_BUILD_TYPE)
  SET(SQISIGN_BUILD_TYPE "ref")
endif()

if (${SQISIGN_BUILD_TYPE} MATCHES "ref")
  list(APPEND GLOBAL_PUBLIC_DEFS "SQISIGN_BUILD_TYPE_REF" "SQISIGN_GF_IMPL_REF")
elseif (${SQISIGN_BUILD_TYPE} MATCHES "opt")
  list(APPEND GLOBAL_PUBLIC_DEFS "SQISIGN_BUILD_TYPE_OPT" "SQISIGN_GF_IMPL_REF")
elseif (${SQISIGN_BUILD_TYPE} MATCHES "broadwell")
  list(APPEND GLOBAL_PUBLIC_DEFS "SQISIGN_BUILD_TYPE_BROADWELL" "SQISIGN_GF_IMPL_BROADWELL")
elseif (${SQISIGN_BUILD_TYPE} MATCHES "arm64crypto")
  list(APPEND GLOBAL_PUBLIC_DEFS "SQISIGN_BUILD_TYPE_ARM64CRYPTO" "SQISIGN_GF_IMPL_REF")
endif()

add_compile_definitions(${GLOBAL_PUBLIC_DEFS})

SET(SVARIANT_S "lvl1;lvl3;lvl5")

include(.cmake/flags.cmake)
include(.cmake/sanitizers.cmake)
include(.cmake/target.cmake)
if(ENABLE_SIGN)
	include(.cmake/gmpconfig.cmake)
	add_compile_definitions(ENABLE_SIGN)
	list(APPEND GLOBAL_PUBLIC_DEFS "ENABLE_SIGN")
endif()

set(BM_BINS "" CACHE INTERNAL "List of benchmark binaries")

set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)
set(SELECT_SQISIGN_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/sqisign_variant.cmake)

set(INC_PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(test)

include(.cmake/bm.cmake)

#if(ENABLE_TESTS)
#    enable_testing()
#    add_subdirectory(test)
#endif()

include(CMakePackageConfigHelpers)

set(CMAKE_INSTALL_BINDIR "bin" CACHE STRING "Install directory for binaries")
set(CMAKE_INSTALL_LIBDIR "lib" CACHE STRING "Install directory for libraries")
set(CMAKE_INSTALL_INCLUDEDIR "include" CACHE STRING "Install directory for headers")
set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/SQIsign" CACHE STRING "Install directory for CMake configs")

set(_INSTALL_TARGETS "")
set(_SQISIGN_VARIANTS "lvl1" "lvl3" "lvl5")

foreach(_VAR ${_SQISIGN_VARIANTS})
    if(TARGET sqisign_${_VAR})
        list(APPEND _INSTALL_TARGETS "sqisign_${_VAR}")
    endif()
    if(TARGET sqisign_${_VAR}_nistapi)
        list(APPEND _INSTALL_TARGETS "sqisign_${_VAR}_nistapi")
    endif()
endforeach()

macro(get_all_targets _result _dir)
    get_property(_subdirs DIRECTORY "${_dir}" PROPERTY SUBDIRECTORIES)
    foreach(_subdir ${_subdirs})
        get_all_targets(${_result} "${_subdir}")
    endforeach()
    get_property(_sub_targets DIRECTORY "${_dir}" PROPERTY BUILDSYSTEM_TARGETS)
    foreach(_tgt ${_sub_targets})
        get_target_property(_tgt_type ${_tgt} TYPE)
        if(_tgt_type MATCHES "STATIC_LIBRARY|SHARED_LIBRARY")
            # Filter out test-only internal modules
            if(NOT _tgt MATCHES "test|benchmark")
                list(APPEND ${_result} ${_tgt})
            endif()
        endif()
    endforeach()
endmacro()

get_all_targets(_INTERNAL_LIBS "${PROJECT_SOURCE_DIR}/src")
if(_INTERNAL_LIBS)
    list(APPEND _INSTALL_TARGETS ${_INTERNAL_LIBS})
    list(REMOVE_DUPLICATES _INSTALL_TARGETS)
endif()

install(DIRECTORY "${INC_PUBLIC}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.h"
)

if(_INSTALL_TARGETS)
    foreach(_tgt ${_INSTALL_TARGETS})
        if(TARGET ${_tgt})
            target_compile_definitions(${_tgt} PUBLIC ${GLOBAL_PUBLIC_DEFS})
        endif()
    endforeach()

    install(TARGETS ${_INSTALL_TARGETS}
        EXPORT SQIsignTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    )

    install(EXPORT SQIsignTargets
        FILE SQIsignTargets.cmake
        NAMESPACE SQIsign::
        DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
    )

    configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/.cmake/SQIsignConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/SQIsignConfig.cmake"
        INSTALL_DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/SQIsignConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/SQIsignConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/SQIsignConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_CMAKEDIR}"
    )
endif()
