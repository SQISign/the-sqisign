if (GMP_LIBRARY STREQUAL "SYSTEM")
  # use system gmp version
  message(STATUS "Using system GMP")

  find_library(GMP gmp)
  find_path(GMP_INCLUDE gmp.h)

  add_library(GMP UNKNOWN IMPORTED)
  set_target_properties(GMP PROPERTIES
    IMPORTED_LOCATION ${GMP}
    INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE}
  )
elseif (GMP_LIBRARY STREQUAL "BUILD")
  # Download and build own libgmp version
  if (POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
  endif()
  SET(GMP_BUILD_CONFIG_ARGS "" CACHE STRING "Some user-specified gmp config options")
  option(ENABLE_GMP_STATIC "Option to statically link. Default is dynamic linking" OFF)

  if (ENABLE_GMP_STATIC)
    set(GMP_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
  else()
    set(GMP_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
  endif()

  cmake_host_system_information(RESULT N QUERY NUMBER_OF_PHYSICAL_CORES)
  if (N EQUAL 0)
    # Choose a "safe" amount
    set(N 8)
  endif()
  set(GMP_PARALLEL_BUILD_ARGS -j${N}) 

  message(STATUS "Building GMP with additional options: ${GMP_BUILD_CONFIG_ARGS}")
  include(ExternalProject)
  find_program(MAKE_EXE NAMES make gmake nmake)
  set(libgmp_INSTALL_DIR "${CMAKE_BINARY_DIR}/libgmp")
  ExternalProject_Add(libgmp_external
    PREFIX ${libgmp_INSTALL_DIR}
    URL               https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
    URL_HASH          SHA256=a3c2b80201b89e68616f4ad30bc66aee4927c3ce50e33929ca819d5c43538898
    CONFIGURE_COMMAND ${libgmp_INSTALL_DIR}/src/libgmp_external/configure --prefix=${libgmp_INSTALL_DIR} ${GMP_BUILD_CONFIG_ARGS}
    BUILD_COMMAND     ${MAKE_EXE} ${GMP_PARALLEL_BUILD_ARGS}
    INSTALL_COMMAND   ${MAKE_EXE} install
    BUILD_BYPRODUCTS  ${libgmp_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${GMP_LIB_SUFFIX}
  )

  # Needed to avoid errors about missing directory when creating the GMP target
  file(MAKE_DIRECTORY ${libgmp_INSTALL_DIR}/include)

  if(ENABLE_GMP_STATIC)
    add_library(GMP STATIC IMPORTED)
    set_target_properties(GMP PROPERTIES
      IMPORTED_LOCATION ${libgmp_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gmp${GMP_LIB_SUFFIX}
      INTERFACE_INCLUDE_DIRECTORIES ${libgmp_INSTALL_DIR}/include
    )
  else()
    add_library(GMP SHARED IMPORTED)
    set_target_properties(GMP PROPERTIES
      IMPORTED_LOCATION ${libgmp_INSTALL_DIR}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}gmp${GMP_LIB_SUFFIX}
      INTERFACE_INCLUDE_DIRECTORIES ${libgmp_INSTALL_DIR}/include
    )
  endif()

  add_dependencies(GMP libgmp_external)
elseif (GMP_LIBRARY STREQUAL "MINI")
  # Use mini-gmp
  message(STATUS "Using mini-GMP")

  include(CheckTypeSize)

  add_library(GMP STATIC
    ${PROJECT_SOURCE_DIR}/src/mini-gmp/mini-gmp.c ${PROJECT_SOURCE_DIR}/src/mini-gmp/mini-gmp-extra.c)
  target_include_directories(GMP PRIVATE ${PROJECT_SOURCE_DIR}/src/common/generic/include) # for tutil.h
  target_include_directories(GMP INTERFACE ${PROJECT_SOURCE_DIR}/src/mini-gmp)
  set_source_files_properties(${PROJECT_SOURCE_DIR}/src/mini-gmp/mini-gmp.c PROPERTIES COMPILE_OPTIONS "-w")

  set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/src/mini-gmp")
  set(CMAKE_EXTRA_INCLUDE_FILES "mini-gmp.h")
  check_type_size("mp_limb_t" MP_LIMB_T_BYTES)

  math(EXPR GMP_LIMB_BITS "${MP_LIMB_T_BYTES} * 8")

  add_compile_definitions(GMP_LIMB_BITS=${GMP_LIMB_BITS})
  add_compile_definitions(MINI_GMP)
else()
  message(FATAL_ERROR "Invalid choice for GMP_LIBRARY: ${GMP_LIBRARY}")
endif()