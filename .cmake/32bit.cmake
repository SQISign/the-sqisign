set(CMAKE_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_NAME})
if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(CMAKE_SYSTEM_PROCESSOR i686)
endif()
set(GMP_LIBRARY "BUILD" CACHE STRING "" FORCE)
set(GMP_BUILD_CONFIG_ARGS "ABI=32" CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS "-m32" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-m32" CACHE STRING "" FORCE)
