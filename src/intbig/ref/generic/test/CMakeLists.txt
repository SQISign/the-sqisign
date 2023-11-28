add_executable(sqisign_test_intbig test_intbig.c)
target_link_libraries(sqisign_test_intbig ${LIB_INTBIG} ${GMP} sqisign_common_sys)
target_include_directories(sqisign_test_intbig PRIVATE ../include ${INC_PUBLIC} ${INC_COMMON})

# MSAN and GMP lead to false positives, see
# https://gmplib.org/list-archives/gmp-bugs/2019-March/004529.html
if(NOT CMAKE_BUILD_TYPE STREQUAL "MSAN")
    add_test(sqisign_test_intbig sqisign_test_intbig ${SQISIGN_TEST_REPS} 256)
endif()
