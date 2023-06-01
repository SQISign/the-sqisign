add_executable(fp2.test_${SVARIANT_LOWER} ${ECX_DIR}/test/fp2-test.c)
	target_include_directories(fp2.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include ../include ${INC_GF_${SVARIANT_UPPER}} ${INC_EC} ${INC_COMMON})
	target_link_libraries(fp2.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}})
		
add_executable(poly-mul.test_${SVARIANT_LOWER} ${ECX_DIR}/test/poly-mul-test.c)
	target_include_directories(poly-mul.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON})
	target_link_libraries(poly-mul.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})
		
add_executable(poly-redc.test_${SVARIANT_LOWER} ${ECX_DIR}/test/poly-redc-test.c)
	target_include_directories(poly-redc.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include  ${INC_EC} ${INC_COMMON})
	target_link_libraries(poly-redc.test_${SVARIANT_LOWER} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})
	
add_executable(mont.test_${SVARIANT_LOWER} ${ECX_DIR}/test/mont-test.c)
	target_include_directories(mont.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} .)
	target_link_libraries(mont.test_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})

add_executable(ec.test_${SVARIANT_LOWER} ${ECX_DIR}/test/ec-test.c ${ECX_DIR}/test/test_extras.c)
	target_include_directories(ec.test_${SVARIANT_LOWER} PUBLIC ${ECX_DIR}/test ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} .)
	target_link_libraries(ec.test_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})

add_executable(velu.test_${SVARIANT_LOWER} ${ECX_DIR}/test/velu-test.c)
	target_include_directories(velu.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} .)
	target_link_libraries(velu.test_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})

add_executable(isog.test_${SVARIANT_LOWER} ${ECX_DIR}/test/isog-test.c)
	target_include_directories(isog.test_${SVARIANT_LOWER} PUBLIC ${INC_GF_${SVARIANT_UPPER}} ${INC_INTBIG} ${INC_PRECOMP_${SVARIANT_UPPER}} ${PROJECT_SOURCE_DIR}/include  ../include ${INC_EC} ${INC_COMMON} .)
	target_link_libraries(isog.test_${SVARIANT_LOWER} ${LIB_PRECOMP_${SVARIANT_UPPER}} ${LIB_INTBIG} ${LIB_GF_${SVARIANT_UPPER}} ${LIB_EC_${SVARIANT_UPPER}})


add_test(ec_fp2.test_${SVARIANT_LOWER} fp2.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_poly-mul.test_${SVARIANT_LOWER} poly-mul.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_poly-redc.test_${SVARIANT_LOWER} poly-redc.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_mont.test_${SVARIANT_LOWER} mont.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_ec.test_${SVARIANT_LOWER} ec.test_${SVARIANT_LOWER} test ${SQISIGN_TEST_REPS})
add_test(ec_velu.test_${SVARIANT_LOWER} velu.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})
add_test(ec_isog.test_${SVARIANT_LOWER} isog.test_${SVARIANT_LOWER} ${SQISIGN_TEST_REPS})