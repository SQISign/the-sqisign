# NIST KAT generation apps
foreach(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    add_executable(PQCgenKAT_sign_${SVARIANT_LOWER} PQCgenKAT_sign.c)
    target_link_libraries(PQCgenKAT_sign_${SVARIANT_LOWER} PRIVATE sqisign_${SVARIANT_LOWER}_test_nistapi)
    target_include_directories(PQCgenKAT_sign_${SVARIANT_LOWER} PRIVATE ../include)
endforeach()

# Examples with NIST API
foreach(SVARIANT ${SVARIANT_S})
    string(TOLOWER ${SVARIANT} SVARIANT_LOWER)
    add_executable(example_nistapi_${SVARIANT_LOWER} example_nistapi.c)
    target_link_libraries(example_nistapi_${SVARIANT_LOWER} PRIVATE sqisign_${SVARIANT_LOWER}_nistapi)
    target_include_directories(example_nistapi_${SVARIANT_LOWER} PRIVATE ../include ../src/${SVARIANT_LOWER})
endforeach()
