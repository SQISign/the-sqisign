set(SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF
    torsion_constants.c
    quaternion_data.c
    endomorphism_action.c
    klpt_constants.c
)

add_library(${LIB_PRECOMP_${SVARIANT_UPPER}} ${SOURCE_FILES_PRECOMP_${SVARIANT_UPPER}_REF})
target_include_directories(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE common ${INC_INTBIG} ${INC_QUATERNION} ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/ec/ref/include ${PROJECT_SOURCE_DIR}/src/ec/ref/${SVARIANT_LOWER}/include ${PROJECT_SOURCE_DIR}/src/gf/ref/${SVARIANT_LOWER}/include ${INC_PRECOMP_${SVARIANT_UPPER}} ${INC_COMMON})
target_compile_options(${LIB_PRECOMP_${SVARIANT_UPPER}} PRIVATE ${C_OPT_FLAGS})

add_custom_command(
    OUTPUT
        "${CMAKE_CURRENT_SOURCE_DIR}/torsion_constants.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/torsion_constants.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/klpt_constants.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/klpt_constants.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/encoded_sizes.h"
        "${PROJECT_SOURCE_DIR}/src/nistapi/lvl5/api.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/quaternion_data.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/quaternion_data.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/endomorphism_action.c"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/endomorphism_action.h"
    COMMAND
        echo "Please run manually: make precomp"
)

find_program(SAGEMATH sage)
add_custom_target(precomp_${SVARIANT_LOWER}
    DEPENDS
        "./sqisign_parameters.txt"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_torsion_constants.sage"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_klpt_constants.sage"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_sizes.sage"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_quaternion_data.sage"
    COMMAND
        "${SAGEMATH}" "${PROJECT_SOURCE_DIR}/scripts/precompute_endomorphism_action.sage"
    WORKING_DIRECTORY
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

set_directory_properties(PROPERTIES CLEAN_NO_CUSTOM true)
set_target_properties(precomp_${SVARIANT_LOWER} PROPERTIES EXCLUDE_FROM_ALL TRUE)

add_dependencies(precomp precomp_lvl5)
