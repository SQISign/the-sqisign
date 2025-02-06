set(SOURCE_FILES_MP_GENERIC_REF
    mp.c
)

add_library(${LIB_MP} STATIC ${SOURCE_FILES_MP_GENERIC_REF})
target_include_directories(${LIB_MP} PRIVATE ${INC_PUBLIC} ${INC_COMMON} ${INC_MP})
target_compile_options(${LIB_MP} PRIVATE ${C_OPT_FLAGS})