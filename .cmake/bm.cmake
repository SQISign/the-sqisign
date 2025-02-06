add_custom_target(bm
    COMMAND ${CMAKE_COMMAND} -E echo "Running all benchmarks..."
)

foreach(bm_bin ${BM_BINS})
    add_custom_command(
        TARGET bm
        POST_BUILD
        COMMAND $<TARGET_FILE:${bm_bin}>
        COMMENT "Running ${bm_bin}"
    )
endforeach()