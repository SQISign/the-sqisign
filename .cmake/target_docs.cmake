# Find Doxygen
find_package(Doxygen)

# Check if Doxygen is found
if(DOXYGEN_FOUND)
  # Set input and output directories
  set(DOXYGEN_CONF ${PROJECT_SOURCE_DIR}/Doxyfile)

  # Add a target to generate Doxygen documentation using 'make doc'
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_CONF}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM )
else()
    message(FATAL_ERROR "Doxygen not found - can't generate docs.")
endif(DOXYGEN_FOUND)
