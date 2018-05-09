# Generate a version file from a template based on the VCS version
# - Determine the current version string
#    => If git is unavailable, fallback to specified string
# - Test if it is the same as the cached version
# - If same, do nothing, otherwise generate output
if(NOT GIT_EXECUTABLE)
  set(VERSION ${VERSION_FALLBACK})
else()
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --long --always --dirty=+dirty
    OUTPUT_VARIABLE VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

if(EXISTS ${VERSION_CACHE_FILE})
  FILE(READ ${VERSION_CACHE_FILE} VERSION_CACHED)
endif()

FILE(WRITE ${VERSION_CACHE_FILE} ${VERSION}${VERSION_EXTRA})

if(NOT VERSION_CACHED STREQUAL VERSION OR ${TEMPLATE} IS_NEWER_THAN ${OUTPUT})
  configure_file ("${TEMPLATE}" "${OUTPUT}")
endif()
