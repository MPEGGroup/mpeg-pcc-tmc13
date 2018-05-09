cmake_minimum_required(VERSION 3.0)

if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:tbb_debug.lib")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:tbb.lib")
endif()

include(CheckSymbolExists)
check_symbol_exists(getrusage sys/resource.h HAVE_GETRUSAGE)

##
# Determine the software version from VCS
# Fallback to descriptive version if VCS unavailable
set(VERSION_FALLBACK "unconfirmed-release-1.0b")
set(VERSION_FILE ${CMAKE_CURRENT_BINARY_DIR}/version.cpp)
set(VERSION_FILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp.in)
set(VERSION_CACHE_FILE ${CMAKE_CURRENT_BINARY_DIR}/version.cache)
find_package(Git)
add_custom_target(
  genversion
  COMMAND ${CMAKE_COMMAND} -D OUTPUT=${VERSION_FILE}
                           -D TEMPLATE=${VERSION_FILE_IN}
                           -D VERSION_FALLBACK=${VERSION_FALLBACK}
                           -D VERSION_EXTRA=${VERSION_EXTRA}
                           -D VERSION_CACHE_FILE=${VERSION_CACHE_FILE}
                           -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                           -P ${CMAKE_SOURCE_DIR}/scripts/genversion.cmake
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

add_custom_command(
  OUTPUT ${VERSION_FILE}
  DEPENDS genversion
)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/TMC3Config.h.in"
  "${PROJECT_BINARY_DIR}/tmc3/TMC3Config.h"
)

file(GLOB PROJECT_IN_FILES "*.in")
file(GLOB PROJECT_INC_FILES
  "*.h"
  "../dependencies/nanoflann/*.hpp"
  "../dependencies/nanoflann/*.h"
  "../dependencies/arithmetic-coding/inc/*.h"
  "../dependencies/program-options-lite/*.h"
)
file(GLOB PROJECT_INL_FILES "*.inl")
file(GLOB PROJECT_CPP_FILES
  "*.cpp"
  "../dependencies/arithmetic-coding/src/*.cpp"
  "../dependencies/program-options-lite/*.cpp"
)
file(GLOB PROJECT_C_FILES "*.c")
file(GLOB PROJECT_CL_FILES "*.cl")
source_group (inc FILES ${PROJECT_INC_FILES})
source_group (input FILES ${PROJECT_IN_FILES})
source_group (inl FILES ${PROJECT_INL_FILES})
source_group (cpp FILES ${PROJECT_CPP_FILES})
source_group (c FILES ${PROJECT_C_FILES})
source_group (cl FILES ${PROJECT_CL_FILES})

include_directories(
  "${PROJECT_BINARY_DIR}/tmc3"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/nanoflann"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/tbb/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/arithmetic-coding/inc"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/program-options-lite"
)

add_executable (tmc3
  ${PROJECT_CPP_FILES}
  ${PROJECT_C_FILES}
  ${PROJECT_CL_FILES}
  ${PROJECT_PUB_FILES}
  ${PROJECT_INC_FILES}
  ${PROJECT_IN_FILES}
  ${PROJECT_INL_FILES}
  ${VERSION_FILE}
)
add_dependencies(tmc3 genversion)
target_link_libraries (tmc3 tbb_static)

install (TARGETS tmc3 DESTINATION bin)
