cmake_minimum_required(VERSION 3.0)

if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:tbb_debug.lib")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:tbb.lib")
endif()

include(CheckSymbolExists)
check_symbol_exists(getrusage sys/resource.h HAVE_GETRUSAGE)

set (TMC3_VERSION_MAJOR 0)
set (TMC3_VERSION_MINOR 1)

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
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/nanoflann"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/tbb/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/arithmetic-coding/inc"
  "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies/program-options-lite"
)

add_executable(tmc3 ${PROJECT_CPP_FILES} ${PROJECT_C_FILES} ${PROJECT_CL_FILES} ${PROJECT_PUB_FILES} ${PROJECT_INC_FILES} ${PROJECT_IN_FILES} ${PROJECT_INL_FILES})
target_link_libraries (tmc3 tbb_static)

install (TARGETS tmc3 DESTINATION bin)
