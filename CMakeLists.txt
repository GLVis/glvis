cmake_minimum_required(VERSION 2.8.11)

project(glvis)

find_package(MFEM REQUIRED)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DGLVIS_DEBUG)
endif()

if(${MFEM_USE_MPI})
  find_package(MPI REQUIRED)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

# TODO: Grab this from MFEM?
# A handy function to add the current source directory to a local
# filename. To be used for creating a list of sources.
function(convert_filenames_to_full_paths NAMES)
  unset(tmp_names)
  foreach(name ${${NAMES}})
    list(APPEND tmp_names ${CMAKE_CURRENT_SOURCE_DIR}/${name})
  endforeach()
  set(${NAMES} ${tmp_names} PARENT_SCOPE)
endfunction()

include_directories(${MFEM_INCLUDE_DIRS})

set(SOURCES)
set(HEADERS)

list(APPEND LIBRARIES pthread)

find_package(PkgConfig)

pkg_check_modules(LIBX11 QUIET x11)
if(${LIBX11_FOUND})
  list(APPEND LIBRARIES ${LIBX11_LIBRARIES})
  include_directories(${LIBX11_INCLUDE_DIRS})
endif()

pkg_check_modules(LIBGLU QUIET glu)
if(${LIBGLU_FOUND})
  list(APPEND LIBRARIES ${LIBGLU_LIBRARIES})
  include_directories(${LIBGLU_INCLUDE_DIRS})
endif()

pkg_check_modules(LIBPNG QUIET libpng)
if(${LIBPNG_FOUND})
  add_definitions(-DGLVIS_USE_LIBPNG)
  list(APPEND LIBRARIES ${LIBPNG_LIBRARIES})
  include_directories(${LIBPNG_INCLUDE_DIRS})
endif()

pkg_check_modules(FONTCONFIG QUIET fontconfig)
if(${FONTCONFIG_FOUND})
  add_definitions(-DGLVIS_USE_FREETYPE)
  list(APPEND LIBRARIES ${FONTCONFIG_LIBRARIES})
  include_directories(${FONTCONFIG_INCLUDE_DIRS})
endif()

add_definitions(
  -DGLVIS_MULTISAMPLE=4
  -DGLVIS_MS_LINEWIDTH=1.4
  )

add_subdirectory(lib)

add_library(glvis ${SOURCES})

target_link_libraries(glvis PUBLIC mfem)

add_executable(glvis-exe glvis.cpp)
set_target_properties(glvis-exe PROPERTIES OUTPUT_NAME glvis)


target_link_libraries(
  glvis-exe
  glvis
  ${LIBRARIES})
