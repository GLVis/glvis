# Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-443271.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see https://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# CMake v3.10 is where "Imported Targets" like "Freetype::Freetype" are
# introduced.
cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

project(glvis NONE)

# Import MFEM. The following variables can be used to help CMake find MFEM:
#  * MFEM_DIR - absolute path to the MFEM build or install prefix.
#  * mfem_DIR - absolute path to where MFEMConfig.cmake is.
message(STATUS "Looking for mfem ...")
set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
if (MFEM_DIR)
   find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
else()
   find_package(mfem REQUIRED NAMES MFEM)
endif()
message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
# Use the same C++ compiler as MFEM. This is needed when MFEM was built using
# an MPI wrapper and we do not have explicitly the MPI compile and link flags.
if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
endif()

enable_language(C)
enable_language(CXX)

# Default options match the Makefile
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

option(GLVIS_USE_LIBTIFF
  "Use libtiff for taking screenshots internally"
  OFF)

option(GLVIS_USE_LIBPNG
  "Use libpng for taking screenshots internally"
  ON)

#
# Handle a few other definitions
#

# Default multisampling mode
if (NOT GLVIS_MULTISAMPLE)
  set(GLVIS_MULTISAMPLE 4)
endif (NOT GLVIS_MULTISAMPLE)

# Default multisampling line-width
if (NOT GLVIS_MS_LINEWIDTH)
  if (NOT APPLE)
    set(GLVIS_MS_LINEWIDTH 1.4)
  else()
    # This value seems to work better on Macs
    set(GLVIS_MS_LINEWIDTH 1.0)
  endif()
endif (NOT GLVIS_MS_LINEWIDTH)

#
# Start finding everything
#

set(_glvis_compile_defs)
set(_glvis_compile_opts)
set(_glvis_include_dirs)
set(_glvis_libraries)

list(APPEND _glvis_compile_defs "GLVIS_MULTISAMPLE=${GLVIS_MULTISAMPLE}")
list(APPEND _glvis_compile_defs "GLVIS_MS_LINEWIDTH=${GLVIS_MS_LINEWIDTH}")
if (NOT WIN32)
  list(APPEND _glvis_compile_defs "GLVIS_USE_LOGO")
else()
  list(APPEND _glvis_compile_defs "_USE_MATH_DEFINES")
endif()

if(NOT WIN32)
  find_program(XXD_FOUND xxd)
  if(NOT XXD_FOUND)
    message(FATAL_ERROR "Required tool not found: xxd.")
  endif()
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND _glvis_compile_defs "GLVIS_DEBUG")
  list(APPEND _glvis_compile_opts "-Wall")
endif()

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
list(APPEND _glvis_include_dirs "${MFEM_INCLUDE_DIRS}")
list(APPEND _glvis_libraries "${MFEM_LIBRARIES}")

if (NOT EMSCRIPTEN)

  # Find OpenGL
  find_package(OpenGL REQUIRED)
  list(APPEND _glvis_libraries OpenGL::GL)

  # Find 'SDL2' which is not part of the CMake standard modules
  list(APPEND SDL2_DIR "../SDL2")
  set(SDL2_DIR_SAVE ${SDL2_DIR})
  find_package(SDL2 QUIET HINTS ${SDL2_DIR})
  # The SDL2 target 'SDL2::SDL2' is not always defined, so instead we use
  # directly 'SDL2_INCLUDE_DIRS' and 'SDL2_LIBRARIES'.
  if (TARGET SDL2::SDL2)
    list(APPEND _glvis_libraries SDL2::SDL2)
  elseif (SDL2_INCLUDE_DIRS AND SDL2_LIBRARIES)
    # SDL2_INCLUDE_DIRS is defined with 'SDL2' at the end of the path we need to
    # strip that.
    set(new_list_)
    foreach (path IN LISTS SDL2_INCLUDE_DIRS)
      string(REGEX REPLACE "/include/SDL2$" "/include" new_path_ "${path}")
      list(APPEND new_list_ "${new_path_}")
    endforeach()
    set(SDL2_INCLUDE_DIRS ${new_list_})
    list(APPEND _glvis_include_dirs "${SDL2_INCLUDE_DIRS}")
    list(APPEND _glvis_libraries "${SDL2_LIBRARIES}")
  else()
    # find_package() did not work; try more direct approach.
    set(SDL2_DIR ${SDL2_DIR_SAVE})
    find_path(SDL2_INCLUDE_DIRS "SDL2/SDL.h"
      HINTS ${SDL2_DIR} ENV SDL2_DIR
      PATH_SUFFIXES "include"
      DOC "SDL2 include path.")
    find_library(SDL2_LIBRARIES "SDL2"
      HINTS ${SDL2_DIR} ENV SDL2_DIR
      PATH_SUFFIXES "lib"
      DOC "SDL2 library.")
    if (NOT (SDL2_INCLUDE_DIRS AND SDL2_LIBRARIES))
      message(FATAL_ERROR "SDL2 library not found. Please set SDL2_DIR.")
    endif()
    list(APPEND _glvis_include_dirs "${SDL2_INCLUDE_DIRS}")
    list(APPEND _glvis_libraries "${SDL2_LIBRARIES}")
  endif()
  message(STATUS "SDL2 found: ${SDL2_LIBRARIES}")
  message(STATUS "SDL2_INCLUDE_DIRS = ${SDL2_INCLUDE_DIRS}")

  # Find GLEW
  # GLEW search path can be set using CMAKE_PREFIX_PATH
  list(APPEND CMAKE_PREFIX_PATH "../glew")
  find_package(GLEW REQUIRED)
  list(APPEND _glvis_libraries GLEW::GLEW)

  # Find 'glm' which is not part of the CMake standard modules
  list(APPEND GLM_DIR "../glm")
  find_package(glm QUIET HINTS ${GLM_DIR})
  if (NOT GLM_INCLUDE_DIRS)
    find_path(GLM_INCLUDE_DIRS "glm/glm.hpp"
      HINTS ${GLM_DIR} ENV GLM_DIR
      PATH_SUFFIXES "include"
      DOC "GLM include path.")
    if (NOT GLM_INCLUDE_DIRS)
      message(FATAL_ERROR "GLM headers not found. Please set GLM_DIR.")
    endif()
  endif()
  message(STATUS "GLM headers found: ${GLM_INCLUDE_DIRS}")
  list(APPEND _glvis_include_dirs "${GLM_INCLUDE_DIRS}")

  # Find TIFF
  if (GLVIS_USE_LIBTIFF)
    find_package(TIFF)
    if (TIFF_FOUND)
      list(APPEND _glvis_compile_defs "GLVIS_USE_LIBTIFF")
      list(APPEND _glvis_libraries TIFF::TIFF)
    else()
      message(WARNING "TIFF library not found. TIFF disabled.")
      set(GLVIS_USE_LIBTIFF OFF)
    endif (TIFF_FOUND)
  endif (GLVIS_USE_LIBTIFF)

  # Find PNG
  if (GLVIS_USE_LIBPNG)
    find_package(PNG)
    if (PNG_FOUND)
      list(APPEND _glvis_compile_defs "GLVIS_USE_LIBPNG")
      list(APPEND _glvis_libraries PNG::PNG)
    else()
      message(WARNING "PNG library not found. PNG disabled.")
      set(GLVIS_USE_LIBPNG OFF)
    endif (PNG_FOUND)
  endif (GLVIS_USE_LIBPNG)

  # Find FreeType.
  find_package(Freetype REQUIRED)
  # Find FontConfig (FindFontconfig.cmake was added in CMake 3.14)
  find_library(FONTCONFIG_LIBRARY fontconfig
    HINTS ${FONTCONFIG_DIR} $ENV{FONTCONFIG_DIR}
    DOC "The fontconfig library for use with FreeType."
    NO_DEFAULT_PATH)
  find_library(FONTCONFIG_LIBRARY fontconfig)
  if (FONTCONFIG_LIBRARY)
    list(APPEND _glvis_libraries Freetype::Freetype)
    # Need FONTCONFIG_INCLUDE_DIRS?
    list(APPEND _glvis_libraries "${FONTCONFIG_LIBRARY}")
    message(STATUS "Found Fontconfig: ${FONTCONFIG_LIBRARY}")
  else()
    message(STATUS "Fontconfig not found. Please set FONTCONFIG_DIR.")
    message(FATAL_ERROR "Fontconfig not found.")
  endif (FONTCONFIG_LIBRARY)

  # Find threading library
  find_package(Threads REQUIRED)
  list(APPEND _glvis_libraries "${CMAKE_THREAD_LIBS_INIT}")
  if (CMAKE_USE_PTHREADS_INIT)
    message(STATUS "System threading library: pthreads")
  elseif (CMAKE_USE_WIN32_THREADS_INIT)
    message(STATUS "System threading library: Win32 threads")
  elseif (Threads_FOUND)
    message(STATUS "System threading library: other")
  endif()

else(NOT EMSCRIPTEN)

  find_package(glm REQUIRED)
  list(APPEND _glvis_include_dirs "${GLM_INCLUDE_DIRS}")
  set(_emscripten_opts)
  # Enable embind
  list(APPEND _emscripten_opts "--bind")

  # OpenGL, SDL2, and GLEW are provided by the Emscripten runtime
  # Enable SDL2
  list(APPEND _emscripten_opts "-s USE_SDL=2")
  list(APPEND _emscripten_opts
      "-s DISABLE_DEPRECATED_FIND_EVENT_TARGET_BEHAVIOR=0")

  # Set WebGL options
  list(APPEND _emscripten_opts
      "-s GL_ASSERTIONS=1 -s GL_DEBUG=1 -s MAX_WEBGL_VERSION=2")

  # Enable Freetype port
  list(APPEND _emscripten_opts "-s USE_FREETYPE=1")

  # Module export options
  list(APPEND _emscripten_opts "-s ENVIRONMENT=web")
  list(APPEND _emscripten_opts "-s MODULARIZE=1")
  list(APPEND _emscripten_opts "-s EXPORT_NAME=glvis")

  # Other emscripten options
  list(APPEND _emscripten_opts "-s ALLOW_MEMORY_GROWTH=1")
  list(APPEND _emscripten_opts "-s SINGLE_FILE=1")
  list(APPEND _emscripten_opts "--no-heap-copy")

  # Since we don't have access to fontconfig, embed a font
  if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/OpenSans.ttf")
    message(FATAL_ERROR "Preload font file \"OpenSans.ttf\" not found.")
  endif()
  list(APPEND _emscripten_opts
       "--embed-file ${CMAKE_CURRENT_SOURCE_DIR}/OpenSans.ttf@OpenSans.ttf")

endif(NOT EMSCRIPTEN)

message(STATUS "GLVis build type: CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "GLVis defines: ${_glvis_compile_defs}")
# message(STATUS "GLVis opts: ${_glvis_compile_opts}")
# message(STATUS "GLVis include dirs: ${_glvis_include_dirs}")
# message(STATUS "GLVis libraries: ${_glvis_libraries}")
message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")


#
# Setup the GLVis library target
#

add_subdirectory(lib)
add_subdirectory(share)

if(NOT EMSCRIPTEN)

  # Setup the GLVis executable
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON CACHE BOOL "")
  set(CMAKE_INSTALL_RPATH "${MFEM_LIBRARY_DIR}" CACHE PATH "")

  add_executable(glvis-exe glvis.cpp)
  set_target_properties(glvis-exe PROPERTIES OUTPUT_NAME glvis)

  target_link_libraries(glvis-exe PRIVATE glvis)

  if (WIN32)
    target_sources(glvis-exe
      PRIVATE
        "share/windows.manifest"
        "${CMAKE_CURRENT_BINARY_DIR}/share/resource.rc")
  else()
    target_link_libraries(glvis-exe PRIVATE glvis_logo)
  endif (WIN32)

  # Install the executable
  install(TARGETS glvis-exe RUNTIME DESTINATION bin)

  # Install the gnutls helper script
  if (MFEM_USE_GNUTLS)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/glvis-keygen.sh
      DESTINATION bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif (MFEM_USE_GNUTLS)

  set(GLVIS_APP_ICON ${CMAKE_CURRENT_SOURCE_DIR}/share/GLVis.icns)
  set(GLVIS_APP_CREDITS ${CMAKE_CURRENT_SOURCE_DIR}/share/Credits.rtf)
  set_source_files_properties(
    ${GLVIS_APP_ICON} ${GLVIS_APP_CREDITS}
    PROPERTIES
      MACOSX_PACKAGE_LOCATION "Resources")

  add_executable(app MACOSX_BUNDLE glvis.cpp ${GLVIS_APP_ICON} ${GLVIS_APP_CREDITS})
  set_target_properties(app
    PROPERTIES
      EXCLUDE_FROM_ALL TRUE
      OUTPUT_NAME GLVis
      MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/share/Info.cmake.plist.in)
  target_link_libraries(app PRIVATE glvis glvis_logo)
  install(TARGETS app
    RUNTIME DESTINATION .
    BUNDLE DESTINATION .)
  install(CODE [[
    include (BundleUtilities)
    fixup_bundle("${CMAKE_INSTALL_PREFIX}/GLVis.app" "" "")
  ]] COMPONENT RUNTIME)

endif(NOT EMSCRIPTEN)

if(ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(ENABLE_TESTS)
