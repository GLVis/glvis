# CMake v3.10 is where "Imported Targets" like "Freetype::Freetype" are
# introduced.
cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 11)
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
  OFF)

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

if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND _glvis_compile_defs "GLVIS_DEBUG")
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
  find_package(SDL2 QUIET)
  if (TARGET SDL2::SDL2)
    list(APPEND _glvis_libraries SDL2::SDL2)
    message(STATUS "SDL2 found: SDL2::SDL2")
  else()
    if (NOT (SDL2_INCLUDE_DIRS AND SDL2_LIBRARIES))
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
    endif()
    list(APPEND _glvis_include_dirs "${SDL2_INCLUDE_DIRS}")
    list(APPEND _glvis_libraries "${SDL2_LIBRARIES}")
    message(STATUS "SDL2 found: ${SDL2_LIBRARIES}")
  endif()

  # Find GLEW
  find_package(GLEW REQUIRED)
  list(APPEND _glvis_libraries GLEW::GLEW)

  # Find 'glm' which is not part of the CMake standard modules
  find_package(glm QUIET)
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

  if (NOT (GLVIS_USE_LIBTIFF OR GLVIS_USE_LIBPNG))
    list(APPEND _glvis_compile_defs "SDL_SCREENSHOTS")
  endif()

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
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads REQUIRED)
  if (NOT CMAKE_USE_PTHREADS_INIT)
    message(FATAL_ERROR "The required pthreads library was not found.")
  else()
    list(APPEND _glvis_libraries "${CMAKE_THREAD_LIBS_INIT}")
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

  # Set WebGL options
  list(APPEND _emscripten_opts
       "-s GL_ASSERTIONS=1 -s GL_DEBUG=1 -s USE_WEBGL2=1")

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

if(NOT EMSCRIPTEN)

  # Setup the GLVis executable
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON CACHE BOOL "")
  set(CMAKE_INSTALL_RPATH "${MFEM_LIBRARY_DIR}" CACHE PATH "")

  add_executable(glvis-exe glvis.cpp)
  set_target_properties(glvis-exe PROPERTIES OUTPUT_NAME glvis)

  target_link_libraries(glvis-exe PRIVATE glvis)

  # Install the executable
  install(TARGETS glvis-exe RUNTIME DESTINATION bin)

  # Install the gnutls helper script
  if (MFEM_USE_GNUTLS)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/glvis-keygen.sh
      DESTINATION bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  endif (MFEM_USE_GNUTLS)

endif(NOT EMSCRIPTEN)
