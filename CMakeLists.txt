cmake_minimum_required(VERSION 3.1)

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
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

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
    set(GLVIS_MS_LINEWIDTH 0.01)
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
# Find X11. Only need libX11.
find_package(X11 REQUIRED)
list(APPEND _glvis_include_dirs "${X11_X11_INCLUDE_PATH}")
list(APPEND _glvis_libraries "${X11_X11_LIB}")

# Find OpenGL
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)
find_package(GLEW REQUIRED)
find_package(glm REQUIRED)

list(APPEND _glvis_include_dirs "${OPENGL_INCLUDE_DIR}")
list(APPEND _glvis_libraries "${OPENGL_LIBRARIES}")
list(APPEND _glvis_include_dirs "${SDL2_INCLUDE_DIRS}")
list(APPEND _glvis_libraries "${SDL2_LIBRARIES}")
list(APPEND _glvis_include_dirs "${GLEW_INCLUDE_DIRS}")
list(APPEND _glvis_libraries "${GLEW_LIBRARIES}")
list(APPEND _glvis_include_dirs "${GLM_INCLUDE_DIRS}")


# Find TIFF
if (GLVIS_USE_LIBTIFF)
  find_package(TIFF)
  if (TIFF_FOUND)
    list(APPEND _glvis_compile_defs "GLVIS_USE_LIBTIFF")
    list(APPEND _glvis_include_dirs "${TIFF_INCLUDE_DIRS}")
    list(APPEND _glvis_libraries "${TIFF_LIBRARIES}")
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
    list(APPEND _glvis_include_dirs "${PNG_INCLUDE_DIRS}")
    list(APPEND _glvis_libraries "${PNG_LIBRARIES}")
  else()
    message(WARNING "PNG library not found. PNG disabled.")
    set(GLVIS_USE_LIBPNG OFF)
  endif (PNG_FOUND)
endif (GLVIS_USE_LIBPNG)


# Find FreeType.
find_package(Freetype)
if (NOT FREETYPE_FOUND)
    message(FATAL_ERROR "FreeType not found.")
endif (NOT FREETYPE_FOUND)
# Find FontConfig
find_library(FONTCONFIG_LIBRARY fontconfig
HINTS ${FONTCONFIG_DIR} $ENV{FONTCONFIG_DIR}
DOC "The fontconfig library for use with FreeType."
NO_DEFAULT_PATH)
find_library(FONTCONFIG_LIBRARY fontconfig)
if (FONTCONFIG_LIBRARY)
  list(APPEND _glvis_compile_defs "GLVIS_USE_FREETYPE")
  list(APPEND _glvis_include_dirs "${FREETYPE_INCLUDE_DIRS}")
  list(APPEND _glvis_libraries "${FREETYPE_LIBRARIES}")
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

else()
  find_package(glm REQUIRED)
  list(APPEND _glvis_include_dirs "${GLM_INCLUDE_DIRS}")
  set(_emscripten_opts)
  # Enable embind
  list(APPEND _emscripten_opts "--bind")

  # OpenGL, SDL2, and GLEW are provided by the Emscripten runtime
  # Enable SDL2
  list(APPEND _emscripten_opts "-s USE_SDL=2")

  # Set WebGL options
  list(APPEND _emscripten_opts "-s GL_ASSERTIONS=1 -s GL_DEBUG=1 -s USE_WEBGL2=1")

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
  list(APPEND _emscripten_opts "--embed-file ${CMAKE_CURRENT_SOURCE_DIR}/OpenSans.ttf@OpenSans.ttf")

endif(NOT EMSCRIPTEN)

message(STATUS "GLVis build type: CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
# message(STATUS "GLVis defines: ${_glvis_compile_defs}")
# message(STATUS "GLVis opts: ${_glvis_compile_opts}")
# message(STATUS "GLVis include dirs: ${_glvis_include_dirs}")
# message(STATUS "GLVis libraries: ${_glvis_libraries}")
message(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")


#
# Setup the GLVis library target
#

add_subdirectory(lib)

if(NOT EMSCRIPTEN)
#
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
