cmake_minimum_required(VERSION 2.8.12)

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

option(GLVIS_USE_FREETYPE
  "Use freetype and fontconfig for rendinering and finding fonts."
  ON)

option(GLVIS_USE_GLX10
  "Use GLX 1.0 calls. Use if X server doesn't support GLX 1.3."
  OFF)

option(GLVIS_USE_CET_COLOR
    "Use additional perceptually uniform color maps."
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

if (GLVIS_USE_GLX10)
  list(APPEND _glvis_compile_defs "GLVIS_GLX10")
endif (GLVIS_USE_GLX10)

if (GLVIS_USE_CET_COLOR)
    list(APPEND _glvis_compile_defs "GLVIS_USE_CET_COLOR")
endif (GLVIS_USE_CET_COLOR)

list(APPEND _glvis_compile_defs "GLVIS_MULTISAMPLE=${GLVIS_MULTISAMPLE}")
list(APPEND _glvis_compile_defs "GLVIS_MS_LINEWIDTH=${GLVIS_MS_LINEWIDTH}")

if (CMAKE_BUILD_TYPE MATCHES "Debug|debug|DEBUG")
  list(APPEND _glvis_compile_defs "GLVIS_DEBUG")
endif()

# Include paths and libraries needed by MFEM
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
list(APPEND _glvis_include_dirs "${MFEM_INCLUDE_DIRS}")
list(APPEND _glvis_libraries "${MFEM_LIBRARIES}")

# Find X11. Only need libX11.
find_package(X11 REQUIRED)
list(APPEND _glvis_include_dirs "${X11_X11_INCLUDE_PATH}")
list(APPEND _glvis_libraries "${X11_X11_LIB}")

# Find OpenGL
if (NOT APPLE)
  find_package(OpenGL REQUIRED)
else()

  find_path(OPENGL_INCLUDE_DIR GL/gl.h
    HINTS ${X11_X11_INCLUDE_PATH} ${OPENGL_DIR} $ENV{OPENGL_DIR}
    DOC "OpenGL include path"
    NO_DEFAULT_PATH)
  find_path(OPENGL_INCLUDE_DIR GL/gl.h)

  if (NOT OPENGL_LIBRARIES)
    # I know this is ugly, but the internet suggests manually editing
    # the CMakeCache.txt file. That seems worse.
    get_filename_component(_x11_base_dir "${X11_X11_LIB}" DIRECTORY)

    find_library(_gl_library GL
      HINTS ${_x11_base_dir} ${OPENGL_DIR} $ENV{OPENGL_DIR}
      DOC "OpenGL GL library."
      NO_DEFAULT_PATH)
    find_library(_gl_library GL)
    if (NOT _gl_library)
      message(FATAL_ERROR "OpenGL not found. Please set OPENGL_DIR.")
    endif()

    find_library(_glu_library GLU
      HINTS ${_x11_base_dir} ${OPENGL_DIR} $ENV{OPENGL_DIR}
      DOC "OpenGL GLU library."
      NO_DEFAULT_PATH)
    find_library(_glu_library GLU)
    if (NOT _glu_library)
      message(FATAL_ERROR "GLU library not found. Please set OPENGL_DIR.")
    endif()

    set(OPENGL_LIBRARIES "${_gl_library}" "${_glu_library}")
    message(STATUS "Found OpenGL: ${OPENGL_LIBRARIES}")
  endif (NOT OPENGL_LIBRARIES)
endif (NOT APPLE)

list(APPEND _glvis_include_dirs "${OPENGL_INCLUDE_DIR}")
list(APPEND _glvis_libraries "${OPENGL_LIBRARIES}")

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

# Find FreeType and Fontconfig.
if (GLVIS_USE_FREETYPE)
  find_package(Freetype)
  if (NOT FREETYPE_FOUND)
    message(WARNING "FreeType not found. Disabling FreeType support.")
  else()
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
      message(WARNING "Fontconfig not found. Disabling FreeType support.")
    endif (FONTCONFIG_LIBRARY)
  endif (NOT FREETYPE_FOUND)
endif (GLVIS_USE_FREETYPE)

# Find threading library
set(CMAKE_THREAD_PREFER_PTHREAD ON)
find_package(Threads REQUIRED)
if (NOT CMAKE_USE_PTHREADS_INIT)
  message(FATAL_ERROR "The required pthreads library was not found.")
else()
  list(APPEND _glvis_libraries "${CMAKE_THREAD_LIBS_INIT}")
endif()

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
