cmake_minimum_required(VERSION 3.0)

project(glvis)

# Default options match the Makefile
option(GLVIS_USE_THREADS
  "Whether to enable threads in GLVis"
  ON)

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

#
# Handle a few other definitions
#

# Default multisampling mode
if (NOT GLVIS_MULTISAMPLE)
  set(GLVIS_MULTISAMPLE 4)
endif (NOT GLVIS_MULTISAMPLE)

# Default multisampling line-width
if (NOT GLVIS_MS_LINEWIDTH)
  set(GLVIS_MS_LINEWIDTH 1.4)
endif (NOT GLVIS_MS_LINEWIDTH)

#
# Start finding everything
#

set(_glvis_compile_defs)
set(_glvis_compile_opts)
set(_glvis_include_dirs)
set(_glvis_libraries)

if (GLVIS_USE_GLX10)
  list(APPEND _glvis_compile_defs "-DGLVIS_GLX10")
endif (GLVIS_USE_GLX10)

list(APPEND _glvis_compile_defs "-DGLVIS_MULTISAMPLE=${GLVIS_MULTISAMPLE}")
list(APPEND _glvis_compile_defs "-DGLVIS_MS_LINEWIDTH=${GLVIS_MS_LINEWIDTH}")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  list(APPEND _glvis_compile_defs "-DGLVIS_DEBUG")
endif()

# Base everything on MFEM
find_package(MFEM REQUIRED)
list(APPEND _glvis_include_dirs "${MFEM_INCLUDE_DIRS}")
list(APPEND _glvis_libraries "${MFEM_LIBRARIES}")

# Find X11. Only appear to need libX11.
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

    find_library(_glu_library GLU
      HINTS ${_x11_base_dir} ${OPENGL_DIR} $ENV{OPENGL_DIR}
      DOC "OpenGL GLU library."
      NO_DEFAULT_PATH)
    find_library(_glu_library GLU)

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
    list(APPEND _glvis_compile_defs "-DGLVIS_USE_LIBTIFF")
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
    list(APPEND _glvis_compile_defs "-DGLVIS_USE_PNG")
    list(APPEND _glvis_include_dirs "${PNG_INCLUDE_DIRS}")
    list(APPEND _glvis_libraries "${PNG_LIBRARIES}")
  else()
    message(WARNING "PNG library not found. PNG disabled.")
    set(GLVIS_USE_LIBPNG OFF)
  endif (PNG_FOUND)
endif (GLVIS_USE_LIBPNG)

# Find fontconfig and freetype.
if (GLVIS_USE_FREETYPE)

  find_package(Freetype)

  if (NOT FREETYPE_FOUND)
    message(WARNING "freetype not found. Cannot use freetype!")
  endif (NOT FREETYPE_FOUND)
  
  # Find FontConfig
  find_library(FONTCONFIG_LIBRARY fontconfig
    HINTS ${FONTCONFIG_DIR} $ENV{FONTCONFIG_DIR}
    DOC "The fontconfig library for use with Freetype."
    NO_DEFAULT_PATH)
  find_library(FONTCONFIG_LIBRARY fontconfig)

  if (FONTCONFIG_LIBRARY)
    message(STATUS "Found Fontconfig: ${FONTCONFIG_LIBRARY}")
  else()
    message(WARNING "fontconfig not found. Cannot use freetype!")
  endif (FONTCONFIG_LIBRARY)
  
  if (FREETYPE_FOUND AND FONTCONFIG_LIBRARY)

    list(APPEND _glvis_compile_defs "-DGLVIS_USE_FREETYPE")

    list(APPEND _glvis_include_dirs "${FREETYPE_INCLUDE_DIRS}")
    
    list(APPEND _glvis_library_dirs "${FREETYPE_LIBRARIES}")
    list(APPEND _glvis_library_dirs "${FONTCONFIG_LIBRARY}")

  else()
    set(GLVIS_USE_FREETYPE OFF)
  endif (FREETYPE_FOUND AND FONTCONFIG_LIBRARY)
endif (GLVIS_USE_FREETYPE)

# FIXME: MFEM exports the MPI include directory in its
# MFEM_INCLUDE_DIRS variable. However, I don't see it exporting the
# MPI libraries, so we need to go find them. For now, I use the
# include path to verify that the MPI library matches.
if (MFEM_USE_MPI)
  find_package(MPI REQUIRED)
  find_package(LAPACK REQUIRED)
  
  list(FIND MFEM_INCLUDE_DIRS ${MPI_INCLUDE_PATH} _found_correct_mpi)
  if (_found_correct_mpi EQUAL -1)
    message(FATAL_ERROR
      "Found a different version of MPI than used to build MFEM.")
  endif (_found_correct_mpi EQUAL -1)

  # Sometimes these cause problems...
  set(_mpi_compile_flags)
  set(_mpi_link_flags)
  if (UNIX)
    separate_arguments(
      _mpi_compile_flags UNIX_COMMAND "${MPI_CXX_COMPILE_FLAGS}")
    separate_arguments(
      _mpi_link_flags UNIX_COMMAND "${MPI_CXX_LINK_FLAGS}")
  else()
    separate_arguments(
      _mpi_compile_flags WINDOWS_COMMAND "${MPI_CXX_COMPILE_FLAGS}")
    separate_arguments(
      _mpi_link_flags WINDOWS_COMMAND "${MPI_CXX_LINK_FLAGS}")
  endif (UNIX)
  
  list(APPEND _glvis_compiler_opts "${_mpi_compile_flags}")
  list(APPEND _glvis_libraries "${_mpi_link_flags}")
  list(APPEND _glvis_libraries "${MPI_CXX_LIBRARIES}")
  list(APPEND _glvis_include_dirs "${MPI_CXX_INCLUDE_PATH}")

  list(APPEND _glvis_libraries "${LAPACK_LINKER_FLAGS}")
  list(APPEND _glvis_libraries "${LAPACK_LIBRARIES}")
endif (MFEM_USE_MPI)

# Find threading library
if (GLVIS_USE_THREADS)
  find_package(Threads)
  if (NOT THREADS_FOUND)
    message(WARNING "No threading model found. Threads disabled.")
    set(GLVIS_USE_THREADS OFF)
  else()
    list(APPEND _glvis_libraries "${CMAKE_THREAD_LIBS_INIT}")
  endif (NOT THREADS_FOUND)
endif (GLVIS_USE_THREADS)


#
# Setup the GLVis library target
#

add_subdirectory(lib)

#
# Setup the GLVis executable

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

