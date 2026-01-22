# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-443271.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see https://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

list(APPEND SOURCES
  gl/renderer.cpp
  gl/renderer_core.cpp
  gl/shader.cpp
  gl/types.cpp
  aux_vis.cpp
  coll_reader.cpp
  data_state.cpp
  file_reader.cpp
  font.cpp
  gltf.cpp
  material.cpp
  openglvis.cpp
  palettes.cpp
  palettes_default.cpp
  palettes_base.cpp
  sdl.cpp
  sdl_helper.cpp
  sdl_main.cpp
  stream_reader.cpp
  vsdata.cpp
  vssolution.cpp
  vssolution3d.cpp
  vsvector.cpp
  vsvector3d.cpp)

list(APPEND HEADERS
  gl/attr_traits.hpp
  gl/platform_gl.hpp
  gl/renderer.hpp
  gl/renderer_core.hpp
  gl/shader.hpp
  gl/types.hpp
  aux_vis.hpp
  coll_reader.hpp
  data_state.hpp
  file_reader.hpp
  font.hpp
  geom_utils.hpp
  gltf.hpp
  logo.hpp
  material.hpp
  openglvis.hpp
  palettes.hpp
  palettes_base.hpp
  sdl.hpp
  sdl_helper.hpp
  sdl_main.hpp
  sdl_mac.hpp
  stream_reader.hpp
  visual.hpp
  vsdata.hpp
  vssolution.hpp
  vssolution3d.hpp
  vsvector.hpp
  vsvector3d.hpp)

if(EMSCRIPTEN)
  # Emscripten build target
  list(APPEND SOURCES aux_js.cpp)
  add_executable(libglvis ${SOURCES} ${HEADERS})

  # Copy font to lib folder
  configure_file(../OpenSans.ttf OpenSans.ttf COPYONLY)

  # Set Emscripten options first, since "-s" emscripten options can't be passed
  # through target_compile_options/target_link_libraries
  string(REPLACE ";" " " _emscripten_opts "${_emscripten_opts}")
  set_target_properties(libglvis PROPERTIES COMPILE_FLAGS "${_emscripten_opts}")
  set_target_properties(libglvis PROPERTIES LINK_FLAGS "${_emscripten_opts}")

  target_compile_definitions(libglvis PUBLIC ${_glvis_compile_defs})
  target_include_directories(libglvis SYSTEM PUBLIC "${_glvis_include_dirs}")
  target_link_libraries(libglvis PUBLIC "${_glvis_libraries}")
else()
  # Desktop build target
  list(APPEND SOURCES gl/renderer_ff.cpp threads.cpp gl2ps.c sdl_x11.cpp)
  list(APPEND HEADERS gl/renderer_ff.hpp threads.hpp gl2ps.h sdl_x11.hpp)
  if (APPLE)
    list(APPEND SOURCES sdl_mac.mm)
    list(APPEND HEADERS sdl_mac.hpp)
  endif()
  if (WIN32)
    list(APPEND SOURCES sdl_windows.cpp)
    list(APPEND HEADERS sdl_windows.hpp)
  endif()
  add_library(glvis ${SOURCES} ${HEADERS})
  target_compile_definitions(glvis PUBLIC "${_glvis_compile_defs}")
  target_compile_options(glvis PUBLIC "${_glvis_compile_opts}")
  target_include_directories(glvis SYSTEM PUBLIC "${_glvis_include_dirs}")
  target_link_libraries(glvis PUBLIC "${_glvis_libraries}")
  if (APPLE)
    target_link_libraries(glvis PUBLIC "-framework Cocoa")
  endif()
endif(EMSCRIPTEN)
