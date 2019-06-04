list(APPEND SOURCES
  aux_vis.cpp
  font.cpp
  gl2ps.c
  material.cpp
  openglvis.cpp
  palettes.cpp
  gl/types.cpp
  gl/renderer.cpp
  gl/renderer_core.cpp
  gl/renderer_ff.cpp
  sdl.cpp 
  threads.cpp
  vsdata.cpp
  vssolution3d.cpp
  vssolution.cpp
  vsvector3d.cpp
  vsvector.cpp)

list(APPEND HEADERS
  aux_vis.hpp
  font.hpp
  gl2ps.h
  material.hpp
  openglvis.hpp
  palettes.hpp
  gl/platform_gl.hpp
  gl/types.hpp
  gl/renderer.hpp
  gl/renderer_core.hpp
  gl/renderer_ff.hpp
  sdl.hpp
  threads.hpp
  visual.hpp
  vsdata.hpp
  vssolution3d.hpp
  vssolution.hpp
  vsvector3d.hpp
  vsvector.hpp)

# Allegedly adding the headers is helpful for IDEs.
add_library(glvis ${SOURCES} ${HEADERS})
target_compile_definitions(glvis PUBLIC "${_glvis_compile_defs}")
target_compile_options(glvis PUBLIC "${_glvis_compile_opts}")
target_include_directories(glvis SYSTEM PUBLIC "${_glvis_include_dirs}")
target_link_libraries(glvis PUBLIC "${_glvis_libraries}")
