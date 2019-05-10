list(APPEND SOURCES
  aux_gl3.cpp
  aux_vis.cpp
  font.cpp
  gl2ps.c
  material.cpp
  openglvis.cpp
  palettes.cpp
  renderer.cpp
  renderer_core.cpp
  renderer_ff.cpp
  sdl.cpp 
  threads.cpp
  vsdata.cpp
  vssolution3d.cpp
  vssolution.cpp
  vsvector3d.cpp
  vsvector.cpp)

list(APPEND HEADERS
  aux_gl3.hpp
  aux_vis.hpp
  font.hpp
  gl2ps.h
  material.hpp
  openglvis.hpp
  palettes.hpp
  platform_gl.hpp
  renderer.hpp
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
