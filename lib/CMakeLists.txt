list(APPEND SOURCES
  aux_gl.cpp
  aux_vis.cpp
  gl2ps.c
  material.cpp
  openglvis.cpp
  threads.cpp
  tk.cpp
  vsdata.cpp
  vssolution3d.cpp
  vssolution.cpp
  vsvector3d.cpp
  vsvector.cpp)

list(APPEND HEADERS
  aux_gl.hpp
  aux_vis.hpp
  gl2ps.h
  material.hpp
  openglvis.hpp
  palettes.hpp
  threads.hpp
  tk.h
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
