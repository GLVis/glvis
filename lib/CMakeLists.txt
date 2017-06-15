list(APPEND SRCS
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

list(APPEND HDRS
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

convert_filenames_to_full_paths(SRCS)
convert_filenames_to_full_paths(HDRS)

set(SOURCES ${SOURCES} ${SRCS} PARENT_SCOPE)
set(HEADERS ${HEADERS} ${HDRS} PARENT_SCOPE)
