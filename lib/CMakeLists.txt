list(APPEND SOURCES
  gl/renderer.cpp
  gl/renderer_core.cpp
  gl/types.cpp
  aux_vis.cpp
  font.cpp
  material.cpp
  openglvis.cpp
  palettes.cpp
  sdl.cpp 
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
  gl/types.hpp
  aux_vis.hpp
  font.hpp
  material.hpp
  openglvis.hpp
  palettes.hpp
  sdl.hpp
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
  list(APPEND SOURCES gl/renderer_ff.cpp threads.cpp gl2ps.c)
  list(APPEND HEADERS gl/renderer_ff.hpp threads.hpp gl2ps.h)
  add_library(glvis ${SOURCES} ${HEADERS})
  target_compile_definitions(glvis PUBLIC "${_glvis_compile_defs}")
  target_compile_options(glvis PUBLIC "${_glvis_compile_opts}>")
  target_include_directories(glvis SYSTEM PUBLIC "${_glvis_include_dirs}")
  target_link_libraries(glvis PUBLIC "${_glvis_libraries}")
endif(EMSCRIPTEN)
