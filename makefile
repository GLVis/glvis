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

define GLVIS_HELP_MSG

GLVis makefile targets:

   make
   make status/info
   make install
   make app
   make clean
   make distclean
   make style

Examples:

make -j 4
   Build GLVis using the current configuration options from MFEM.
   (GLVis requires the MFEM finite element library, and uses its compiler and
    linker options in its build process.)
make status
   Display information about the current configuration.
make install PREFIX=<dir>
   Install the glvis executable in <dir>.
make app
   Build a Mac OS application bundle.
make clean
   Clean the glvis executable, library and object files.
make distclean
   In addition to "make clean", remove the local installation directory and some
   run-time generated files.
make style
   Format the GLVis C++ source files using the Artistic Style (astyle) settings
   from MFEM.

endef

# Custom configuration flags
GLVIS_CONFIG_MK ?=
-include $(GLVIS_CONFIG_MK)

# Default installation location
PREFIX ?= ./bin
INSTALL ?= /usr/bin/install

# Archiver: AR and ARFLAGS are defined by default, RANLIB is not.
# The default value of AR is 'ar' and we do not want to change that.
# The default value of ARFLAGS is 'rv', however, we want to set a different
# default, so we modify ARFLAGS, unless it was already changed on the command
# line or in the configuration file $(GLVIS_CONFIG_MK).
ifeq ($(origin ARFLAGS),default)
   ARFLAGS = cruv
endif
RANLIB ?= ranlib

# Use the MFEM build directory
MFEM_DIR ?= ../mfem
CONFIG_MK ?= $(MFEM_DIR)/config/config.mk
# Use the MFEM install directory
# MFEM_DIR = ../mfem/mfem
# CONFIG_MK = $(MFEM_DIR)/config.mk

# Use the compiler used by MFEM. Get the compiler and the options for compiling
# and linking from MFEM's config.mk. (Skip this if the target does not require
# building.)
ifeq (,$(filter help clean distclean style,$(MAKECMDGOALS)))
   -include $(CONFIG_MK)
endif

# GLVis requires c++11 which is also required by MFEM version >= 4.0
CXX = $(MFEM_CXX)
CPPFLAGS = $(MFEM_CPPFLAGS)
CXXFLAGS = $(MFEM_CXXFLAGS)

# MFEM config does not define C compiler
CC     ?= gcc
CFLAGS ?= -O3

# GLVis uses xxd to convert the logo  to a compilable file
XXD_FOUND := $(shell command -v xxd 2> /dev/null)

# Optional compile/link flags
GLVIS_OPTS ?=
GLVIS_LDFLAGS ?=

# emcc is used when building the wasm/js version
EMCC      ?= emcc -std=c++11
FONT_FILE ?= OpenSans.ttf
EMCC_OPTS ?= -s USE_SDL=2 -s USE_FREETYPE=1 -s USE_LIBPNG=1
# NOTE: we don't want to have DISABLE_DEPRECATED_FIND_EVENT_TARGET_BEHAVIOR=0
# longterm but until the SDL layer supports non-default canvas ids we need this,
# which will probably be a while.
# NOTE: we don't want to have --minify 0 longterm but we need to patch
# `_JSEvents_requestFullscreen' to be a noop for the live page and that is
# not easy to do in a minified build.
EMCC_LIBS ?= -s USE_SDL=2 --bind -s ALLOW_MEMORY_GROWTH=1 -s SINGLE_FILE=1 \
 --no-heap-copy -s ENVIRONMENT=web -s MODULARIZE=1 -s EXPORT_NAME=glvis \
 -s GL_ASSERTIONS=1 -s GL_DEBUG=1 -s USE_FREETYPE=1 -s MAX_WEBGL_VERSION=2 \
 -s DISABLE_DEPRECATED_FIND_EVENT_TARGET_BEHAVIOR=0 --minify 0

# Flags used when $(GLVIS_DEBUG) is not the same as $(MFEM_DEBUG)
CXX11FLAG ?= -std=c++11
OPTIM_OPTS ?= $(CXX11FLAG) -O3
DEBUG_OPTS ?= $(CXX11FLAG) -g -Wall
GLVIS_DEBUG ?= $(MFEM_DEBUG)
ifneq ($(GLVIS_DEBUG),$(MFEM_DEBUG))
   ifeq ($(GLVIS_DEBUG),YES)
      CXXFLAGS = $(DEBUG_OPTS)
   else
      CXXFLAGS = $(OPTIM_OPTS)
   endif
endif

GLVIS_FLAGS = $(CPPFLAGS) $(CXXFLAGS) $(MFEM_INCFLAGS) $(GLVIS_OPTS)
GLVIS_LIBS = $(MFEM_LIBS)

ifeq ($(GLVIS_DEBUG),YES)
   GLVIS_FLAGS += -DGLVIS_DEBUG
endif

NOTMAC := $(subst Darwin,,$(shell uname -s))

# Default multisampling mode and multisampling line-width
GLVIS_MULTISAMPLE  ?= 4
GLVIS_MS_LINEWIDTH ?= $(if $(NOTMAC),1.4,1.0)
GLVIS_FONT_SIZE ?= 12
DEFINES = -DGLVIS_MULTISAMPLE=$(GLVIS_MULTISAMPLE)\
 -DGLVIS_MS_LINEWIDTH=$(GLVIS_MS_LINEWIDTH)\
 -DGLVIS_FONT_SIZE=$(GLVIS_FONT_SIZE)\
 -DGLVIS_OGL3

# Enable logo setting via SDL (disabled on Windows/CMake build)
GLVIS_USE_LOGO ?= YES
ifeq ($(GLVIS_USE_LOGO),YES)
  DEFINES += -DGLVIS_USE_LOGO
endif

GLVIS_FLAGS += $(DEFINES)

# We don't want most of the stuff below because Emscripten handles that for us
EMCC_OPTS += $(CPPFLAGS) $(CXXFLAGS) $(MFEM_INCFLAGS) $(DEFINES)
EMCC_LIBS += $(MFEM_LIBS)

# Macro that searches for a file in a list of directories returning the first
# directory that contains the file.
# $(1) - the file to search for
# $(2) - list of directories to search
define find_dir
$(patsubst %/$(1),%,$(firstword $(wildcard $(foreach d,$(2),$(d)/$(1)))))
endef

# Macro to find the proper library sub-directory, 'lib64' or 'lib', given a
# tentative prefix and a library name. Returns empty path if prefix is empty,
# '/usr', or the library is not found.
# $(1) - the prefix to search, e.g. $(SDL_DIR)
# $(2) - library name without 'lib' prefix, e.g. 'SDL2'
define dir2lib
$(if $(filter-out /usr,$(1)),$(patsubst %/,%,$(dir $(firstword $(wildcard\
 $(1)/lib64/lib$(2).* $(1)/lib/lib$(2).*)))))
endef

BREW_PREFIX := $(if $(NOTMAC),,$(shell brew --prefix 2> /dev/null))

FREETYPE_SEARCH_PATHS = $(BREW_PREFIX) /usr /opt/X11
FREETYPE_SEARCH_FILE = include/freetype2/ft2build.h
FREETYPE_DIR = $(call find_dir,$(FREETYPE_SEARCH_FILE),$(FREETYPE_SEARCH_PATHS))
FREETYPE_LIB_DIR = $(call dir2lib,$(FREETYPE_DIR),freetype)
FREETYPE_LIBS = -lfreetype -lfontconfig

# If GLEW is in /usr, there's no need to add search paths
GLEW_SEARCH_PATHS = /usr/local $(BREW_PREFIX) $(abspath ../glew)
GLEW_SEARCH_FILE = include/GL/glew.h
GLEW_DIR ?= $(call find_dir,$(GLEW_SEARCH_FILE),$(GLEW_SEARCH_PATHS))
GLEW_LIB_DIR = $(call dir2lib,$(GLEW_DIR),GLEW)
GLEW_LIBS = -lGLEW

# If SDL is in /usr, there's no need to add search paths
SDL_SEARCH_PATHS := /usr/local $(BREW_PREFIX) $(abspath ../SDL2)
SDL_SEARCH_FILE = include/SDL2/SDL.h
SDL_DIR ?= $(call find_dir,$(SDL_SEARCH_FILE),$(SDL_SEARCH_PATHS))
SDL_LIB_DIR = $(call dir2lib,$(SDL_DIR),SDL2)
SDL_LIBS = -lSDL2

# If GLM is in /usr/include, there's no need to add search paths
GLM_SEARCH_PATHS = /usr/local/include \
 $(if $(BREW_PREFIX),$(BREW_PREFIX)/include) $(abspath ../glm)
GLM_SEARCH_FILE = glm/glm.hpp
GLM_DIR ?= $(call find_dir,$(GLM_SEARCH_FILE),$(GLM_SEARCH_PATHS))

# If OpenGL is in /usr, there's no need to add search paths
OPENGL_SEARCH_PATHS = /usr/local /opt/local
OPENGL_SEARCH_FILE = include/GL/gl.h
OPENGL_DIR ?= $(call find_dir,$(OPENGL_SEARCH_FILE),$(OPENGL_SEARCH_PATHS))
OPENGL_LIB_DIR = $(if $(NOTMAC),$(call dir2lib,$(OPENGL_DIR),GL))
OPENGL_LIBS = $(if $(NOTMAC),-lGL,-framework OpenGL -framework Cocoa)


# Regarding -DGLEW_NO_GLU, see https://github.com/nigels-com/glew/issues/192
GL_OPTS ?= $(if $(FREETYPE_DIR),-I$(FREETYPE_DIR)/include/freetype2) \
 $(if $(SDL_DIR),-I$(SDL_DIR)/include) \
 $(if $(GLEW_DIR),-I$(GLEW_DIR)/include) -DGLEW_NO_GLU \
 $(if $(GLM_DIR),-I$(GLM_DIR)) \
 $(if $(OPENGL_DIR),-I$(OPENGL_DIR)/include)

rpath=-Wl,-rpath,
GL_LIBS ?= $(if $(FREETYPE_LIB_DIR),-L$(FREETYPE_LIB_DIR)) \
 $(if $(SDL_LIB_DIR),-L$(SDL_LIB_DIR) $(rpath)$(SDL_LIB_DIR)) \
 $(if $(NOTMAC),$(if $(OPENGL_LIB_DIR),-L$(OPENGL_LIB_DIR) \
   $(rpath)$(OPENGL_LIB_DIR))) \
 $(if $(GLEW_LIB_DIR),-L$(GLEW_LIB_DIR) $(rpath)$(GLEW_LIB_DIR)) \
 $(FREETYPE_LIBS) $(SDL_LIBS) $(GLEW_LIBS) $(OPENGL_LIBS)

EMCC_OPTS += $(if $(GLM_DIR),-I$(GLM_DIR))

GLVIS_FLAGS += $(GL_OPTS)
GLVIS_LIBS  += $(GL_LIBS)

# Take screenshots internally with libtiff, libpng, or sdl2?
GLVIS_USE_LIBTIFF ?= NO
GLVIS_USE_LIBPNG  ?= YES
TIFF_OPTS = -DGLVIS_USE_LIBTIFF -I/sw/include
TIFF_LIBS = -L/sw/lib -ltiff
PNG_OPTS = -DGLVIS_USE_LIBPNG
PNG_LIBS = -lpng
ifeq ($(GLVIS_USE_LIBTIFF),YES)
   GLVIS_FLAGS += $(TIFF_OPTS)
   GLVIS_LIBS  += $(TIFF_LIBS)
else ifeq ($(GLVIS_USE_LIBPNG),YES)
   GLVIS_FLAGS += $(PNG_OPTS)
   GLVIS_LIBS  += $(PNG_LIBS)
else
   # no flag --> SDL screenshots
endif

PTHREAD_LIB = -lpthread
GLVIS_LIBS += $(PTHREAD_LIB)

LIBS = $(strip $(GLVIS_LIBS) $(GLVIS_LDFLAGS))
CCC  = $(strip $(CXX) $(GLVIS_FLAGS))
Ccc  = $(strip $(CC) $(CFLAGS) $(GL_OPTS))

# generated with 'echo lib/gl/*.c* lib/*.c*', does not include lib/*.m (Obj-C)
ALL_SOURCE_FILES = \
 lib/gl/renderer.cpp lib/gl/renderer_core.cpp lib/gl/renderer_ff.cpp    \
 lib/gl/shader.cpp lib/gl/types.cpp lib/aux_js.cpp lib/aux_vis.cpp      \
 lib/coll_reader.cpp lib/data_state.cpp lib/file_reader.cpp             \
 lib/font.cpp lib/gl2ps.c lib/gltf.cpp lib/material.cpp                 \
 lib/openglvis.cpp lib/palettes_base.cpp lib/palettes.cpp lib/sdl.cpp   \
 lib/script_controller.cpp lib/sdl_helper.cpp lib/sdl_main.cpp          \
 lib/sdl_windows.cpp lib/sdl_x11.cpp lib/stream_reader.cpp              \
 lib/threads.cpp lib/vsdata.cpp lib/vssolution.cpp lib/vssolution3d.cpp \
 lib/vsvector.cpp lib/vsvector3d.cpp lib/window.cpp
OBJC_SOURCE_FILES = $(if $(NOTMAC),,lib/sdl_mac.mm)
DESKTOP_ONLY_SOURCE_FILES = \
 lib/gl/renderer_ff.cpp lib/gl2ps.c lib/script_controller.cpp           \
 lib/sdl_x11.cpp
WEB_ONLY_SOURCE_FILES = lib/aux_js.cpp
LOGO_FILE = share/logo.rgba
LOGO_FILE_CPP = $(LOGO_FILE).bin.cpp
COMMON_SOURCE_FILES = $(filter-out \
 $(DESKTOP_ONLY_SOURCE_FILES) $(WEB_ONLY_SOURCE_FILES),$(ALL_SOURCE_FILES))

# generated with 'echo lib/gl/*.h* lib/*.h*'
HEADER_FILES = \
 lib/gl/attr_traits.hpp lib/gl/platform_gl.hpp lib/gl/renderer.hpp      \
 lib/gl/shader.hpp lib/gl/renderer_core.hpp lib/gl/renderer_ff.hpp      \
 lib/gl/types.hpp lib/aux_vis.hpp lib/coll_reader.hpp                   \
 lib/data_state.hpp lib/file_reader.hpp lib/font.hpp                    \
 lib/geom_utils.hpp lib/gl2ps.h lib/gltf.hpp lib/logo.hpp               \
 lib/material.hpp lib/openglvis.hpp lib/palettes_base.hpp               \
 lib/palettes.hpp lib/script_controller.hpp lib/sdl_helper.hpp          \
 lib/sdl.hpp lib/sdl_mac.hpp lib/sdl_main.hpp lib/sdl_windows.hpp       \
 lib/sdl_x11.hpp lib/stream_reader.hpp lib/threads.hpp lib/visual.hpp   \
 lib/vsdata.hpp lib/vssolution.hpp lib/vssolution3d.hpp                 \
 lib/vsvector.hpp lib/vsvector3d.hpp lib/window.hpp

DESKTOP_SOURCE_FILES = $(COMMON_SOURCE_FILES) $(DESKTOP_ONLY_SOURCE_FILES) $(LOGO_FILE_CPP)
WEB_SOURCE_FILES     = $(COMMON_SOURCE_FILES) $(WEB_ONLY_SOURCE_FILES)
OBJECT_FILES1        = $(DESKTOP_SOURCE_FILES:.cpp=.o)
OBJECT_FILES         = $(OBJECT_FILES1:.c=.o) $(OBJC_SOURCE_FILES:.mm=.o)
BYTECODE_FILES       = $(WEB_SOURCE_FILES:.cpp=.bc)

# Targets
.PHONY: clean distclean install status info opt debug style js

%.o: %.cpp
	$(CCC) -o $@ -c $<

%.o: %.c %.h
	$(Ccc) -o $@ -c $<

%.o: %.mm
	$(CCC) -o $@ -c $<

%.bc: %.cpp
	$(EMCC) $(EMCC_OPTS) $(GLVIS_FLAGS) -c $< -o $@

glvis:	glvis.cpp lib/libglvis.a $(CONFIG_MK) $(MFEM_LIB_FILE)
	$(CCC) -o glvis glvis.cpp -Llib -lglvis $(LIBS)

$(LOGO_FILE_CPP): $(LOGO_FILE)
ifndef XXD_FOUND
	$(error Required tool not found: xxd)
endif
	cd $(dir $(LOGO_FILE)) && xxd -i $(notdir $(LOGO_FILE)) > \
		$(notdir $(LOGO_FILE_CPP))

# Generate an error message if the MFEM library is not built and exit
$(CONFIG_MK) $(MFEM_LIB_FILE):
ifeq (,$(and $(findstring B,$(MAKEFLAGS)),$(wildcard $(CONFIG_MK))))
	$(error The MFEM library is not built)
endif

opt:
	$(MAKE) "GLVIS_DEBUG=NO"

debug:
	$(MAKE) "GLVIS_DEBUG=YES"

$(OBJECT_FILES): $(HEADER_FILES) $(CONFIG_MK)

lib/libglvis.a: $(OBJECT_FILES)
	$(AR) $(ARFLAGS) $@ $^; $(RANLIB) $@

js: lib/libglvis.js
lib/libglvis.js: $(BYTECODE_FILES) $(CONFIG_MK) $(MFEM_LIB_FILE)
	$(EMCC) $(EMCC_OPTS) -o $@ $(BYTECODE_FILES) $(EMCC_LIBS) --embed-file $(FONT_FILE)

clean:
	rm -rf lib/*.o lib/*.bc lib/gl/*.o lib/gl/*.bc lib/*~ *~ glvis
	rm -rf $(LOGO_FILE_CPP) share/*.o
	rm -rf lib/libglvis.a lib/libglvis.js *.dSYM
	rm -rf GLVis.app

distclean: clean
	rm -rf bin/
	rm -f GLVis_coloring.gf

install: glvis
	mkdir -p $(PREFIX)
	$(INSTALL) -m 750 glvis $(PREFIX)
ifeq ($(MFEM_USE_GNUTLS),YES)
	$(INSTALL) -m 750 glvis-keygen.sh $(PREFIX)
endif

app: glvis
	mkdir -p GLVis.app/Contents/MacOS
	mkdir -p GLVis.app/Contents/Resources
	cp share/Info.plist GLVis.app/Contents
	cp glvis GLVis.app/Contents/MacOS
	cp share/GLVis.icns GLVis.app/Contents/Resources
	cp share/Credits.rtf GLVis.app/Contents/Resources

help:
	$(info $(value GLVIS_HELP_MSG))
	@true

status info:
	$(info MFEM_DIR    = $(MFEM_DIR))
	$(info GLVIS_FLAGS = $(GLVIS_FLAGS))
	$(info GLVIS_LIBS  = $(value GLVIS_LIBS))
	$(info GLVIS_LIBS  = $(GLVIS_LIBS))
	$(info PREFIX      = $(PREFIX))
	@true

# Print the contents of a makefile variable, e.g.: 'make print-MFEM_LIBS'.
print-%:
	$(info [ variable name]: $*)
	$(info [        origin]: $(origin $*))
	$(info [         value]: $(value $*))
	$(info [expanded value]: $($*))
	$(info )
	@true

ASTYLE_BIN ?= astyle
ASTYLE = $(ASTYLE_BIN) --options=$(MFEM_DIR)/config/mfem.astylerc
ALL_FILES = ./glvis.cpp $(ALL_SOURCE_FILES) $(HEADER_FILES)
EXT_FILES = lib/gl2ps.c lib/gl2ps.h
FORMAT_FILES := $(filter-out $(EXT_FILES), $(ALL_FILES))

style:
	@if ! $(ASTYLE) $(FORMAT_FILES) | grep Formatted; then\
	   echo "No source files were changed.";\
	fi
