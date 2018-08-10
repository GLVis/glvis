# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see http://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

define GLVIS_HELP_MSG

GLVis makefile targets:

   make
   make status/info
   make install
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
make clean
   Clean the glvis executable, library and object files.
make distclean
   In addition to "make clean", remove the local installation directory and some
   run-time generated files.
make style
   Format the GLVis C++ source files using the Artistic Style (astyle) settings
   from MFEM.

endef

# Default installation location
PREFIX = ./bin
INSTALL = /usr/bin/install

# Use the MFEM build directory
MFEM_DIR = ../mfem
CONFIG_MK = $(MFEM_DIR)/config/config.mk
# Use the MFEM install directory
# MFEM_DIR = ../mfem/mfem
# CONFIG_MK = $(MFEM_DIR)/config.mk

# Use two relative paths to MFEM: first one for compilation in '.' and second
# one for compilation in 'lib'.
MFEM_DIR1 := $(MFEM_DIR)
MFEM_DIR2 := $(realpath $(MFEM_DIR))

# Use the compiler used by MFEM. Get the compiler and the options for compiling
# and linking from MFEM's config.mk. (Skip this if the target does not require
# building.)
ifeq (,$(filter help clean distclean style,$(MAKECMDGOALS)))
   -include $(CONFIG_MK)
endif

CXX = $(MFEM_CXX)
CPPFLAGS = $(MFEM_CPPFLAGS)
CXXFLAGS = $(MFEM_CXXFLAGS)

# MFEM config does not define C compiler
CC     = emcc 
CFLAGS = -O3 

# Optional link flags
LDFLAGS =

OPTIM_OPTS = -O3
DEBUG_OPTS = -g -Wall
GLVIS_DEBUG = $(MFEM_DEBUG)
ifneq ($(GLVIS_DEBUG),$(MFEM_DEBUG))
   ifeq ($(GLVIS_DEBUG),YES)
      CXXFLAGS = $(DEBUG_OPTS)
   else
      CXXFLAGS = $(OPTIM_OPTS)
   endif
endif

GLVIS_FLAGS = $(CPPFLAGS) $(CXXFLAGS) $(MFEM_INCFLAGS)
GLVIS_LIBS = $(MFEM_LIBS)

ifeq ($(GLVIS_DEBUG),YES)
   GLVIS_FLAGS += -DGLVIS_DEBUG
endif

NOTMAC := $(subst Darwin,,$(shell uname -s))
SO_EXT = $(if $(NOTMAC),so,dylib)

# Default multisampling mode and multisampling line-width
GLVIS_MULTISAMPLE  = 4
GLVIS_MS_LINEWIDTH = $(if $(NOTMAC),1.4,0.01)
GLVIS_FLAGS += -DGLVIS_MULTISAMPLE=$(GLVIS_MULTISAMPLE)\
 -DGLVIS_MS_LINEWIDTH=$(GLVIS_MS_LINEWIDTH)\
 -DGLVIS_OGL3

# Macro that searches for a file in a list of directories returning the first
# directory that contains the file.
# $(1) - the file to search for
# $(2) - list of directories to search
define find_dir
$(patsubst %/$(1),%,$(firstword $(wildcard $(foreach d,$(2),$(d)/$(1)))))
endef

# The X11 and OpenGL libraries

BREW_PATH = /Users/yang39/usr/brew

X11_SEARCH_PATHS = /usr /usr/X11 /opt/X11 /usr/X11R6
X11_SEARCH_FILE = include/X11/Xlib.h
X11_DIR = $(call find_dir,$(X11_SEARCH_FILE),$(X11_SEARCH_PATHS))
X11_LIB_DIR = $(call find_dir,libX11.$(SO_EXT),$(X11_DIR)/lib64 $(X11_DIR)/lib)
GL_OPTS = -Ithirdparty/glm -s USE_SDL=2 -s USE_SDL_TTF=2
# for servers not supporting GLX 1.3:
# GL_OPTS = -I$(X11_DIR)/include -DGLVIS_GLX10

GL_LIBS = -L$(BREW_PATH)/lib -s USE_SDL=2 -lGLEW -lSDL2

GLVIS_FLAGS += $(GL_OPTS)
GLVIS_LIBS  += $(GL_LIBS)

# Take screenshots internally with libtiff, libpng, or externally with xwd?
USE_LIBTIFF = NO
USE_LIBPNG  = NO
TIFF_OPTS = -DGLVIS_USE_LIBTIFF -I/sw/include
TIFF_LIBS = -L/sw/lib -ltiff
PNG_OPTS = -DGLVIS_USE_LIBPNG
PNG_LIBS = -lpng
ifeq ($(USE_LIBTIFF),YES)
   GLVIS_FLAGS += $(TIFF_OPTS)
   GLVIS_LIBS  += $(TIFF_LIBS)
endif
ifeq ($(USE_LIBPNG),YES)
   GLVIS_FLAGS += $(PNG_OPTS)
   GLVIS_LIBS  += $(PNG_LIBS)
endif

# Render fonts using the freetype library and use the fontconfig library to
# find font files.
USE_FREETYPE = YES
# libfreetype + libfontconfig
# get cflags with: freetype-config --cflags  or  pkg-config freetype2 --cflags
# get libs with:   freetype-config --libs    or  pkg-config freetype2 --libs
# libfontconfig:   pkg-config fontconfig --cflags
#                  pkg-config fontconfig --libs
FT_OPTS = -DGLVIS_USE_FREETYPE -s USE_FREETYPE=1 #-I$(X11_DIR)/include/freetype2
FT_LIBS = -lfreetype -lfontconfig
ifeq ($(USE_FREETYPE),YES)
   GLVIS_FLAGS += $(FT_OPTS)
   GLVIS_LIBS  += $(FT_LIBS)
endif

PTHREAD_LIB = -lpthread
GLVIS_LIBS += $(PTHREAD_LIB)

LIBS = $(strip $(GLVIS_LIBS) $(LDFLAGS))
CCC  = $(strip $(CXX) $(GLVIS_FLAGS))
Ccc  = $(strip $(CC) $(CFLAGS) $(GL_OPTS))

# generated with 'echo lib/*.c*'
SOURCE_FILES = lib/aux_vis.cpp lib/aux_gl3.cpp lib/font.cpp lib/sdl.cpp \
 lib/material.cpp lib/openglvis.cpp lib/palettes.cpp lib/threads.cpp lib/vsdata.cpp \
 lib/vssolution.cpp lib/vssolution3d.cpp lib/vsvector.cpp lib/vsvector3d.cpp lib/glstate.cpp
OBJECT_FILES1 = $(SOURCE_FILES:.cpp=.o)
OBJECT_FILES = $(OBJECT_FILES1:.c=.o)
# generated with 'echo lib/*.h*'
HEADER_FILES = lib/aux_vis.hpp lib/aux_gl3.hpp lib/font.hpp lib/sdl.hpp lib/material.hpp \
 lib/openglvis.hpp lib/palettes.hpp lib/threads.hpp lib/visual.hpp \
 lib/vsdata.hpp lib/vssolution.hpp lib/vssolution3d.hpp lib/vsvector.hpp lib/vsvector3d.hpp lib/glstate.hpp

# Targets

.PHONY: clean distclean install status info opt debug style

.SUFFIXES: .c .cpp .o
.cpp.o:
	$(CCC) -g -c lib/$(<F) -o lib/$(*F).bc
.c.o:
	$(Ccc) -g -c lib/$(<F) -o lib/$(*F).bc

glvis: override MFEM_DIR = $(MFEM_DIR1)
glvis:	glvis.cpp lib/libglvis.a $(CONFIG_MK) $(MFEM_LIB_FILE)
	$(CCC) -o glvis glvis.cpp -Llib -lglvis $(LIBS)

# Generate an error message if the MFEM library is not built and exit
$(CONFIG_MK) $(MFEM_LIB_FILE):
ifeq (,$(and $(findstring B,$(MAKEFLAGS)),$(wildcard $(CONFIG_MK))))
	$(error The MFEM library is not built)
endif

opt:
	$(MAKE) "GLVIS_DEBUG=NO"

debug:
	$(MAKE) "GLVIS_DEBUG=YES"

$(OBJECT_FILES): override MFEM_DIR = $(MFEM_DIR2)
$(OBJECT_FILES): $(HEADER_FILES) $(CONFIG_MK)

lib/libglvis.a: $(OBJECT_FILES)
	cd lib;	ar cruv libglvis.a *.o;	ranlib libglvis.a

clean:
	rm -rf lib/*.o lib/*~ *~ glvis lib/libglvis.a *.dSYM

distclean: clean
	rm -rf bin/
	rm -f GLVis_coloring.gf

install: glvis
	mkdir -p $(PREFIX)
	$(INSTALL) -m 750 glvis $(PREFIX)
ifeq ($(MFEM_USE_GNUTLS),YES)
	$(INSTALL) -m 750 glvis-keygen.sh $(PREFIX)
endif

help:
	$(info $(value GLVIS_HELP_MSG))
	@true

status info:
	$(info MFEM_DIR    = $(MFEM_DIR))
	$(info GLVIS_FLAGS = $(GLVIS_FLAGS))
	$(info GLVIS_LIBS  = $(value GLVIS_LIBS))
	$(info PREFIX      = $(PREFIX))
	@true

ASTYLE = astyle --options=$(MFEM_DIR1)/config/mfem.astylerc
ALL_FILES = ./glvis.cpp $(SOURCE_FILES) $(HEADER_FILES)
EXT_FILES = lib/aux_gl.cpp lib/aux_gl.hpp lib/gl2ps.c lib/gl2ps.h \
  lib/tk.cpp lib/tk.h
FORMAT_FILES := $(filter-out $(EXT_FILES), $(ALL_FILES))

style:
	@if ! $(ASTYLE) $(FORMAT_FILES) | grep Formatted; then\
	   echo "No source files were changed.";\
	fi
