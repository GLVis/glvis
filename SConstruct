# Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see http://glvis.googlecode.com.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

import os
import sys

Help("""
       Type: 'scons' to build the production library,
             'scons -c' to clean the build,
             'scons debug=1' to build the debug version.
       """)

env = Environment()

CC_OPTS    = '-O3'
DEBUG_OPTS = '-g -DGLVIS_DEBUG -Wall'

# GLVis-specific options
env.Append(CPPDEFINES = ['GLVIS_MULTISAMPLE=4'])

# Debug options
debug = ARGUMENTS.get('debug', 0)
if int(debug):
   env.Append(CCFLAGS = DEBUG_OPTS)
else:
   env.Append(CCFLAGS = CC_OPTS)

# Mac-specific options
if (sys.platform == "darwin"):
   env.Append(LIBPATH = ["/sw/lib", "/usr/local/lib"])
   env.Append(CPPPATH = ["/sw/include", "/usr/local/include"])
   env.Append(LINKFLAGS = """-Wl,-dylib_file,/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib:/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries/libGL.dylib""")

conf = Configure(env)

# Check for LAPACK
if conf.CheckLib('lapack', 'dsyevr_'):
   env.Append(LIBS = ['blas'])
   print 'Using LAPACK'
else:
   print 'Did not find LAPACK, continuing without it'

# Check for libtiff
if conf.CheckLibWithHeader('tiff', 'tiff.h', 'c++'):
   env.Append(CPPDEFINES = ['GLVIS_USE_LIBTIFF'])
   print 'Using libtiff'
else:
   print 'Did not find libtiff, continuing without it'

env = conf.Finish()

env.Append(CPPPATH = ['../mfem','lib','/usr/X11R6/include'])
env.Append(LIBS = ['glvis','mfem','X11','GL','GLU'])
env.Append(LIBPATH = ['lib','../mfem',os.environ['HOME']+'/lib','/usr/X11R6/lib'])

# libglvis.a library
lib_src = [Glob('lib/*.cpp'),Glob('lib/*.c')]
env.Library('lib/glvis',lib_src)

# glvis binary
env.Program('glvis','glvis.cpp')
