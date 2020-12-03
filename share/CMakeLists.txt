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

add_custom_command(
    OUTPUT logo.rgba.bin.cpp
    COMMAND ${CMAKE_COMMAND} -E copy
	    ${CMAKE_CURRENT_SOURCE_DIR}/logo.rgba
	    ${CMAKE_CURRENT_BINARY_DIR}/logo.rgba
    COMMAND xxd -i logo.rgba > logo.rgba.bin.cpp
    DEPENDS logo.rgba
    COMMENT "Generating logo.rgba.bin.o")

add_library(glvis_logo OBJECT
    ${CMAKE_CURRENT_BINARY_DIR}/logo.rgba.bin.cpp)
