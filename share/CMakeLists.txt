# Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-443271.
#
# This file is part of the GLVis visualization tool and library. For more
# information and source code availability see https://glvis.org.
#
# GLVis is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

add_custom_command(
    OUTPUT logo.rgba.bin.cpp
    COMMAND ${CMAKE_COMMAND} -E copy
	    ${CMAKE_CURRENT_SOURCE_DIR}/logo.rgba
	    ${CMAKE_CURRENT_BINARY_DIR}/logo.rgba
    COMMAND xxd -i logo.rgba > logo.rgba.bin.cpp
    DEPENDS logo.rgba
    COMMENT "Generating logo.rgba.bin.o")

add_library(glvis_logo STATIC
    ${CMAKE_CURRENT_BINARY_DIR}/logo.rgba.bin.cpp)
