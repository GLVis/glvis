// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_WINDOW_HPP
#define GLVIS_WINDOW_HPP

#include <string>

#include "data_state.hpp"
#include "stream_reader.hpp"

class VisualizationSceneScalarData;
class communication_thread;

extern const char *string_none;
extern const char *string_default;

struct Window
{
   DataState data_state;
   VisualizationSceneScalarData *vs = NULL;
   communication_thread *comm_thread = NULL;

   int         window_x        = 0; // not a command line option
   int         window_y        = 0; // not a command line option
   int         window_w        = 400;
   int         window_h        = 350;
   const char *window_title    = string_default;
   std::string plot_caption;
   std::string extra_caption;

   /// Visualize the data in the global variables mesh, sol/grid_f, etc
   bool GLVisInitVis(StreamCollection input_streams);
   void GLVisStartVis();

private:
   /// Thread-local singleton for key handlers
   static thread_local Window *locwin;

   /// Switch representation of the quadrature function
   static void SwitchQuadSolution();
};

#endif // GLVIS_WINDOW_HPP
