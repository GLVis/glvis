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

class SdlWindow;
class VisualizationSceneScalarData;
class communication_thread;
class GLVisCommand;

struct Window
{
private:
   struct
   {
      std::unique_ptr<SdlWindow> wnd;
      std::unique_ptr<VisualizationSceneScalarData> vs;
      std::unique_ptr<communication_thread> comm_thread;
      std::unique_ptr<GLVisCommand> glvis_command;
   } internal;

public:
   DataState data_state;
   const std::unique_ptr<SdlWindow> &wnd{internal.wnd};
   const std::unique_ptr<VisualizationSceneScalarData> &vs{internal.vs};
   const std::unique_ptr<communication_thread> &comm_thread {internal.comm_thread};
   const std::unique_ptr<GLVisCommand> &glvis_command{internal.glvis_command};

   int         window_x        = 0; // not a command line option
   int         window_y        = 0; // not a command line option
   int         window_w        = 400;
   int         window_h        = 350;
   const char *window_title    = nullptr;
   std::string plot_caption;
   std::string extra_caption;

   Window() = default;
   Window(Window &&w) { *this = std::move(w); }
   Window& operator=(Window &&w);

   /// Visualize the data in the global variables mesh, sol/grid_f, etc
   bool GLVisInitVis(StreamCollection input_streams);
   void GLVisStartVis();

   /// Switch the quadrature function representation and update the visualization
   void SwitchQuadSolution(DataState::QuadSolution type);

   /// Sets a new mesh and solution from another DataState object, and
   /// updates the VisualizationScene with the new data.
   ///
   /// Mesh space and grid function dimensions must both match the original
   /// dimensions of the current DataState. If there is a mismatch in either
   /// value, the function will return false, and the mesh/solution will not be
   /// updated.
   bool SetNewMeshAndSolution(DataState new_state);

   /// Updates the VisualizationScene with the new data of the given DataState object.
   /// @note: Use with caution when the update is compatible
   /// @see SetNewMeshAndSolution()
   void ResetMeshAndSolution(DataState &ss);

private:
   /// Thread-local singleton for key handlers
   static thread_local Window *locwin;

   /// Switch representation of the quadrature function
   static void SwitchQuadSolution();
};

#endif // GLVIS_WINDOW_HPP
