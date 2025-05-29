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

#include "window.hpp"
#include "visual.hpp"

Window &Window::operator=(Window &&w)
{
   internal = std::move(w.internal);

   data_state = std::move(w.data_state);

   window_x = w.window_x;
   window_y = w.window_y;
   window_w = w.window_w;
   window_h = w.window_h;
   window_title = w.window_title;
   headless = w.headless;
   plot_caption = std::move(w.plot_caption);
   extra_caption = std::move(w.extra_caption);

   return *this;
}

// Visualize the data in the global variables mesh, sol/grid_f, etc
bool Window::GLVisInitVis(StreamCollection input_streams)
{
   DataState::FieldType field_type = data_state.GetType();

   if (field_type <= DataState::FieldType::MIN
       || field_type >= DataState::FieldType::MAX)
   {
      return false;
   }

   static const char *window_titles[] = { "GLVis [mesh]",
                                          "GLVis [scalar data]",
                                          "GLVis [vector data]",
                                        };

   const char *win_title = (window_title == nullptr) ?
                           window_titles[(int)field_type] : window_title;

   GLWindow *new_wnd = InitVisualization(win_title, window_x, window_y, window_w,
                                         window_h, headless);
   if (new_wnd != wnd.get()) { internal.wnd.reset(new_wnd); }
   if (!wnd)
   {
      std::cerr << "Initializing the visualization failed." << std::endl;
      return false;
   }

   if (input_streams.size() > 0)
   {
      if (!headless)
      {
         wnd->setOnKeyDown(SDLK_SPACE, ThreadsPauseFunc);
      }
      internal.glvis_command.reset(new GLVisCommand(*this));
      SetGLVisCommand(glvis_command.get());
      internal.comm_thread.reset(new communication_thread(std::move(input_streams),
                                                          glvis_command.get(), headless));
   }

   locwin = this;

   if (data_state.quad_f)
   {
      wnd->setOnKeyDown('Q', SwitchQuadSolution);
   }

   double mesh_range = -1.0;
   if (field_type == DataState::FieldType::SCALAR
       || field_type == DataState::FieldType::MESH)
   {
      if (data_state.mesh->SpaceDimension() == 2)
      {
         internal.vs.reset(new VisualizationSceneSolution(*this));

         if (field_type == DataState::FieldType::MESH)
         {
            vs->OrthogonalProjection = 1;
            vs->SetLight(false);
            vs->Zoom(1.8);
            // Use the 'bone' palette when visualizing a 2D mesh only (otherwise
            // the 'jet-like' palette is used in 2D, see vssolution.cpp).
            vs->palette.SetFallbackIndex(4);
         }
      }
      else if (data_state.mesh->SpaceDimension() == 3)
      {
         VisualizationSceneSolution3d *vss;
         vss = new VisualizationSceneSolution3d(*this);
         internal.vs.reset(vss);

         if (field_type == DataState::FieldType::MESH)
         {
            if (data_state.mesh->Dimension() == 3)
            {
               // Use the 'white' palette when visualizing a 3D volume mesh only
               vss->palette.SetFallbackIndex(11);
               vss->SetLightMatIdx(4);
            }
            else
            {
               // Use the 'bone' palette when visualizing a surface mesh only
               vss->palette.SetFallbackIndex(4);
            }
            // Otherwise, the 'vivid' palette is used in 3D see vssolution3d.cpp
            vss->ToggleDrawAxes();
            vss->ToggleDrawMesh();
         }
      }
      if (field_type == DataState::FieldType::MESH)
      {
         if (data_state.grid_f)
         {
            mesh_range = data_state.grid_f->Max() + 1.0;
         }
         else
         {
            mesh_range = data_state.sol->Max() + 1.0;
         }
      }
   }
   else if (field_type == DataState::FieldType::VECTOR)
   {
      if (data_state.mesh->SpaceDimension() == 2)
      {
         internal.vs.reset(new VisualizationSceneVector(*this));
      }
      else if (data_state.mesh->SpaceDimension() == 3)
      {
         if (data_state.grid_f)
         {
            data_state.ProjectVectorFEGridFunction();
         }
         internal.vs.reset(new VisualizationSceneVector3d(*this));
      }
   }

   if (vs)
   {
      // increase the refinement factors if visualizing a GridFunction
      if (data_state.grid_f)
      {
         vs->AutoRefine();
         vs->SetShading(VisualizationSceneScalarData::Shading::Noncomforming, true);
      }
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
         vs->SetAutoscale(0);
      }
      if (data_state.mesh->SpaceDimension() == 2
          && field_type == DataState::FieldType::MESH)
      {
         SetVisualizationScene(vs.get(), 2, data_state.keys.c_str());
      }
      else
      {
         SetVisualizationScene(vs.get(), 3, data_state.keys.c_str());
      }
   }
   return true;
}

void Window::GLVisStartVis()
{
   RunVisualization();
   internal.vs.reset();
   internal.wnd.reset();
   if (glvis_command)
   {
      glvis_command->Terminate();
      internal.comm_thread.reset();
      internal.glvis_command.reset();
   }
   std::cout << "GLVis window closed." << std::endl;
}

void Window::SwitchQuadSolution(DataState::QuadSolution quad_type)
{
   data_state.SwitchQuadSolution(quad_type);
   ResetMeshAndSolution(data_state);
}

bool Window::SetNewMeshAndSolution(DataState new_state)
{
   if (new_state.mesh->SpaceDimension() == data_state.mesh->SpaceDimension() &&
       new_state.GetType() == data_state.GetType() &&
       (((new_state.grid_f && data_state.grid_f) &&
         (new_state.grid_f->VectorDim() == data_state.grid_f->VectorDim()))
        ||(!new_state.grid_f && !data_state.grid_f)))
   {
      ResetMeshAndSolution(new_state);

      data_state = std::move(new_state);

      return true;
   }
   else
   {
      return false;
   }
}

void Window::ResetMeshAndSolution(DataState &ss)
{
   if (ss.mesh->SpaceDimension() == 3 &&
       ss.GetType() == DataState::FieldType::VECTOR)
   {
      ss.ProjectVectorFEGridFunction();
   }

   vs->NewMeshAndSolution(ss);
}


thread_local Window *Window::locwin = NULL;

void Window::SwitchQuadSolution()
{
   int iqs = ((int)locwin->data_state.GetQuadSolution()+1)
             % ((int)DataState::QuadSolution::MAX);
   locwin->SwitchQuadSolution((DataState::QuadSolution)iqs);
   SendExposeEvent();
}
