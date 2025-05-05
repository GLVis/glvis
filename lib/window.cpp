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

extern thread_local GLVisCommand* glvis_command;

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

   const char *win_title = (window_title == string_default) ?
                           window_titles[(int)field_type] : window_title;

   if (InitVisualization(win_title, window_x, window_y, window_w, window_h))
   {
      cerr << "Initializing the visualization failed." << endl;
      return false;
   }

   if (input_streams.size() > 0)
   {
      GetAppWindow()->setOnKeyDown(SDLK_SPACE, ThreadsPauseFunc);
      glvis_command = new GLVisCommand(*this);
      comm_thread = new communication_thread(std::move(input_streams), glvis_command);
   }

   locwin = this;

   if (data_state.quad_f)
   {
      GetAppWindow()->setOnKeyDown('Q', SwitchQuadSolution);
   }

   double mesh_range = -1.0;
   if (field_type == DataState::FieldType::SCALAR
       || field_type == DataState::FieldType::MESH)
   {
      if (data_state.grid_f)
      {
         data_state.grid_f->GetNodalValues(data_state.sol);
      }
      if (data_state.mesh->SpaceDimension() == 2)
      {
         vs = new VisualizationSceneSolution(*this);

         if (field_type == DataState::FieldType::MESH)
         {
            vs->OrthogonalProjection = 1;
            vs->SetLight(false);
            vs->Zoom(1.8);
            // Use the 'bone' palette when visualizing a 2D mesh only (otherwise
            // the 'jet-like' palette is used in 2D, see vssolution.cpp).
            vs->palette.SetIndex(4);
         }
      }
      else if (data_state.mesh->SpaceDimension() == 3)
      {
         VisualizationSceneSolution3d *vss;
         vs = vss = new VisualizationSceneSolution3d(*this);

         if (field_type == DataState::FieldType::MESH)
         {
            if (data_state.mesh->Dimension() == 3)
            {
               // Use the 'white' palette when visualizing a 3D volume mesh only
               vss->palette.SetIndex(11);
               vss->SetLightMatIdx(4);
            }
            else
            {
               // Use the 'bone' palette when visualizing a surface mesh only
               vss->palette.SetIndex(4);
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
            mesh_range = data_state.sol.Max() + 1.0;
         }
      }
   }
   else if (field_type == DataState::FieldType::VECTOR)
   {
      if (data_state.mesh->SpaceDimension() == 2)
      {
         vs = new VisualizationSceneVector(*this);
      }
      else if (data_state.mesh->SpaceDimension() == 3)
      {
         if (data_state.grid_f)
         {
            data_state.ProjectVectorFEGridFunction();
         }
         vs = new VisualizationSceneVector3d(*this);
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
         SetVisualizationScene(vs, 2, data_state.keys.c_str());
      }
      else
      {
         SetVisualizationScene(vs, 3, data_state.keys.c_str());
      }
   }
   return true;
}

void Window::GLVisStartVis()
{
   RunVisualization(); // deletes vs
   vs = NULL;
   if (glvis_command)
   {
      glvis_command->Terminate();
      delete comm_thread;
      delete glvis_command;
      glvis_command = NULL;
   }
   cout << "GLVis window closed." << endl;
}

thread_local Window *Window::locwin = NULL;

void Window::SwitchQuadSolution()
{
   int iqs = ((int)locwin->data_state.GetQuadSolution()+1)
             % ((int)DataState::QuadSolution::MAX);
   locwin->data_state.SwitchQuadSolution((DataState::QuadSolution)iqs, locwin->vs);
   SendExposeEvent();
}
