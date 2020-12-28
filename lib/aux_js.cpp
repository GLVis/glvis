// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "visual.hpp"
#include "palettes.hpp"
#include "stream_reader.hpp"
#include <SDL2/SDL_hints.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>

std::string plot_caption;
std::string extra_caption; // used in extern context
mfem::GeometryRefiner GLVisGeometryRefiner; // used in extern context

namespace js
{

using namespace mfem;

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
GridFunction *ProjectVectorFEGridFunction(GridFunction *gf)
{
   if ((gf->VectorDim() == 3) && (gf->FESpace()->GetVDim() == 1))
   {
      int p = gf->FESpace()->GetOrder(0);
      cout << "Switching to order " << p
           << " discontinuous vector grid function..." << endl;
      int dim = gf->FESpace()->GetMesh()->Dimension();
      FiniteElementCollection *d_fec = new L2_FECollection(p, dim, 1);
      FiniteElementSpace *d_fespace =
         new FiniteElementSpace(gf->FESpace()->GetMesh(), d_fec, 3);
      GridFunction *d_gf = new GridFunction(d_fespace);
      d_gf->MakeOwner(d_fec);
      gf->ProjectVectorFieldOn(*d_gf);
      delete gf;
      return d_gf;
   }
   return gf;
}

bool startVisualization(const std::string input, const std::string data_type,
                        int w, int h)
{
   std::stringstream ss(input);
   const int field_type = ReadStream(ss, data_type);

   std::string line;
   while (ss >> line)
   {
      if (line == "keys")
      {
         std::cout << "parsing 'keys'" << std::endl;
         ss >> stream_state.keys;
      }
      else
      {
         std::cout << "unknown line '" << line << "'" << std::endl;
      }
   }

   if (field_type < 0 || field_type > 2)
   {
      return false;
   }

   if (InitVisualization("glvis", 0, 0, w, h))
   {
      return false;
   }

   VisualizationSceneScalarData * vs;

   double mesh_range = -1.0;
   if (field_type == 0 || field_type == 2)
   {
      if (stream_state.grid_f)
      {
         stream_state.grid_f->GetNodalValues(stream_state.sol);
      }
      if (stream_state.mesh->SpaceDimension() == 2)
      {
         VisualizationSceneSolution * vss;
         if (field_type == 2)
         {
            paletteSet(4);
         }
         if (stream_state.normals.Size() > 0)
         {
            vs = vss = new VisualizationSceneSolution(*stream_state.mesh, stream_state.sol,
                                                      &stream_state.normals);
         }
         else
         {
            vs = vss = new VisualizationSceneSolution(*stream_state.mesh, stream_state.sol);
         }
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(*stream_state.grid_f);
         }
         if (field_type == 2)
         {
            vs->OrthogonalProjection = 1;
            vs->SetLight(0);
            vs->Zoom(1.8);
         }
      }
      else if (stream_state.mesh->SpaceDimension() == 3)
      {
         VisualizationSceneSolution3d * vss;
         vs = vss = new VisualizationSceneSolution3d(*stream_state.mesh,
                                                     stream_state.sol);
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(stream_state.grid_f);
         }
         if (field_type == 2)
         {
            if (stream_state.mesh->Dimension() == 3)
            {
               paletteSet(4);
               // paletteSet(11);
               // Set_Material_And_Light(4,3);
            }
            else
            {
               paletteSet(4);
            }
            vss->ToggleDrawAxes();
            vss->ToggleDrawMesh();
         }
      }
      if (field_type == 2)
      {
         if (stream_state.grid_f)
         {
            mesh_range = stream_state.grid_f->Max() + 1.0;
         }
         else
         {
            mesh_range = stream_state.sol.Max() + 1.0;
         }
      }
   }
   else if (field_type == 1)
   {
      if (stream_state.mesh->SpaceDimension() == 2)
      {
         if (stream_state.grid_f)
         {
            vs = new VisualizationSceneVector(*stream_state.grid_f);
         }
         else
         {
            vs = new VisualizationSceneVector(*stream_state.mesh, stream_state.solu,
                                              stream_state.solv);
         }
      }
      else if (stream_state.mesh->SpaceDimension() == 3)
      {
         if (stream_state.grid_f)
         {
            stream_state.grid_f = ProjectVectorFEGridFunction(stream_state.grid_f);
            vs = new VisualizationSceneVector3d(*stream_state.grid_f);
         }
         else
         {
            vs = new VisualizationSceneVector3d(*stream_state.mesh, stream_state.solu,
                                                stream_state.solv, stream_state.solw);
         }
      }
   }

   if (vs)
   {
      // increase the refinement factors if visualizing a GridFunction
      if (stream_state.grid_f)
      {
         vs->AutoRefine();
         vs->SetShading(2, true);
      }
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
         vs->SetAutoscale(0);
      }
      if (stream_state.mesh->SpaceDimension() == 2 && field_type == 2)
      {
         SetVisualizationScene(vs, 2);
      }
      else
      {
         SetVisualizationScene(vs, 3);
      }
   }

   CallKeySequence(stream_state.keys.c_str());

   SendExposeEvent();
   return true;
}

void iterVisualization()
{
   GetAppWindow()->mainIter();
}

void setCanvasId(const std::string & id)
{
   std::cout << "glvis: setting canvas id to " << id << std::endl;
   GetAppWindow()->setCanvasId(id);
}

void disableKeyHandling()
{
   SDL_EventState(SDL_KEYDOWN, SDL_DISABLE);
   SDL_EventState(SDL_KEYUP, SDL_DISABLE);
}

void enableKeyHandling()
{
   SDL_EventState(SDL_KEYDOWN, SDL_ENABLE);
   SDL_EventState(SDL_KEYUP, SDL_ENABLE);
}

void setKeyboardListeningElementId(const std::string & id)
{
   SDL_SetHint(SDL_HINT_EMSCRIPTEN_KEYBOARD_ELEMENT, id.c_str());
}

void setupResizeEventCallback(const std::string & id)
{
   // typedef EM_BOOL (*em_ui_callback_func)(int eventType, const EmscriptenUiEvent *uiEvent, void *userData);
   std::cout << "glvis: adding resize callback for " << id << std::endl;
   auto err = emscripten_set_resize_callback(id.c_str(), nullptr,
                                             true, [](int eventType, const EmscriptenUiEvent *uiEvent,
                                                      void *userData) -> EM_BOOL
   {
      std::cout << "got resize event" << std::endl;
      return true;
   });
   // TODO: macro to wrap this
   if (err != EMSCRIPTEN_RESULT_SUCCESS)
   {
      std::cerr << "error (emscripten_set_resize_callback): " << err << std::endl;
   }
}
} // namespace js

namespace em = emscripten;
EMSCRIPTEN_BINDINGS(js_funcs)
{
   em::function("startVisualization", &js::startVisualization);
   em::function("iterVisualization", &js::iterVisualization);
   em::function("sendExposeEvent", &SendExposeEvent);
   em::function("disableKeyHanding", &js::disableKeyHandling);
   em::function("enableKeyHandling", &js::enableKeyHandling);
   em::function("setKeyboardListeningElementId",
                js::setKeyboardListeningElementId);
   em::function("getTextureMode", &GetUseTexture);
   em::function("setTextureMode", &SetUseTexture);
   em::function("resizeWindow", &ResizeWindow);
   em::function("setCanvasId", &js::setCanvasId);
   em::function("setupResizeEventCallback", &js::setupResizeEventCallback);
}
