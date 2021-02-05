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

static VisualizationSceneScalarData * vs = nullptr;
StreamState stream_state;
GLVisWindow * mainWindow = nullptr;

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

   // 0 - scalar data, 1 - vector data, 2 - mesh only, (-1) - unknown
   const int field_type = stream_state.ReadStream(ss, data_type);

   std::string line;
   double minv = 0.0, maxv = 0.0;
   while (ss >> line)
   {
      if (line == "keys")
      {
         std::cout << "parsing 'keys'" << std::endl;
         ss >> stream_state.keys;
      }
      else if (line == "valuerange")
      {
         std::cout << "parsing 'valuerange'" << std::endl;
         ss >> minv >> maxv;
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

   try
   {
       mainWindow = new GLVisWindow("glvis", 0, 0, w, h, false);
   }
   catch (std::runtime_error& ex)
   {
      cerr << "Initializing the visualization failed: " << endl
           << ex.what() << endl;
      return false;
   }
   catch(...)
   {
       cerr << "Initializing the visualization failed - unknown error."
            << endl;
      return false;
   }

   delete vs;
   vs = nullptr;

   double mesh_range = -1.0;
   mainWindow->InitVisualization(field_type, std::move(stream_state));

   if (minv || maxv)
   {
      vs->SetValueRange(minv, maxv);
   }
   mainWindow->SendExposeEvent();
   return true;
}


int updateVisualization(std::string data_type, std::string stream)
{
   std::stringstream ss(stream);

   if (data_type != "solution")
   {
      std::cerr << "unsupported data type '" << data_type << "' for stream update" <<
                std::endl;
      return 1;
   }

   StreamState new_state;
   new_state.ReadStream(ss, data_type);


   auto * new_m = new Mesh(ss, 1, 0, stream_state.fix_elem_orient);
   auto * new_g = new GridFunction(new_m, ss);
   double mesh_range = -1.0;

   if (stream_state.SetNewMeshAndSolution(std::move(new_state), vs))
   {
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
      }
      mainWindow->SendExposeEvent();
      return 0;
   }
   else
   {
      cout << "Stream: field type does not match!" << endl;
      return 1;
   }
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

std::string getHelpString()
{
   VisualizationSceneScalarData* vss
      = dynamic_cast<VisualizationSceneScalarData*>(GetVisualizationScene());
   return vss->GetHelpString();
}

void ResizeWindow(int w, int h)
{
    mainWindow->ResizeWindow(w, h);
}

int GetUseTexture()
{
   return vs->GetPalette().GetSmoothSetting();
}

void SetUseTexture(int ut)
{
   if (ut == 0)
   {
      vs->GetPalette().UseDiscrete();
   }
   else
   {
      vs->GetPalette().UseSmooth();
   }
}

} // namespace js

namespace em = emscripten;
EMSCRIPTEN_BINDINGS(js_funcs)
{
   em::function("startVisualization", &js::startVisualization);
   em::function("updateVisualization", &js::updateVisualization);
   em::function("iterVisualization", &js::iterVisualization);
   em::function("sendExposeEvent", &SendExposeEvent);
   em::function("disableKeyHanding", &js::disableKeyHandling);
   em::function("enableKeyHandling", &js::enableKeyHandling);
   em::function("setKeyboardListeningElementId",
                js::setKeyboardListeningElementId);
   em::function("getTextureMode", &js::GetUseTexture);
   em::function("setTextureMode", &js::SetUseTexture);
   em::function("resizeWindow", &js::ResizeWindow);
   em::function("setCanvasId", &js::setCanvasId);
   em::function("setupResizeEventCallback", &js::setupResizeEventCallback);
   em::function("getHelpString", &js::getHelpString);
}
