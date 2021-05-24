// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
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

namespace js
{

using namespace mfem;

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

   if (InitVisualization("glvis", 0, 0, w, h))
   {
      return false;
   }

   delete vs;
   vs = nullptr;

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
            // Use the 'bone' palette when visualizing a 2D mesh only (otherwise
            // the 'jet-like' palette is used in 2D, see vssolution.cpp).
            vs->palette.SetIndex(4);
         }
      }
      else if (stream_state.mesh->SpaceDimension() == 3)
      {
         VisualizationSceneSolution3d * vss;
         vs = vss = new VisualizationSceneSolution3d(*stream_state.mesh,
                                                     stream_state.sol);
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(stream_state.grid_f.get());
         }
         if (field_type == 2)
         {
            if (stream_state.mesh->Dimension() == 3)
            {
               // Use the 'white' palette when visualizing a 3D volume mesh only
               // vs->palette.SetIndex(4);
               vs->palette.SetIndex(11);
               vss->SetLightMatIdx(4);
            }
            else
            {
               // Use the 'bone' palette when visualizing a surface mesh only
               // (the same as when visualizing a 2D mesh only)
               vs->palette.SetIndex(4);
            }
            // Otherwise, the 'vivid' palette is used in 3D see vssolution3d.cpp

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
            stream_state.grid_f
               = ProjectVectorFEGridFunction(std::move(stream_state.grid_f));
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

   if (minv || maxv)
   {
      vs->SetValueRange(minv, maxv);
   }

   SendExposeEvent();
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
   double mesh_range = -1.0;

   if (stream_state.SetNewMeshAndSolution(std::move(new_state), vs))
   {
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
      }

      SendExposeEvent();
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

void processKeys(const std::string & keys)
{
   CallKeySequence(keys.c_str());
}

void processKey(int sym, bool ctrl=false, bool shift=false, bool alt=false)
{
   Uint16 mod = 0;
   mod |= ctrl ? KMOD_CTRL : 0;
   mod |= shift ? KMOD_SHIFT : 0;
   mod |= alt ? KMOD_ALT : 0;
   GetAppWindow()->callKeyDown(sym, mod);
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
} // namespace js

// Info on type conversion:
// https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html#built-in-type-conversions
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
   em::function("getTextureMode", &GetUseTexture);
   em::function("setTextureMode", &SetUseTexture);
   em::function("resizeWindow", &ResizeWindow);
   em::function("setCanvasId", &js::setCanvasId);
   em::function("setupResizeEventCallback", &js::setupResizeEventCallback);
   em::function("getHelpString", &js::getHelpString);
   em::function("processKeys", &js::processKeys);
   em::function("processKey", &js::processKey);
}
