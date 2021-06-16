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

static GLVisWindow * mainWindow = nullptr;
static VisualizationSceneScalarData * vssd = nullptr;
static bool keep_attr = false;

namespace js
{

using namespace mfem;

bool startVisualization(const std::string input, const std::string data_type,
                        int w, int h)
{
   std::stringstream ss(input);

   StreamState stream_state;
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

   if (mainWindow == nullptr)
   {
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
      catch (...)
      {
         cerr << "Initializing the visualization failed - unknown error."
              << endl;
         return false;
      }
   }

   double mesh_range = -1.0;
   mainWindow->InitVisualization(field_type, std::move(stream_state), keep_attr);
   vssd = dynamic_cast<VisualizationSceneScalarData*>(mainWindow->getScene());

   if (minv || maxv)
   {
      vssd->SetValueRange(minv, maxv);
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
   double mesh_range = -1.0;

   StreamState & stream_state = mainWindow->getStreamState();
   if (stream_state.SetNewMeshAndSolution(std::move(new_state), vssd))
   {
      if (mesh_range > 0.0)
      {
         vssd->SetValueRange(-mesh_range, mesh_range);
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
   mainWindow->getSdl()->mainIter();
}

void setCanvasId(const std::string & id)
{
   std::cout << "glvis: setting canvas id to " << id << std::endl;
   mainWindow->getSdl()->setCanvasId(id);
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
   mainWindow->CallKeySequence(keys.c_str());
}

void processKey(int sym, bool ctrl=false, bool shift=false, bool alt=false)
{
   Uint16 mod = 0;
   mod |= ctrl ? KMOD_CTRL : 0;
   mod |= shift ? KMOD_SHIFT : 0;
   mod |= alt ? KMOD_ALT : 0;
   mainWindow->getSdl()->callKeyDown(sym, mod);
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
   return vssd->GetHelpString();
}

void ResizeWindow(int w, int h)
{
   mainWindow->ResizeWindow(w, h);
}

int GetUseTexture()
{
   return vssd->palette.GetSmoothSetting();
}

void SetUseTexture(int ut)
{
   if (ut == 0)
   {
      vssd->palette.UseDiscrete();
   }
   else
   {
      vssd->palette.UseSmooth();
   }
}

void SendExposeEvent()
{
   mainWindow->SendExposeEvent();
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
   em::function("sendExposeEvent", &js::SendExposeEvent);
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
   em::function("processKeys", &js::processKeys);
   em::function("processKey", &js::processKey);
}
