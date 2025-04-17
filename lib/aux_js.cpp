// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "aux_vis.hpp"
#include "palettes.hpp"
#include "stream_reader.hpp"
#include "visual.hpp"

#ifdef GLVIS_USE_LIBPNG
#include <png.h>
#endif

#include <SDL2/SDL_hints.h>
#include <emscripten/bind.h>
#include <emscripten/html5.h>
#include <emscripten/val.h>

// used in extern context
thread_local std::string plot_caption;
thread_local std::string extra_caption;
thread_local mfem::GeometryRefiner GLVisGeometryRefiner;

static VisualizationSceneScalarData * vs = nullptr;

// either bitmap data or png bytes
std::vector<unsigned char> * screen_state = nullptr;

DataState stream_state;

int last_stream_nproc = 1;

namespace em = emscripten;

namespace js
{

using namespace mfem;

/// Switch representation of the quadrature function
void SwitchQuadSolution();

/// Display a new stream
void display(std::stringstream & commands, const int w, const int h)
{
   // reset antialiasing
   GetAppWindow()->getRenderer().setAntialiasing(0);

   std::string word;
   double minv = 0.0, maxv = 0.0;
   while (commands >> word)
   {
      if (word == "keys")
      {
         std::cout << "parsing 'keys'" << std::endl;
         commands >> stream_state.keys;
      }
      else if (word == "valuerange")
      {
         std::cout << "parsing 'valuerange'" << std::endl;
         commands >> minv >> maxv;
      }
      else
      {
         std::cout << "unknown command '" << word << "'" << std::endl;
      }
   }

   DataState::FieldType field_type = stream_state.GetType();

   if (field_type <= DataState::FieldType::MIN
       || field_type >= DataState::FieldType::MAX)
   {
      return;
   }

   if (InitVisualization("glvis", 0, 0, w, h))
   {
      return;
   }

   delete vs;
   vs = nullptr;

   if (stream_state.quad_f)
   {
      GetAppWindow()->setOnKeyDown('Q', SwitchQuadSolution);
   }

   double mesh_range = -1.0;
   if (field_type == DataState::FieldType::SCALAR
       || field_type == DataState::FieldType::MESH)
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
                                                      stream_state.mesh_quad.get(), &stream_state.normals);
         }
         else
         {
            vs = vss = new VisualizationSceneSolution(*stream_state.mesh, stream_state.sol,
                                                      stream_state.mesh_quad.get());
         }
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(*stream_state.grid_f);
         }
         if (field_type == DataState::FieldType::MESH)
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
                                                     stream_state.sol, stream_state.mesh_quad.get());
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(stream_state.grid_f.get());
         }
         if (field_type == DataState::FieldType::MESH)
         {
            if (stream_state.mesh->Dimension() == 3)
            {
               // Use the 'white' palette when visualizing a 3D volume mesh only
               vs->palette.SetIndex(11);
               vss->SetLightMatIdx(4);
            }
            else
            {
               // Use the 'bone' palette when visualizing a surface mesh only
               vs->palette.SetIndex(4);
            }
            // Otherwise, the 'vivid' palette is used in 3D see vssolution3d.cpp

            vss->ToggleDrawAxes();
            vss->ToggleDrawMesh();
         }
      }
      if (field_type == DataState::FieldType::MESH)
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
   else if (field_type == DataState::FieldType::VECTOR)
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
                                              stream_state.solv, stream_state.mesh_quad.get());
         }
      }
      else if (stream_state.mesh->SpaceDimension() == 3)
      {
         if (stream_state.grid_f)
         {
            stream_state.ProjectVectorFEGridFunction();
            vs = new VisualizationSceneVector3d(*stream_state.grid_f,
                                                stream_state.mesh_quad.get());
         }
         else
         {
            vs = new VisualizationSceneVector3d(*stream_state.mesh, stream_state.solu,
                                                stream_state.solv, stream_state.solw,
                                                stream_state.mesh_quad.get());
         }
      }
   }

   if (vs)
   {
      // increase the refinement factors if visualizing a GridFunction
      if (stream_state.grid_f)
      {
         vs->AutoRefine();
         vs->SetShading(VisualizationSceneScalarData::Shading::Noncomforming, true);
      }
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
         vs->SetAutoscale(0);
      }
      if (stream_state.mesh->SpaceDimension() == 2 &&
          field_type == DataState::FieldType::MESH)
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
}

//
// StreamReader::ReadStream requires a list of unique_ptr to istream and since
// we cannot directly pass a list of string we need to repack the strings into
// a new list.
//
// each string in streams must start with `parallel <nproc> <rank>'
//
using StringArray = std::vector<std::string>;
void processParallelStreams(DataState & state,
                            const StringArray & streams,
                            std::stringstream * commands = nullptr)
{
   // std::cerr << "got " << streams.size() << " streams" << std::endl;
   // HACK: match unique_ptr<istream> interface for ReadStreams:
   std::vector<std::stringstream> sstreams(streams.size());
   StreamCollection istreams(streams.size());
   for (int i = 0; i < streams.size(); ++i)
   {
      sstreams[i] = std::stringstream(streams[i]);
      // pull off the first list
      std::string word;
      int nproc, rank;
      sstreams[i] >> word >> nproc >> rank;
      // std::cerr << "packing " << rank+1 << "/" << nproc << std::endl;
      istreams[i] = std::unique_ptr<std::istream>(&sstreams[i]);
   }

   StreamReader reader(state);
   reader.ReadStreams(istreams);

   if (commands)
   {
      commands->seekg(istreams[0]->tellg());
   }

   // HACK: don't let unique_ptr free the data
   for (int i = 0; i < streams.size(); ++i)
   {
      istreams[i].release();
   }

   last_stream_nproc = streams.size();
}

void displayParallelStreams(const StringArray & streams, const int w,
                            const int h)
{
   std::stringstream commands(streams[0]);
   processParallelStreams(stream_state, streams, &commands);

   display(commands, w, h);
}

void displayStream(const std::string & stream, const int w, const int h)
{
   std::stringstream ss(stream);
   std::string data_type;
   ss >> data_type;

   StreamReader reader(stream_state);
   reader.ReadStream(ss, data_type);

   display(ss, w, h);
}

//
// update the existing stream
//
int update(DataState & new_state)
{
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
      std::cerr << "update: field type does not match!" << std::endl;
      return 1;
   }
}


int updateStream(std::string stream)
{
   std::stringstream ss(stream);
   std::string data_type;
   ss >> data_type;

   DataState new_state;
   StreamReader reader(new_state);
   reader.ReadStream(ss, data_type);

   return update(new_state);
}


int updateParallelStreams(const StringArray & streams)
{
   if (streams.size() != last_stream_nproc)
   {
      std::cout << "current stream-list size does not match the previous: " <<
                streams.size() << " != " << last_stream_nproc << std::endl;
      return 1;
   }
   DataState new_state;
   processParallelStreams(new_state, streams);

   return update(new_state);
}

//
// other methods
//
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

em::val getScreenBuffer(bool flip_y=false)
{
   MyExpose();
   auto * wnd = GetAppWindow();
   int w, h;
   wnd->getGLDrawSize(w, h);

   // 4 bytes for rgba
   const size_t buffer_size = w*h*4;
   if (screen_state == nullptr)
   {
      screen_state = new std::vector<unsigned char>(buffer_size);
   }
   else if (buffer_size > screen_state->size())
   {
      screen_state->resize(buffer_size);
   }

   glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, screen_state->data());

   // `canvas` starts at top left but `glReadPixels' starts at bottom left
   if (flip_y)
   {
      auto * orig = screen_state;
      auto * flip = new std::vector<unsigned char>(buffer_size);
      // from `glReadPixels` man page: Pixels are returned in row order from
      // the lowest to the highest row, left to right in each row.
      for (int j = 0; j < h; ++j)
      {
         for (int i = 0; i < w*4; ++i)
         {
            (*flip)[4*w*j + i] = (*orig)[4*w*(h-j-1) + i];
         }
      }
      screen_state = flip;
      delete orig;
   }

   return em::val(em::typed_memory_view(screen_state->size(),
                                        screen_state->data()));
}

void SwitchQuadSolution()
{
   int iqs = ((int)stream_state.GetQuadSolution()+1)
             % ((int)DataState::QuadSolution::MAX);
   stream_state.SwitchQuadSolution((DataState::QuadSolution)iqs, vs);
   SendExposeEvent();
}

#ifdef GLVIS_USE_LIBPNG
em::val getPNGByteArray()
{
   constexpr const char * filename = "im.png";
   auto * wnd = GetAppWindow();
   int w, h;
   wnd->getGLDrawSize(w, h);

   MyExpose();
   // save to in-memory file
   int status = SaveAsPNG(filename, w, h, wnd->isHighDpi(), true);
   if (status != 0)
   {
      fprintf(stderr, "unable to generate png\n");
      return em::val::null();
   }

   // load in-memory file to byte array
   std::ifstream ifs(filename, std::ios::binary);
   if (!ifs)
   {
      fprintf(stderr, "unable to load png\n");
      return em::val::null();
   }

   if (!screen_state)
   {
      screen_state = new std::vector<unsigned char>(std::istreambuf_iterator<char>
                                                    (ifs), {});
   }
   else
   {
      screen_state->assign(std::istreambuf_iterator<char>
                           (ifs), {});
   }

   return em::val(em::typed_memory_view(screen_state->size(),
                                        screen_state->data()));
}
#endif // GLVIS_USE_LIBPNG
} // namespace js

// Info on type conversion:
// https://emscripten.org/docs/porting/connecting_cpp_and_javascript/embind.html#built-in-type-conversions
EMSCRIPTEN_BINDINGS(js_funcs)
{
   em::enum_<DataState::FieldType>("FieldType")
   .value("UNKNOWN", DataState::FieldType::UNKNOWN)
   .value("MIN", DataState::FieldType::MIN)
   .value("MESH", DataState::FieldType::MESH)
   .value("SCALAR", DataState::FieldType::SCALAR)
   .value("VECTOR", DataState::FieldType::VECTOR)
   .value("MAX", DataState::FieldType::MAX)
   ;
   em::function("displayStream", &js::displayStream);
   em::function("displayParallelStreams", &js::displayParallelStreams);
   em::function("updateStream", &js::updateStream);
   em::function("updateParallelStreams", &js::updateParallelStreams);
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
   em::function("getScreenBuffer", &js::getScreenBuffer);
#ifdef GLVIS_USE_LIBPNG
   em::function("getPNGByteArray", &js::getPNGByteArray);
#endif
   em::register_vector<std::string>("StringArray");
}
