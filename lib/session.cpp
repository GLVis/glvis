// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "threads.hpp" // IWYU pragma: keep for GLVisCommand
#include "aux_vis.hpp" // SetUseHiDPI

#include "session.hpp"

#include <fem/geom.hpp>

#ifdef NVTX_DBG_HPP
#undef NVTX_COLOR
#define NVTX_COLOR ::nvtx::kGold
#include NVTX_DBG_HPP
#else
#define dbg(...)
#endif

extern thread_local mfem::GeometryRefiner GLVisGeometryRefiner;

void Session::StartSession()
{
   auto funcThread = [](Window w, StreamCollection is)
   {
      if (w.GLVisInitVis(std::move(is)))
      {
         w.GLVisStartVis();
      }
   };
   handler = std::thread {funcThread,
                          std::move(win), std::move(input_streams)};
   handler.detach();
}

bool Session::StartSavedSession(std::string stream_file)
{
   std::unique_ptr<std::ifstream> ifs(new std::ifstream(stream_file));
   if (!(*ifs))
   {
      std::cout << "Can not open stream file: " << stream_file << std::endl;
      return false;
   }
   std::string data_type;
   *ifs >> data_type >> std::ws;
   StreamReader reader(win.data_state);
   reader.ReadStream(*ifs, data_type);
   input_streams.emplace_back(std::move(ifs));

   StartSession();
   return true;
}

int Session::StartStreamSession(std::unique_ptr<mfem::socketstream> &&stream,
                                const std::string &data_type)
{
   StreamReader reader(win.data_state);
   int ierr = reader.ReadStream(*stream, data_type);
   if (ierr) { return ierr; }
   input_streams.emplace_back(std::move(stream));

   StartSession();
   return 0;
}

int Session::StartStreamSession(StreamCollection &&streams)
{
   dbg();
   StreamReader reader(win.data_state);
   int ierr = reader.ReadStreams(streams);
   if (ierr) { return ierr; }
   input_streams = std::move(streams);

   StartSession();
   return 0;
}

int Session::StartSerialStreamSession(std::unique_ptr<std::istream> &&stream,
                                      const std::string &data_type)
{
   dbg("data_type: '{}'", data_type);
   StreamReader reader(win.data_state);
   int ierr = reader.ReadStream(*stream, data_type);
   if (ierr) { dbg("❌ ERROR ❌"); return ierr; }

   input_streams.emplace_back(std::move(stream));
   StartSession();
   return 0;
}

///////////////////////////////////////////////////////////////////////////////
int GLVisLibWindow(//void *win_ptr,
   bool fix_elem_orient,
   bool save_coloring, bool headless,
   const std::string &plot_caption,
   std::unique_ptr<std::istream> &&stream,
   const std::string &data_type)
{
   const int geom_ref_type = mfem::Quadrature1D::ClosedUniform;
   const bool enable_hidpi = true;

#ifdef _WIN32
   // Call needed to avoid SDL_Init failure when not substituting main() for
   // SDL_main().
   SDL_SetMainReady();
#endif

   Window win;
   dbg();
   // std::this_thread::sleep_for(std::chrono::milliseconds(500));

   std::cout << std::endl
             << "       _/_/_/  _/      _/      _/  _/"          << std::endl
             << "    _/        _/      _/      _/        _/_/_/" << std::endl
             << "   _/  _/_/  _/      _/      _/  _/  _/_/"      << std::endl
             << "  _/    _/  _/        _/  _/    _/      _/_/"   << std::endl
             << "   _/_/_/  _/_/_/_/    _/      _/  _/_/_/"      << std::endl
             << std::endl;

   SetUseHiDPI(enable_hidpi);
   GLVisGeometryRefiner.SetType(geom_ref_type);

   dbg("Main Window structure");
   //    auto *win = static_cast<Window*>(win_ptr);

   std::vector<Session> current_sessions;

   dbg("Get main thread");
   GetMainThread(win.headless);


   dbg("StartStreamSession");
   Session new_session(fix_elem_orient, save_coloring, plot_caption, headless);

   new_session.StartSerialStreamSession(std::move(stream), data_type);
   current_sessions.emplace_back(std::move(new_session));


   dbg("Starting message loop in main thread");
   MainThreadLoop();

   dbg("EXIT_SUCCESS");
   return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
void *GLVisLibGetWindow()
{
   static Window win;
   return (void*) &win;
}