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
#pragma once

#include "stream_reader.hpp"
#include "window.hpp"

#include <thread>

class Session
{
   StreamCollection input_streams;
   Window win;
   std::thread handler;

public:
   Session(bool fix_elem_orient,
           bool save_coloring,
           std::string plot_caption,
           bool headless)
   {
      win.data_state.fix_elem_orient = fix_elem_orient;
      win.data_state.save_coloring = save_coloring;
      win.plot_caption = plot_caption;
      win.headless = headless;
   }

   Session(Window other_win): win(std::move(other_win)) { }

   ~Session() = default;

   Session(Session&& from) = default;
   Session& operator= (Session&& from) = default;

   inline DataState& GetState() { return win.data_state; }
   inline const DataState& GetState() const { return win.data_state; }

   void StartSession();

   bool StartSavedSession(std::string stream_file);

   int StartStreamSession(std::unique_ptr<mfem::socketstream> &&stream,
                          const std::string &data_type);

   int StartStreamSession(std::unique_ptr<std::istream> &&stream,
                          const std::string &data_type);

   int StartStreamSession(StreamCollection &&streams);

   int StartSerialStreamSession(std::istream &stream,
                                const std::string &data_type);

};