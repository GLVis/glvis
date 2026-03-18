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

#include <memory>

#include "session.hpp"

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
   StreamReader reader(win.data_state);
   int ierr = reader.ReadStreams(streams);
   if (ierr) { return ierr; }
   input_streams = std::move(streams);

   StartSession();
   return 0;
}

int Session::StartStreamSession(std::unique_ptr<std::istream> &&stream,
                                const std::string &data_type)
{
   StreamReader reader(win.data_state);
   int ierr = reader.ReadStream(*stream, data_type);
   if (ierr) { return ierr; }
   input_streams.emplace_back(std::move(stream));

   StartSession();
   return 0;
}
