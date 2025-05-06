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

#ifndef GLVIS_STREAM_READER_HPP
#define GLVIS_STREAM_READER_HPP

#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include "mfem.hpp"
#include "data_state.hpp"

using StreamCollection = std::vector<std::unique_ptr<std::istream>>;

class StreamReader
{
   DataState &data;

public:

   StreamReader(DataState &data_): data(data_) { }

   /// Prints available commands
   static void PrintCommands();

   /// Tests if the data type is supported
   static bool SupportsDataType(const std::string &data_type);

   /// Read the content of an input stream (e.g. from socket/file)
   int ReadStream(std::istream &is, const std::string &data_type);

   /// Read the content of an input streams (e.g. from sockets/files)
   int ReadStreams(const StreamCollection& input_streams);

   /// Write the state to a stream (e.g. to socket/file)
   void WriteStream(std::ostream &os);
};

#endif // GLVIS_STREAM_READER_HPP
