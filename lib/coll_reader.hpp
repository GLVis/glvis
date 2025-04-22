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

#ifndef GLVIS_COLL_READER_HPP
#define GLVIS_COLL_READER_HPP

#include <iostream>
#include "mfem.hpp"
#include "data_state.hpp"

class DataCollectionReader
{
   DataState &data;
   std::string protocol;

public:
   enum class CollType
   {
      VISIT,
      PARAVIEW,
      SIDRE,
      FMS,
      CONDUIT,
      ADIOS2,
   };

   DataCollectionReader(DataState &data_) : data(data_) { }

   /// Set the I/O protocol to for loading of collections
   void SetProtocol(const char *protocol_) { protocol = protocol_; }

   /// Read the mesh from a file and solution from a collection
   int ReadSerial(CollType ct, const char *collection, int ti,
                  const char *field = NULL, bool quad = false, int component = -1);
};

#endif // GLVIS_COLL_READER_HPP
