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

#ifndef GLVIS_COLL_READER_HPP
#define GLVIS_COLL_READER_HPP

#include <iostream>
#include "mfem.hpp"
#include "data_state.hpp"

class DataCollectionReader
{
   DataState &data;
   int pad_digits = 6;
   std::string protocol;

public:
   enum class CollType
   {
      VISIT,
      PARAVIEW, // not currently implemented
      SIDRE,
      FMS,
      CONDUIT,
      ADIOS2, // not currently implemented
   };

   DataCollectionReader(DataState &data_) : data(data_) { }

   /// Set the padding digits
   void SetPadDigits(int digits) { pad_digits = digits; }

   /// Set the I/O protocol to for loading of collections
   void SetProtocol(const char *protocol_) { protocol = protocol_; }

   /// Load the mesh and solution from a collection
   /** Loads the mesh and optionally a grid or quadrature function from the
       provided data collection.
       @param ct        collection type
       @param dc        data collection
       @param ti        cycle to load
       @param field     name of the (Q-)field to load (NULL for mesh only)
       @param quad      if true, Q-field is loaded, otherwise a regular field
       @param component component of the field (-1 means all components) */
   int ReadSerial(CollType ct, const char *dc, int ti,
                  const char *field = NULL, bool quad = false, int component = -1);
};

#endif // GLVIS_COLL_READER_HPP
