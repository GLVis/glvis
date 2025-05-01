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

#include "coll_reader.hpp"

#include <vector>

using namespace std;
using namespace mfem;

int DataCollectionReader::ReadSerial(CollType ct, const char *coll_file,
                                     int ti, const char *field, bool quad, int component)
{
   switch (ct)
   {
      case CollType::VISIT:
      {
         auto dc = new VisItDataCollection(coll_file);
         dc->SetPadDigits(pad_digits);
         data.SetDataCollectionField(dc, ti, field, quad, component);
      }
      break;
#ifdef MFEM_USE_SIDRE
      case CollType::SIDRE:
      {
         auto dc = new SidreDataCollection(coll_file);
         dc->SetPadDigits(pad_digits);
         data.SetDataCollectionField(dc, ti, field, quad, component);
      }
      break;
#endif // MFEM_USE_SIDRE
#ifdef MFEM_USE_FMS
      case CollType::FMS:
      {
         auto dc = new FMSDataCollection(coll_file);
         dc->SetPadDigits(pad_digits);
         dc->SetProtocol(protocol);
         data.SetDataCollectionField(dc, ti, field, quad, component);
      }
      break;
#endif // MFEM_USE_FMS
#ifdef MFEM_USE_CONDUIT
      case CollType::CONDUIT:
      {
         auto dc = new ConduitDataCollection(coll_file);
         dc->SetPadDigits(pad_digits);
         dc->SetProtocol(protocol);
         data.SetDataCollectionField(dc, ti, field, quad, component);
      }
      break;
#endif // MFEM_USE_CONDUIT
      default:
         cerr << "Unknown collection type. Exit" << endl;
         exit(1);
   }

   data.ExtrudeMeshAndSolution();
   return 0;
}
