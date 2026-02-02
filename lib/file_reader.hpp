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

#ifndef GLVIS_FILE_READER_HPP
#define GLVIS_FILE_READER_HPP

#include "data_state.hpp"

class FileReader
{
   DataState &data;
   int pad_digits;

   static bool CheckStreamIsComplex(std::istream &sol, bool parallel = false);

   int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                                  const char *sol_prefix, int component = -1);
   int ReadParMeshAndQuadFunction(int np, const char *mesh_prefix,
                                  const char *sol_prefix, int component = -1);

public:
   enum class FileType
   {
      MESH,
      SCALAR_SOL,
      VECTOR_SOL,
      GRID_FUNC,
      QUAD_FUNC,
   };

   FileReader(DataState &data_, int pad_digits_ = 6)
      : data(data_), pad_digits(pad_digits_) { }

   /// Read the mesh and the solution from a file
   int ReadSerial(FileType ft, const char *mesh_file, const char *sol_file,
                  int component = -1);

   /// Read the mesh and the solution from multiple files
   int ReadParallel(int np, FileType ft, const char *mesh_file,
                    const char *sol_file, int component = -1);
};

#endif // GLVIS_FILE_READER_HPP
