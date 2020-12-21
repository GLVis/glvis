// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
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
#include "mfem.hpp"

struct StreamState
{
   mfem::Vector sol, solu, solv, solw, normals;
   std::string keys;
   mfem::Mesh *mesh{nullptr};
   mfem::GridFunction *grid_f{nullptr};
   int is_gf{0};
   bool fix_elem_orient{false};
   bool save_coloring{false};
};
extern StreamState stream_state;

void Extrude1DMeshAndSolution(mfem::Mesh **mesh_p,
                              mfem::GridFunction **grid_f_p,
                              mfem::Vector *sol);
void SetMeshSolution(mfem::Mesh *mesh, mfem::GridFunction *&grid_f,
                     bool save_coloring);
int ReadStream(std::istream &is, const std::string &data_type);

#endif // GLVIS_STREAM_READER_HPP
