// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
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
#include <memory>
#include "mfem.hpp"
#include "openglvis.hpp"

struct StreamState
{
   mfem::Vector sol, solu, solv, solw, normals;
   std::string keys;
   std::unique_ptr<mfem::Mesh> mesh;
   std::unique_ptr<mfem::GridFunction> grid_f;
   int is_gf{0};
   bool fix_elem_orient{false};
   bool save_coloring{false};

   /// Helper function for visualizing 1D data
   void Extrude1DMeshAndSolution();

   /// Set a (checkerboard) solution when only the mesh is given
   void SetMeshSolution();

   int ReadStream(std::istream &is, const std::string &data_type);

   /// Sets a new mesh and solution from another StreamState object, and
   /// updates the given VisualizationScene pointer with the new data.
   ///
   /// Mesh space and grid function dimensions must both match the original
   /// dimensions of the current StreamState. If there is a mismatch in either
   /// value, the function will return false, and the mesh/solution will not be
   /// updated.
   bool SetNewMeshAndSolution(StreamState new_state,
                              VisualizationScene* vs);
};

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
std::unique_ptr<mfem::GridFunction>
ProjectVectorFEGridFunction(std::unique_ptr<mfem::GridFunction> gf);

#endif // GLVIS_STREAM_READER_HPP
