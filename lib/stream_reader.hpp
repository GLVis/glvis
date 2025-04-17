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

#ifndef GLVIS_STREAM_READER_HPP
#define GLVIS_STREAM_READER_HPP

#include <string>
#include <vector>
#include <iostream>
#include <memory>
#include "mfem.hpp"
#include "openglvis.hpp"

using StreamCollection = std::vector<std::unique_ptr<std::istream>>;

struct StreamState
{
private:
   struct
   {
      std::unique_ptr<mfem::Mesh> mesh;
      std::unique_ptr<mfem::Mesh> mesh_quad;
      std::unique_ptr<mfem::GridFunction> grid_f;
      std::unique_ptr<mfem::QuadratureFunction> quad_f;
   } internal;

public:
   mfem::Vector sol, solu, solv, solw, normals;
   std::string keys;
   const std::unique_ptr<mfem::Mesh> &mesh{internal.mesh};
   const std::unique_ptr<mfem::Mesh> &mesh_quad{internal.mesh_quad};
   const std::unique_ptr<mfem::GridFunction> &grid_f{internal.grid_f};
   const std::unique_ptr<mfem::QuadratureFunction> &quad_f{internal.quad_f};
   bool fix_elem_orient{false};
   bool save_coloring{false};
   bool keep_attr{false};

   enum class FieldType
   {
      UNKNOWN = -1,
      MIN = -1,
      //----------
      SCALAR,
      VECTOR,
      MESH,
      //----------
      MAX
   };

   enum class QuadSolution
   {
      NONE = -1,
      MIN = -1,
      //----------
      LOR_ClosedGL,
      HO_L2_collocated,
      HO_L2_projected,
      //----------
      MAX
   } quad_sol {QuadSolution::NONE};

   StreamState() = default;
   StreamState(StreamState &&ss) { *this = std::move(ss); }
   StreamState& operator=(StreamState &&ss);

   void SetMesh(mfem::Mesh *mesh);
   void SetMesh(std::unique_ptr<mfem::Mesh> &&pmesh);

   void SetGridFunction(mfem::GridFunction *gf);
   void SetGridFunction(std::unique_ptr<mfem::GridFunction> &&pgf);

   void SetQuadFunction(mfem::QuadratureFunction *qf);
   void SetQuadFunction(std::unique_ptr<mfem::QuadratureFunction> &&pqf);

   /// Helper function for visualizing 1D or 2D3V data
   void ExtrudeMeshAndSolution();

   /// Helper function for visualizing 1D data
   void Extrude1DMeshAndSolution();

   /// Helper function for visualization of 2D3V data
   void Extrude2D3VMeshAndSolution();

   /// Helper function to build the quadrature function from pieces
   void CollectQuadratures(mfem::QuadratureFunction *qf_array[], int npieces);

   /// Set a (checkerboard) solution when only the mesh is given
   void SetMeshSolution();

   /// Set the quadrature function representation producing a proxy grid function
   void SetQuadSolution(QuadSolution type = QuadSolution::LOR_ClosedGL);

   /// Switch the quadrature function representation and update the visualization
   void SwitchQuadSolution(QuadSolution type, VisualizationScene* vs);

   /// Get the current representation of quadrature solution
   inline QuadSolution GetQuadSolution() const { return quad_sol; }

   FieldType ReadStream(std::istream &is, const std::string &data_type);

   FieldType ReadStreams(const StreamCollection& input_streams);

   void WriteStream(std::ostream &os);

   // Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
   // or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
   // product vector grid function of the same order.
   static std::unique_ptr<mfem::GridFunction>
   ProjectVectorFEGridFunction(std::unique_ptr<mfem::GridFunction> gf);

   void ProjectVectorFEGridFunction()
   { internal.grid_f = ProjectVectorFEGridFunction(std::move(internal.grid_f)); }

   /// Sets a new mesh and solution from another StreamState object, and
   /// updates the given VisualizationScene pointer with the new data.
   ///
   /// Mesh space and grid function dimensions must both match the original
   /// dimensions of the current StreamState. If there is a mismatch in either
   /// value, the function will return false, and the mesh/solution will not be
   /// updated.
   bool SetNewMeshAndSolution(StreamState new_state,
                              VisualizationScene* vs);

   /// Updates the given VisualizationScene pointer with the new data
   /// of the given StreamState object.
   /// @note: Use with caution when the update is compatible
   /// @see SetNewMeshAndSolution()
   static void ResetMeshAndSolution(StreamState &ss, VisualizationScene* vs);
};

#endif // GLVIS_STREAM_READER_HPP
