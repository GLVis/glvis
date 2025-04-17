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

#ifndef GLVIS_DATA_STATE_HPP
#define GLVIS_DATA_STATE_HPP

#include <string>
#include <vector>
#include <memory>
#include "mfem.hpp"
#include "openglvis.hpp"

struct DataState
{
   enum class FieldType
   {
      UNKNOWN = -1,
      MIN = -1,
      //----------
      MESH,
      SCALAR,
      VECTOR,
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
   };

private:
   friend class StreamReader;
   friend class FileReader;
   struct
   {
      std::unique_ptr<mfem::Mesh> mesh;
      std::unique_ptr<mfem::Mesh> mesh_quad;
      std::unique_ptr<mfem::GridFunction> grid_f;
      std::unique_ptr<mfem::QuadratureFunction> quad_f;
   } internal;

   FieldType type {FieldType::UNKNOWN};
   QuadSolution quad_sol {QuadSolution::NONE};

   void SetGridFunctionSolution(int component = -1);
   void SetQuadFunctionSolution(int component = -1);

public:
   mfem::Vector sol, solu, solv, solw, normals;
   const std::unique_ptr<mfem::Mesh> &mesh{internal.mesh};
   const std::unique_ptr<mfem::Mesh> &mesh_quad{internal.mesh_quad};
   const std::unique_ptr<mfem::GridFunction> &grid_f{internal.grid_f};
   const std::unique_ptr<mfem::QuadratureFunction> &quad_f{internal.quad_f};

   std::string keys;
   bool fix_elem_orient{false};
   bool save_coloring{false};
   bool keep_attr{false};

   DataState() = default;
   DataState(DataState &&ss) { *this = std::move(ss); }
   DataState& operator=(DataState &&ss);

   inline FieldType GetType() const { return type; }

   void SetMesh(mfem::Mesh *mesh);
   void SetMesh(std::unique_ptr<mfem::Mesh> &&pmesh);

   void SetGridFunction(mfem::GridFunction *gf, int component = -1);
   void SetGridFunction(std::unique_ptr<mfem::GridFunction> &&pgf,
                        int component = -1);

   void SetQuadFunction(mfem::QuadratureFunction *qf, int component = -1);
   void SetQuadFunction(std::unique_ptr<mfem::QuadratureFunction> &&pqf,
                        int component = -1);
   void SetQuadFunction(const std::vector<mfem::QuadratureFunction*> &qf_array,
                        int component = -1);

   /// Helper function for visualizing 1D or 2D3V data
   void ExtrudeMeshAndSolution();

   /// Helper function for visualizing 1D data
   void Extrude1DMeshAndSolution();

   /// Helper function for visualization of 2D3V data
   void Extrude2D3VMeshAndSolution();

   /// Helper function to build the quadrature function from pieces

   /// Set a (checkerboard) solution when only the mesh is given
   void SetMeshSolution();

   /// Set the quadrature function representation producing a proxy grid function
   void SetQuadSolution(QuadSolution type = QuadSolution::LOR_ClosedGL);

   /// Switch the quadrature function representation and update the visualization
   void SwitchQuadSolution(QuadSolution type, VisualizationScene* vs);

   /// Get the current representation of quadrature solution
   inline QuadSolution GetQuadSolution() const { return quad_sol; }

   // Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
   // or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
   // product vector grid function of the same order.
   static std::unique_ptr<mfem::GridFunction>
   ProjectVectorFEGridFunction(std::unique_ptr<mfem::GridFunction> gf);

   void ProjectVectorFEGridFunction()
   { internal.grid_f = ProjectVectorFEGridFunction(std::move(internal.grid_f)); }

   /// Sets a new mesh and solution from another DataState object, and
   /// updates the given VisualizationScene pointer with the new data.
   ///
   /// Mesh space and grid function dimensions must both match the original
   /// dimensions of the current DataState. If there is a mismatch in either
   /// value, the function will return false, and the mesh/solution will not be
   /// updated.
   bool SetNewMeshAndSolution(DataState new_state,
                              VisualizationScene* vs);

   /// Updates the given VisualizationScene pointer with the new data
   /// of the given DataState object.
   /// @note: Use with caution when the update is compatible
   /// @see SetNewMeshAndSolution()
   static void ResetMeshAndSolution(DataState &ss, VisualizationScene* vs);
};

#endif // GLVIS_DATA_STATE_HPP
