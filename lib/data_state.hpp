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

#ifndef GLVIS_DATA_STATE_HPP
#define GLVIS_DATA_STATE_HPP

#include <string>
#include <memory>
#include <vector>

#include <mfem.hpp>

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
      std::unique_ptr<mfem::DataCollection> data_coll;
   } internal;

   FieldType type {FieldType::UNKNOWN};
   QuadSolution quad_sol {QuadSolution::NONE};

   void SetGridFunctionSolution(int component = -1);
   void SetQuadFunctionSolution(int component = -1);

   /// Updates the given VisualizationScene pointer with the new data
   /// of the given DataState object.
   /// @note: Use with caution when the update is compatible
   /// @see SetNewMeshAndSolution()
   void ResetMeshAndSolution(DataState &ss, VisualizationScene* vs);

public:
   mfem::Vector sol, solu, solv, solw, normals;
   const std::unique_ptr<mfem::Mesh> &mesh{internal.mesh};
   const std::unique_ptr<mfem::Mesh> &mesh_quad{internal.mesh_quad};
   const std::unique_ptr<mfem::GridFunction> &grid_f{internal.grid_f};
   const std::unique_ptr<mfem::QuadratureFunction> &quad_f{internal.quad_f};
   const std::unique_ptr<mfem::DataCollection> &data_coll{internal.data_coll};

   std::string keys;
   bool fix_elem_orient{false};
   bool save_coloring{false};
   bool keep_attr{false};

   DataState() = default;
   DataState(DataState &&ss) { *this = std::move(ss); }
   DataState& operator=(DataState &&ss);

   /// Get type of the contained data
   inline FieldType GetType() const { return type; }

   /// Set a mesh (plain pointer version)
   /** Note that ownership is passed from the caller.
       @see SetMesh(std::unique_ptr<mfem::Mesh> &&pmesh) */
   void SetMesh(mfem::Mesh *mesh);

   /// Set a mesh (unique pointer version)
   /** Sets the mesh and resets grid/quadrature functions if they do not use
       the same one. */
   void SetMesh(std::unique_ptr<mfem::Mesh> &&pmesh);

   /// Set a grid function (plain pointer version)
   /** Note that ownership is passed from the caller.
       @see SetGridFunction(std::unique_ptr<mfem::GridFunction> &&, int ) */
   void SetGridFunction(mfem::GridFunction *gf, int component = -1);

   /// Set a grid function (unique pointer version)
   /** Sets the grid function or its component (-1 means all components). */
   void SetGridFunction(std::unique_ptr<mfem::GridFunction> &&pgf,
                        int component = -1);

   /// Set a quadrature function (plain pointer version)
   /** Note that ownership is passed from the caller.
       @see SetQuadFunction(std::unique_ptr<mfem::QuadFunction> &&, int ) */
   void SetQuadFunction(mfem::QuadratureFunction *qf, int component = -1);

   /// Set a quadrature function (unique pointer version)
   /** Sets the quadrature function or its component (-1 means all components). */
   void SetQuadFunction(std::unique_ptr<mfem::QuadratureFunction> &&pqf,
                        int component = -1);

   /// Set a quadrature function from pieces
   /** Serializes the pieces of a quadrature function and sets it or its
       component (-1 means all components) */
   void SetQuadFunction(const std::vector<mfem::QuadratureFunction*> &qf_array,
                        int component = -1);

   /// Set a data collection field
   /** Sets the mesh and optionally a grid or quadrature function from the
       provided data collection.
       @param dc        data collection
       @param ti        cycle to load
       @param field     name of the (Q-)field to load (NULL for mesh only)
       @param quad      if true, Q-field is loaded, otherwise a regular field
       @param component component of the field (-1 means all components) */
   void SetDataCollectionField(mfem::DataCollection *dc, int ti,
                               const char *field = NULL, bool quad = false, int component = -1);

   /// Helper function for visualizing 1D or 2D3V data
   void ExtrudeMeshAndSolution();

   /// Helper function for visualizing 1D data
   void Extrude1DMeshAndSolution();

   /// Helper function for visualization of 2D3V data
   void Extrude2D3VMeshAndSolution();

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
};

#endif // GLVIS_DATA_STATE_HPP
