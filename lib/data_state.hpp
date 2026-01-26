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

#include <map>
#include <string>
#include <memory>
#include <vector>
#include <utility>

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

   enum class ComplexSolution
   {
      NONE = -1,
      MIN = -1,
      //----------
      Magnitude,
      Phase,
      Real,
      Imag,
      //----------
      MAX
   };

   // Class used for storing offsets and map of DOFs for each rank
   class Offset
   {
      std::map<std::pair<int, int>, int> dof;
   public:
      int nelems, nedges, nverts;
#ifdef GLVIS_DEBUG
      // in debug mode, we store the element centers
      // to be able to compare them with the ones of the global mesh,
      // as it could depend on the way the global mesh is constructed
      // from the array of 'local' ones.
      struct xy {double x,y;};
      std::map<std::pair<int, int>, xy> exy_map;
#endif
      Offset() = default;
      int& operator[](const std::pair<int, int> &key) { return dof[key]; }
      const int& operator[](const std::pair<int, int> &key) const { return dof.at(key); }
   };
   using Offsets = std::vector<Offset>;

private:
   struct
   {
      std::unique_ptr<mfem::Vector> sol, solx, soly, solz;
      std::unique_ptr<mfem::Vector> normals;
      std::unique_ptr<mfem::Mesh> mesh;
      std::unique_ptr<mfem::Mesh> mesh_quad;
      std::unique_ptr<mfem::GridFunction> grid_f;
      std::unique_ptr<mfem::ComplexGridFunction> cgrid_f;
      std::unique_ptr<mfem::QuadratureFunction> quad_f;
      std::unique_ptr<mfem::DataCollection> data_coll;
      std::unique_ptr<Offsets> offsets;
   } internal;

   FieldType type {FieldType::UNKNOWN};
   ComplexSolution cmplx_sol {ComplexSolution::NONE};
   QuadSolution quad_sol {QuadSolution::NONE};

   void SetGridFunctionSolution(int component = -1);
   void SetComplexFunctionSolution(int component = -1);
   void SetQuadFunctionSolution(int component = -1);

   /// Compute the dofs offsets from the grid function vector
   void ComputeDofsOffsets(std::vector<mfem::GridFunction*> &gf_array);

public:
   const std::unique_ptr<mfem::Vector> &sol{internal.sol};
   const std::unique_ptr<mfem::Vector> &solx{internal.solx};
   const std::unique_ptr<mfem::Vector> &soly{internal.soly};
   const std::unique_ptr<mfem::Vector> &solz{internal.solz};
   const std::unique_ptr<mfem::Vector> &normals{internal.normals};
   const std::unique_ptr<mfem::Mesh> &mesh{internal.mesh};
   const std::unique_ptr<mfem::Mesh> &mesh_quad{internal.mesh_quad};
   const std::unique_ptr<mfem::GridFunction> &grid_f{internal.grid_f};
   const std::unique_ptr<mfem::ComplexGridFunction> &cgrid_f{internal.cgrid_f};
   const std::unique_ptr<mfem::QuadratureFunction> &quad_f{internal.quad_f};
   const std::unique_ptr<mfem::DataCollection> &data_coll{internal.data_coll};
   const std::unique_ptr<Offsets> &offsets{internal.offsets};

   std::string keys;
   bool fix_elem_orient{false};
   bool save_coloring{false};
   bool keep_attr{false};
   double cmplx_phase{};
   double cmplx_mag_max{0.};

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

   /// Set scalar data
   void SetScalarData(mfem::Vector sol);

   /// Set normals
   void SetNormals(mfem::Vector normals);

   /// Set 2D vector data
   void SetVectorData(mfem::Vector solx, mfem::Vector soly);

   /// Set 3D vector data
   void SetVectorData(mfem::Vector solx, mfem::Vector soly, mfem::Vector solz);

   /// Set a grid function (plain pointer version)
   /** Note that ownership is passed from the caller.
       @see SetGridFunction(std::unique_ptr<mfem::GridFunction> &&, int ) */
   void SetGridFunction(mfem::GridFunction *gf, int component = -1);

   /// Set a grid function (unique pointer version)
   /** Sets the grid function or its component (-1 means all components). */
   void SetGridFunction(std::unique_ptr<mfem::GridFunction> &&pgf,
                        int component = -1);

   /// Set a grid function from pieces
   /** Serializes the pieces of a grid function and sets it or its
       component (-1 means all components) */
   void SetGridFunction(std::vector<mfem::GridFunction*> &gf_array,
                        int num_pieces, int component = -1);

   /// Set a complex grid function (plain pointer version)
   /** Note that ownership is passed from the caller.
       @see SetCmplxGridFunction(std::unique_ptr<mfem::ComplexGridFunction> &&, int ) */
   void SetCmplxGridFunction(mfem::ComplexGridFunction *gf, int component = -1);

   /// Set a complex grid function (unique pointer version)
   /** Sets the complex grid function or its component (-1 means all
       components). */
   void SetCmplxGridFunction(std::unique_ptr<mfem::ComplexGridFunction> &&pgf,
                             int component = -1);

   /// Set a complex grid function from pieces
   /** Serializes the pieces of a complex grid function and sets it or its
       component (-1 means all components) */
   void SetCmplxGridFunction(const std::vector<mfem::ComplexGridFunction*>
                             &cgf_array, int component = -1);

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

   /// Set the complex function representation producing a proxy grid function
   void SetComplexSolution(ComplexSolution type = ComplexSolution::Magnitude,
                           bool print = true);

   /// Switch the complex function representation
   void SwitchComplexSolution(ComplexSolution type, bool print = true);

   /// Get the current representation of complex solution
   inline ComplexSolution GetComplexSolution() const { return cmplx_sol; }

   /// Set the quadrature function representation producing a proxy grid function
   void SetQuadSolution(QuadSolution type = QuadSolution::LOR_ClosedGL);

   /// Switch the quadrature function representation
   void SwitchQuadSolution(QuadSolution type);

   /// Get the current representation of quadrature solution
   inline QuadSolution GetQuadSolution() const { return quad_sol; }

   // Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
   // or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
   // product vector grid function of the same order.
   static std::unique_ptr<mfem::GridFunction>
   ProjectVectorFEGridFunction(std::unique_ptr<mfem::GridFunction> gf);

   static std::unique_ptr<mfem::ComplexGridFunction>
   ProjectVectorFEGridFunction(std::unique_ptr<mfem::ComplexGridFunction> gf);

   void ProjectVectorFEGridFunction();
};

#endif // GLVIS_DATA_STATE_HPP
