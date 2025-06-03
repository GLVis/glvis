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

#ifndef GLVIS_VSSOLUTION_HPP
#define GLVIS_VSSOLUTION_HPP

#include "gl/types.hpp"
#include "vsdata.hpp"

#include <mfem.hpp>

// Visualization header file

class VisualizationSceneSolution : public VisualizationSceneScalarData
{
protected:
   mfem::Vector *v_normals;
   mfem::GridFunction *rsol;

   int drawmesh, drawelems, draworder;
   enum class Numbering { NONE, ELEMENTS, EDGES, VERTICES, DOFS, MAX };
   Numbering drawnums;
   int drawbdr, draw_cp;

   int refine_func = 0;

   double minv_sol, maxv_sol;
   bool have_sol_range = false; // true when minv_sol and maxv_sol are set

   int TimesToRefine, EdgeRefineFactor;

   gl3::GlDrawable disp_buf;

   bool e_nums_buf_ready = false;
   bool v_nums_buf_ready = false;
   bool f_nums_buf_ready = false;
   bool d_nums_buf_ready = false;
   gl3::GlDrawable e_nums_buf, v_nums_buf, f_nums_buf, d_nums_buf;

   gl3::GlDrawable lcurve_buf;
   gl3::GlDrawable line_buf;
   gl3::GlDrawable bdr_buf;
   gl3::GlDrawable cp_buf;

   gl3::GlDrawable order_buf;
   gl3::GlDrawable order_noarrow_buf;

   void Init();

   void FindNewBox(double rx[], double ry[], double rval[]);

   void DrawCPLine(gl3::GlBuilder& bld,
                   mfem::DenseMatrix &pointmat, mfem::Vector &values, mfem::Array<int> &ind);

   void GetRefinedDetJ(int i, const mfem::IntegrationRule &ir,
                       mfem::Vector &vals, mfem::DenseMatrix &tr);

   // redefined for vector solution
   virtual void GetRefinedValues(int i, const mfem::IntegrationRule &ir,
                                 mfem::Vector &vals, mfem::DenseMatrix &tr);
   virtual int GetRefinedValuesAndNormals(int i, const mfem::IntegrationRule &ir,
                                          mfem::Vector &vals, mfem::DenseMatrix &tr,
                                          mfem::DenseMatrix &normals);

   void DrawLevelCurves(gl3::GlBuilder& buf, mfem::Array<int> &RG,
                        mfem::DenseMatrix &pointmat,
                        mfem::Vector &values, int sides, mfem::Array<double> &lvl,
                        int flat = 0);

   int GetFunctionAutoRefineFactor() override;

   // Used for drawing markers for element and vertex numbering
   double GetElementLengthScale(int k);

public:
   int attr_to_show, bdr_attr_to_show;
   mfem::Array<int> el_attr_to_show, bdr_el_attr_to_show;

   VisualizationSceneSolution();
   VisualizationSceneSolution(mfem::Mesh &m, mfem::Vector &s,
                              mfem::Mesh *mc = nullptr,
                              mfem::Vector *normals = nullptr);

   virtual ~VisualizationSceneSolution();

   std::string GetHelpString() const override;

   void SetGridFunction(mfem::GridFunction & u);

   void NewMeshAndSolution(mfem::Mesh *new_m, mfem::Mesh *new_mc,
                           mfem::Vector *new_sol,
                           mfem::GridFunction *new_u = nullptr);

   void SetNewScalingFromBox() override;
   void FindNewBox(bool prepare) override;
   void FindNewValueRange(bool prepare) override;
   void FindNewBoxAndValueRange(bool prepare) override
   { FindNewBox(prepare); }
   void FindMeshBox(bool prepare) override;

   void ToggleLogscale(bool print) override;
   void EventUpdateBackground() override;
   void EventUpdateColors() override;
   void UpdateLevelLines()  override { PrepareLevelCurves(); }
   void UpdateValueRange(bool prepare) override;

   void PrepareWithNormals();
   void PrepareFlat();
   void PrepareFlat2();

   void PrepareLines() override;
   void PrepareLines2();
   void PrepareLines3();

   void Prepare() override;
   void PrepareLevelCurves();
   void PrepareLevelCurves2();

   void PrepareBoundary();

   void PrepareOrderingCurve();
   void PrepareOrderingCurve1(gl3::GlDrawable& buf, bool arrows, bool color);

   void PrepareNumbering(bool invalidate = true);
   void PrepareElementNumbering();
   void PrepareElementNumbering1();
   void PrepareElementNumbering2();
   void PrepareVertexNumbering();
   void PrepareVertexNumbering1();
   void PrepareVertexNumbering2();
   void PrepareEdgeNumbering();
   void PrepareDofNumbering();

   void PrepareCP();

   gl3::SceneInfo GetSceneObjs() override;

   void glTF_ExportBoundary(glTF_Builder &bld,
                            glTF_Builder::buffer_id buffer,
                            glTF_Builder::material_id black_mat);
   void glTF_Export() override;

   void ToggleDrawBdr();

   virtual void ToggleDrawElems();

   void ToggleDrawMesh() { drawmesh = (drawmesh+1)%3; }

   // 0 - none, 1 - no arrows (color), 2 - with arrows (color),
   //           3 - no arrows (black), 4 - with arrows (black)
   void ToggleDrawOrdering() { draworder = (draworder+1)%5; }

   // 0 - none, 1 - elements, 2 - edges, 3 - vertices, 4 - DOFs
   void ToggleDrawNumberings()
   {
      drawnums = (Numbering) (((int)drawnums + 1) % (int)Numbering::MAX);
      PrepareNumbering(false);
   }

   void SetShading(Shading, bool) override;
   void ToggleShading() override;

   void ToggleDrawCP() { draw_cp = !draw_cp; PrepareCP(); }

   void ToggleRefinements();
   void ToggleRefinementFunction();

   void SetRefineFactors(int, int) override;
   void AutoRefine() override;
   void ToggleAttributes(mfem::Array<int> &attr_list) override;

   virtual void SetDrawMesh(int i) { drawmesh = i % 3; }
   virtual int GetDrawMesh() { return drawmesh; }
};

void DrawNumberedMarker(gl3::GlDrawable& buff, const double x[3], double dx,
                        int n);

#endif // GLVIS_VSSOLUTION_HPP
