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

#ifndef GLVIS_VSSOLUTION_3D_HPP
#define GLVIS_VSSOLUTION_3D_HPP

#include "mfem.hpp"
#include "gl/types.hpp"
#include "vsdata.hpp"
#include <map>
using namespace mfem;

class VisualizationSceneSolution3d : public VisualizationSceneScalarData
{
protected:

   int drawmesh, drawelems, draworder;
   int cplane;
   int cp_drawmesh, cp_drawelems, drawlsurf;
   // Algorithm used to draw the cutting plane when shading is 2 and cplane is 1
   // 0 - slower, more accurate algorithm for curved meshes (default)
   // 1 - faster algorithm suitable for meshes with planar faces
   int cp_algo;

   gl3::GlDrawable disp_buf;
   gl3::GlDrawable line_buf;
   gl3::GlDrawable cplane_buf;
   gl3::GlDrawable cplines_buf;
   gl3::GlDrawable lsurf_buf;
   gl3::GlDrawable other_buf;
   gl3::GlDrawable order_buf, order_noarrow_buf;

   double *node_pos;

   int nlevels;
   Array<double> levels;

   GridFunction *GridF{};

   void Init();

   void NewMeshAndSolution(Mesh *new_m, Mesh *new_mc, Vector *new_sol,
                           GridFunction *new_u);

   void GetFaceNormals(const int FaceNo, const int side,
                       const IntegrationRule &ir, DenseMatrix &normals);

   void DrawRefinedSurf (int n, double *points, int elem, int func,
                         int part = -1);
   void DrawRefinedSurf (int n, DenseMatrix &pointmat,
                         Vector &values, Array<int> &RefGeoms);
   void DrawRefinedSurfLevelLines (int n, DenseMatrix &pointmat,
                                   Vector &values, Array<int> &RefGeoms);
   void DrawRefinedSurfEdges (int n, DenseMatrix &pointmat,
                              Vector &values, Array<int> &RefEdges,
                              int part = -1);
   void DrawBdrElCoarseSurfEdges(gl3::GlBuilder &line, int be,
                                 DenseMatrix &pointmat, const IntegrationRule *ir = NULL,
                                 Array<int> *idxs = NULL);
   void DrawFaceCoarseSurfEdges(gl3::GlBuilder &line, int f, DenseMatrix &pointmat,
                                const IntegrationRule *ir = NULL, Array<int> *idxs = NULL);
   void DrawCoarseSurfEdges(gl3::GlBuilder &line, int f, int e1, int e2,
                            DenseMatrix &pointmat, const IntegrationRule *ir = NULL,
                            Array<int> *idxs = NULL);
   void LiftRefinedSurf (int n, DenseMatrix &pointmat,
                         Vector &values, int *RG);
   void DrawTetLevelSurf(gl3::GlDrawable& target, const DenseMatrix &verts,
                         const Vector &vals,
                         const int *ind, const Array<double> &levels,
                         const DenseMatrix *grad = NULL);

   static int GetPyramidFaceSplits(const Array<bool> &quad_diag,
                                   const Array<int> &faces,
                                   const Array<int> &ofaces);
   void DrawRefinedPyramidLevelSurf(gl3::GlDrawable& target,
                                    const DenseMatrix &verts,
                                    const Vector &vals, const int *RG,
                                    const int np, const int face_splits,
                                    const DenseMatrix *grad = NULL);

   static int GetWedgeFaceSplits(const Array<bool> &quad_diag,
                                 const Array<int> &faces,
                                 const Array<int> &ofaces);
   void DrawRefinedWedgeLevelSurf(gl3::GlDrawable& target,
                                  const DenseMatrix &verts,
                                  const Vector &vals, const int *RG,
                                  const int np, const int face_splits,
                                  const DenseMatrix *grad = NULL);

   static int GetHexFaceSplits(const Array<bool> &quad_diag,
                               const Array<int> &faces,
                               const Array<int> &ofaces);
   void DrawRefinedHexLevelSurf(gl3::GlDrawable& target,
                                const DenseMatrix &verts,
                                const Vector &vals, const int *RG,
                                const int nh, const int face_splits,
                                const DenseMatrix *grad = NULL);

   int GetFunctionAutoRefineFactor() override;

   bool CheckPositions(Array<int> &vertices) const
   {
      int n = 0;
      for (int j = 0; j < vertices.Size(); j++)
      {
         if (node_pos[vertices[j]] >= 0.0) { n++; }
      }
      return (n < vertices.Size());
   }

public:
   int TimesToRefine;
   double FaceShiftScale;

   Array<int> bdr_attr_to_show;

   VisualizationSceneSolution3d(Window &win, bool init = true);

   void NewMeshAndSolution(Mesh *new_m, Mesh *new_mc, Vector *new_sol)
   { NewMeshAndSolution(new_m, new_mc, new_sol, NULL); }

   void NewMeshAndSolution(GridFunction *new_u, Mesh *new_mc);

   virtual ~VisualizationSceneSolution3d();

   std::string GetHelpString() const override;

   void FindNewBox(bool prepare) override;
   void FindNewValueRange(bool prepare) override;

   void PrepareRuler() override
   { VisualizationSceneScalarData::PrepareRuler(false); }
   virtual void PrepareFlat();
   void PrepareLines() override;
   void Prepare() override;
   virtual void PrepareOrderingCurve();
   virtual void PrepareOrderingCurve1(gl3::GlDrawable& buf, bool arrows,
                                      bool color);
   gl3::SceneInfo GetSceneObjs() override;

   void glTF_Export() override;

   void ToggleDrawElems()
   { drawelems = !drawelems; Prepare(); }

   void ToggleDrawMesh();

   // 0 - none, 1 - no arrows (color), 2 - with arrows (color),
   //           3 - no arrows (black), 4 - with arrows (black)
   void ToggleDrawOrdering() { draworder = (draworder+1)%5; }

   void SetShading(Shading, bool) override;
   void ToggleShading() override;
   void SetRefineFactors(int, int) override;
   void AutoRefine() override;
   void ToggleAttributes(Array<int> &attr_list) override;

   void FindNodePos();

   void CuttingPlaneFunc (int type);
   // func: 0 - draw surface, 1 - draw level lines
   void CutRefinedElement(gl3::GlDrawable& target,
                          const DenseMatrix &verts, const Vector &vert_dist,
                          const Vector &vals, const Geometry::Type geom,
                          const int *elems, int num_elems, int func);
   void CutRefinedFace(gl3::GlDrawable& target,
                       const DenseMatrix &verts, const Vector &vert_dist,
                       const Vector &vals, const Geometry::Type geom,
                       const int *faces, int num_faces);
   void CPPrepare();
   void CPMoved();
   void PrepareFlat2();
   void PrepareLines2();
   virtual void PrepareCuttingPlane();
   void PrepareCuttingPlane2();
   void PrepareCuttingPlaneLines();
   void PrepareCuttingPlaneLines2();
   void PrepareLevelSurf();
   void ToggleCuttingPlane();
   void ToggleCPDrawElems();
   void ToggleCPDrawMesh();
   void ToggleCPAlgorithm();
   void MoveLevelSurf(int);
   void NumberOfLevelSurf(int);
   void EventUpdateColors() override;
   void UpdateLevelLines() override
   { PrepareLines(); PrepareCuttingPlaneLines(); }
   void UpdateValueRange(bool prepare) override;

   virtual void SetDrawMesh(int i)
   {
      if (drawmesh != i % 3)
      {
         drawmesh = i % 3;
         PrepareLines();
      }
   }
   virtual int GetDrawMesh() { return drawmesh; }
};

#endif
