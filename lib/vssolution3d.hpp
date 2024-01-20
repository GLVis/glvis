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

   int drawmesh, drawelems, shading, draworder;
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

   GridFunction *GridF;

   void Init();

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

   int GetAutoRefineFactor();

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

   VisualizationSceneSolution3d();
   VisualizationSceneSolution3d(Mesh & m, Vector & s);

   void SetGridFunction (GridFunction *gf) { GridF = gf; }

   void NewMeshAndSolution(Mesh *new_m, Vector *new_sol,
                           GridFunction *new_u = NULL);

   virtual ~VisualizationSceneSolution3d();

   virtual std::string GetHelpString() const;

   virtual void FindNewBox(bool prepare);
   virtual void FindNewValueRange(bool prepare);

   virtual void PrepareRuler()
   { VisualizationSceneScalarData::PrepareRuler(false); }
   virtual void PrepareFlat();
   virtual void PrepareLines();
   virtual void Prepare();
   virtual void PrepareOrderingCurve();
   virtual void PrepareOrderingCurve1(gl3::GlDrawable& buf, bool arrows,
                                      bool color);
   virtual gl3::SceneInfo GetSceneObjs();

   virtual void glTF_Export();

   void ToggleDrawElems()
   { drawelems = !drawelems; Prepare(); }

   void ToggleDrawMesh();

   // 0 - none, 1 - no arrows (color), 2 - with arrows (color),
   //           3 - no arrows (black), 4 - with arrows (black)
   void ToggleDrawOrdering() { draworder = (draworder+1)%5; }

   void ToggleShading();
   int GetShading() { return shading; };
   virtual void SetShading(int, bool);
   virtual void SetRefineFactors(int, int);
   virtual void AutoRefine();
   virtual void ToggleAttributes(Array<int> &attr_list);

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
   virtual void EventUpdateColors();
   virtual void UpdateLevelLines()
   { PrepareLines(); PrepareCuttingPlaneLines(); }
   virtual void UpdateValueRange(bool prepare);

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
