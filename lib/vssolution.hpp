// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef GLVIS_VSSOLUTION
#define GLVIS_VSSOLUTION

#include "mfem.hpp"
using namespace mfem;

// Visualization header file

class VisualizationSceneSolution : public VisualizationSceneScalarData
{
protected:
   Vector *v_normals;
   GridFunction *rsol;

   int drawmesh, drawelems, drawnums, draworder;
   int displlist, linelist, lcurvelist;
   int bdrlist, drawbdr, draw_cp, cp_list;
   int e_nums_list, v_nums_list;
   int order_list, order_list_noarrow;

   void Init();

   void FindNewBox(double rx[], double ry[], double rval[]);

   void DrawCPLine(DenseMatrix &pointmat, Vector &values, Array<int> &ind);

   void GetRefinedDetJ(int i, const IntegrationRule &ir,
                       Vector &vals, DenseMatrix &tr);

   // redefined for vector solution
   virtual void GetRefinedValues(int i, const IntegrationRule &ir,
                                 Vector &vals, DenseMatrix &tr);
   virtual int GetRefinedValuesAndNormals(int i, const IntegrationRule &ir,
                                          Vector &vals, DenseMatrix &tr,
                                          DenseMatrix &normals);

   void DrawLevelCurves(Array<int> &RG, DenseMatrix &pointmat,
                        Vector &values, int sides, Array<double> &lvl,
                        int flat = 0);

   int GetAutoRefineFactor();

   // Used for drawing markers for element and vertex numbering
   double GetElementLengthScale(int k);

   // Rendering large numbers of text objects for element or vertex numbering is
   // slow.  Turn it off above some entity count.
   static const int MAX_RENDER_NUMBERING = 1000;

public:
   int shading, TimesToRefine, EdgeRefineFactor;

   int attr_to_show, bdr_attr_to_show;
   Array<int> el_attr_to_show, bdr_el_attr_to_show;

   VisualizationSceneSolution();
   VisualizationSceneSolution(Mesh &m, Vector &s, Vector *normals = NULL);

   virtual ~VisualizationSceneSolution();

   void SetGridFunction(GridFunction & u) { rsol = &u; }

   void NewMeshAndSolution(Mesh *new_m, Vector *new_sol,
                           GridFunction *new_u = NULL);

   virtual void SetNewScalingFromBox();
   virtual void FindNewBox(bool prepare);
   virtual void FindNewValueRange(bool prepare);
   virtual void FindNewBoxAndValueRange(bool prepare)
   { FindNewBox(prepare); }
   virtual void FindMeshBox(bool prepare);

   virtual void ToggleLogscale(bool print);
   virtual void UpdateLevelLines() { PrepareLevelCurves(); }
   virtual void UpdateValueRange(bool prepare);

   void PrepareWithNormals();
   void PrepareFlat();
   void PrepareFlat2();

   virtual void PrepareLines();
   void PrepareLines2();
   void PrepareLines3();

   virtual void Prepare();
   void PrepareLevelCurves();
   void PrepareLevelCurves2();

   void PrepareBoundary();

   void PrepareOrderingCurve();
   void PrepareOrderingCurve1(int list, bool arrows);

   void PrepareNumbering();
   void PrepareElementNumbering();
   void PrepareElementNumbering1();
   void PrepareElementNumbering2();
   void PrepareVertexNumbering();
   void PrepareVertexNumbering1();
   void PrepareVertexNumbering2();

   void PrepareCP();

   virtual void Draw();

   void ToggleDrawBdr()
   { drawbdr = !drawbdr; }

   virtual void ToggleDrawElems();

   void ToggleDrawMesh() { drawmesh = (drawmesh+1)%3; }

   // 0 - none, 1 - with arrows, 2 - no arrows
   void ToggleDrawOrdering() {  draworder = (draworder+1)%3; }

   // 0 - none, 1 - elements, 2 - vertices
   void ToggleDrawNumberings() { drawnums = (drawnums+1)%3; }

   virtual void SetShading(int, bool);
   void ToggleShading();

   void ToggleDrawCP() { draw_cp = !draw_cp; PrepareCP(); }

   virtual void SetRefineFactors(int, int);
   virtual void AutoRefine();
   virtual void ToggleAttributes(Array<int> &attr_list);
};

void DrawNumberedMarker(const double x[3], double dx, int n);

void DrawTriangle(const double pts[][3], const double cv[],
                  const double minv, const double maxv);

void DrawQuad(const double pts[][3], const double cv[],
              const double minv, const double maxv);

void DrawPatch(const DenseMatrix &pts, Vector &vals, DenseMatrix &normals,
               const int n, const Array<int> &ind, const double minv,
               const double maxv, const int normals_opt = 0);

#endif
