// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.googlecode.com.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#ifndef GLVIS_VSSOLUTION
#define GLVIS_VSSOLUTION

// Visualization header file

class VisualizationSceneSolution : public VisualizationSceneScalarData
{
protected:
   GridFunction * rsol;

   int drawmesh, drawelems;
   int displlist, linelist, lcurvelist;
   int bdrlist, drawbdr, draw_cp, cp_list;

   void Init();

   void FindNewBox(double rx[], double ry[], double rval[]);
   void FixValueRange();
   void NewZRange(double z1, double z2);

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

public:
   int shading, TimesToRefine, EdgeRefineFactor;

   int attr_to_show;
   Array<int> el_attr_to_show;

   VisualizationSceneSolution();
   VisualizationSceneSolution (Mesh & m, Vector & s);

   virtual ~VisualizationSceneSolution();

   void SetGridFunction(GridFunction & u) { rsol = &u; };

   void NewMeshAndSolution(Mesh *new_m, Vector *new_sol,
                           GridFunction *new_u = NULL, int rescale = 0);

   virtual void SetNewScalingFromBox();
   virtual void FindNewBox();

   virtual void UpdateLevelLines() { PrepareLevelCurves(); };
   virtual void UpdateValueRange();

   void PrepareFlat();
   void PrepareFlat2();

   virtual void PrepareLines();
   void PrepareLines2();
   void PrepareLines3();

   virtual void Prepare();
   void PrepareLevelCurves();
   void PrepareLevelCurves2();
   void DefaultLevelLines() { SetLevelLines (minv, maxv, 15); };

   void PrepareBoundary();

   void PrepareCP();

   virtual void Draw();

   void ToggleDrawBdr()
   { drawbdr = !drawbdr; }

   void ToggleDrawElems()
   {
      drawelems = (drawelems+3)%4;
      if (drawelems != 0 && shading == 2)
      {
         FindNewBox();
         PrepareAxes();
         PrepareLines();
         Prepare();
         DefaultLevelLines();
         PrepareLevelCurves();
         PrepareCP();
      }
   }

   void ToggleDrawMesh() { drawmesh = (drawmesh+1)%3; }

   virtual void SetShading(int);
   void ToggleShading();

   void ToggleDrawCP() { draw_cp = !draw_cp; PrepareCP(); }

   virtual void SetRefineFactors(int, int);
   virtual void ToggleAttributes(Array<int> &attr_list);
};


void DrawTriangle(const double pts[][3], const double cv[],
                  const double minv, const double maxv);

void DrawQuad(const double pts[][3], const double cv[],
              const double minv, const double maxv);

void DrawPatch(const DenseMatrix &pts, Vector &vals, DenseMatrix &normals,
               const int n, const Array<int> &ind, const double minv,
               const double maxv, const int normals_opt = 0);

#endif
