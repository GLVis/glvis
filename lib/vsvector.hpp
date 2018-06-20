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

#ifndef GLVIS_VSVECTOR
#define GLVIS_VSVECTOR

#include "mfem.hpp"
#include "aux_gl3.hpp"
using namespace mfem;

class VisualizationSceneVector : public VisualizationSceneSolution
{
protected:

   Vector *solx, *soly;
   int vectorlist, displinelist, drawdisp, drawvector;

   gl3::GlDrawable displine_buf;
   GridFunction *VecGridF;

   void Init();

   virtual void GetRefinedValues(int i, const IntegrationRule &ir,
                                 Vector &vals, DenseMatrix &tr);
   virtual int GetRefinedValuesAndNormals(int i, const IntegrationRule &ir,
                                          Vector &vals, DenseMatrix &tr,
                                          DenseMatrix &normals);

   double (*Vec2Scalar)(double, double);

   void DrawVector(double, double, double, double, double);

   double maxlen;

   Vector vc0;
   IsoparametricTransformation T0;

public:
   VisualizationSceneVector(Mesh &m, Vector &sx, Vector &sy);
   VisualizationSceneVector(GridFunction &vgf);

   void NewMeshAndSolution(GridFunction &vgf);

   virtual ~VisualizationSceneVector();

   void NPressed();
   void PrepareDisplacedMesh();
   virtual void PrepareLines()
   { VisualizationSceneSolution::PrepareLines(); PrepareDisplacedMesh(); }

   virtual void ToggleDrawElems();

   virtual void PrepareVectorField();
   void ToggleVectorField();

   void ToggleDisplacements()
   {
      drawdisp = (drawdisp+1)%4;
      if (drawdisp != 1)
      {
         PrepareDisplacedMesh();
      }
   }

   virtual void Draw();

   virtual void EventUpdateColors() { Prepare(); PrepareVectorField(); }

   // refinement factor for the vectors
   int RefineFactor;

   double ArrowScale;

   void CycleVec2Scalar(int print = 0);
};

#endif
