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

#ifndef GLVIS_VSVECTOR
#define GLVIS_VSVECTOR

class VisualizationSceneVector : public VisualizationSceneSolution
{
protected:

   Vector *solx, *soly;
   int vectorlist, displinelist, drawdisp, drawvector;

   GridFunction *VecGridF;

   void Init();

   virtual void GetRefinedValues(int i, const IntegrationRule &ir,
                                 Vector &vals, DenseMatrix &tr);

   double (*Vec2Scalar)(double, double);

   void DrawVector(double, double, double, double, double);

   double maxlen;

   Vector vc0;
   IsoparametricTransformation T0;

public:
   VisualizationSceneVector(Mesh &m, Vector &sx, Vector &sy);
   VisualizationSceneVector(GridFunction &vgf);

   void NewMeshAndSolution(GridFunction &vgf, int rescale = 0);

   virtual ~VisualizationSceneVector();

   void NPressed();
   void PrepareDisplacedMesh();
   virtual void PrepareLines()
   { VisualizationSceneSolution::PrepareLines(); PrepareDisplacedMesh(); }

   virtual void PrepareVectorField();
   void ToggleVectorField();

   void ToggleDisplacements()
   {
      drawdisp = (drawdisp+1)%4;
      if (drawdisp != 1)
         PrepareDisplacedMesh();
   }

   virtual void Draw();

   virtual void EventUpdateColors() { Prepare(); PrepareVectorField(); }

   int RefineFactor;

   double ArrowScale;

   void CycleVec2Scalar();
};

#endif
