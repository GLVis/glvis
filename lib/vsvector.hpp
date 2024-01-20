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

#ifndef GLVIS_VSVECTOR_HPP
#define GLVIS_VSVECTOR_HPP

#include "mfem.hpp"
#include "gl/types.hpp"
#include "vssolution.hpp"
using namespace mfem;

class VisualizationSceneVector : public VisualizationSceneSolution
{
protected:

   Vector *solx, *soly;
   int drawdisp, drawvector;

   gl3::GlDrawable vector_buf;
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

   virtual std::string GetHelpString() const;

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

   virtual gl3::SceneInfo GetSceneObjs();

   virtual void glTF_Export();

   virtual void EventUpdateColors() { Prepare(); PrepareVectorField(); }

   // refinement factor for the vectors
   int RefineFactor;

   double ArrowScale;

   void CycleVec2Scalar(int print = 0);
};

#endif
