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

   void GetRefinedValues(int i, const IntegrationRule &ir,
                         Vector &vals, DenseMatrix &tr) override;
   int GetRefinedValuesAndNormals(int i, const IntegrationRule &ir,
                                  Vector &vals, DenseMatrix &tr,
                                  DenseMatrix &normals) override;

   double (*Vec2Scalar)(double, double);

   void DrawVector(double, double, double, double, double);

   double maxlen;

   Vector vc0;
   IsoparametricTransformation T0;

   int GetFunctionAutoRefineFactor() override;

public:
   VisualizationSceneVector(Mesh &m, Vector &sx, Vector &sy, Mesh *mc = NULL);
   VisualizationSceneVector(GridFunction &vgf);

   void NewMeshAndSolution(GridFunction &vgf, Mesh *mc = NULL);

   virtual ~VisualizationSceneVector();

   std::string GetHelpString() const override;

   void NPressed();
   void PrepareDisplacedMesh();
   void PrepareLines() override
   { VisualizationSceneSolution::PrepareLines(); PrepareDisplacedMesh(); }

   void ToggleDrawElems() override;

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

   gl3::SceneInfo GetSceneObjs() override;

   void glTF_Export() override;

   void EventUpdateColors() override { Prepare(); PrepareVectorField(); }

   // refinement factor for the vectors
   int RefineFactor;

   double ArrowScale;

   void CycleVec2Scalar(int print = 0);
};

#endif
