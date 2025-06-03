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

#ifndef GLVIS_VSVECTOR_3D_HPP
#define GLVIS_VSVECTOR_3D_HPP

#include <mfem.hpp>
#include "gl/types.hpp"
#include "vssolution3d.hpp"

class VisualizationSceneVector3d : public VisualizationSceneSolution3d
{
protected:

   mfem::Vector *solx, *soly, *solz;
   int drawvector, scal_func;
   double mesh_volume;
   gl3::GlDrawable vector_buf;
   gl3::GlDrawable displine_buf;

   mfem::GridFunction *VecGridF;
   mfem::FiniteElementSpace *sfes;

   void Init();

   mfem::Array<int> vflevel;
   mfem::Array<double> dvflevel;

   int GetFunctionAutoRefineFactor() override;

public:
   int ianim, ianimd, ianimmax, drawdisp;

   VisualizationSceneVector3d(mfem::Mesh & m, mfem::Vector & sx, mfem::Vector & sy,
                              mfem::Vector & sz,
                              mfem::Mesh *mc = nullptr);
   VisualizationSceneVector3d(mfem::GridFunction &vgf, mfem::Mesh *mc = nullptr);

   void NewMeshAndSolution(mfem::Mesh *new_m, mfem::Mesh *new_mc,
                           mfem::GridFunction *new_v);

   virtual ~VisualizationSceneVector3d();

   std::string GetHelpString() const override;

   void NPressed();
   void PrepareFlat() override;
   void Prepare() override;
   void PrepareLines() override;

   void PrepareFlat2();
   void PrepareLines2();

   void DrawVector (gl3::GlDrawable& buf,
                    int type, double v0, double v1, double v2,
                    double sx, double sy, double sz, double s);
   virtual void PrepareVectorField();
   void PrepareDisplacedMesh();
   void ToggleVectorField(int i);

   void SetScalarFunction();
   void ToggleScalarFunction();

   void PrepareCuttingPlane() override;

   void ToggleDisplacements() {drawdisp = (drawdisp+1)%2;};

   gl3::SceneInfo GetSceneObjs() override;

   void EventUpdateColors() override
   { Prepare(); PrepareVectorField(); PrepareCuttingPlane(); };

   void ToggleVectorFieldLevel(int v);
   void AddVectorFieldLevel();
   void RemoveVectorFieldLevel();
};

#endif
