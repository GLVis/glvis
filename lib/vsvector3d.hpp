// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
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
   int drawvector;

   enum class ScalarFunction
   {
      Invalid = -1,
      Min = -1,
      //----------
      Magnitude,
      Component_X,
      Component_Y,
      Component_Z,
      //----------
      Max
   } scal_func;
   static const char *scal_func_name[];

   double mesh_volume;
   double vector_h, vector_hh;
   int arrows_nl;
   gl3::GlDrawable vector_buf;
   gl3::GlDrawable displine_buf;

   mfem::GridFunction *VecGridF{};
   mfem::FiniteElementSpace *sfes{};

   void Init();

   void NewMeshAndSolution(mfem::Mesh *new_m, mfem::Mesh *new_mc,
                           mfem::Vector *new_sol_x, mfem::Vector *new_sol_y, mfem::Vector *new_sol_z,
                           mfem::GridFunction *new_u = nullptr);

   mfem::Array<int> vflevel;
   mfem::Array<double> dvflevel;

   int GetFunctionAutoRefineFactor() override;

   void ArrowsDrawOrNot(mfem::Array<int> l[], int nv, mfem::Vector & sol, int nl,
                        mfem::Array<double> & level);
   int ArrowDrawOrNot(double v, int nl, mfem::Array<double> & level);

   // key handlers
   static thread_local VisualizationSceneVector3d  *vsvector3d;
   int ianim, ianimd, ianimmax, drawdisp;

   static void KeyDPressed();
   static void KeyNPressed();
   static void KeyBPressed();
   static void KeyrPressed();
   static void KeyRPressed();
   static void KeyuPressed();
   static void KeyUPressed();
   static void KeywPressed();
   static void KeyWPressed();
   static void KeyvPressed();
   static void KeyVPressed();
   static void VectorKeyFPressed();

   void NPressed();

public:

   VisualizationSceneVector3d(Window &win);

   void NewMeshAndSolution(const DataState &s) override;

   void FindNewValueRange(bool prepare) override;

   virtual ~VisualizationSceneVector3d();

   std::string GetHelpString() const override;

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
