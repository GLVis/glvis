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

#ifndef GLVIS_VSVECTOR_3D
#define GLVIS_VSVECTOR_3D

#include "mfem.hpp"
#include "aux_gl3.hpp"
using namespace mfem;

class VisualizationSceneVector3d : public VisualizationSceneSolution3d
{
protected:

   Vector *solx, *soly, *solz;
   int vectorlist, displinelist, drawvector, scal_func;
   gl3::GlDrawable vector_buf;
   gl3::GlDrawable displine_buf;

   GridFunction *VecGridF;
   FiniteElementSpace *sfes;

   void Init();

   Array<int> vflevel;
   Array<double> dvflevel;

public:
   int ianim, ianimd, ianimmax, drawdisp;

   VisualizationSceneVector3d(Mesh & m, Vector & sx, Vector & sy, Vector & sz);
   VisualizationSceneVector3d (GridFunction &vgf);

   void NewMeshAndSolution(Mesh *new_m, GridFunction *new_v);

   virtual ~VisualizationSceneVector3d();

   void NPressed();
   virtual void PrepareFlat();
   virtual void Prepare();
   virtual void PrepareLines();

   void PrepareFlat2();
   void PrepareLines2();

   void DrawVector (gl3::GlBuilder& builder,
                    int type, double v0, double v1, double v2,
                    double sx, double sy, double sz, double s);
   virtual void PrepareVectorField();
   void PrepareDisplacedMesh();
   void ToggleVectorField(int i);

   void SetScalarFunction();
   void ToggleScalarFunction();

   virtual void PrepareCuttingPlane();

   void ToggleDisplacements() {drawdisp = (drawdisp+1)%2;};

   virtual void Draw();

   virtual void EventUpdateColors()
   { Prepare(); PrepareVectorField(); PrepareCuttingPlane(); };

   void ToggleVectorFieldLevel(int v);
   void AddVectorFieldLevel();
   void RemoveVectorFieldLevel();
};

#endif
