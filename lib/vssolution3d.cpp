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

#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <X11/keysym.h>

#include "mfem.hpp"
#include "visual.hpp"


VisualizationSceneSolution3d *vssol3d;

// Definitions of some more keys

static void Solution3dKeyHPressed()
{
   cout << endl
        << "+------------------------------------+" << endl
        << "| Keys                               |" << endl
        << "+------------------------------------+" << endl
        << "| a -  Displays/Hides the axes       |" << endl
        << "| A -  Turns antialiasing on/off     |" << endl
        << "| c -  Displays/Hides the colorbar   |" << endl
        << "| e -  Displays/Hides the elements   |" << endl
        << "| E -  Toggle the elements in the CP |" << endl
        << "| f -  Smooth/Flat/discont. shading  |" << endl
        << "| g -  Toggle background             |" << endl
        << "| h -  Displays help menu            |" << endl
        << "| i -  Toggle cutting plane          |" << endl
        << "| j -  Turn on/off perspective       |" << endl
        << "| k/K  Adjust the transparency level |" << endl
        << "| ,/<  Adjust color transparency     |" << endl
        << "| l -  Turns on/off the light        |" << endl
        << "| L -  Toggle logarithmic scale      |" << endl
        << "| m -  Displays/Hides the mesh       |" << endl
        << "| M -  Toggle the mesh in the CP     |" << endl
        << "| o/O  (De)refine elem, disc shading |" << endl
        << "| p -  Cycle through color palettes  |" << endl
        << "| P -  Print to PostScript file      |" << endl
        << "| q -  Quits                         |" << endl
        << "| r -  Reset the plot to 3D view     |" << endl
        << "| R -  Reset the plot to 2D view     |" << endl
        << "| s -  Turn on/off unit cube scaling |" << endl
        << "| S -  Take snapshot/Record a movie  |" << endl
        << "| t -  Cycle materials and lights    |" << endl
        << "| \\ -  Set light source position     |" << endl
        << "| u/U  Move the level surface        |" << endl
        << "| v/V  Add/Delete a level surface    |" << endl
        << "| w/W  Move bdr elements up/down     |" << endl
        << "| x/X  Rotate clipping plane (phi)   |" << endl
        << "| y/Y  Rotate clipping plane (theta) |" << endl
        << "| z/Z  Translate clipping plane      |" << endl
        << "+------------------------------------+" << endl
        << "| Function keys                      |" << endl
        << "+------------------------------------+" << endl
        << "| F1 - X window info and keystrokes  |" << endl
        << "| F2 - Update colors, etc.           |" << endl
        << "| F3/F4 - Shrink/Zoom bdr elements   |" << endl
        << "| F5 - Set level lines               |" << endl
        << "| F6 - Palete options                |" << endl
        << "| F7 - Manually set min/max value    |" << endl
        << "| F8 - List of subdomains to show    |" << endl
        << "| F9/F10 - Walk through subdomains   |" << endl
        << "| F11/F12 - Shrink/Zoom subdomains   |" << endl
        << "+------------------------------------+" << endl
        << "| Keypad                             |" << endl
        << "+------------------------------------+" << endl
        << "| 1-9  Small rotation, reset with 5  |" << endl
        << "| *,/  Scale up/down                 |" << endl
        << "| +/-  Change z-scaling              |" << endl
        << "| . -  Start/stop spinning           |" << endl
        << "| 0/Enter - Spinning speed and dir.  |" << endl
        << "+------------------------------------+" << endl
        << "| Mouse                              |" << endl
        << "+------------------------------------+" << endl
        << "| left   btn    - Rotation           |" << endl
        << "| middle btn    - Translation        |" << endl
        << "| right  btn    - Scaling            |" << endl
        << "| left  + Alt   - Tilt               |" << endl
        << "| left  + Shift - Spinning           |" << endl
        << "| right + Shift - Change light pos.  |" << endl
        << "| left  + Ctrl  - Spherical rotation |" << endl
        << "| middle+ Ctrl  - Object translation |" << endl
        << "| right + Ctrl  - Object scaling     |" << endl
        << "| left  + Ctrl + Shift - z-Spinning  |" << endl
        << "+------------------------------------+" << endl;
}

static void KeyIPressed()
{
   vssol3d -> ToggleCuttingPlane();
   SendExposeEvent();
}

void VisualizationSceneSolution3d::CPPrepare()
{
   PrepareCuttingPlane();
   PrepareCuttingPlaneLines();
}

void VisualizationSceneSolution3d::CPMoved()
{
   CPPrepare();
   if (cplane == 2)
   {
      Prepare();
      PrepareLines();
   }
}

static void KeyxPressed()
{
   vssol3d -> CuttingPlane -> IncreasePhi();
   vssol3d -> FindNodePos();
   vssol3d -> CPMoved();
   SendExposeEvent();
}

static void KeyXPressed()
{
   vssol3d -> CuttingPlane -> DecreasePhi();
   vssol3d -> FindNodePos();
   vssol3d -> CPMoved();
   SendExposeEvent();
}

static void KeyyPressed()
{
   vssol3d -> CuttingPlane -> IncreaseTheta();
   vssol3d -> FindNodePos();
   vssol3d -> CPMoved();
   SendExposeEvent();
}

static void KeyYPressed()
{
   vssol3d -> CuttingPlane -> DecreaseTheta();
   vssol3d -> FindNodePos();
   vssol3d -> CPMoved();
   SendExposeEvent();
}

static void KeyzPressed()
{
   vssol3d -> CuttingPlane -> IncreaseDistance();
   vssol3d -> FindNodePos();
   vssol3d -> CPMoved();
   SendExposeEvent();
}

static void KeyZPressed()
{
   vssol3d -> CuttingPlane -> DecreaseDistance();
   vssol3d -> FindNodePos();
   vssol3d -> CPMoved();
   SendExposeEvent();
}

static void KeymPressed()
{
   vssol3d -> ToggleDrawMesh();
   SendExposeEvent();
}

static void KeyePressed()
{
   vssol3d -> ToggleDrawElems();
   SendExposeEvent();
}

static void KeyMPressed()
{
   vssol3d -> ToggleCPDrawMesh();
   SendExposeEvent();
}

static void KeyEPressed()
{
   vssol3d -> ToggleCPDrawElems();
   SendExposeEvent();
}

static void KeyFPressed()
{
   vssol3d -> ToggleShading();
   SendExposeEvent();
}

static void KeyoPressed()
{
   if (vssol3d -> TimesToRefine < 32)
   {
      cout << "Refinement times: " << ++vssol3d -> TimesToRefine
           << endl;
      if (vssol3d -> GetShading() == 2)
      {
         vssol3d -> Prepare();
         vssol3d -> PrepareLines();
         vssol3d -> CPPrepare();
         SendExposeEvent();
      }
   }
}

static void KeyOPressed()
{
   if (vssol3d -> TimesToRefine > 1)
   {
      cout << "Refinement times: " << --vssol3d -> TimesToRefine
           << endl;
      if (vssol3d -> GetShading() == 2)
      {
         vssol3d -> Prepare();
         vssol3d -> PrepareLines();
         vssol3d -> CPPrepare();
         SendExposeEvent();
      }
   }
}

static void KeywPressed()
{
   if (vssol3d -> GetShading() == 2)
   {
      if ( fabs(vssol3d -> FaceShiftScale += 0.01) < 0.001 )
         vssol3d -> FaceShiftScale = 0.0;
      cout << "New Shift Scale: " << vssol3d -> FaceShiftScale
           << endl;
      vssol3d -> Prepare();
      vssol3d -> PrepareLines();
      vssol3d -> CPPrepare();
      SendExposeEvent();
   }
}

static void KeyWPressed()
{
   if (vssol3d -> GetShading() == 2)
   {
      if ( fabs(vssol3d -> FaceShiftScale -= 0.01) < 0.001 )
         vssol3d -> FaceShiftScale = 0.0;
      cout << "New Shift Scale: " << vssol3d -> FaceShiftScale
           << endl;
      vssol3d -> Prepare();
      vssol3d -> PrepareLines();
      vssol3d -> CPPrepare();
      SendExposeEvent();
   }
}

static void KeyuPressed()
{
   vssol3d -> MoveLevelSurf(+1);
   SendExposeEvent();
}

static void KeyUPressed()
{
   vssol3d -> MoveLevelSurf(-1);
   SendExposeEvent();
}

static void KeyvPressed()
{
   vssol3d -> NumberOfLevelSurf(+1);
   SendExposeEvent();
}

static void KeyVPressed()
{
   vssol3d -> NumberOfLevelSurf(-1);
   SendExposeEvent();
}

static int magic_key_pressed = 0;
void ToggleMagicKey()
{
   magic_key_pressed = 1-magic_key_pressed;
}

void KeyF3Pressed()
{
   if (vssol3d->GetShading() == 2)
   {
      if (vssol3d->bdrc.Width() == 0)
         vssol3d->ComputeBdrAttrCenter();
      vssol3d->shrink *= 0.9;
      if (magic_key_pressed)
         vssol3d -> Scale(1.11111111111111111111111);
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

void KeyF4Pressed()
{
   if (vssol3d->GetShading() == 2)
   {
      if (vssol3d->bdrc.Width() == 0)
         vssol3d->ComputeBdrAttrCenter();
      vssol3d->shrink *= 1.11111111111111111111111;
      if (magic_key_pressed)
         vssol3d -> Scale(0.9);
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

void KeyF11Pressed()
{
   if (vssol3d->GetShading() == 2)
   {
      if (vssol3d->matc.Width() == 0)
         vssol3d->ComputeElemAttrCenter();
      vssol3d->shrinkmat *= 0.9;
      if (magic_key_pressed)
         vssol3d -> Scale(1.11111111111111111111111);
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

void KeyF12Pressed()
{
   if (vssol3d->GetShading() == 2)
   {
      if (vssol3d->matc.Width() == 0)
         vssol3d->ComputeElemAttrCenter();
      vssol3d->shrinkmat *= 1.11111111111111111111111;
      if (magic_key_pressed)
         vssol3d -> Scale(0.9);
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

static void KeyF8Pressed()
{
   const Array<int> &attr_list = vssol3d->GetMesh()->bdr_attributes;
   Array<int> &attr_marker = vssol3d->bdr_attr_to_show;
   int attr;

   cout << "Bdr attributes ON: ";
   for (int i = 0; i < attr_list.Size(); i++)
      if (attr_marker[attr_list[i]-1])
         cout << " " << attr_list[i];
   cout << endl;

   cout << "Bdr attribute to toggle : " << flush;
   cin >> attr;
   if (attr < 1)
   {
      cout << "Hiding all bdr attributes." << endl;
      attr_marker = 0;
   }
   else if (attr > attr_marker.Size())
   {
      cout << "Showing all bdr attributes." << endl;
      attr_marker = 1;
   }
   else
   {
      attr_marker[attr-1] = !attr_marker[attr-1];
   }
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

static void KeyF9Pressed()
{
   const Array<int> &attr_list = vssol3d->GetMesh()->bdr_attributes;
   Array<int> &attr_marker = vssol3d->bdr_attr_to_show;
   int attr, n, j;

   if (attr_list.Size() == 0)
      return;
   n = 0;
   for (int i = 0; i < attr_list.Size(); i++)
      if (attr_marker[attr_list[i]-1])
         j = i, n++;
   if (n == 1)
      j = (j + 1) % attr_list.Size();
   else
      j = 0;
   attr = attr_list[j];
   attr_marker = 0;
   attr_marker[attr-1] = 1;
   cout << "Showing bdr attribute " << attr << endl;
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

static void KeyF10Pressed()
{
   const Array<int> &attr_list = vssol3d->GetMesh()->bdr_attributes;
   Array<int> &attr_marker = vssol3d->bdr_attr_to_show;
   int attr, n, j;

   if (attr_list.Size() == 0)
      return;
   n = 0;
   for (int i = 0; i < attr_list.Size(); i++)
      if (attr_marker[attr_list[i]-1])
         j = i, n++;
   if (n == 1)
      j = (j + attr_list.Size() - 1) % attr_list.Size();
   else
      j = attr_list.Size() - 1;
   attr = attr_list[j];
   attr_marker = 0;
   attr_marker[attr-1] = 1;
   cout << "Showing bdr attribute " << attr << endl;
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

// 'v' must define three vectors that add up to the zero vector
// that is v[0], v[1], and v[2] are the sides of a triangle
int UnitCrossProd(double v[][3], double nor[])
{
   // normalize the three vectors
   for (int i = 0; i < 3; i++)
      if (Normalize(v[i]))
         return 1;

   // find the pair that forms an angle closest to pi/2, i.e. having
   // the longest cross product:
   double cp[3][3], max_a = 0.;
   int k = 0;
   for (int i = 0; i < 3; i++)
   {
      CrossProd(v[(i+1)%3], v[(i+2)%3], cp[i]);
      double a = sqrt(InnerProd(cp[i], cp[i]));
      if (max_a < a)
         max_a = a, k = i;
   }
   if (max_a == 0.)
      return 1;
   for (int i = 0; i < 3; i++)
      nor[i] = cp[k][i] / max_a;

   return 0;
}

int Compute3DUnitNormal(const double p1[], const double p2[],
                        const double p3[], double nor[])
{
   double v[3][3];

   for (int i = 0; i < 3; i++)
   {
      v[0][i] = p2[i] - p1[i];
      v[1][i] = p3[i] - p2[i];
      v[2][i] = p1[i] - p3[i];
   }

   return UnitCrossProd(v, nor);
}

int Compute3DUnitNormal (const double p1[], const double p2[],
                         const double p3[], const double p4[], double nor[])
{
   double v[3][3];

   for (int i = 0; i < 3; i++)
   {
      // cross product of the two diagonals:
      /*
        v[0][i] = p3[i] - p1[i];
        v[1][i] = p4[i] - p2[i];
        v[2][i] = (p1[i] + p2[i]) - (p3[i] + p4[i]);
      */

      // cross product of the two vectors connecting the midpoints of the
      // two pairs of opposing sides; this gives the normal vector in the
      // midpoint of the quad:
      v[0][i] = 0.5 * ((p2[i] + p3[i]) - (p1[i] + p4[i]));
      v[1][i] = 0.5 * ((p4[i] + p3[i]) - (p1[i] + p2[i]));
      v[2][i] = p1[i] - p3[i];
   }

   return UnitCrossProd(v, nor);
}

VisualizationSceneSolution3d::VisualizationSceneSolution3d()
{}

VisualizationSceneSolution3d::VisualizationSceneSolution3d(Mesh &m, Vector &s)
{
   mesh = &m;
   sol = &s;
   GridF = NULL;

   Init();

   auxKeyFunc (AUX_h, Solution3dKeyHPressed);
   auxKeyFunc (AUX_H, Solution3dKeyHPressed);
}

void Set_Palette(int);

void VisualizationSceneSolution3d::Init()
{
   vssol3d = this;

   cplane = 0;
   cp_drawmesh = 0; cp_drawelems = 1;
   drawlsurf = 0;

   drawelems = shading = 1;
   drawmesh = 0;
   scaling = 0;

   shrink = 1.0;
   shrinkmat = 1.0;
   bdrc.SetSize(3,0);
   matc.SetSize(3,0);

   TimesToRefine = 1;
   FaceShiftScale = 0.0;

   minv = sol->Min();
   maxv = sol->Max();

   bdr_attr_to_show.SetSize (mesh->bdr_attributes.Max());
   bdr_attr_to_show = 1;

   VisualizationSceneScalarData::Init();

   scaling = 0;
   SetNewScalingFromBox(); // No scaling for 3D

   node_pos = new double[mesh->GetNV()];

   Set_Palette(12); // use the 'vivid' palette in 3D
   SetUseTexture(1);

   CuttingPlane = new Plane(-1.0,0.0,0.0,0.5*(x[0]+x[1]));

   // Question: How do we want to determine the levels?
   // a) from the solution in the whole domain or
   // b) only from the values in the CP
   // This implements a)
   SetLevelLines (minv, maxv, 15);  //  15 level curves

   nlevels = 1;

   FindNodePos();

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      auxKeyFunc (AUX_m, KeymPressed);
      auxKeyFunc (AUX_M, KeyMPressed);

      auxKeyFunc (AUX_e, KeyePressed);
      auxKeyFunc (AUX_E, KeyEPressed);

      auxKeyFunc (AUX_f, KeyFPressed);
      auxKeyFunc (AUX_F, KeyFPressed);

      auxKeyFunc (AUX_i, KeyIPressed);
      auxKeyFunc (AUX_I, KeyIPressed);

      auxKeyFunc (AUX_o, KeyoPressed);
      auxKeyFunc (AUX_O, KeyOPressed);

      auxKeyFunc (AUX_w, KeywPressed);
      auxKeyFunc (AUX_W, KeyWPressed);

      auxKeyFunc (AUX_x, KeyxPressed);
      auxKeyFunc (AUX_X, KeyXPressed);

      auxKeyFunc (AUX_y, KeyyPressed);
      auxKeyFunc (AUX_Y, KeyYPressed);

      auxKeyFunc (AUX_z, KeyzPressed);
      auxKeyFunc (AUX_Z, KeyZPressed);

      auxKeyFunc (AUX_u, KeyuPressed);
      auxKeyFunc (AUX_U, KeyUPressed);

      auxKeyFunc (AUX_v, KeyvPressed);
      auxKeyFunc (AUX_V, KeyVPressed);

      auxKeyFunc (XK_F3, KeyF3Pressed);
      auxKeyFunc (XK_F4, KeyF4Pressed);
      auxKeyFunc (XK_Num_Lock, ToggleMagicKey);

      auxKeyFunc (XK_F8, KeyF8Pressed);
      auxKeyFunc (XK_F9, KeyF9Pressed);
      auxKeyFunc (XK_F10, KeyF10Pressed);

      auxKeyFunc (XK_F11, KeyF11Pressed);
      auxKeyFunc (XK_F12, KeyF12Pressed);
   }
   displlist  = glGenLists (1);
   linelist   = glGenLists (1);
   cplanelist = glGenLists (1);
   cplanelineslist = glGenLists (1);
   lsurflist = glGenLists (1);

   Prepare();
   PrepareLines();
   CPPrepare();
   PrepareLevelSurf();
}

VisualizationSceneSolution3d::~VisualizationSceneSolution3d()
{
   glDeleteLists (displlist, 1);
   glDeleteLists (linelist, 1);
   glDeleteLists (cplanelist, 1);
   glDeleteLists (cplanelineslist, 1);
   glDeleteLists (lsurflist, 1);
   delete [] node_pos;
}

void VisualizationSceneSolution3d::NewMeshAndSolution(
   Mesh *new_m, Vector *new_sol, GridFunction *new_u, int rescale)
{
   if (mesh->GetNV() != new_m->GetNV())
   {
      delete [] node_pos;
      node_pos = new double[new_m->GetNV()];
   }
   mesh = new_m;
   sol = new_sol;
   GridF = new_u;
   FindNodePos();

   if (rescale == 1)
   {
      FindNewBox();
      PrepareAxes();
      minv = sol->Min();
      maxv = sol->Max();
   }
   else if (rescale == 2)
   {
      minv = sol->Min();
      maxv = sol->Max();
   }

   Prepare();
   PrepareLines();
   CPPrepare();
   PrepareLevelSurf();
}

void VisualizationSceneSolution3d::SetShading(int s)
{
   if (shading == s || s < 0)
      return;

   if (s > 2 || (GridF == NULL && s > 1))
      return;
   int os = shading;
   shading = s;
   Prepare();
   if (GridF != NULL && (s == 2 || os == 2))
   {
      PrepareLines();
      CPPrepare();
   }

   static const char *shading_type[3] =
      {"flat", "smooth", "non-conforming (with subdivision)"};
   cout << "Shading type : " << shading_type[shading] << endl;
}

void VisualizationSceneSolution3d::ToggleShading()
{
   if (GridF)
      SetShading((shading+1)%3);
   else
      SetShading(1-shading);
}

void VisualizationSceneSolution3d::SetRefineFactors(int f, int ignored)
{
   if (TimesToRefine == f || f < 1)
      return;

   TimesToRefine = f;

   if (shading == 2)
   {
      Prepare();
      PrepareLines();
      CPPrepare();
   }
}

void VisualizationSceneSolution3d::FindNewBox()
{
   int i;

   int nv = mesh -> GetNV();

   double *coord = mesh->GetVertex(0);

   x[0] = x[1] = coord[0];
   y[0] = y[1] = coord[1];
   z[0] = z[1] = coord[2];

   for(i = 1; i < nv; i++)
   {
      coord = mesh->GetVertex(i);
      if (coord[0] < x[0]) x[0] = coord[0];
      if (coord[1] < y[0]) y[0] = coord[1];
      if (coord[2] < z[0]) z[0] = coord[2];
      if (coord[0] > x[1]) x[1] = coord[0];
      if (coord[1] > y[1]) y[1] = coord[1];
      if (coord[2] > z[1]) z[1] = coord[2];
   }

   SetNewScalingFromBox();
}

void VisualizationSceneSolution3d::FindNodePos()
{
   int i, nnodes = mesh -> GetNV();

   for (i = 0; i < nnodes; i++)
      node_pos[i] = CuttingPlane -> Transform (mesh -> GetVertex (i));
}

void VisualizationSceneSolution3d::ToggleDrawMesh()
{
   drawmesh = (drawmesh+1)%3;
   if (drawmesh)
      PrepareLines();
}

void VisualizationSceneSolution3d::ToggleCuttingPlane()
{
   if (cplane == 2 && cp_drawmesh == 3)
      cp_drawmesh = 2;

   cplane = (cplane+1)%3;
   if (cplane)
      CPPrepare();
   if (cplane == 0 || cplane == 2)
   {
      Prepare();
      PrepareLines();
   }
}

void VisualizationSceneSolution3d::ToggleCPDrawElems()
{
   cp_drawelems = 1-cp_drawelems;
   if (cp_drawelems)
      PrepareCuttingPlane();
}

void VisualizationSceneSolution3d::ToggleCPDrawMesh()
{
   if (cplane == 1)
      cp_drawmesh = (cp_drawmesh+1)%3;
   else if (cplane == 2)
      cp_drawmesh = (cp_drawmesh+1)%4;
   if (cp_drawmesh)
      PrepareCuttingPlaneLines();
}

void VisualizationSceneSolution3d::MoveLevelSurf(int move)
{
   drawlsurf += move;
   if (drawlsurf < 0)
      drawlsurf = 0;
   if (drawlsurf > 49)
      drawlsurf = 49;
   PrepareLevelSurf();
}

void VisualizationSceneSolution3d::NumberOfLevelSurf(int c)
{
   nlevels += c;
   if (nlevels < 1)
      nlevels = 1;
   PrepareLevelSurf();
}

int Normalize(DenseMatrix &normals)
{
   int err = 0;
   for (int i = 0; i < normals.Width(); i++)
      err += Normalize(&normals(0, i));
   return err;
}

void VisualizationSceneSolution3d::GetFaceNormals(
   const int FaceNo, const int side, const IntegrationRule &ir,
   DenseMatrix &normals)
{
   // match the way GridFunction::GetFaceValues works
   double aJInv[9], alnor[3];
   DenseMatrix JInv(aJInv, 3, 3);
   Vector lnor(alnor, 3), nr;
   IntegrationRule eir(ir.GetNPoints());
   FaceElementTransformations *Tr =
      mesh->GetFaceElementTransformations(FaceNo);
   normals.SetSize(3, ir.GetNPoints());
   ElementTransformation *ETr;
   IntegrationPointTransformation *LTr;
   if (side == 0)
   {
      ETr = Tr->Elem1;
      LTr = &Tr->Loc1;
   }
   else
   {
      ETr = Tr->Elem2;
      LTr = &Tr->Loc2;
   }
   LTr->Transform(ir, eir);
   for (int i = 0; i < normals.Width(); i++)
   {
      LTr->Transf.SetIntPoint(&ir.IntPoint(i));
      const DenseMatrix &LJac = LTr->Transf.Jacobian();
      lnor(0) = LJac(1,0) * LJac(2,1) - LJac(2,0) * LJac(1,1);
      lnor(1) = LJac(2,0) * LJac(0,1) - LJac(0,0) * LJac(2,1);
      lnor(2) = LJac(0,0) * LJac(1,1) - LJac(1,0) * LJac(0,1);
      ETr->SetIntPoint(&eir.IntPoint(i));
      const DenseMatrix &Jac = ETr->Jacobian();
      CalcInverse(Jac, JInv);
      normals.GetColumnReference(i, nr);
      JInv.MultTranspose(lnor, nr);
   }
   if (side)
      normals *= -1.;
   JInv.ClearExternalData();
}

void VisualizationSceneSolution3d::DrawRefinedSurf(
   int n, double *points, int elem, int func, int part)
{
   int i, j;
   RefinedGeometry *RefG;
   IntegrationPointTransformation ip_transf;

   switch (n)
   {
   case 3:
      RefG = GlobGeometryRefiner.Refine (Geometry::TRIANGLE, TimesToRefine);
      ip_transf.Transf.SetFE (&TriangleFE);
      break;
   case 4:
      RefG = GlobGeometryRefiner.Refine (Geometry::SQUARE, TimesToRefine);
      ip_transf.Transf.SetFE (&QuadrilateralFE);
      break;
   case 5:
      DrawRefinedSurf (3, points, elem, func, 0); // draw (0,1,2)
      for (i = 0; i < 3; i++)
         points[4+i] = points[i]; // move point 0 to point 1
      DrawRefinedSurf (4, points+4, elem, func, 1); // draw (0,2,3,4)
      return;
   case 6:
      DrawRefinedSurf (4, points, elem, func, 0); // draw (0,1,2,3)
      for (i = 0; i < 3; i++)
         points[8+i] = points[i]; // move point 0 to point 2
      DrawRefinedSurf (4, points+8, elem, func, 1); // draw (0,3,4,5)
      return;
   default:
      return;
   }
   DenseMatrix &pm = ip_transf.Transf.GetPointMat();
   pm.SetSize (3, n);
   for (i = 0; i < n; i++)
      for (j = 0; j < 3; j++)
         pm(j,i) = points[4*i+j];
   IntegrationRule tir (RefG -> RefPts.GetNPoints());
   ip_transf.Transform (RefG -> RefPts, tir);
   DenseMatrix pointmat;
   Vector values;
   GridF -> GetValues (elem, tir, values, pointmat);

   LiftRefinedSurf (n, pointmat, values, NULL);

   switch (func)
   {
   case 1:
      DrawRefinedSurf (n, pointmat, values, RefG->RefGeoms);
      break;

   case 2:
      DrawRefinedSurfEdges (n, pointmat, values, RefG->RefEdges, part);
      break;

   case 3:
      DrawRefinedSurfLevelLines (n, pointmat, values, RefG->RefGeoms);
      break;
   }
}

void VisualizationSceneSolution3d::LiftRefinedSurf(
   int n, DenseMatrix &pointmat, Vector &values, int *RG)
{
   int i, j;

   if (FaceShiftScale == 0.0)
      return;

   double norm[3];

   if (RG == NULL) // use the normal vector of the cutting plane
   {
      double *eqn = CuttingPlane -> Equation();
      norm[0] = -eqn[0]; norm[1] = -eqn[1]; norm[2] = -eqn[2];
   }
   else
   {
      // use the normal vector to the polygon defined by
      // the points with indexes given by the first n integers in RG
      double pts[4][3];
      for (i = 0; i < n; i++)
         for (j = 0; j < 3; j++)
            pts[i][j] = pointmat(j,RG[i]);
      if (n > 3)
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], pts[3], norm);
      else
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], norm);
      if (j)
      {
         // could not compute normal
         cerr << "WARNING: VisualizationSceneSolution3d::LiftRefinedSurf"
              << endl;
         return;
      }
   }

   double bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                             (y[1]-y[0])*(y[1]-y[0]) +
                             (z[1]-z[0])*(z[1]-z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   for (i = 0; i < pointmat.Width(); i++)
   {
      double val = sc * (values(i) - minv) / (maxv - minv);
      for (j = 0; j < 3; j++)
         pointmat(j, i) +=  val * norm[j];
   }
}

void VisualizationSceneSolution3d::DrawRefinedSurf(
   int n, DenseMatrix &pointmat, Vector &values, Array<int> &RefGeoms)
{
   double norm[3], pts[4][3];

   for (int i = 0; i < RefGeoms.Size()/n; i++)
   {
      int *RG = &(RefGeoms[i*n]);
      int j;
      for (j = 0; j < n; j++)
         for (int l = 0; l < 3; l++)
            pts[j][l] = pointmat(l, RG[j]);
      if (n > 3)
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], pts[3], norm);
      else
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], norm);
      if (!j)
      {
         glBegin (GL_POLYGON);
         glNormal3dv (norm);
         for (j = 0; j < n; j++)
         {
            MySetColor (values(RG[j]), minv, maxv);
            glVertex3dv (pts[j]);
         }
         glEnd();
      }
      /*
        else
        cerr << "WARNING: VisualizationSceneSolution3d::DrawRefinedSurf"
        << endl;
      */
   }
}

void VisualizationSceneSolution3d::DrawRefinedSurfEdges(
   int n, DenseMatrix &pointmat, Vector &values, Array<int> &RefEdges,
   int part)
{
   int k, k_start, k_end;

   k_start = 0;
   k_end   = RefEdges.Size();
   if (part == 0)
      k_end = (k_end/n) * (n-1);
   if (part == 1)
      k_start = (k_end/n);

   if (part != 1)
      glBegin (GL_LINES);
   for (k = k_start; k < k_end; k++)
   {
      int RE = RefEdges[k];

      glVertex3d (pointmat(0, RE), pointmat(1, RE),
                  pointmat(2, RE));
   }
   if (part != 0)
      glEnd();
}

void VisualizationSceneSolution3d::DrawRefinedSurfLevelLines(
   int n, DenseMatrix &pointmat, Vector &values, Array<int> &RefGeoms)
{
   int j, k;
   int *RG;
   double point[4][4];

   for (k = 0; k < RefGeoms.Size()/n; k++)
   {
      RG = &(RefGeoms[k*n]);

      for (j = 0; j < n; j++)
      {
         for (int i = 0; i < 3; i++)
            point[j][i] = pointmat(i, RG[j]);
         point[j][3] = values(RG[j]);
      }
      DrawPolygonLevelLines (point[0], n, level);
   }
}

void VisualizationSceneSolution3d::PrepareFlat()
{
   int i, j;

   glNewList (displlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material();

   int ne = mesh -> GetNBE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double p[4][3], c[4];

   for (i = 0; i < ne; i++)
   {
      if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) continue;

      mesh->GetBdrElementVertices (i, vertices);

      if (cplane == 2)
      {
         int n = 0;
         for (j = 0; j < vertices.Size(); j ++)
            if (node_pos[vertices[j]] >= 0.0)
               n++;
         if (n < vertices.Size())
            continue;  // with the next boundary element
      }

      mesh->GetBdrPointMatrix (i, pointmat);

      for (j = 0; j < pointmat.Width(); j++)
      {
         p[j][0] = pointmat(0, j);
         p[j][1] = pointmat(1, j);
         p[j][2] = pointmat(2, j);
         c[j] = (*sol)(vertices[j]);
      }
      if (j == 3)
         DrawTriangle(p, c, minv, maxv);
      else
         DrawQuad(p, c, minv, maxv);
   }
   glEndList();
}

void VisualizationSceneSolution3d::PrepareFlat2()
{
   int i, j, k, fn, fo, ft, di, have_normals;
   double bbox_diam, vmin, vmax, mm;

   glNewList (displlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material();

   int nbe = mesh -> GetNBE();
   DenseMatrix pointmat, normals;
   Vector values;
   RefinedGeometry * RefG;
   Array<int> vertices;
   double norm[3];

   bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                      (y[1]-y[0])*(y[1]-y[0]) +
                      (z[1]-z[0])*(z[1]-z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   ft = 1;
   for (i = 0; i < nbe; i++)
   {
      if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) continue;

      mesh->GetBdrElementVertices (i, vertices);

      if (cplane == 2)
      {
         int n = 0;
         for (j = 0; j < vertices.Size(); j ++)
            if (node_pos[vertices[j]] >= 0.0)
               n++;
         if (n < vertices.Size())
            continue;  // with the next boundary element
      }

      mesh -> GetBdrElementFace (i, &fn, &fo);
      RefG = GlobGeometryRefiner.Refine (mesh -> GetFaceBaseGeometry (fn),
                                         TimesToRefine);
      // di = GridF -> GetFaceValues (fn, 2, RefG->RefPts, values, pointmat);

      // this assumes the interior boundary faces are properly oriented ...
      di = fo % 2;
      GridF -> GetFaceValues (fn, di, RefG->RefPts, values, pointmat);
      GetFaceNormals(fn, di, RefG->RefPts, normals);
      have_normals = 1;
      ShrinkPoints3D(pointmat, i, fn, fo);
      if (ft)
      {
         vmin = values.Min();
         vmax = values.Max();
         ft = 0;
      }
      else
      {
         mm = values.Min();
         if (mm < vmin)  vmin = mm;
         mm = values.Max();
         if (mm > vmax)  vmax = mm;
      }

      int sides;
      switch (mesh -> GetBdrElementType (i))
      {
      case Element::TRIANGLE:
         sides = 3;
         break;

      case Element::QUADRILATERAL:
         sides = 4;
         break;
      }

      // compute an average normal direction for the current face
      if (sc != 0.0)
      {
         for (int i = 0; i < 3; i++)
            norm[i] = 0.0;
         Normalize(normals);
         for (k = 0; k < normals.Width(); k++)
            for (int j = 0; j < 3; j++)
               norm[j] += normals(j, k);
         Normalize(norm);
         for (int i = 0; i < pointmat.Width(); i++)
         {
            double val = sc * (values(i) - minv) / (maxv - minv);
            for (int j = 0; j < 3; j++)
               pointmat(j, i) += val * norm[j];
         }
         have_normals = 0;
      }

      have_normals = have_normals ? 2 : 0;
      if (di)
         have_normals = -1 - have_normals;
      DrawPatch(pointmat, values, normals, sides, RefG->RefGeoms,
                minv, maxv, have_normals);
   }
   glEndList();
   cout << "VisualizationSceneSolution3d::PrepareFlat2() : [min,max] = ["
        << vmin << "," << vmax << "]" << endl;
}

void VisualizationSceneSolution3d::Prepare()
{
   int i,j;

   switch (shading)
   {
   case 0:
      PrepareFlat();
      return;
   case 2:
      PrepareFlat2();
      return;
   default:
      break;
   }

   glNewList (displlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
   Set_Material();

   int ne = mesh -> GetNBE();
   int nv = mesh -> GetNV();
   DenseMatrix pointmat;
   Array<int> vertices;
   double nor[3];

   Vector nx(nv);
   Vector ny(nv);
   Vector nz(nv);

   Table ba_to_be; // boundary_attribute--to--boundary_element
   {
      Table be_to_ba;
      be_to_ba.MakeI(ne);
      for (int i = 0; i < ne; i++)
         be_to_ba.AddAColumnInRow(i);
      be_to_ba.MakeJ();
      for (int i = 0; i < ne; i++)
         be_to_ba.AddConnection(i, mesh->GetBdrAttribute(i)-1);
      be_to_ba.ShiftUpI();

      Transpose(be_to_ba, ba_to_be);
   }

   for (int d = 0; d < mesh -> bdr_attributes.Size(); d++)
   {
      const int attr = mesh->bdr_attributes[d]-1;

      if (!bdr_attr_to_show[attr]) continue;

      const int nelem = ba_to_be.RowSize(attr);
      const int *elem = ba_to_be.GetRow(attr);

      for (i = 0; i < nelem; i++)
      {
         mesh->GetBdrElementVertices(elem[i], vertices);
         for (j = 0; j < vertices.Size(); j++)
            nx(vertices[j]) = ny(vertices[j]) = nz(vertices[j]) = 0.;
      }
      for (i = 0; i < nelem; i++)
      {
         mesh->GetBdrPointMatrix(elem[i], pointmat);
         mesh->GetBdrElementVertices(elem[i], vertices);

         if (pointmat.Width() == 3)
            j = Compute3DUnitNormal(&pointmat(0,0), &pointmat(0,1),
                                    &pointmat(0,2), nor);
         else
            j = Compute3DUnitNormal(&pointmat(0,0), &pointmat(0,1),
                                    &pointmat(0,2), &pointmat(0,3), nor);
         if (j == 0)
            for (j = 0; j < pointmat.Size(); j++)
            {
               nx(vertices[j]) += nor[0];
               ny(vertices[j]) += nor[1];
               nz(vertices[j]) += nor[2];
            }
      }

      for (i = 0; i < nelem; i++)
      {
         mesh->GetBdrElementVertices(elem[i], vertices);
         if (cplane == 2)
         {
            int n = 0;
            for (j = 0; j < vertices.Size(); j ++)
               if (node_pos[vertices[j]] >= 0.0)
                  n++;
            if (n < vertices.Size())
               continue;  // with the next boundary element
         }
         switch (mesh->GetBdrElementType(elem[i]))
         {
         case Element::TRIANGLE:
            glBegin (GL_TRIANGLES);
            break;

         case Element::QUADRILATERAL:
            glBegin (GL_QUADS);
            break;
         }
         mesh->GetBdrPointMatrix(elem[i], pointmat);

         for (j = 0; j < pointmat.Size(); j++)
         {
            MySetColor((*sol)(vertices[j]), minv ,maxv);
            glNormal3d(nx(vertices[j]), ny(vertices[j]), nz(vertices[j]));
            glVertex3dv(&pointmat(0, j));
         }
         glEnd();
      }

   }
   glEndList();
}

void VisualizationSceneSolution3d::PrepareLines()
{
   if (!drawmesh) return;

   if (shading == 2)
   {
      PrepareLines2();
      return;
   }

   int i, j, k, ne = mesh -> GetNBE();
   DenseMatrix pointmat;

   glNewList(linelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   Array<int> vertices;

   for (i = 0; i < ne; i++)
   {
      if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) continue;

      mesh -> GetBdrElementVertices (i, vertices);
      if (cplane == 2)
      {
         int n = 0;
         for (j = 0; j < vertices.Size(); j ++)
            if (node_pos[vertices[j]] >= 0.0)
               n++;
         if (n < vertices.Size())
            continue;  // with the next boundary element
      }

      double point[4][4];
      mesh->GetBdrPointMatrix (i, pointmat);

      switch (drawmesh)
      {
      case 1:
         switch (mesh->GetBdrElementType(i))
         {
         case Element::TRIANGLE:
            glBegin (GL_TRIANGLES);
            break;

         case Element::QUADRILATERAL:
            glBegin (GL_QUADS);
            break;
         }

         for (j = 0; j < pointmat.Size(); j++)
            glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j));
         glEnd();
         break;

      case 2:
         for (j = 0; j < pointmat.Size(); j++)
         {
            for (k = 0; k < 3; k++)
               point[j][k] = pointmat(k,j);
            point[j][3] = (*sol)(vertices[j]);
         }
         DrawPolygonLevelLines (point[0], pointmat.Size(), level);
         break;
      }
   }

   glEndList();
}

void VisualizationSceneSolution3d::PrepareLines2()
{
   int i, j, k, fn, fo, di;
   double bbox_diam;

   glNewList (linelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   int nbe = mesh -> GetNBE();
   DenseMatrix pointmat, normals;
   Vector values;
   RefinedGeometry * RefG;
   Array<int> vertices;

   bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                      (y[1]-y[0])*(y[1]-y[0]) +
                      (z[1]-z[0])*(z[1]-z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   for (i = 0; i < nbe; i++)
   {
      if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) continue;

      mesh->GetBdrElementVertices (i, vertices);

      if (cplane == 2)
      {
         int n = 0;
         for (j = 0; j < vertices.Size(); j ++)
            if (node_pos[vertices[j]] >= 0.0)
               n++;
         if (n < vertices.Size())
            continue;  // with the next boundary element
      }

      mesh -> GetBdrElementFace (i, &fn, &fo);
      RefG = GlobGeometryRefiner.Refine (mesh -> GetFaceBaseGeometry (fn),
                                         TimesToRefine);
      // di = GridF -> GetFaceValues (fn, 2, RefG->RefPts, values, pointmat);
      di = fo % 2;
      GridF -> GetFaceValues (fn, di, RefG->RefPts, values, pointmat);
      ShrinkPoints3D(pointmat, i, fn, fo);

      if (sc != 0.0)
      {
         GetFaceNormals(fn, di, RefG->RefPts, normals);
         double norm[3];
         for (int i = 0; i < 3; i++)
            norm[i] = 0.0;
         for (k = 0; k < normals.Width(); k++)
            for (int j = 0; j < 3; j++)
               norm[j] += normals(j, k);
         double len = sqrt(InnerProd(norm, norm));
         if (len > 0.0)
            len = 1.0 / len;
         for (int i = 0; i < 3; i++)
            norm[i] *= len;
         for (int i = 0; i < pointmat.Width(); i++)
         {
            double val = sc * (values(i) - minv) / (maxv - minv);
            for (int j = 0; j < 3; j++)
               pointmat(j, i) += val * norm[j];
         }
      }

      if (drawmesh == 1)
      {
         Array<int> &REdges = RefG->RefEdges;

         glBegin(GL_LINES);
         for (k = 0; k < REdges.Size(); k++)
            glVertex3dv(&pointmat(0, REdges[k]));
         glEnd();
      }
      else if (drawmesh == 2)
      {
         double point[4][4];
         int sides;
         switch (mesh -> GetBdrElementType (i))
         {
         case Element::TRIANGLE:
            sides = 3;
            break;

         case Element::QUADRILATERAL:
            sides = 4;
            break;
         }
         for (k = 0; k < RefG->RefGeoms.Size()/sides; k++)
         {
            int *RG = &(RefG->RefGeoms[k*sides]);

            for (j = 0; j < sides; j++)
            {
               for (int ii = 0; ii < 3; ii++)
                  point[j][ii] = pointmat(ii, RG[j]);
               point[j][3] = values(RG[j]);
            }
            DrawPolygonLevelLines (point[0], sides, level);
         }
      }
   }
   glEndList();
}

void VisualizationSceneSolution3d::CuttingPlaneFunc(int func)
{
   int i, j, k, m, n, n2;
   int flag[8], oedges[6];
   static const int tet_edges[12]={0,3, 0,2, 0,1, 1,2, 1,3, 2,3};
   static const int hex_edges[24] =
      { 0,1, 1,2, 3,2, 0,3, 4,5, 5,6, 7,6, 4,7, 0,4, 1,5, 2,6, 3,7 };
   static const int hex_cutting[24][3] =
      {{ 3,  4,  6}, {17,  9, 18}, { 4,  6,  1}, {19, 11, 20},
       {21, 12, 22}, { 6,  1,  3}, {23, 14, 16}, { 1,  3,  4},
       {18,  0, 17}, {15, 13, 10}, {20,  2, 19}, { 8, 15, 13},
       {10,  8, 15}, {22,  5, 21}, {13, 10,  8}, {16,  7, 23},
       { 9, 18,  0}, { 7, 23, 14}, {11, 20,  2}, { 0, 17,  9},
       {12, 22,  5}, { 2, 19, 11}, {14, 16,  7}, { 5, 21, 12}};
   const int *ev;
   double t, point[6][4], norm[3];

   DenseMatrix pointmat;

   Array<int> nodes;
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = n2 = 0; // n will be the number of intersection points
      mesh -> GetElementVertices(i, nodes);
      for (j = 0; j < nodes.Size(); j++)
         if (node_pos[nodes[j]] >= 0.0)
            flag[j] = 1;
         else
            flag[j] = -1;
      switch (mesh -> GetElementType(i))
      {
      case Element::TETRAHEDRON:
         ev = tet_edges;
         for (j = 0; j < 6; j++, ev += 2)
            if (flag[ev[0]] != flag[ev[1]])
               oedges[n++] = j;
         ev = tet_edges;
         break;
      case Element::HEXAHEDRON:
         int emark[12];
         ev = hex_edges;
         for (j = 0; j < 12; j++, ev += 2)
            emark[j] = flag[ev[1]] - flag[ev[0]];
         do
         {
            for (j = 0; j < 12; j++)
               if (emark[j])
                  break;
            if (j == 12)
               break;
            k = 2 * j;
            if (emark[j] > 0)
               k++;
            do
            {
               for (j = 0; j < 3; j++)
               {
                  m = hex_cutting[k][j];
                  ev = hex_edges + 2 * (m / 2);
                  if ((m % 2 == 0 && flag[ev[0]] > flag[ev[1]]) ||
                      (m % 2 == 1 && flag[ev[1]] > flag[ev[0]]))
                     break;
               }
               oedges[n2++] = k/2;
               emark[k/2] = 0;
               k = m;
            }
            while (k/2 != oedges[n]);
            if (n == 0)
            {
               n = n2;
            }
            else
            {
               break;
            }
         }
         while (1);
         n2 -= n;
         ev = hex_edges;
         break;
      default:
         break;
      }

      while (n > 2)
      {
         if (shading != 2)
            mesh -> GetPointMatrix (i, pointmat);
         else
         {
            const IntegrationRule *ir;
            ir = Geometries.GetVertices (mesh -> GetElementBaseGeometry(i));
            pointmat.SetSize (3, ir -> GetNPoints());
            for (j = 0; j < ir -> GetNPoints(); j++)
            {
               const IntegrationPoint &ip = ir -> IntPoint (j);
               pointmat(0,j) = ip.x;
               pointmat(1,j) = ip.y;
               pointmat(2,j) = ip.z;
            }
         }
         for (j = 0; j < n; j++)
         {
            const int *en = ev + 2*oedges[j];
            t =
               node_pos[ nodes[en[1]] ] /
               ( node_pos[ nodes[en[1]] ] -
                 node_pos[ nodes[en[0]] ] );
            for (k = 0; k < 3; k++)
               point[j][k] =
                  t*(pointmat(k,en[0])-pointmat(k,en[1])) +
                  pointmat(k,en[1]);
            point[j][3] =
               t * ((*sol)(nodes[en[0]]) -
                    (*sol)(nodes[en[1]])) +
               (*sol)(nodes[en[1]]);
         }

         switch (func)
         {
         case 1:  // PrepareCuttingPlane()
         {
            if (shading == 2)
            {
               // changes point for n > 4
               DrawRefinedSurf (n, point[0], i, 1);
            }
            else
            {
               while (1)
               {
                  if (n > 3)
                  {
                     j = Compute3DUnitNormal(point[0], point[1], point[2],
                                             point[3], norm);
                     if (j && n > 4)
                     {
                        for (int j = 3; j < n; j++)
                           for (int i = 0; i < 4; i++)
                              point[j-2][i] = point[j][i];
                        n -= 2;
                        continue;
                     }
                  }
                  else
                     j = Compute3DUnitNormal(point[0], point[1], point[2],
                                             norm);
                  break;
               }

               if (!j)
               {
                  glBegin(GL_POLYGON);
                  glNormal3dv (norm);
                  for(j = 0; j < n; j++)
                  {
                     MySetColor(point[j][3], minv, maxv);
                     glVertex3dv(point[j]);
                  }
                  glEnd();
               }
            }
         }
         break;

         case 2:  // PrepareCuttingPlaneLines() with mesh
         {
            if (shading == 2)
            {
               // changes point for n > 4
               DrawRefinedSurf (n, point[0], i, 2);
            }
            else
            {
               glBegin (GL_POLYGON);
               for(j = 0; j < n; j++)
                  glVertex3dv (point[j]);
               glEnd();
            }
         }
         break;

         case 3:  // PrepareCuttingPlaneLines() with level lines
         {
            if (shading == 2)
            {
               // changes point for n > 4
               DrawRefinedSurf (n, point[0], i, 3);
            }
            else
               DrawPolygonLevelLines(point[0], n, level);
         }
         break;
         }

         for (j = 0; j < n2; j++)
            oedges[j] = oedges[j+n];
         n = n2;
         n2 = 0;
      }
   }
}

void VisualizationSceneSolution3d::PrepareCuttingPlane()
{
   if (cp_drawelems == 0) return;
   if (cplane == 0) return;
   if (cplane == 2)
   {
      PrepareCuttingPlane2();
      return;
   }

   glNewList(cplanelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material();

   CuttingPlaneFunc(1);

   glEndList();
}

void VisualizationSceneSolution3d::PrepareCuttingPlane2()
{
   int i, j, n = 0;
   double p[4][3], c[4];
   DenseMatrix pointmat, normals;
   Vector values;
   RefinedGeometry * RefG;

   glNewList(cplanelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material();

   double * coord;

   Array<int> nodes;
   Array<int> partition (mesh -> GetNE());
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = 0; // n will be the number of nodes behind the cutting plane
      mesh -> GetElementVertices(i, nodes);
      for (j = 0; j < nodes.Size(); j++)
         if (node_pos[nodes[j]] >= 0.0)
            n++;

      if (n == nodes.Size())
         partition[i] = 0;
      else
         partition[i] = 1;
   }

   for (i = 0; i < mesh -> GetNFaces(); i++)
   {
      int e1, e2;
      mesh -> GetFaceElements (i, &e1, &e2);
      if (e2 >= 0 && partition[e1] != partition[e2])
         if (shading != 2)
         {
            mesh -> GetFaceVertices (i, nodes);
            for (j = 0; j < nodes.Size(); j++)
            {
               coord = mesh -> GetVertex(nodes[j]);
               p[j][0] = coord[0];
               p[j][1] = coord[1];
               p[j][2] = coord[2];
               c[j] = (*sol)(nodes[j]);
            }

            if (nodes.Size() == 3)
               DrawTriangle(p, c, minv, maxv);
            else
               DrawQuad(p, c, minv, maxv);
         }
         else // shading == 2
         {
            RefG = GlobGeometryRefiner.Refine (mesh -> GetFaceBaseGeometry (i),
                                               TimesToRefine);
            // partition[e1] is 0 if e1 is behind the cutting plane
            // and 1 otherwise
            int di = partition[e1];
            GridF -> GetFaceValues (i, di, RefG->RefPts, values, pointmat);
            GetFaceNormals(i, di, RefG->RefPts, normals);
            switch (mesh -> GetFaceBaseGeometry (i))
            {
            case Geometry::TRIANGLE:  n = 3; break;
            case Geometry::SQUARE:    n = 4; break;
            }
            // DrawRefinedSurf (n, pointmat, values, RefG->RefGeoms);
            DrawPatch(pointmat, values, normals, n, RefG->RefGeoms,
                      minv, maxv, di ? -3 : 2);
         } // end shading == 2
   }
   glEndList();
}

void VisualizationSceneSolution3d::PrepareCuttingPlaneLines()
{
   if (cp_drawmesh == 0) return;
   if (cplane == 0) return;
   if (cplane == 2 && cp_drawmesh != 3)
   {
      PrepareCuttingPlaneLines2();
      return;
   }

   int func;

   glNewList(cplanelineslist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   switch (cp_drawmesh)
   {
   case 1: func = 2; break;
   case 2:
   case 3: func = 3; break;
   }

   CuttingPlaneFunc(func);

   glEndList();
}

void VisualizationSceneSolution3d::PrepareCuttingPlaneLines2()
{
   int i, j, n = 0;
   double point[4][4];
   DenseMatrix pointmat;
   Vector values;
   RefinedGeometry * RefG;

   glNewList(cplanelineslist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   double * coord;

   Array<int> nodes;
   Array<int> partition (mesh -> GetNE());
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = 0;  // n will be the number of nodes behind the cutting plane
      mesh -> GetElementVertices(i,nodes);
      for(j=0; j<nodes.Size(); j++)
         if (node_pos[nodes[j]] >= 0.0)
            n++;

      if (n == nodes.Size())
         partition[i] = 0;
      else
         partition[i] = 1;
   }

   for (i = 0; i < mesh -> GetNFaces(); i++)
   {
      int e1, e2;
      mesh -> GetFaceElements (i, &e1, &e2);
      if (e2 >= 0 && partition[e1] != partition[e2])
      {
         if (shading != 2)
         {
            mesh -> GetFaceVertices (i, nodes);
            for(j = 0; j < nodes.Size(); j++)
            {
               coord = mesh -> GetVertex(nodes[j]);
               point[j][0] = coord[0];
               point[j][1] = coord[1];
               point[j][2] = coord[2];
               point[j][3] = (*sol)(nodes[j]);
            }
            switch (cp_drawmesh)
            {
            case 1:
               glBegin (GL_POLYGON);
               for (j = 0; j < nodes.Size(); j++)
                  glVertex3dv (point[j]);
               glEnd();
               break;
            case 2:
               DrawPolygonLevelLines (point[0], nodes.Size(), level);
               break;
            }
         }
         else // shading == 2
         {
            RefG = GlobGeometryRefiner.Refine (mesh -> GetFaceBaseGeometry (i),
                                               TimesToRefine);
            // partition[e1] is 0 if e1 is behind the cutting plane
            // and 1 otherwise
            int di = partition[e1];
            GridF -> GetFaceValues (i, di, RefG->RefPts, values, pointmat);
            switch (mesh -> GetFaceBaseGeometry (i))
            {
            case Geometry::TRIANGLE:  n = 3; break;
            case Geometry::SQUARE:    n = 4; break;
            }
            switch (cp_drawmesh)
            {
            case 1:
               DrawRefinedSurfEdges (n, pointmat, values, RefG->RefEdges);
               break;
            case 2:
               DrawRefinedSurfLevelLines (n, pointmat, values,
                                          RefG->RefGeoms);
               break;
            }
         } // end shading == 2
      }
   }
   glEndList();
}

void VisualizationSceneSolution3d::PrepareLevelSurf()
{
   double t, lvl, pts[4][4], norm[3];
   int ie, i, j, pos[4], ne, l;
   DenseMatrix pointmat;
   Array<int> vertices;

   if (drawlsurf == 0)
   {
      //  Create empty list
      glNewList(lsurflist, GL_COMPILE);
      glEndList();
      return;
   }

   ne = mesh -> GetNE();

   glNewList(lsurflist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material();

   levels.SetSize(nlevels);
   for (l = 0; l < nlevels; l++)
   {
      lvl = ((double)(50*l+drawlsurf) / (nlevels*50));
      levels[l] = (1.0-lvl) * minv + lvl * maxv;
   }

   for (ie = 0; ie < ne; ie++)
   {
      // Ignore element that is not tetrahedron
      if (mesh -> GetElementType(ie) != Element::TETRAHEDRON)
         continue;

      mesh -> GetPointMatrix (ie, pointmat);
      mesh -> GetElementVertices (ie, vertices);
      pts[0][0] = pointmat(0,0);
      pts[0][1] = pointmat(1,0);
      pts[0][2] = pointmat(2,0);
      pts[0][3] = (*sol)(vertices[0]);
      pts[1][0] = pointmat(0,1);
      pts[1][1] = pointmat(1,1);
      pts[1][2] = pointmat(2,1);
      pts[1][3] = (*sol)(vertices[1]);
      pts[2][0] = pointmat(0,2);
      pts[2][1] = pointmat(1,2);
      pts[2][2] = pointmat(2,2);
      pts[2][3] = (*sol)(vertices[2]);
      pts[3][0] = pointmat(0,3);
      pts[3][1] = pointmat(1,3);
      pts[3][2] = pointmat(2,3);
      pts[3][3] = (*sol)(vertices[3]);

      for (l = 0; l < levels.Size(); l++)
      {
         lvl = levels[l];
         MySetColor (lvl,minv,maxv);

         j = 0;
         for (i = 0; i < 4; i++)
            if (pts[i][3] < lvl)
               pos[j++] = i;
            else
               pos[3-i+j] = i;
         if (j == 3)
            i = pos[3], pos[3] = pos[0], pos[0] = i, j = 1;
         if (j == 1)
         {
            double vert[3][3];

            t = (lvl-pts[pos[0]][3])/(pts[pos[1]][3]-pts[pos[0]][3]);
            vert[0][0] = (1.0-t) * pts[pos[0]][0] + t * pts[pos[1]][0];
            vert[0][1] = (1.0-t) * pts[pos[0]][1] + t * pts[pos[1]][1];
            vert[0][2] = (1.0-t) * pts[pos[0]][2] + t * pts[pos[1]][2];

            t = (lvl-pts[pos[0]][3])/(pts[pos[2]][3]-pts[pos[0]][3]);
            vert[1][0] = (1.0-t) * pts[pos[0]][0] + t * pts[pos[2]][0];
            vert[1][1] = (1.0-t) * pts[pos[0]][1] + t * pts[pos[2]][1];
            vert[1][2] = (1.0-t) * pts[pos[0]][2] + t * pts[pos[2]][2];

            t = (lvl-pts[pos[0]][3])/(pts[pos[3]][3]-pts[pos[0]][3]);
            vert[2][0] = (1.0-t) * pts[pos[0]][0] + t * pts[pos[3]][0];
            vert[2][1] = (1.0-t) * pts[pos[0]][1] + t * pts[pos[3]][1];
            vert[2][2] = (1.0-t) * pts[pos[0]][2] + t * pts[pos[3]][2];

            if (!Compute3DUnitNormal (vert[0], vert[1], vert[2], norm))
            {
               glNormal3dv (norm);

               glBegin (GL_TRIANGLES);
               glVertex3dv (vert[0]);
               glVertex3dv (vert[1]);
               glVertex3dv (vert[2]);
               glEnd();
            }
         }
         else if (j == 2)
         {
            double vert[4][3];

            t = (lvl-pts[pos[0]][3])/(pts[pos[2]][3]-pts[pos[0]][3]);
            vert[0][0] = (1.0-t) * pts[pos[0]][0] + t * pts[pos[2]][0];
            vert[0][1] = (1.0-t) * pts[pos[0]][1] + t * pts[pos[2]][1];
            vert[0][2] = (1.0-t) * pts[pos[0]][2] + t * pts[pos[2]][2];

            t = (lvl-pts[pos[0]][3])/(pts[pos[3]][3]-pts[pos[0]][3]);
            vert[1][0] = (1.0-t) * pts[pos[0]][0] + t * pts[pos[3]][0];
            vert[1][1] = (1.0-t) * pts[pos[0]][1] + t * pts[pos[3]][1];
            vert[1][2] = (1.0-t) * pts[pos[0]][2] + t * pts[pos[3]][2];

            t = (lvl-pts[pos[1]][3])/(pts[pos[2]][3]-pts[pos[1]][3]);
            vert[3][0] = (1.0-t) * pts[pos[1]][0] + t * pts[pos[2]][0];
            vert[3][1] = (1.0-t) * pts[pos[1]][1] + t * pts[pos[2]][1];
            vert[3][2] = (1.0-t) * pts[pos[1]][2] + t * pts[pos[2]][2];

            t = (lvl-pts[pos[1]][3])/(pts[pos[3]][3]-pts[pos[1]][3]);
            vert[2][0] = (1.0-t) * pts[pos[1]][0] + t * pts[pos[3]][0];
            vert[2][1] = (1.0-t) * pts[pos[1]][1] + t * pts[pos[3]][1];
            vert[2][2] = (1.0-t) * pts[pos[1]][2] + t * pts[pos[3]][2];

            if (!Compute3DUnitNormal (vert[0], vert[1], vert[2], vert[3],
                                      norm))
            {
               glNormal3dv (norm);

               glBegin (GL_QUADS);
               glVertex3dv (vert[0]);
               glVertex3dv (vert[1]);
               glVertex3dv (vert[2]);
               glVertex3dv (vert[3]);
               glEnd();
            }
         }
      }
   }

   glEndList();
}

void VisualizationSceneSolution3d::ShrinkPoints3D(DenseMatrix &pointmat,
                                                  int i, int fn, int fo)
{
   if (shrink != 1.0)
   {
      int attr = mesh->GetBdrAttribute(i);
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < 3; d++)
            pointmat(d,k) = shrink*pointmat(d,k) + (1-shrink)*bdrc(d,attr-1);
   }

   if (shrinkmat != 1.0)
   {
      int attr, elem1, elem2;
      mesh->GetFaceElements(fn, &elem1, &elem2);
      if (fo % 2 == 0)
         attr = mesh->GetAttribute(elem1);
      else
         attr = mesh->GetAttribute(elem2);
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < 3; d++)
            pointmat(d,k) = shrinkmat*pointmat(d,k) + (1-shrinkmat)*matc(d,attr-1);
   }
}

void VisualizationSceneSolution3d::ComputeBdrAttrCenter()
{
   DenseMatrix pointmat;
   Vector nbdrc(mesh->bdr_attributes.Max());

   bdrc.SetSize(3,mesh->bdr_attributes.Max());
   bdrc = 0.0;
   nbdrc = 0.0;

   for (int i = 0; i < mesh -> GetNBE(); i++)
   {
      mesh->GetBdrPointMatrix(i, pointmat);
      nbdrc(mesh->GetBdrAttribute(i)-1) += pointmat.Width();
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < 3; d++)
            bdrc(d,mesh->GetBdrAttribute(i)-1) += pointmat(d,k);
   }

   for (int i = 0; i < mesh->bdr_attributes.Max(); i++)
      if (nbdrc(i) != 0)
         for (int d = 0; d < 3; d++)
            bdrc(d,i) /= nbdrc(i);
}

void VisualizationSceneSolution3d::ComputeElemAttrCenter()
{
   DenseMatrix pointmat;
   Vector nmatc(mesh->attributes.Max());

   matc.SetSize(3,mesh->attributes.Max());
   matc = 0.0;
   nmatc = 0.0;

   for (int i = 0; i < mesh -> GetNE(); i++)
   {
      mesh->GetPointMatrix(i, pointmat);
      nmatc(mesh->GetAttribute(i)-1) += pointmat.Width();
      for (int k = 0; k < pointmat.Width(); k++)
         for (int d = 0; d < 3; d++)
            matc(d,mesh->GetAttribute(i)-1) += pointmat(d,k);
   }

   for (int i = 0; i < mesh->attributes.Max(); i++)
      if (nmatc(i) != 0)
         for (int d = 0; d < 3; d++)
            matc(d,i) /= nmatc(i);
}

void VisualizationSceneSolution3d::Draw()
{
   glEnable(GL_DEPTH_TEST);

   Set_Background();
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // draw colored faces
   glPolygonOffset (1, 1);
   glEnable (GL_POLYGON_OFFSET_FILL);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   // model transformation
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity();

   glMultMatrixd (rotmat);
   glMultMatrixd (translmat);
   glScaled(xscale, yscale, zscale);
   glTranslated(-(x[0]+x[1])/2, -(y[0]+y[1])/2, -(z[0]+z[1])/2);

   glDisable(GL_CLIP_PLANE0);
   // draw colorbar
   glDisable(GL_LIGHTING);
   if (colorbar)
      if (drawmesh == 2 || cp_drawmesh >= 2)
         if (drawlsurf)
            DrawColorBar(minv,maxv,&level,&levels);
         else
            DrawColorBar(minv,maxv,&level);
      else
         if (drawlsurf)
            DrawColorBar(minv,maxv,NULL,&levels);
         else
            DrawColorBar(minv,maxv);

   Set_Black_Material();
   // draw axes
   if (drawaxes){
      glCallList(axeslist);
      DrawCoordinateCross();
   }
   DrawRuler();
   if (light)
      glEnable(GL_LIGHTING);

   // define and enable the clipping plane
   if (cplane)
   {
      //  double *eqn
      //  eqn = CuttingPlane->Equation();
      //  double tr_eqn[4];
      //  tr_eqn[0] = eqn[0];
      //  tr_eqn[1] = eqn[1];
      //  tr_eqn[2] = eqn[2];
      //  tr_eqn[3] = eqn[3];
      //  tr_eqn[3] = eqn[3] + 1e-2;
      //  glClipPlane(GL_CLIP_PLANE0,tr_eqn);
      glClipPlane (GL_CLIP_PLANE0, CuttingPlane->Equation());
      Set_Black_Material();
      glDisable(GL_CLIP_PLANE0);
      if ( cp_drawmesh )
         glCallList(cplanelineslist);
      glEnable(GL_CLIP_PLANE0);
   }

   Set_Black_Material();

   glDisable(GL_LIGHTING);
   // draw lines
   if (drawmesh)
      glCallList(linelist);
   if (light)
      glEnable(GL_LIGHTING);

   if (MatAlpha < 1.0)
      Set_Transparency();

   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
   if (GetUseTexture())
   {
      glEnable (GL_TEXTURE_1D);
      glColor4d(1, 1, 1, 1);
   }

   if (drawlsurf)
   {
      // glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      glCallList (lsurflist);
      // Set_Black_Material();
      // glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      // glCallList (lsurflist);
   }

   // draw elements
   if (drawelems)
      glCallList(displlist);

   if (cplane && cp_drawelems)
   {
      glDisable(GL_CLIP_PLANE0);
      glCallList(cplanelist);
      glEnable(GL_CLIP_PLANE0);
   }

   if (GetUseTexture())
      glDisable (GL_TEXTURE_1D);

   if (MatAlpha < 1.0)
      Remove_Transparency();

   glFlush();
   glXSwapBuffers (auxXDisplay(), auxXWindow());
}
