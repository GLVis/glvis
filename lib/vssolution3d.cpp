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

void Set_Material();
void Set_Black_Material();

VisualizationSceneSolution3d * vssol3d;
extern VisualizationScene   * locscene;
extern VisualizationSceneScalarData * vsdata;

// Definitions of some more keys

static void Solution3dKeyHPressed (){
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
        << "| F3/F4 - Shrink/Zoom subdomains     |" << endl
        << "| F5 - Set level lines               |" << endl
        << "| F6 - Palete options                |" << endl
        << "| F7 - Manually set min/max value    |" << endl
        << "| F8 - List of subdomains to show    |" << endl
        << "| F9/F10 - Walk through subdomains   |" << endl
        << "| F11/F12 - Shrink/Zoom elements     |" << endl
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

static void KeyIPressed ()
{
   vssol3d -> ToggleCuttingPlane();
   SendExposeEvent();
}

void VisualizationSceneSolution3d::CPPrepare()
{
   PrepareCuttingPlane();
   PrepareCuttingPlaneLines();
   if (cplane == 2)
   {
      Prepare();
      PrepareLines();
   }
   SendExposeEvent();
}

static void KeyxPressed()
{
   vssol3d -> CuttingPlane -> IncreasePhi();
   vssol3d -> FindNodePos();
   vssol3d -> CPPrepare();
   SendExposeEvent();
}

static void KeyXPressed()
{
   vssol3d -> CuttingPlane -> DecreasePhi();
   vssol3d -> FindNodePos();
   vssol3d -> CPPrepare();
   SendExposeEvent();
}

static void KeyyPressed()
{
   vssol3d -> CuttingPlane -> IncreaseTheta();
   vssol3d -> FindNodePos();
   vssol3d -> CPPrepare();
   SendExposeEvent();
}

static void KeyYPressed()
{
   vssol3d -> CuttingPlane -> DecreaseTheta();
   vssol3d -> FindNodePos();
   vssol3d -> CPPrepare();
   SendExposeEvent();
}

static void KeyzPressed()
{
   vssol3d -> CuttingPlane -> IncreaseDistance();
   vssol3d -> FindNodePos();
   vssol3d -> CPPrepare();
   SendExposeEvent();
}

static void KeyZPressed()
{
   vssol3d -> CuttingPlane -> DecreaseDistance();
   vssol3d -> FindNodePos();
   vssol3d -> CPPrepare();
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
         locscene -> Prepare();
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
         locscene -> Prepare();
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
      locscene -> Prepare();
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
      locscene -> Prepare();
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
   (vsdata -> GetMesh()) -> ScaleSubdomains(0.9);
   if (magic_key_pressed)
      locscene -> Scale(1.11111111111111111111111);
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

void KeyF4Pressed()
{
   (vsdata -> GetMesh()) -> ScaleSubdomains(1.11111111111111111111111);
   if (magic_key_pressed)
      locscene -> Scale(0.9);
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

void KeyF11Pressed()
{
   (vsdata -> GetMesh()) -> ScaleElements(0.9);
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

void KeyF12Pressed()
{
   (vsdata -> GetMesh()) -> ScaleElements(1.11111111111111111111111);
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

static void KeyF8Pressed()
{
   int attr;

   cout << "Bdr attributes ON: ";
   for (attr = 0; attr < vssol3d -> bdr_attr_to_show.Size(); attr++)
      if (vssol3d -> bdr_attr_to_show[attr])
         cout << " " << attr;
   cout << endl;

   cout << "Bdr attribute to toggle : " << flush;
   cin >> attr;
   if (attr < 0 )
   {
      for (attr = 0; attr < vssol3d -> bdr_attr_to_show.Size(); attr++)
         vssol3d -> bdr_attr_to_show[attr] = 0;
   }
   else if (attr >= vssol3d -> bdr_attr_to_show.Size())
   {
      for (attr = 0; attr < vssol3d -> bdr_attr_to_show.Size(); attr++)
         vssol3d -> bdr_attr_to_show[attr] = 1;
   }
   else
   {
      vssol3d -> bdr_attr_to_show[attr] = !vssol3d -> bdr_attr_to_show[attr];
   }
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

static void KeyF9Pressed()
{
   int attr;

   if (vssol3d -> attr_to_show == -1) {
      for (attr = 0; attr < vssol3d -> bdr_attr_to_show.Size(); attr++)
         vssol3d -> bdr_attr_to_show[attr] = 0;
      vssol3d -> attr_to_show = 0;
   } else
      vssol3d -> bdr_attr_to_show[vssol3d -> attr_to_show++] = 0;

   if (vssol3d -> attr_to_show >= vssol3d -> bdr_attr_to_show.Size())
      vssol3d -> attr_to_show = 0;
   vssol3d -> bdr_attr_to_show[vssol3d -> attr_to_show] = 1;
   cout << "Showing bdr attribute " << vssol3d -> attr_to_show << endl;
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

static void KeyF10Pressed()
{
   int attr;

   if (vssol3d -> attr_to_show == -1) {
      for (attr = 0; attr < vssol3d -> bdr_attr_to_show.Size(); attr++)
         vssol3d -> bdr_attr_to_show[attr] = 0;
   } else
      vssol3d -> bdr_attr_to_show[vssol3d -> attr_to_show--] = 0;

   if (vssol3d -> attr_to_show < 0)
      vssol3d -> attr_to_show = vssol3d -> bdr_attr_to_show.Size()-1;
   vssol3d -> bdr_attr_to_show[vssol3d -> attr_to_show] = 1;
   cout << "Showing bdr attribute " << vssol3d -> attr_to_show << endl;
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

double Distance (double p1[], double p2[])
{
   double v[] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };

   return sqrt (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double InnerProduct (double p1[], double p2[], double p3[])
{
   double v1[] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };
   double v2[] = { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] };

   return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

int Compute3DUnitNormal (double p1[], double p2[], double p3[], double nor[])
{
   int err = 0;
   double d1 = Distance (p2, p3);
   double d2 = Distance (p1, p3);
   double d3 = Distance (p1, p2);
   double *pp[3], pr, dmin;
   if (d1 < d2)
      if (d1 < d3)
         pp[0] = p1, pp[1] = p2, pp[2] = p3, pr = d2*d3, dmin = d1;
      else
         pp[0] = p3, pp[1] = p1, pp[2] = p2, pr = d1*d2, dmin = d3;
   else
      if (d2 < d3)
         pp[0] = p2, pp[1] = p3, pp[2] = p1, pr = d1*d3, dmin = d2;
      else
         pp[0] = p3, pp[1] = p1, pp[2] = p2, pr = d1*d2, dmin = d3;
   double inpr = InnerProduct (pp[0], pp[1], pp[2]);
   if (dmin == 0.0 || inpr >= (1.-1.e-9) * pr)
      err += 2;
   double v1[] = { pp[1][0]-pp[0][0], pp[1][1]-pp[0][1], pp[1][2]-pp[0][2] };
   double v2[] = { pp[2][0]-pp[0][0], pp[2][1]-pp[0][1], pp[2][2]-pp[0][2] };
   double n[] = {  v1[1] * v2[2] - v1[2] * v2[1],
                   v1[2] * v2[0] - v1[0] * v2[2],
                   v1[0] * v2[1] - v1[1] * v2[0]  };
   double rlen = 1.0 / sqrt (n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

   if (!finite(rlen))
      rlen = 1.0, err += 1;

   nor[0] = rlen * n[0];
   nor[1] = rlen * n[1];
   nor[2] = rlen * n[2];

   return err;
}

int Compute3DUnitNormal (double p1[], double p2[], double p3[], double p4[],
                         double nor[])
{
   if (Compute3DUnitNormal(p1, p2, p3, nor) == 0)
      return 0;
   if (Compute3DUnitNormal(p1, p2, p4, nor) == 0)
      return 0;
   if (Compute3DUnitNormal(p1, p3, p4, nor) == 0)
      return 0;
   return Compute3DUnitNormal(p2, p3, p4, nor);
}

VisualizationSceneSolution3d::VisualizationSceneSolution3d(){}

VisualizationSceneSolution3d
::VisualizationSceneSolution3d ( Mesh & m, Vector & s )
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

   TimesToRefine = 1;
   FaceShiftScale = 0.0;

   attr_to_show = -1;

   minv = sol->Min();
   maxv = sol->Max();

   int i;
   bdr_attr_to_show.SetSize (mesh->bdr_attributes.Max());
   for (i = 0; i < bdr_attr_to_show.Size(); i++)
      bdr_attr_to_show[i] = 1;

   VisualizationSceneScalarData :: Init();

   scaling = 0;             //
   SetNewScalingFromBox (); // No scaling for 3D

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

VisualizationSceneSolution3d::~VisualizationSceneSolution3d (){
   glDeleteLists (displlist, 1);
   glDeleteLists (linelist, 1);
   glDeleteLists (cplanelist, 1);
   glDeleteLists (cplanelineslist, 1);
   glDeleteLists (lsurflist, 1);
   delete [] node_pos;
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

   SetNewScalingFromBox ();
}

void VisualizationSceneSolution3d :: FindNodePos()
{
   int i, nnodes = mesh -> GetNV();

   for(i = 0; i < nnodes; i++)
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
      CPPrepare ();
   else
   {
      Prepare();
      PrepareLines();
   }
}

void VisualizationSceneSolution3d::ToggleCPDrawElems ()
{
   cp_drawelems = 1-cp_drawelems;
   if (cp_drawelems)
      PrepareCuttingPlane();
}

void VisualizationSceneSolution3d::ToggleCPDrawMesh ()
{
   if (cplane == 1)
      cp_drawmesh = (cp_drawmesh+1)%3;
   else if (cplane == 2)
      cp_drawmesh = (cp_drawmesh+1)%4;
   if (cp_drawmesh)
      PrepareCuttingPlaneLines();
}

void VisualizationSceneSolution3d::MoveLevelSurf (int move)
{
   drawlsurf += move;
   if (drawlsurf < 0)
      drawlsurf = 0;
   if (drawlsurf > 49)
      drawlsurf = 49;
   PrepareLevelSurf ();
}

void VisualizationSceneSolution3d::NumberOfLevelSurf (int c)
{
   nlevels += c;
   if (nlevels < 1)
      nlevels = 1;
   PrepareLevelSurf ();
}

void VisualizationSceneSolution3d::DrawRefinedSurf (
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

void VisualizationSceneSolution3d::LiftRefinedSurf (
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

void VisualizationSceneSolution3d::DrawRefinedSurf (
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
         glEnd ();
      }
      /*
        else
        cerr << "WARNING: VisualizationSceneSolution3d::DrawRefinedSurf"
        << endl;
      */
   }
}

void VisualizationSceneSolution3d::DrawRefinedSurfEdges (
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
      glEnd ();
}

void VisualizationSceneSolution3d::DrawRefinedSurfLevelLines (
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

void VisualizationSceneSolution3d :: PrepareFlat (){
   int i, j;

   glNewList (displlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

   int ne = mesh -> GetNBE();
   DenseMatrix pointmat;
   Array<int> vertices;

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

      switch (mesh->GetBdrElementType(i))
      {
      case Element::TRIANGLE:
         glBegin (GL_TRIANGLES);
         break;

      case Element::QUADRILATERAL:
         glBegin (GL_QUADS);
         break;
      }
      mesh->GetBdrPointMatrix (i, pointmat);

      double v10[] = { pointmat(0,1)-pointmat(0,0),
                       pointmat(1,1)-pointmat(1,0),
                       pointmat(2,1)-pointmat(2,0) };
      double v21[] = { pointmat(0,2)-pointmat(0,1),
                       pointmat(1,2)-pointmat(1,1),
                       pointmat(2,2)-pointmat(2,1) };

      double norm[] = { v10[1]*v21[2]-v10[2]*v21[1],
                        v10[2]*v21[0]-v10[0]*v21[2],
                        v10[0]*v21[1]-v10[1]*v21[0] };
      double rlen = 1.0/sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);

      glNormal3d (norm[0]*rlen, norm[1]*rlen, norm[2]*rlen);

      for (j = 0; j < pointmat.Size(); j++) {
         MySetColor ( (*sol)(vertices[j]) , minv, maxv );
         glVertex3d (pointmat(0, j),
                     pointmat(1, j),
                     pointmat(2, j));
      }
      glEnd ();
   }
   glEndList ();
}

void VisualizationSceneSolution3d::PrepareFlat2()
{
   int i, j, k, fn, fo, ft, di;
   double bbox_diam, vmin, vmax, mm;

   glNewList (displlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

   int nbe = mesh -> GetNBE();
   DenseMatrix pointmat;
   Vector values;
   RefinedGeometry * RefG;
   Array<int> vertices;

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
      int *RG = &(RefG->RefGeoms[0]);
      double pts[][3] = { { pointmat(0,RG[0]), pointmat(1,RG[0]),
                            pointmat(2,RG[0]) },
                          { pointmat(0,RG[1]), pointmat(1,RG[1]),
                            pointmat(2,RG[1]) },
                          { pointmat(0,RG[2]), pointmat(1,RG[2]),
                            pointmat(2,RG[2]) } };
      double norm[3];
      if (Compute3DUnitNormal (pts[0], pts[1], pts[2], norm))
         cerr << "WARNING: VisualizationSceneSolution3d::PrepareFlat2()"
              << endl;
      if (di != 0 && sc != 0.0)
      {
         norm[0] = -norm[0];
         norm[1] = -norm[1];
         norm[2] = -norm[2];
      }

      for (k = 0; k < RefG->RefGeoms.Size()/sides; k++)
      {
         RG = &(RefG->RefGeoms[k*sides]);
         switch (sides)
         {
         case 3:
            glBegin (GL_TRIANGLES);
            break;

         case 4:
            glBegin (GL_QUADS);
            break;
         }

         if (sc == 0.0)
         {
            glNormal3dv (norm);
            for (j = 0; j < sides; j++)
            {
               MySetColor ( values(RG[j]) , minv , maxv );
               glVertex3d (pointmat(0, RG[j]), pointmat(1, RG[j]),
                           pointmat(2, RG[j]));
            }
         }
         else
         {
            double nnorm[3];
            for (j = 0; j < 3; j++)
            {
               double val = sc * (values(RG[j]) - minv) / (maxv - minv);
               for (int l = 0; l < 3; l++)
                  pts[j][l] = pointmat(l, RG[j]) + val * norm[l];
            }
            if (Compute3DUnitNormal (pts[0], pts[1], pts[2], nnorm))
               cerr << "WARNING: VisualizationSceneSolution3d::PrepareFlat2()"
                    << endl;
            glNormal3dv (nnorm);
            for (j = 0; j < 3; j++)
            {
               MySetColor ( values(RG[j]) , minv , maxv );
               glVertex3dv (pts[j]);
            }
            for ( ; j < sides; j++)
            {
               double val = (values(RG[j]) - minv) / (maxv - minv);
               MySetColor ( values(RG[j]) , minv , maxv );
               val *= sc;
               glVertex3d (pointmat(0, RG[j])+val*norm[0],
                           pointmat(1, RG[j])+val*norm[1],
                           pointmat(2, RG[j])+val*norm[2]);
            }
         }
         glEnd ();
      }
   }
   glEndList ();
   cout << "VisualizationSceneSolution3d::PrepareFlat2() : [min,max] = ["
        << vmin << "," << vmax << "]" << endl;
}

void VisualizationSceneSolution3d :: Prepare ()
{
   int i,j;

   switch (shading){
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
   Set_Material ();

   int ne = mesh -> GetNBE();
   int nv = mesh -> GetNV();
   DenseMatrix pointmat;
   Array<int> vertices;

   Vector nx(nv);
   Vector ny(nv);
   Vector nz(nv);

   for (int d = 0; d < mesh -> bdr_attributes.Size(); d++)
   {
      if (!bdr_attr_to_show[mesh -> bdr_attributes[d]-1]) continue;

      nx=0.;
      ny=0.;
      nz=0.;

      for (i = 0; i < ne; i++)
         if (mesh -> GetBdrAttribute(i) == mesh -> bdr_attributes[d])
         {
            mesh->GetBdrPointMatrix (i, pointmat);
            mesh->GetBdrElementVertices (i, vertices);

            double v10[] = { pointmat(0,1)-pointmat(0,0),
                             pointmat(1,1)-pointmat(1,0),
                             pointmat(2,1)-pointmat(2,0)};
            double v21[] = { pointmat(0,2)-pointmat(0,1),
                             pointmat(1,2)-pointmat(1,1),
                             pointmat(2,2)-pointmat(2,1) };
            double norm[] = { v10[1]*v21[2]-v10[2]*v21[1],
                              v10[2]*v21[0]-v10[0]*v21[2],
                              v10[0]*v21[1]-v10[1]*v21[0] };
            double rlen = 1.0/sqrt(norm[0]*norm[0]+norm[1]*norm[1]+
                                   norm[2]*norm[2]);

            for (j = 0; j < pointmat.Size(); j++)
            {
               nx(vertices[j]) += norm[0]*rlen;
               ny(vertices[j]) += norm[1]*rlen;
               nz(vertices[j]) += norm[2]*rlen;
            }
         }

      for (i = 0; i < ne; i++)
         if (mesh -> GetBdrAttribute(i) == mesh -> bdr_attributes[d])
         {
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
            switch (mesh->GetBdrElementType(i))
            {
            case Element::TRIANGLE:
               glBegin (GL_TRIANGLES);
               break;

            case Element::QUADRILATERAL:
               glBegin (GL_QUADS);
               break;
            }
            mesh->GetBdrPointMatrix (i, pointmat);

            for (j = 0; j < pointmat.Size(); j++)
            {
               MySetColor ( (*sol)(vertices[j]) , minv , maxv);
               glNormal3d (nx(vertices[j]),ny(vertices[j]),nz(vertices[j]));
               glVertex3d (pointmat.Elem(0, j),
                           pointmat.Elem(1, j),
                           pointmat.Elem(2, j));
            }
            glEnd ();
         }
   }
   glEndList ();
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

   for (i = 0; i < ne; i++){
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
         glEnd ();
         break;

      case 2:
         for (j = 0; j < pointmat.Size(); j++) {
            for (k = 0; k < 3; k++)
               point[j][k] = pointmat(k,j);
            point[j][3] = (*sol)(vertices[j]);
         }
         DrawPolygonLevelLines (point[0], pointmat.Size(), level);
         break;
      }
   }

   glEndList ();
}

void VisualizationSceneSolution3d::PrepareLines2()
{
   int i, j, k, fn, fo, di;
   double bbox_diam;

   glNewList (linelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   int nbe = mesh -> GetNBE();
   DenseMatrix pointmat;
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

      int *RG = &(RefG->RefGeoms[0]);
      double pts[][3] = { { pointmat(0,RG[0]), pointmat(1,RG[0]),
                            pointmat(2,RG[0]) },
                          { pointmat(0,RG[1]), pointmat(1,RG[1]),
                            pointmat(2,RG[1]) },
                          { pointmat(0,RG[2]), pointmat(1,RG[2]),
                            pointmat(2,RG[2]) } };
      double norm[3];
      if (Compute3DUnitNormal (pts[0], pts[1], pts[2], norm))
         cerr << "WARNING: VisualizationSceneSolution3d::PrepareLines2()"
              << endl;
      if (di != 0 && sc != 0.0)
      {
         norm[0] = -norm[0];
         norm[1] = -norm[1];
         norm[2] = -norm[2];
      }

      if (drawmesh == 1)
      {
         Array<int> &REdges = RefG->RefEdges;

         glBegin (GL_LINES);
         for (k = 0; k < REdges.Size(); k++)
         {
            int *RE = &(REdges[k]);

            if (sc == 0.0)
            {
               glVertex3d (pointmat(0, RE[0]), pointmat(1, RE[0]),
                           pointmat(2, RE[0]));
            }
            else
            {
               double val = sc * (values(RE[0]) - minv) / (maxv - minv);
               glVertex3d (pointmat(0, RE[0])+val*norm[0],
                           pointmat(1, RE[0])+val*norm[1],
                           pointmat(2, RE[0])+val*norm[2]);
            }
         }
         glEnd ();
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
            RG = &(RefG->RefGeoms[k*sides]);

            if (sc == 0.0)
            {
               for (j = 0; j < sides; j++)
               {
                  for (int ii = 0; ii < 3; ii++)
                     point[j][ii] = pointmat(ii, RG[j]);
                  point[j][3] = values(RG[j]);
               }
            }
            else
            {
               for (j = 0; j < sides; j++)
               {
                  double val = (values(RG[j]) - minv) / (maxv - minv);
                  val *= sc;
                  for (int ii = 0; ii < 3; ii++)
                     point[j][ii] = pointmat(ii, RG[j])+val*norm[ii];
                  point[j][3] = values(RG[j]);
               }
            }
            DrawPolygonLevelLines (point[0], sides, level);
         }
      }
   }
   glEndList ();
}

void VisualizationSceneSolution3d::CuttingPlaneFunc (int func)
{
   int i, j, k, l, m, n = 0;
   int flag[8], oedges[6];
   static const int tet_edges[12]={0,3, 0,2, 0,1, 1,2, 1,3, 2,3};
   static const int hex_edges[24] =
      { 0,1, 1,2, 3,2, 0,3, 4,5, 5,6, 7,6, 4,7, 0,4, 1,5, 2,6, 3,7 };
   static const int hex_cutting[24][3] =
      { { 3,  5,  7}, { 8, 16, 19},  { 1,  5,  7}, {10, 18, 21},
        { 1,  3,  7}, {12, 20, 23},  { 1,  3,  5}, {14, 17, 22},
        {11, 13, 15}, { 0, 16, 19},  { 9, 13, 15}, { 2, 18, 21},
        { 9, 11, 15}, { 4, 20, 23},  { 9, 11, 13}, { 6, 17, 22},
        { 6, 14, 22}, { 0,  8, 19},  { 0,  8, 16}, { 2, 10, 21},
        { 2, 10, 18}, { 4, 12, 23},  { 4, 12, 20}, { 6, 14, 17} };
   const int *ev;
   double t, point[6][4], norm[3];

   DenseMatrix pointmat;

   Array<int> nodes;
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = 0;             // n will be the number of intersection points
      mesh -> GetElementVertices(i,nodes);
      for (j = 0; j < nodes.Size(); j++)
         if (node_pos[nodes[j]] >= 0.0)
            flag[j] = 1;
         else
            flag[j] = -1;
      switch (mesh -> GetElementType(i))
      {
      case Element::TETRAHEDRON:
         ev = tet_edges;
         for(j=0; j<6; j++, ev += 2)
            if (flag[ev[0]] != flag[ev[1]])
               oedges[n++] = j;
         ev = tet_edges;
         break;
      case Element::HEXAHEDRON:
         ev = hex_edges;
         for (j = 0; j < 12; j++, ev += 2)
            if (flag[ev[0]] != flag[ev[1]])
               break;
         if (j < 12)
         {
            k = 2*j;
            do
            {
               for (l = j = 0; j < 3; j++)
               {
                  ev = hex_edges + 2 * (hex_cutting[k][j] / 2);
                  if (flag[ev[0]] != flag[ev[1]])
                  {
                     m = hex_cutting[k][j];
                     l++;
                     break;
                  }
               }
               for (j++; j < 3; j++)
               {
                  ev = hex_edges + 2 * (hex_cutting[k][j] / 2);
                  if (flag[ev[0]] != flag[ev[1]])
                     l++;
               }
               if (l != 1 || n == 6)
               {
                  cerr << "VisualizationSceneSolution3d::"
                       << "CuttingPlaneFunc ("<< func << ")" << endl;
                  *((int *)NULL) = 0; // force segmentation fault
               }
               oedges[n++] = k/2;
               k = m;
            }
            while (k/2 != oedges[0]);
         }
         ev = hex_edges;
         break;
      default:
         break;
      }

      if (n > 2)
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
            for(k = 0; k < 3; k++)
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
         case 1:  //  PrepareCuttingPlane()
         {
            if (shading == 2)
            {
               // changes point for n > 4
               DrawRefinedSurf (n, point[0], i, 1);
            }
            else
            {
               if (n > 3)
                  j = Compute3DUnitNormal (point[0], point[1], point[2],
                                           point[3], norm);
               else
                  j = Compute3DUnitNormal (point[0], point[1], point[2],
                                           norm);
               if (!j)
               {
                  glBegin (GL_POLYGON);
                  glNormal3dv (norm);
                  for(j = 0; j < n; j++)
                  {
                     MySetColor ( point[j][3] , minv , maxv);
                     glVertex3dv (point[j]);
                  }
                  glEnd();
               }
            }
         }
         break;

         case 2:  //  PrepareCuttingPlaneLines() with mesh
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

         case 3:  //  PrepareCuttingPlaneLines() with level lines
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
      }
   }
}

void VisualizationSceneSolution3d::PrepareCuttingPlane(){

   if (cp_drawelems == 0) return;
   if (cplane == 0) return;
   if (cplane == 2)
   {
      PrepareCuttingPlane2();
      return;
   }

   glNewList(cplanelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

   CuttingPlaneFunc (1);

   glEndList();
}

void VisualizationSceneSolution3d::PrepareCuttingPlane2()
{
   int i, j, n = 0;
   double point[4][4], norm[3];
   DenseMatrix pointmat;
   Vector values;
   RefinedGeometry * RefG;

   glNewList(cplanelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

   double * coord;

   Array<int> nodes;
   Array<int> partition (mesh -> GetNE());
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = 0;   // n will be the number of nodes behind the cutting plane
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
            if (Compute3DUnitNormal (point[0], point[1], point[2], norm))
               cerr << "WARNING: "
                    << "VisualizationSceneSolution3d::PrepareCuttingPlane2()"
                    << endl;

            glBegin (GL_POLYGON);
            for(j = 0; j < nodes.Size(); j++)
            {
               MySetColor (point[j][3], minv, maxv);
               glNormal3dv ( norm );
               glVertex3dv (point[j]);
            }
            glEnd();
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
            DrawRefinedSurf (n, pointmat, values, RefG->RefGeoms);
         } // end shading == 2
   }
   glEndList();
}

void VisualizationSceneSolution3d::PrepareCuttingPlaneLines(){

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

   CuttingPlaneFunc (func);

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
      n = 0;   // n will be the number of nodes behind the cutting plane
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

   ne = mesh -> GetNE ();

   glNewList(lsurflist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

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
               glEnd ();
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
               glEnd ();
            }
         }
      }
   }

   glEndList();
}

void VisualizationSceneSolution3d :: Draw (){
   glEnable(GL_DEPTH_TEST);

   Set_Background ();
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

   Set_Black_Material ();
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
      Set_Black_Material ();
      glDisable(GL_CLIP_PLANE0);
      if ( cp_drawmesh )
         glCallList(cplanelineslist);
      glEnable(GL_CLIP_PLANE0);
   }

   Set_Black_Material ();

   glDisable(GL_LIGHTING);
   // draw lines
   if (drawmesh)
      glCallList(linelist);
   if (light)
      glEnable(GL_LIGHTING);

   if (MatAlpha < 1.0)
      Set_Transparency ();

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
      // Set_Black_Material ();
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
      Remove_Transparency ();

   glFlush();
   glXSwapBuffers (auxXDisplay(), auxXWindow());
}
