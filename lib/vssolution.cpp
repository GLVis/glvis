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
#include <limits>
#include <math.h>

#include <X11/keysym.h>

#include "mfem.hpp"
#include "visual.hpp"


VisualizationSceneSolution *vssol;
extern VisualizationScene  *locscene;

// Definitions of some more keys

static void SolutionKeyHPressed()
{
   cout << endl
        << "+------------------------------------+" << endl
        << "| Keys                               |" << endl
        << "+------------------------------------+" << endl
        << "| a -  Displays/Hides the axes       |" << endl
        << "| A -  Turns antialiasing on/off     |" << endl
        << "| b -  Displays/Hides the boundary   |" << endl
        << "| c -  Displays/Hides the colorbar   |" << endl
        << "| e -  Displays/Hides the elements   |" << endl
        << "| f -  Smooth/Nonconf/Flat shading   |" << endl
        << "| g -  Toggle background             |" << endl
        << "| h -  Displays help menu            |" << endl
        << "| i -  (De)refine elem. (NC shading) |" << endl
        << "| I -  Switch 'i' func. (NC shading) |" << endl
        << "| j -  Turn on/off perspective       |" << endl
        << "| k/K  Adjust the transparency level |" << endl
        << "| ,/<  Adjust color transparency     |" << endl
        << "| l -  Turns on/off the light        |" << endl
        << "| L -  Toggle logarithmic scale      |" << endl
        << "| m -  Displays/Hides the mesh       |" << endl
        << "| p -  Cycle through color palettes  |" << endl
        << "| P -  Print to PostScript file      |" << endl
        << "| q -  Quits                         |" << endl
        << "| r -  Reset the plot to 3D view     |" << endl
        << "| R -  Reset the plot to 2D view     |" << endl
        << "| s -  Turn on/off unit cube scaling |" << endl
        << "| S -  Take snapshot/Record a movie  |" << endl
        << "| t -  Cycle materials and lights    |" << endl
        << "| w -  Toggle the clipping plane     |" << endl
        << "| y/Y  Rotate the clipping plane     |" << endl
        << "| z/Z  Move the clipping plane       |" << endl
        << "| \\ -  Set light source position     |" << endl
        << "+------------------------------------+" << endl
        << "| Function keys                      |" << endl
        << "+------------------------------------+" << endl
        << "| F1 - X window info and keystrokes  |" << endl
        << "| F2 - Update colors, etc.           |" << endl
        << "| F3/F4 - Shrink/Zoom elements       |" << endl
        << "| F5 - Set level lines               |" << endl
        << "| F6 - Palete options                |" << endl
        << "| F7 - Manually set min/max value    |" << endl
        << "| F8 - List of subdomains to show    |" << endl
        << "| F9/F10 - Walk through subdomains   |" << endl
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

static void KeyF8Pressed()
{
   int attr;

   cout << "El attributes ON: ";
   for (attr = 0; attr < vssol -> el_attr_to_show.Size(); attr++)
      if (vssol -> el_attr_to_show[attr])
         cout << " " << attr;
   cout << endl;

   cout << "El attribute to toggle : " << flush;
   cin >> attr;
   if (attr == -1)
   {
      for (attr = 0; attr < vssol -> el_attr_to_show.Size(); attr++)
         vssol -> el_attr_to_show[attr] = 0;
   }
   else if (attr >= vssol -> el_attr_to_show.Size() || attr < 0)
   {
      for (attr = 0; attr < vssol -> el_attr_to_show.Size(); attr++)
         vssol -> el_attr_to_show[attr] = 1;
   }
   else
   {
      vssol -> el_attr_to_show[attr] = !vssol -> el_attr_to_show[attr];
   }
   vssol -> PrepareLines();
   vssol -> Prepare();
   SendExposeEvent();
}

static void KeyF9Pressed()
{
   int attr;

   if (vssol -> attr_to_show == -1) {
      for (attr = 0; attr < vssol -> el_attr_to_show.Size(); attr++)
         vssol -> el_attr_to_show[attr] = 0;
      vssol -> attr_to_show = 0;
   } else
      vssol -> el_attr_to_show[vssol -> attr_to_show++] = 0;

   if (vssol -> attr_to_show >= vssol -> el_attr_to_show.Size())
      vssol -> attr_to_show = 0;
   vssol -> el_attr_to_show[vssol -> attr_to_show] = 1;
   cout << "Showing el attribute " << vssol -> attr_to_show << endl;
   vssol -> PrepareLines();
   vssol -> Prepare();
   SendExposeEvent();
}

static void KeyF10Pressed()
{
   int attr;

   if (vssol -> attr_to_show == -1) {
      for (attr = 0; attr < vssol -> el_attr_to_show.Size(); attr++)
         vssol -> el_attr_to_show[attr] = 0;
   } else
      vssol -> el_attr_to_show[vssol -> attr_to_show--] = 0;

   if (vssol -> attr_to_show < 0)
      vssol -> attr_to_show = vssol -> el_attr_to_show.Size()-1;
   vssol -> el_attr_to_show[vssol -> attr_to_show] = 1;
   cout << "Showing el attribute " << vssol -> attr_to_show << endl;
   vssol -> PrepareLines();
   vssol -> Prepare();
   SendExposeEvent();
}

static void KeyBPressed()
{
   vssol -> ToggleDrawBdr();
   SendExposeEvent();
}

static void KeyMPressed()
{
   vssol -> ToggleDrawMesh();
   SendExposeEvent();
}

static void KeyEPressed()
{
   vssol -> ToggleDrawElems();
   SendExposeEvent();
}

static void KeyFPressed()
{
   vssol -> ToggleShading();
   SendExposeEvent();
}

int refine_func = 0;

void KeyiPressed()
{
   int update = 1;
   switch (refine_func)
   {
   case 0:
      vssol -> TimesToRefine += vssol -> EdgeRefineFactor;
      break;
   case 1:
      if (vssol -> TimesToRefine > vssol -> EdgeRefineFactor)
         vssol -> TimesToRefine -= vssol -> EdgeRefineFactor;
      else
         update = 0;
      break;
   case 2:
      vssol -> TimesToRefine /= vssol -> EdgeRefineFactor;
      vssol -> EdgeRefineFactor ++;
      vssol -> TimesToRefine *= vssol -> EdgeRefineFactor;
      break;
   case 3:
      if (vssol -> EdgeRefineFactor > 1)
      {
         vssol -> TimesToRefine /= vssol -> EdgeRefineFactor;
         vssol -> EdgeRefineFactor --;
         vssol -> TimesToRefine *= vssol -> EdgeRefineFactor;
      }
      else
         update = 0;
      break;
   }
   if (update && vssol -> shading == 2)
   {
      vssol -> FindNewBox();
      locscene -> PrepareAxes();
      vssol -> PrepareLines();
      locscene -> Prepare();
      vssol -> DefaultLevelLines();
      vssol -> PrepareLevelCurves();
      vssol -> PrepareCP();
      SendExposeEvent();
   }
   cout << "Subdivision factors = " << vssol -> TimesToRefine
        << ", " << vssol -> EdgeRefineFactor << endl;
}

void KeyIPressed()
{
   refine_func = (refine_func+1)%4;
   cout << "Key 'i' will: ";
   switch (refine_func)
   {
   case 0:
      cout << "Increase subdivision factor" << endl;
      break;
   case 1:
      cout << "Decrease subdivision factor" << endl;
      break;
   case 2:
      cout << "Increase bdr subdivision factor" << endl;
      break;
   case 3:
      cout << "Decrease bdr subdivision factor" << endl;
      break;
   }
}

static void KeywPressed()
{
   vssol->ToggleDrawCP();
   SendExposeEvent();
}

static void KeyyPressed()
{
   vssol->CuttingPlane->IncreaseTheta();
   vssol->PrepareCP();
   SendExposeEvent();
}

static void KeyYPressed()
{
   vssol->CuttingPlane->DecreaseTheta();
   vssol->PrepareCP();
   SendExposeEvent();
}

static void KeyzPressed()
{
   vssol->CuttingPlane->IncreaseDistance();
   vssol->PrepareCP();
   SendExposeEvent();
}

static void KeyZPressed()
{
   vssol->CuttingPlane->DecreaseDistance();
   vssol->PrepareCP();
   SendExposeEvent();
}

static void KeyF3Pressed()
{
   if (vssol->shading == 2)
   {
      vssol->shrink *= 0.9;
      vssol->Prepare();
      vssol->PrepareLines();
      SendExposeEvent();
   }
}

static void KeyF4Pressed()
{
   if (vssol->shading == 2)
   {
      vssol->shrink *= 1.11111111111111111111111;
      vssol->Prepare();
      vssol->PrepareLines();
      SendExposeEvent();
   }
}

VisualizationSceneSolution::VisualizationSceneSolution()
{}

VisualizationSceneSolution::VisualizationSceneSolution(Mesh &m, Vector &s)
{
   mesh = &m;
   sol = &s;

   Init();

   auxKeyFunc (AUX_h, SolutionKeyHPressed);
   auxKeyFunc (AUX_H, SolutionKeyHPressed);
}

void VisualizationSceneSolution::Init()
{
   rsol  = NULL;
   vssol = this;

   drawelems = shading = 1;
   drawmesh  = 0;

   shrink = 1.0;

   TimesToRefine = EdgeRefineFactor = 1;

   attr_to_show = -1;
   el_attr_to_show.SetSize (mesh->attributes.Max());
   for (int i = 0; i < el_attr_to_show.Size(); i++)
      el_attr_to_show[i] = 1;
   drawbdr = 0;

   VisualizationSceneScalarData::Init();  //  Calls FindNewBox() !!!

   SetUseTexture(1);

   CuttingPlane = new Plane(-1.0,0.0,0.0,0.5*(x[0]+x[1]));
   draw_cp = 0;

   //  static int init = 0;
   //  if (!init)
   {
      //      init = 1;

      auxKeyFunc (AUX_b, KeyBPressed);
      auxKeyFunc (AUX_B, KeyBPressed);

      auxKeyFunc (AUX_m, KeyMPressed);
      auxKeyFunc (AUX_M, KeyMPressed);

      auxKeyFunc (AUX_e, KeyEPressed);
      auxKeyFunc (AUX_E, KeyEPressed);

      auxKeyFunc (AUX_f, KeyFPressed);
      auxKeyFunc (AUX_F, KeyFPressed);

      auxKeyFunc (AUX_i, KeyiPressed);
      auxKeyFunc (AUX_I, KeyIPressed);

      auxKeyFunc (AUX_w, KeywPressed);
      auxKeyFunc (AUX_y, KeyyPressed);
      auxKeyFunc (AUX_Y, KeyYPressed);
      auxKeyFunc (AUX_z, KeyzPressed);
      auxKeyFunc (AUX_Z, KeyZPressed);

      auxKeyFunc (XK_F3, KeyF3Pressed);
      auxKeyFunc (XK_F4, KeyF4Pressed);
      auxKeyFunc (XK_F8, KeyF8Pressed);
      auxKeyFunc (XK_F9, KeyF9Pressed);
      auxKeyFunc (XK_F10, KeyF10Pressed);
   }

   displlist  = glGenLists (1);
   linelist   = glGenLists (1);
   lcurvelist = glGenLists (1);
   bdrlist    = glGenLists (1);
   cp_list    = glGenLists (1);

   Prepare();
   PrepareLines();
   DefaultLevelLines();
   PrepareLevelCurves();
   PrepareBoundary();
}

VisualizationSceneSolution::~VisualizationSceneSolution()
{
   glDeleteLists (displlist, 1);
   glDeleteLists (linelist, 1);
   glDeleteLists (lcurvelist, 1);
   glDeleteLists (bdrlist, 1);
   glDeleteLists (cp_list, 1);
}


void VisualizationSceneSolution::NewMeshAndSolution(
   Mesh *new_m, Vector *new_sol, GridFunction *new_u, int rescale)
{
   mesh = new_m;
   sol = new_sol;
   rsol = new_u;

   if (rescale == 1)
   {
      FindNewBox();
      PrepareAxes();
   }
   else if (rescale == 2)
   {
      double rx[2], ry[2], rv[2];
      FindNewBox(rx, ry, rv);
      minv = rv[0];
      maxv = rv[1];
      FixValueRange();
      NewZRange(minv, maxv);
      PrepareAxes();
   }

   Prepare();
   PrepareLines();
   PrepareLevelCurves();
   PrepareBoundary();
   PrepareCP();
}


void VisualizationSceneSolution::GetRefinedDetJ(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr)
{
   ElementTransformation *T = mesh->GetElementTransformation(i);

   T->Transform(ir, tr);

   vals.SetSize(ir.GetNPoints());
   for (int j = 0; j < ir.GetNPoints(); j++)
   {
      T->SetIntPoint(&ir.IntPoint(j));
      vals(j) = T->Weight();
   }

   if (drawelems == 2)
   {
      for (int j = 0; j < vals.Size(); j++)
      {
         if (vals(j) <= 0.0)
         {
            vals = 0.0;
            break;
         }
         vals(j) = 1.0 / vals(j);
      }
   }
}

void VisualizationSceneSolution::GetRefinedValues(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr)
{
   if (drawelems < 2)
      rsol->GetValues(i, ir, vals, tr);
   else
      GetRefinedDetJ(i, ir, vals, tr);

   if (shrink != 1.0)
      ShrinkPoints(tr);
}

int VisualizationSceneSolution::GetRefinedValuesAndNormals(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr,
   DenseMatrix &normals)
{
   int have_normals = 0;

   if (drawelems < 2)
   {
      rsol->GetGradients(i, ir, tr);
      normals.SetSize(3, tr.Width());
      for (int j = 0; j < tr.Width(); j++)
      {
         normals(0, j) = -tr(0, j);
         normals(1, j) = -tr(1, j);
         normals(2, j) = 1.;
      }
      have_normals = 1;
      rsol->GetValues(i, ir, vals, tr);
   }
   else
      GetRefinedDetJ(i, ir, vals, tr);

   if (shrink != 1.0)
   {
      ShrinkPoints(tr);
      if (have_normals)
      {
         for (int j = 0; j < tr.Width(); j++)
         {
            normals(0, j) /= shrink;
            normals(1, j) /= shrink;
         }
      }
   }

   return have_normals;
}

void VisualizationSceneSolution::SetShading(int s)
{
   if (shading == s || s < 0)
      return;

   if (rsol)
   {
      if (s > 2)
         return;

      if (s == 2 || shading == 2)
      {
         shading = s;
         PrepareLines();
         PrepareLevelCurves();
         PrepareCP();
      }
      else
         shading = s;
   }
   else
   {
      if (s > 1)
         return;
      shading = s;
   }

   Prepare();
}

void VisualizationSceneSolution::ToggleShading()
{
   if (rsol)
   {
      shading = (shading + 1) % 3;
      if (shading != 1)
      {
         FindNewBox();
         PrepareAxes();
         PrepareLines();
         PrepareLevelCurves();
         PrepareCP();
      }
      Prepare();
   }
   else
      SetShading(1 - shading);

   static const char *shading_type[3] =
      {"flat", "smooth", "non-conforming (with subdivision)"};
   cout << "Shading type : " << shading_type[shading] << endl;
}

void VisualizationSceneSolution::SetRefineFactors(int tot, int bdr)
{
   if ((tot == TimesToRefine && bdr == EdgeRefineFactor) || tot < 1 || bdr < 1)
      return;

   if (tot % bdr)
      tot += bdr - tot % bdr;

   TimesToRefine = tot;
   EdgeRefineFactor = bdr;

   if (shading == 2)
   {
      PrepareAxes();
      PrepareLines();
      Prepare();
      PrepareLevelCurves();
      PrepareCP();
   }
}

void VisualizationSceneSolution::NewZRange(double z1, double z2)
{
   // preserve the current box z-size
   zscale *= (z[1]-z[0])/(z2-z1);
   z[0] = z1;
   z[1] = z2;
}

void VisualizationSceneSolution::SetNewScalingFromBox()
{
   if (scaling)
      VisualizationSceneScalarData::SetNewScalingFromBox();
   else
   {
      xscale = x[1]-x[0];
      yscale = y[1]-y[0];
      zscale = z[1]-z[0];
      if (xscale < yscale)
      {
         if (xscale > 0.0)
            zscale *= (yscale/xscale);
         xscale = yscale;
      }
      else
      {
         if (yscale > 0.0)
            zscale *= (xscale/yscale);
      }
      xscale = (xscale > 0.0) ? ( 1.0 / xscale ) : 1.0;
      yscale = xscale;
      zscale = (zscale > 0.0) ? ( 1.0 / zscale ) : 1.0;
   }
   zscale /= ((1. + sqrt(5.)) / 2.);
}

void VisualizationSceneSolution::FindNewBox(double rx[], double ry[],
                                            double rval[])
{
   int i, j;

   if (shading != 2)
   {
      int nv = mesh -> GetNV();

      double *coord = mesh->GetVertex(0);

      rval[0] = rval[1] = (*sol)(0);
      for (i = 1; i < sol->Size(); i++)
      {
         if ((*sol)(i) < rval[0]) rval[0] = (*sol)(i);
         if ((*sol)(i) > rval[1]) rval[1] = (*sol)(i);
      }
      rx[0] = rx[1] = coord[0];
      ry[0] = ry[1] = coord[1];

      for(i = 1; i < nv; i++)
      {
         coord = mesh->GetVertex(i);
         if (coord[0] < rx[0]) rx[0] = coord[0];
         if (coord[1] < ry[0]) ry[0] = coord[1];
         if (coord[0] > rx[1]) rx[1] = coord[0];
         if (coord[1] > ry[1]) ry[1] = coord[1];
      }
   }
   else
   {
      int ne = mesh -> GetNE();
      DenseMatrix pointmat;
      Vector values;
      RefinedGeometry *RefG;

      for (i = 0; i < ne; i++)
      {
         RefG = GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                            TimesToRefine, EdgeRefineFactor);
         GetRefinedValues (i, RefG->RefPts, values, pointmat);
         if (i == 0)
         {
            rx[0] = rx[1] = pointmat(0,0);
            ry[0] = ry[1] = pointmat(1,0);
            rval[0] = rval[1] = values(0);
         }
         for (j = 0; j < values.Size(); j++)
         {
            if (pointmat(0,j) < rx[0]) rx[0] = pointmat(0,j);
            if (pointmat(0,j) > rx[1]) rx[1] = pointmat(0,j);
            if (pointmat(1,j) < ry[0]) ry[0] = pointmat(1,j);
            if (pointmat(1,j) > ry[1]) ry[1] = pointmat(1,j);
            if (values(j) < rval[0]) rval[0] = values(j);
            if (values(j) > rval[1]) rval[1] = values(j);
         }
      }
   }
}

void VisualizationSceneSolution::FixValueRange()
{
   double am = fabs(minv);
   if (am < fabs(maxv)) am = fabs(maxv);
   if (float(am) < 100*numeric_limits<float>::min()) am = 1e-3;
   if ((maxv-minv) < 1e-5*am)
   {
      // Shading quality may be bad since OpenGL uses single precision.
      // We should probably pre-scale the solution before feeding it
      // to OpenGL
      int old_prec = cout.precision(12);
      cout << "[minv,maxv] = " << "[" << minv << "," << maxv
           << "] (maxv-minv = " << maxv-minv << ")\n --> ";
      minv -= 0.49999e-5*am;
      maxv += 0.50001e-5*am;
      cout << "[" << minv << "," << maxv << "]" << endl;
      cout.precision(old_prec);
   }
}

void VisualizationSceneSolution::FindNewBox()
{
   FindNewBox(x, y, z);

   minv = z[0];
   maxv = z[1];

   FixValueRange();

   z[0] = minv;
   z[1] = maxv;

   SetNewScalingFromBox();
}

void DrawTriangle(const double pts[][3], const double cv[],
                  const double minv, const double maxv)
{
   double nor[3];
   if (Compute3DUnitNormal(pts[0], pts[1], pts[2], nor))
      return;
   glBegin(GL_TRIANGLES);
   glNormal3dv(nor);
   for (int j = 0; j < 3; j++)
   {
      MySetColor(cv[j], minv, maxv);
      glVertex3dv(pts[j]);
   }
   glEnd();
}

void DrawQuad(const double pts[][3], const double cv[],
              const double minv, const double maxv)
{
   double nor[3];
   if (Compute3DUnitNormal(pts[0], pts[1], pts[2], pts[3], nor))
      return;
   glBegin(GL_QUADS);
   glNormal3dv(nor);
   for (int j = 0; j < 4; j++)
   {
      MySetColor(cv[j], minv, maxv);
      glVertex3dv(pts[j]);
   }
   glEnd();
}

void DrawPatch(const DenseMatrix &pts, Vector &vals, DenseMatrix &normals,
               const int n, const Array<int> &ind, const double minv,
               const double maxv, const int normals_opt)
{
   double na[3];

   if (normals_opt == 1 || normals_opt == -2)
   {
      normals.SetSize(3, pts.Width());
      normals = 0.;
      for (int i = 0; i < ind.Size(); i += n)
      {
         int j;
         if (n == 3)
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), na);
         else
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), &pts(0, ind[i+3]), na);
         if (j == 0)
            for ( ; j < n; j++)
               for (int k = 0; k < 3; k++)
                  normals(k, ind[i+j]) += na[k];
      }
   }

   if (n == 3)
      glBegin(GL_TRIANGLES);
   else
      glBegin(GL_QUADS);
   if (normals_opt != 0 && normals_opt != -1)
   {
      if (normals_opt > 0)
      {
         for (int i = 0; i < ind.Size(); i++)
         {
            glNormal3dv(&normals(0, ind[i]));
            MySetColor(vals(ind[i]), minv, maxv);
            glVertex3dv(&pts(0, ind[i]));
         }
      }
      else
      {
         for (int i = ind.Size()-1; i >= 0; i--)
         {
            glNormal3dv(&normals(0, ind[i]));
            MySetColor(vals(ind[i]), minv, maxv);
            glVertex3dv(&pts(0, ind[i]));
         }
      }
   }
   else
   {
      for (int i = 0; i < ind.Size(); i += n)
      {
         int j;
         if (n == 3)
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), na);
         else
            j = Compute3DUnitNormal(&pts(0, ind[i]), &pts(0, ind[i+1]),
                                    &pts(0, ind[i+2]), &pts(0, ind[i+3]), na);
         if (j == 0)
         {
            if (normals_opt == 0)
            {
               glNormal3dv(na);
               for ( ; j < n; j++)
               {
                  MySetColor(vals(ind[i+j]), minv, maxv);
                  glVertex3dv(&pts(0, ind[i+j]));
               }
            }
            else
            {
               glNormal3d(-na[0], -na[1], -na[2]);
               for (j = n-1; j >= 0; j--)
               {
                  MySetColor(vals(ind[i+j]), minv, maxv);
                  glVertex3dv(&pts(0, ind[i+j]));
               }
            }
         }
      }
   }
   glEnd();
}

void VisualizationSceneSolution::PrepareFlat()
{
   int i, j;

   glNewList (displlist, GL_COMPILE);

   Set_Material();

   int ne = mesh -> GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double pts[4][3], col[4];

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) continue;

      mesh->GetPointMatrix (i, pointmat);
      mesh->GetElementVertices (i, vertices);

      for (j = 0; j < pointmat.Width(); j++)
      {
         pts[j][0] = pointmat(0, j);
         pts[j][1] = pointmat(1, j);
         pts[j][2] = col[j] = (*sol)(vertices[j]);
      }
      if (j == 3)
         DrawTriangle(pts, col, minv, maxv);
      else
         DrawQuad(pts, col, minv, maxv);
   }

   glEndList();
}

// determines how quads and their level lines are drawn
// when using subdivision:
// 0 - draw a quad
// 1 - draw 2 triangles (split using the '0-2' diagonal)
// 2 - draw 4 triangles (split using both diagonals)
const int split_quads = 0;

void VisualizationSceneSolution::PrepareFlat2()
{
   int i, j, k;

   glNewList (displlist, GL_COMPILE);

   Set_Material();

   int ne = mesh -> GetNE();
   DenseMatrix pointmat, pts3d, normals;
   Vector values;
   RefinedGeometry *RefG;

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) continue;

      RefG = GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      j = GetRefinedValuesAndNormals(i, RefG->RefPts, values, pointmat,
                                     normals);
      Array<int> &RG = RefG->RefGeoms;
      int sides = mesh->GetElement(i)->GetNVertices();

#if 1
      pts3d.SetSize(3, pointmat.Width());
      for (k = 0; k < pointmat.Width(); k++)
      {
         pts3d(0, k) = pointmat(0, k);
         pts3d(1, k) = pointmat(1, k);
         pts3d(2, k) = values(k);
      }
      j = (j != 0) ? 2 : 0;
      DrawPatch(pts3d, values, normals, sides, RG, minv, maxv, j);
#else
      for (k = 0; k < RG.Size()/sides; k++)
      {
         int *ind = &RG[sides*k];
         if (split_quads == 0 || sides == 3)
         {
            for (j = 0; j < sides; j++)
            {
               pts[j][0] = pointmat(0, ind[j]);
               pts[j][1] = pointmat(1, ind[j]);
               pts[j][2] = col[j] = values(ind[j]);
            }
            if (sides == 3)
               DrawTriangle(pts, col, minv, maxv);
            else
               DrawQuad(pts, col, minv, maxv);
         }
         else if (split_quads == 1)
         {
            // draw 2 triangles for each quad
            // (split with the 0-2 diagonal)
            const int vt[2][3] = {{ 0, 1, 2 }, { 2, 3, 0 }};
            for (int it = 0; it < 2; it++)
            {
               for (j = 0; j < 3; j++)
               {
                  pts[j][0] = pointmat(0, ind[vt[it][j]]);
                  pts[j][1] = pointmat(1, ind[vt[it][j]]);
                  pts[j][2] = col[j] = values(ind[vt[it][j]]);
               }
               DrawTriangle(pts, col, minv, maxv);
            }
         }
         else
         {
            // draw 4 triangles for each quad
            // (split with both diagonals)
            pts[2][0] = pts[2][1] = pts[2][2] = 0.0;
            for (j = 0; j < 4; j++)
            {
               pts[2][0] += pointmat(0, ind[j]);
               pts[2][1] += pointmat(1, ind[j]);
               pts[2][2] += values(ind[j]);
            }
            pts[2][0] *= 0.25;
            pts[2][1] *= 0.25;
            pts[2][2] *= 0.25;
            col[2] = pts[2][2];
            for (j = 0; j < 4; j++)
            {
               pts[0][0] = pointmat(0, ind[j]);
               pts[0][1] = pointmat(1, ind[j]);
               pts[0][2] = col[0] = values(ind[j]);
               int l = (j+1)%4;
               pts[1][0] = pointmat(0, ind[l]);
               pts[1][1] = pointmat(1, ind[l]);
               pts[1][2] = col[1] = values(ind[l]);
               DrawTriangle(pts, col, minv, maxv);
            }
         }
      }
#endif
   }

   glEndList();
}

void VisualizationSceneSolution::Prepare()
{
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

   int i, j;

   glNewList (displlist, GL_COMPILE);

   Set_Material();

   int ne = mesh -> GetNE();
   int nv = mesh -> GetNV();
   DenseMatrix pointmat;
   Array<int> vertices;
   double p[4][3], nor[3];

   Vector nx(nv);
   Vector ny(nv);
   Vector nz(nv);

   for (int d = 0; d < mesh -> attributes.Size(); d++)
   {
      if (!el_attr_to_show[mesh -> attributes[d]-1]) continue;

      nx = 0.;
      ny = 0.;
      nz = 0.;

      for (i = 0; i < ne; i++)
         if (mesh -> GetAttribute(i) == mesh -> attributes[d])
         {
            mesh->GetPointMatrix (i, pointmat);
            mesh->GetElementVertices (i, vertices);

            for (j = 0; j < pointmat.Size(); j++)
            {
               p[j][0] = pointmat(0, j);
               p[j][1] = pointmat(1, j);
               p[j][2] = (*sol)(vertices[j]);
            }

            if (pointmat.Width() == 3)
               j = Compute3DUnitNormal(p[0], p[1], p[2], nor);
            else
               j = Compute3DUnitNormal(p[0], p[1], p[2], p[3], nor);

            if (j == 0)
               for (j = 0; j < pointmat.Size(); j++)
               {
                  nx(vertices[j]) += nor[0];
                  ny(vertices[j]) += nor[1];
                  nz(vertices[j]) += nor[2];
               }
         }

      for (i = 0; i < ne; i++)
         if (mesh -> GetAttribute(i) == mesh -> attributes[d])
         {
            switch (mesh->GetElementType(i))
            {
            case Element::TRIANGLE:
               glBegin (GL_TRIANGLES);
               break;

            case Element::QUADRILATERAL:
               glBegin (GL_QUADS);
               break;
            }
            mesh->GetPointMatrix (i, pointmat);
            mesh->GetElementVertices (i, vertices);

            for (j = 0; j < pointmat.Size(); j++)
            {
               MySetColor((*sol)(vertices[j]), minv, maxv);
               glNormal3d(nx(vertices[j]), ny(vertices[j]), nz(vertices[j]));
               glVertex3d(pointmat(0, j), pointmat(1, j), (*sol)(vertices[j]));
            }
            glEnd();
         }
   }

   glEndList();
}

void VisualizationSceneSolution::PrepareLevelCurves()
{
   if (shading == 2)
   {
      PrepareLevelCurves2();
      return;
   }

   int i, j, ne = mesh -> GetNE();
   Array<int> vertices;
   double point[4][4];

   glNewList(lcurvelist, GL_COMPILE);

   for (i = 0; i < ne; i++)
   {
      mesh -> GetElementVertices (i, vertices);

      for (j = 0; j < vertices.Size(); j++)
      {
         double *v = mesh -> GetVertex (vertices[j]);
         point[j][0] = v[0];
         point[j][1] = v[1];
         point[j][2] = point[j][3] = (*sol)(vertices[j]);
      }
      DrawPolygonLevelLines (point[0], vertices.Size(), level);
   }

   glEndList();
}

void VisualizationSceneSolution::DrawLevelCurves(
   Array<int> &RG, DenseMatrix &pointmat, Vector &values,
   int sides, Array<double> &lvl, int flat)
{
   double point[4][4];
   // double zc = 0.5*(z[0]+z[1]);
   double zc = z[1];

   for (int k = 0; k < RG.Size()/sides; k++)
   {
      if (split_quads == 0 || sides == 3)
      {
         for (int j = 0; j < sides; j++)
         {
            int vv = RG[sides*k+j];
            point[j][0] = pointmat(0, vv);
            point[j][1] = pointmat(1, vv);
            point[j][3] = values(vv);
            point[j][2] = (flat) ? zc : point[j][3];
         }
         DrawPolygonLevelLines(point[0], sides, lvl);
      }
      else if (split_quads == 1)
      {
         // split the quad into 2 triangles
         // (with the 0-2 diagonal)
         int *ind = &RG[sides*k];
         const int vt[2][3] = {{ 0, 1, 2 }, { 2, 3, 0 }};
         for (int it = 0; it < 2; it++)
         {
            for (int j = 0; j < 3; j++)
            {
               point[j][0] = pointmat(0, ind[vt[it][j]]);
               point[j][1] = pointmat(1, ind[vt[it][j]]);
               point[j][3] = values(ind[vt[it][j]]);
               point[j][2] = (flat) ? zc : point[j][3];
            }
            DrawPolygonLevelLines(point[0], 3, lvl);
         }
      }
      else
      {
         // split the quad into 4 triangles
         // (with the two diagonals)
         int *ind = &RG[sides*k];
         point[2][0] = point[2][1] = point[2][2] = 0.0;
         for (int j = 0; j < 4; j++)
         {
            point[2][0] += pointmat(0, ind[j]);
            point[2][1] += pointmat(1, ind[j]);
            point[2][2] += values(ind[j]);
         }
         point[2][0] *= 0.25;
         point[2][1] *= 0.25;
         point[2][2] *= 0.25;
         point[2][3] = point[2][2];
         if (flat)
            point[2][2] = zc;

         for (int j = 0; j < 4; j++)
         {
            point[0][0] = pointmat(0, ind[j]);
            point[0][1] = pointmat(1, ind[j]);
            point[0][3] = values(ind[j]);
            point[0][2] = (flat) ? zc : point[0][3];
            int l = (j+1)%4;
            point[1][0] = pointmat(0, ind[l]);
            point[1][1] = pointmat(1, ind[l]);
            point[1][3] = values(ind[l]);
            point[1][2] = (flat) ? zc : point[1][3];

            DrawPolygonLevelLines(point[0], 3, lvl);
         }
      }
   }
}

void VisualizationSceneSolution::PrepareLevelCurves2()
{
   int i, ne = mesh -> GetNE();
   Vector values;
   DenseMatrix pointmat;
   RefinedGeometry *RefG;

   glNewList(lcurvelist, GL_COMPILE);

   for (i = 0; i < ne; i++)
   {
      RefG = GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      GetRefinedValues (i, RefG->RefPts, values, pointmat);
      Array<int> &RG = RefG->RefGeoms;
      int sides = mesh->GetElement(i)->GetNVertices();

      DrawLevelCurves(RG, pointmat, values, sides, level);
   }

   glEndList();
}

void VisualizationSceneSolution::PrepareLines()
{
   if (shading == 2)
   {
      // PrepareLines2();
      PrepareLines3();
      return;
   }

   int i, j, ne = mesh -> GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;

   glNewList(linelist, GL_COMPILE);

   // Set_Black_Material();
   // glColor3f (0, 0, 0);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   for (i = 0; i < ne; i++){
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) continue;

      switch (mesh->GetElementType(i))
      {
      case Element::TRIANGLE:
         glBegin (GL_TRIANGLES);
         break;

      case Element::QUADRILATERAL:
         glBegin (GL_QUADS);
         break;
      }
      mesh->GetPointMatrix (i, pointmat);
      mesh->GetElementVertices (i, vertices);

      for (j = 0; j < pointmat.Size(); j++)
         glVertex3d (pointmat(0, j), pointmat(1, j),
                     (*sol)(vertices[j]));
      glEnd();
   }

   glEndList();
}

void VisualizationSceneSolution::PrepareLines2()
{
   int i, j, k, ne = mesh -> GetNE();
   Vector values;
   DenseMatrix pointmat;
   RefinedGeometry *RefG;

   glNewList(linelist, GL_COMPILE);

   // Set_Black_Material();
   // glColor3f (0, 0, 0);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) continue;

      RefG = GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      GetRefinedValues (i, RefG->RefPts, values, pointmat);
      Array<int> &RG = RefG->RefGeoms;
      int sides = mesh->GetElement(i)->GetNVertices();

      for (k = 0; k < RG.Size()/sides; k++)
      {
         switch (sides)
         {
         case 3:
            glBegin (GL_TRIANGLES);
            break;

         case 4:
            glBegin (GL_QUADS);
            break;
         }

         for (j = 0; j < sides; j++)
            glVertex3d (pointmat(0, RG[sides*k+j]),
                        pointmat(1, RG[sides*k+j]),
                        values(RG[sides*k+j]));
         glEnd();
      }
   }

   glEndList();
}

void VisualizationSceneSolution::PrepareLines3()
{
   int i, k, ne = mesh -> GetNE();
   Vector values;
   DenseMatrix pointmat;
   RefinedGeometry *RefG;

   glNewList(linelist, GL_COMPILE);

   // Set_Black_Material();
   // glColor3f (0, 0, 0);
   // glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) continue;
      RefG = GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      GetRefinedValues (i, RefG->RefPts, values, pointmat);
      Array<int> &RE = RefG->RefEdges;

      glBegin (GL_LINES);
      for (k = 0; k < RE.Size()/2; k++)
      {
         glVertex3d (pointmat(0, RE[2*k]),
                     pointmat(1, RE[2*k]),
                     values(RE[2*k]));
         glVertex3d (pointmat(0, RE[2*k+1]),
                     pointmat(1, RE[2*k+1]),
                     values(RE[2*k+1]));
      }
      glEnd();
   }

   glEndList();
}

void VisualizationSceneSolution::UpdateValueRange()
{
   SetLevelLines (minv, maxv, nl);
   UpdateLevelLines();
   NewZRange(minv, maxv);
   PrepareAxes();
}

void VisualizationSceneSolution::PrepareBoundary()
{
   int i, j, ne = mesh -> GetNBE();
   DenseMatrix pointmat;

   glNewList(bdrlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   Array<int> vertices;

   for (i = 0; i < ne; i++){
      mesh -> GetBdrElementVertices (i, vertices);
      mesh -> GetBdrPointMatrix (i, pointmat);
      glBegin (GL_LINE_LOOP);
      for (j = 0; j < pointmat.Size(); j++)
         glVertex3d (pointmat(0, j), pointmat(1, j), (*sol)(vertices[j]));
      glEnd();
   }

   glEndList();
}

void VisualizationSceneSolution::PrepareCP()
{
   Vector values;
   DenseMatrix pointmat;
   Array<int> ind;

   if (draw_cp == 0)
      return;

   glNewList(cp_list, GL_COMPILE);
   glBegin(GL_LINES);

   if (shading != 2)
   {
      Array<int> vertices;

      for (int i = 0; i < mesh->GetNE(); i++)
      {
         mesh->GetPointMatrix(i, pointmat);
         int n = 0;
         for (int j = 0; j < pointmat.Width(); j++)
         {
            const double s =
               CuttingPlane->Transform(pointmat(0, j),
                                       pointmat(1, j), 0.0);
            if (s >= 0.0)
               n++;
         }
         if (n == 0 || n == pointmat.Width())
            continue;

         mesh->GetElementVertices(i, vertices);
         values.SetSize(vertices.Size());
         ind.SetSize(vertices.Size());
         for (int j = 0; j < values.Size(); j++)
         {
            values(j) = (*sol)(vertices[j]);
            ind[j] = j;
         }

         DrawCPLine(pointmat, values, ind);
      }
   }
   else
   {
      RefinedGeometry *RefG;

      for (int i = 0; i < mesh->GetNE(); i++)
      {
         RefG = GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                            TimesToRefine, EdgeRefineFactor);
         GetRefinedValues (i, RefG->RefPts, values, pointmat);
         Array<int> &RG = RefG->RefGeoms;
         int sides = mesh->GetElement(i)->GetNVertices();

         ind.SetSize(sides);
         for (int k = 0; k < RG.Size()/sides; k++)
         {
            for (int j = 0; j < sides; j++)
               ind[j] = RG[k*sides+j];

            int n = 0;
            for (int j = 0; j < sides; j++)
            {
               const double s =
                  CuttingPlane->Transform(pointmat(0, ind[j]),
                                          pointmat(1, ind[j]), 0.0);
               if (s >= 0.0)
                  n++;
            }
            if (n == 0 || n == sides)
               continue;

            DrawCPLine(pointmat, values, ind);
         }
      }
   }

   glEnd();
   glEndList();
}

void VisualizationSceneSolution::DrawCPLine(
   DenseMatrix &pointmat, Vector &values, Array<int> &ind)
{
   int n, js, nv = ind.Size();
   double s, xs, ys;

   js = nv-1;
   xs = pointmat(0, ind[js]);
   ys = pointmat(1, ind[js]);
   s = CuttingPlane->Transform(xs, ys, 0.0);
   n = 0;
   for (int j = 0; j < nv; j++)
   {
      const double xt = pointmat(0, ind[j]);
      const double yt = pointmat(1, ind[j]);
      const double t = CuttingPlane->Transform(xt, yt, 0.0);
      if ((s >= 0.0 && t < 0.0) || (s < 0.0 && t >= 0.0))
      {
         double a = fabs(s) / (fabs(s) + fabs(t));

         glVertex3d((1.-a) * xs + a * xt,
                    (1.-a) * ys + a * yt,
                    (1.-a) * values(ind[js]) + a * values(ind[j]));
         n++;
      }
      s = t;
      js = j;
      xs = xt;
      ys = yt;
   }
   if (n != 2 && n != 4)
   {
      cerr << "n = " << n << endl;
      mfem_error("VisualizationSceneSolution::DrawCPLine");
   }
}

void VisualizationSceneSolution::Draw()
{
   glEnable(GL_DEPTH_TEST);

   Set_Background();
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // model transformation
   glMatrixMode (GL_MODELVIEW);
   glLoadIdentity();

   glMultMatrixd (translmat);
   glMultMatrixd (rotmat);
   glScaled(xscale, yscale, zscale);
   glTranslated(-(x[0]+x[1])/2, -(y[0]+y[1])/2, -(z[0]+z[1])/2);

   glPolygonOffset (1, 1);
   glEnable (GL_POLYGON_OFFSET_FILL);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   glDisable(GL_CLIP_PLANE0);
   // draw colorbar
   glDisable(GL_LIGHTING);
   if (colorbar)
      if (drawmesh == 2)
         DrawColorBar(minv,maxv,&level);
      else
         DrawColorBar(minv,maxv);

   Set_Black_Material();

   // draw axes
   if (drawaxes)
   {
      glCallList(axeslist);
      DrawCoordinateCross();
   }

   DrawRuler();

   if (draw_cp)
   {
      glClipPlane (GL_CLIP_PLANE0, CuttingPlane->Equation());
      glDisable(GL_CLIP_PLANE0);
      glCallList(cp_list);
      glEnable(GL_CLIP_PLANE0);
   }

   if (drawbdr)
      glCallList(bdrlist);

   // draw lines
   if (drawmesh == 1)
      glCallList(linelist);
   else
      if (drawmesh == 2)
      {
         // glLineWidth(1.0);
         glCallList(lcurvelist);
      }
   if (light)
      glEnable(GL_LIGHTING);

   if (MatAlpha < 1.0)
      Set_Transparency();

   // draw elements
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
   if (drawelems)
   {
      if (GetUseTexture())
      {
         glEnable (GL_TEXTURE_1D);
         glColor4d(1, 1, 1, 1);
      }
      glCallList(displlist);
      if (GetUseTexture())
         glDisable (GL_TEXTURE_1D);
   }

   if (MatAlpha < 1.0)
      Remove_Transparency();

   glFlush();
   glXSwapBuffers (auxXDisplay(), auxXWindow());
}
