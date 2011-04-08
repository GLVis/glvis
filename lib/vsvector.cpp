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

#include "mfem.hpp"
#include "visual.hpp"


static void VectorKeyHPressed()
{
   cout << endl
        << "+------------------------------------+" << endl
        << "| Keys                               |" << endl
        << "+------------------------------------+" << endl
        << "| a -  Displays/Hides the axes       |" << endl
        << "| A -  Turns antialiasing on/off     |" << endl
        << "| b -  Displacements step back       |" << endl
        << "| c -  Displays/Hides the colorbar   |" << endl
        << "| d -  Displays/Hides displacements  |" << endl
        << "| e -  Displays/Hides the elements   |" << endl
        << "| f -  Smooth/Flat shading           |" << endl
        << "| g -  Toggle background             |" << endl
        << "| h -  Displays help menu            |" << endl
        << "| j -  Turn on/off perspective       |" << endl
        << "| k/K  Adjust the transparency level |" << endl
        << "| ,/<  Adjust color transparency     |" << endl
        << "| l -  Turns on/off the light        |" << endl
        << "| L -  Toggle logarithmic scale      |" << endl
        << "| m -  Displays/Hides the mesh       |" << endl
        << "| n -  Displacements step forward    |" << endl
        << "| p -  Cycle through color palettes  |" << endl
        << "| P -  Print to PostScript file      |" << endl
        << "| q -  Quits                         |" << endl
        << "| r -  Reset the plot to 3D view     |" << endl
        << "| R -  Reset the plot to 2D view     |" << endl
        << "| s -  Turn on/off unit cube scaling |" << endl
        << "| S -  Take snapshot/Record a movie  |" << endl
        << "| t -  Cycle materials and lights    |" << endl
        << "| \\ -  Set light source position     |" << endl
        << "| v -  Cycle through vector fields   |" << endl
        << "| V -  Change the arrows scaling     |" << endl
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

VisualizationSceneVector  * vsvector;
extern VisualizationScene * locscene;
extern VisualizationSceneSolution * vssol;

static int ianim = 0;
static int ianimmax = 10;

void KeyDPressed()
{
   vsvector -> ToggleDisplacements();
   SendExposeEvent();
}

void KeyNPressed()
{
   ianim = (ianim + 1) % (ianimmax + 1);
   vsvector -> NPressed();
}

void KeyBPressed()
{
   ianim = (ianim + ianimmax) % (ianimmax + 1);
   vsvector -> NPressed();
}

void VisualizationSceneVector::NPressed()
{
   if (drawdisp == 0)
   {
      drawdisp = 1;
      ianim = 0;
   }

   PrepareDisplacedMesh();

   SendExposeEvent();
}

void KeyvPressed()
{
   vsvector -> ToggleVectorField();
   SendExposeEvent();
}

void KeyVPressed()
{
   cout << "New arrow scale: " << flush;
   cin >> vsvector -> ArrowScale;
   cout << "New arrow scale = " << vsvector -> ArrowScale << endl;
   vsvector -> PrepareVectorField();
   SendExposeEvent();
}

int key_u_func = 0;

void KeyuPressed()
{
   int update = 1;

   switch (key_u_func)
   {
   case 0:
      vsvector->RefineFactor++;
      break;

   case 1:
      if (vsvector->RefineFactor > 1)
         vsvector->RefineFactor--;
      else
         update = 0;
      break;

   case 2:
      vsvector->CycleVec2Scalar();
      SendExposeEvent();
      break;
   }

   switch (key_u_func)
   {
   case 0:
   case 1:
      if (update && vsvector->shading == 2)
      {
         vsvector->PrepareVectorField();
         SendExposeEvent();
      }
      cout << "Vector subdivision factor = "
           << vsvector->RefineFactor << endl;
      break;

   case 2:
      break;
   }
}

void KeyUPressed()
{
   key_u_func = (key_u_func+1)%3;
   cout << "Key 'u' will: ";
   switch(key_u_func)
   {
   case 0:
      cout << "Increase vector subdivision factor" << endl;
      break;

   case 1:
      cout << "Decrease vector subdivision factor" << endl;
      break;

   case 2:
      cout << "Cycle through vector-to-scalar functions" << endl;
      break;
   }
}

void VisualizationSceneVector::ToggleVectorField()
{
   drawvector = (drawvector+1)%4;
   PrepareVectorField();
}

VisualizationSceneVector::VisualizationSceneVector(Mesh & m,
                                                   Vector & sx, Vector & sy)
{
   mesh = &m;
   solx = &sx;
   soly = &sy;

   sol  = new Vector(mesh -> GetNV());

   VecGridF = NULL;

   Init();
}

VisualizationSceneVector::VisualizationSceneVector(GridFunction &vgf)
{
   FiniteElementSpace *fes = vgf.FESpace();
   if (fes == NULL || fes->GetVDim() != 2)
   {
      cout << "VisualizationSceneVector::VisualizationSceneVector" << endl;
      exit(1);
   }

   VecGridF = &vgf;

   mesh = fes->GetMesh();

   solx = new Vector(mesh -> GetNV());
   soly = new Vector(mesh -> GetNV());

   vgf.GetNodalValues (*solx, 1);
   vgf.GetNodalValues (*soly, 2);

   sol = new Vector(mesh -> GetNV());

   //  VisualizationSceneSolution::Init()  sets rsol = NULL !
   {
      Init();
      SetGridFunction(vgf);
   }

   mesh->GetNodes(vc0);
   if (vc0.Size() != vgf.Size())
      vc0.Destroy();
   else
      vc0 += vgf;
}

double VecLength(double x, double y)
{
   return sqrt(x*x+y*y);
}

double VecDirection(double x, double y)
{
   return atan2(y, x);
}

double VecDotNx(double x, double y)
{
   return x;
}

double VecDotNy(double x, double y)
{
   return y;
}

double VecDivSubst(double x, double y)
{
   return 0.0;
}

double VecCurlSubst(double x, double y)
{
   return 0.0;
}

double VecAnisotrSubst(double x, double y)
{
   return 0.0;
}

double (*Vec2ScalarFunctions[7])(double, double) =
{ VecLength, VecDirection, VecDotNx, VecDotNy, VecDivSubst, VecCurlSubst,
  VecAnisotrSubst };

void VisualizationSceneVector::CycleVec2Scalar()
{
   int i;

   for (i = 0; Vec2Scalar != Vec2ScalarFunctions[i]; i++)
      ;

   i = (i + 1) % 7;

   Vec2Scalar = Vec2ScalarFunctions[i];

   for (i = 0; i < mesh->GetNV(); i++)
      (*sol)(i) = Vec2Scalar((*solx)(i), (*soly)(i));

   // update scalar stuff
   FindNewBox();
   PrepareAxes();
   PrepareLines();
   DefaultLevelLines();
   PrepareLevelCurves();
   PrepareCP();
   Prepare();

   if (i == 0)
      maxlen = maxv;

   // update the colors of the vectors
   if (drawvector > 1)
      PrepareVectorField();
}

void VisualizationSceneVector::NewMeshAndSolution(GridFunction &vgf,
                                                  int rescale)
{
   delete sol;

   if (VecGridF)
   {
      delete soly;
      delete solx;
   }

   VecGridF = &vgf;

   mesh = vgf.FESpace()->GetMesh();

   solx = new Vector(mesh->GetNV());
   soly = new Vector(mesh->GetNV());

   vgf.GetNodalValues(*solx, 1);
   vgf.GetNodalValues(*soly, 2);

   mesh->GetNodes(vc0);
   if (vc0.Size() != vgf.Size())
      vc0.Destroy();
   else
      vc0 += vgf;

   sol = new Vector(mesh->GetNV());
   for (int i = 0; i < mesh->GetNV(); i++)
      (*sol)(i) = Vec2Scalar((*solx)(i), (*soly)(i));

   VisualizationSceneSolution::NewMeshAndSolution(mesh, sol, &vgf, rescale);

   if (rescale)
      if (Vec2Scalar == VecLength)
      {
         maxlen = maxv;
      }
      else
      {
         cout << "VisualizationSceneVector::NewMeshAndSolution() : "
            " maxlen not updated!" << endl;
      }

   PrepareVectorField();
}

void VisualizationSceneVector::Init()
{
   drawdisp = 0;
   drawvector = 0;
   ArrowScale = 1.0;
   RefineFactor = 1;
   Vec2Scalar = VecLength;

   for (int i = 0; i < mesh->GetNV(); i++)
      (*sol)(i) = Vec2Scalar((*solx)(i), (*soly)(i));

   vectorlist = glGenLists(1);
   displinelist = glGenLists(1);

   VisualizationSceneSolution::Init();

   PrepareVectorField();
   // PrepareDisplacedMesh(); // called by PrepareLines()

   vsvector = this;

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      auxKeyFunc (AUX_h, VectorKeyHPressed);
      auxKeyFunc (AUX_H, VectorKeyHPressed);

      auxKeyFunc (AUX_d, KeyDPressed);
      auxKeyFunc (AUX_D, KeyDPressed);
      auxKeyFunc (AUX_n, KeyNPressed);
      auxKeyFunc (AUX_N, KeyNPressed);
      auxKeyFunc (AUX_b, KeyBPressed);
      auxKeyFunc (AUX_B, KeyBPressed);
      auxKeyFunc (AUX_v, KeyvPressed);
      auxKeyFunc (AUX_V, KeyVPressed);
      auxKeyFunc (AUX_u, KeyuPressed);
      auxKeyFunc (AUX_U, KeyUPressed);
   }

   // Vec2Scalar is VecLength
   maxlen = maxv;
}

VisualizationSceneVector::~VisualizationSceneVector()
{
   glDeleteLists (displinelist, 1);
   glDeleteLists (vectorlist, 1);

   delete sol;

   if (VecGridF)
   {
      delete soly;
      delete solx;
   }
}

void VisualizationSceneVector::GetRefinedValues(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr)
{
   if (drawelems < 2)
   {
      DenseMatrix vec_vals;
      ElementTransformation *T;
      DenseMatrix gv;
      Vector ev;

      VecGridF->GetVectorValues(i, ir, vec_vals, tr);
      vals.SetSize(vec_vals.Width());
      for (int j = 0; j < vec_vals.Width(); j++)
      {
         if (Vec2Scalar == VecDivSubst || Vec2Scalar == VecCurlSubst ||
             Vec2Scalar == VecAnisotrSubst)
         {
            T = mesh->GetElementTransformation(i);
            T->SetIntPoint(&ir.IntPoint(j));
            if (Vec2Scalar == VecDivSubst)
            {
               vals(j) = VecGridF->GetDivergence(*T);
            }
            else
            {
               VecGridF->GetVectorGradient(*T, gv);
               if (Vec2Scalar == VecCurlSubst)
               {
                  vals(j) = gv(1, 0) - gv(0, 1); // curl v in 2D
               }
               else
               {
                  gv.Symmetrize();
                  gv.Eigenvalues(ev);
                  vals(j) = ev(1) - ev(0);
                  // vals(j) = ev(1);
                  // vals(j) = (ev(0) <= 0.0) ? ev(0) : ev(1);
                  // vals(j) = (fabs(ev(0)) >= fabs(ev(1))) ? ev(0) : ev(1);
               }
            }
         }
         else
            vals(j) = Vec2Scalar(vec_vals(0, j), vec_vals(1, j));
      }
   }
   else
   {
      ElementTransformation *T = mesh->GetElementTransformation(i);

      T->Transform(ir, tr);

      vals.SetSize(ir.GetNPoints());
      if (vc0.Size() == 0)
      {
         for (int j = 0; j < ir.GetNPoints(); j++)
         {
            T->SetIntPoint(&ir.IntPoint(j));
            vals(j) = T->Weight();
         }
      }
      else
      {
         mesh->GetElementTransformation(i, vc0, &T0);
         for (int j = 0; j < ir.GetNPoints(); j++)
         {
            T->SetIntPoint(&ir.IntPoint(j));
            T0.SetIntPoint(&ir.IntPoint(j));
            vals(j) = T->Weight() / T0.Weight();
         }
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

   if (shrink != 1.0)
      ShrinkPoints(tr, i, 0, 0);
}

int VisualizationSceneVector::GetRefinedValuesAndNormals(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr,
   DenseMatrix &normals)
{
   int have_normals = 0;

   GetRefinedValues(i, ir, vals, tr);

   return have_normals;
}

void VisualizationSceneVector::PrepareDisplacedMesh()
{
   int i, j, ne = mesh -> GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double zc = 0.5*(z[0]+z[1]);

   // prepare the displaced mesh
   glNewList(displinelist, GL_COMPILE);

   // Set_Material();
   // glColor3f (1, 0, 0);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   if (shading != 2)
   {
      for (i = 0; i < ne; i++)
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
            glVertex3d (pointmat.Elem(0, j)+
                        (*solx)(vertices[j])*(ianim)/ianimmax,
                        pointmat.Elem(1, j)+
                        (*soly)(vertices[j])*(ianim)/ianimmax,
                        zc);
         glEnd();
      }
   }
   else if (drawdisp < 2)
   {
      double sc = double(ianim)/ianimmax;
      DenseMatrix vvals, pm;

      for (i = 0; i < ne; i++)
      {
         RefinedGeometry *RefG =
            GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                        TimesToRefine, EdgeRefineFactor);
         VecGridF->GetVectorValues(i, RefG->RefPts, vvals, pm);

         Array<int> &RE = RefG->RefEdges;

         glBegin (GL_LINES);
         for (int k = 0; k+1 < RE.Size(); k++)
         {
            glVertex3d (pm(0, RE[k]) + sc * vvals(0, RE[k]),
                        pm(1, RE[k]) + sc * vvals(1, RE[k]), zc);
            k++;
            glVertex3d (pm(0, RE[k]) + sc * vvals(0, RE[k]),
                        pm(1, RE[k]) + sc * vvals(1, RE[k]), zc);
         }
         glEnd();
      }
   }
   else
   {
      Vector vals;
      DenseMatrix vvals, pm;

      double x_min, x_max, y_min, y_max;
      Array<double> levels_x, levels_y;

      x_min = y_min = numeric_limits<double>::infinity();
      x_max = y_max = -x_min;

      for (i = 0; i < ne; i++)
      {
         RefinedGeometry *RefG =
            GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                        TimesToRefine, EdgeRefineFactor);
         VecGridF->GetVectorValues(i, RefG->RefPts, vvals, pm);

         vvals += pm;

         for (int j = 0; j < vvals.Width(); j++)
         {
            double x = vvals(0, j), y = vvals(1, j);

            if (drawdisp == 3)
            {
               double phi, rho;
               phi = atan2(y, x);
               rho = sqrt(x * x + y * y);

               x = phi;
               y = rho;
            }

            if (x_min > x)  x_min = x;
            if (x_max < x)  x_max = x;
            if (y_min > y)  y_min = y;
            if (y_max < y)  y_max = y;
         }
      }

      // define levels_x, levels_y
      {
         int n = level.Size()-1;
         int nx, ny;
         if (drawdisp == 2)
         {
            if ((x_max - x_min) < (y_max - y_min))
            {
               nx = n;
               ny = (int)ceil(n * (y_max - y_min) / (x_max - x_min));
            }
            else
            {
               nx = (int)ceil(n * (x_max - x_min) / (y_max - y_min));
               ny = n;
            }
         }
         else
         {
            nx = ny = n;
         }
         levels_x.SetSize(nx+1);
         levels_y.SetSize(ny+1);
         for (i = 0; i <= nx; i++)
         {
            double a = double(i) / nx;
            levels_x[i] = (1. - a) * x_min + a * x_max;
         }
         double offs = 1e-3 / nx;
         levels_x[0]  = (1. - offs) * x_min + offs * x_max;
         levels_x[nx] = offs * x_min + (1. - offs) * x_max;
         for (i = 0; i <= ny; i++)
         {
            double a = double(i) / ny;
            levels_y[i] = (1. - a) * y_min + a * y_max;
         }
         offs = 1e-3 / ny;
         levels_y[0]  = (1. - offs) * y_min + offs * y_max;
         levels_y[ny] = offs * y_min + (1. - offs) * y_max;
      }

      for (i = 0; i < ne; i++)
      {
         RefinedGeometry *RefG =
            GlobGeometryRefiner.Refine (mesh->GetElementBaseGeometry(i),
                                        TimesToRefine, EdgeRefineFactor);
         VecGridF->GetVectorValues(i, RefG->RefPts, vvals, pm);

         {
            double sc = double(ianim)/ianimmax;
            // pm = pm + sc * vvals
            pm.Add(sc, vvals);
            // vvals = pm + (1-sc) * vvals
            Add(pm, vvals, 1.-sc, vvals);
         }

         if (drawdisp == 3)
         {
            for (int j = 0; j < vvals.Width(); j++)
            {
               double x = vvals(0, j), y = vvals(1, j);
               {
                  double phi, rho;
                  phi = atan2(y, x);
                  rho = sqrt(x * x + y * y);

                  vvals(0, j) = phi;
                  vvals(1, j) = rho;
               }
            }
         }

         Array<int> &RG = RefG->RefGeoms;
         int sides = mesh->GetElement(i)->GetNVertices();

         vals.SetSize(vvals.Width());
         for (int j = 0; j < vvals.Width(); j++)
            vals(j) = vvals(0, j);
         DrawLevelCurves(RG, pm, vals, sides, levels_x, 1);
         for (int j = 0; j < vvals.Width(); j++)
            vals(j) = vvals(1, j);
         DrawLevelCurves(RG, pm, vals, sides, levels_y, 1);
      }
   }

   glEndList();
}

double new_maxlen;

void VisualizationSceneVector::DrawVector(double px, double py, double vx,
                                          double vy, double cval)
{
   double zc = 0.5*(z[0]+z[1]);

   if (drawvector == 1)
   {
      arrow_type = 0;
      arrow_scaling_type = 0;
      Arrow(px, py, zc, vx, vy, 0.0, 1.0, 1./4./3.);
   }
   else if (drawvector > 0)
   {
      double area = (x[1]-x[0])*(y[1]-y[0]);
      double h = sqrt(area/mesh->GetNV()) * ArrowScale;

      arrow_type = 1;
      arrow_scaling_type = 1;

      MySetColor(cval, minv, maxv);

      if (drawvector == 2)
         Arrow(px, py, zc, vx, vy, 0.0, h, 0.125);
      else
      {
         double len = VecLength(vx, vy);
         Arrow(px, py, zc, vx, vy, 0.0,
               h*max(0.01, len/maxlen), 0.125);
         if (len > new_maxlen)
            new_maxlen = len;
      }
   }
}

void VisualizationSceneVector::PrepareVectorField()
{
   int rerun;
   do
   {
      rerun = 0;

      glNewList(vectorlist, GL_COMPILE);

      if (drawvector > 0)
      {
         int i;

         if (drawvector == 3)
            new_maxlen = 0.0;
         for (i = 0; i < mesh->GetNV(); i++)
         {
            double *v = mesh->GetVertex(i);
            DrawVector(v[0], v[1], (*solx)(i), (*soly)(i), (*sol)(i));
         }

         if (shading == 2 && RefineFactor > 1)
         {
            DenseMatrix vvals, pm;
            for (i = 0; i < mesh->GetNE(); i++)
            {
               const IntegrationRule *ir =
                  GlobGeometryRefiner.RefineInterior(
                     mesh->GetElementBaseGeometry(i), RefineFactor);
               if (ir == NULL)
                  continue;
               VecGridF->GetVectorValues(i, *ir, vvals, pm);
               for (int j = 0; j < vvals.Width(); j++)
               {
                  DrawVector(pm(0, j), pm(1,j), vvals(0, j), vvals(1, j),
                             Vec2Scalar(vvals(0, j), vvals(1, j)));
               }
            }
            for (i = 0; i < mesh->GetNEdges(); i++)
            {
               const IntegrationRule *ir =
                  GlobGeometryRefiner.RefineInterior(
                     mesh->GetFaceBaseGeometry(i), RefineFactor);
               if (ir == NULL)
                  continue;
               VecGridF->GetFaceVectorValues(i, 0, *ir, vvals, pm);
               for (int j = 0; j < vvals.Width(); j++)
               {
                  DrawVector(pm(0, j), pm(1,j), vvals(0, j), vvals(1, j),
                             Vec2Scalar(vvals(0, j), vvals(1, j)));
               }
            }
         }

         if (drawvector == 3 && new_maxlen != maxlen)
         {
            maxlen = new_maxlen;
            rerun = 1;
         }
      }

      glEndList();

   }
   while(rerun);
}

void VisualizationSceneVector::Draw()
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

   // draw colored faces
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

   // draw lines
   if ( drawmesh == 1 )
      glCallList(linelist);
   else if (drawmesh == 2 )
      glCallList(lcurvelist);

   if (drawvector == 1)
      glCallList(vectorlist);

   if (drawdisp > 0)
   {
      if (drawmesh == 1)
         glColor3d(1., 0., 0.);
      glCallList(displinelist);
   }

   if (light)
      glEnable(GL_LIGHTING);

   if (GetUseTexture())
   {
      glEnable (GL_TEXTURE_1D);
      glColor4d(1, 1, 1, 1);
   }

   // draw vector field
   if (drawvector > 1)
      glCallList(vectorlist);

   if (MatAlpha < 1.0)
      Set_Transparency();

   // draw elements
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
   if (drawelems)
      glCallList(displlist);

   if (MatAlpha < 1.0)
      Remove_Transparency();

   if (GetUseTexture())
      glDisable (GL_TEXTURE_1D);

   glFlush();
   glXSwapBuffers (auxXDisplay(), auxXWindow());
}
