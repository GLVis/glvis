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

#include <cstdlib>
#include <iostream>
#include <limits>
#include <cmath>

#include <mfem.hpp>
#include "vsvector.hpp"

using namespace mfem;
using namespace std;

std::string VisualizationSceneVector::GetHelpString() const
{
   std::stringstream os;
   os << endl
      << "+------------------------------------+" << endl
      << "| Keys                               |" << endl
      << "+------------------------------------+" << endl
      << "| a -  Displays/Hides the axes       |" << endl
      << "| A -  Turns antialiasing on/off     |" << endl
      << "| b -  Displacements step back       |" << endl
      << "| B -  Toggle 2D boundary            |" << endl
      << "| c -  Toggle colorbar and caption   |" << endl
      << "| C -  Change the main plot caption  |" << endl
      << "| d -  Displays/Hides displacements  |" << endl
      << "| e -  Displays/Hides the elements   |" << endl
      << "| f -  Smooth/Flat shading           |" << endl
      << "| g -  Toggle background             |" << endl
      << "| h -  Displays help menu            |" << endl
      << "| i -  Toggle the cutting plane      |" << endl
      << "| j -  Turn on/off perspective       |" << endl
      << "| k/K  Adjust the transparency level |" << endl
      << "| ,/<  Adjust color transparency     |" << endl
      << "| l -  Turns on/off the light        |" << endl
      << "| L -  Toggle logarithmic scale      |" << endl
      << "| m -  Displays/Hides the mesh       |" << endl
      << "| n -  Displacements step forward    |" << endl
      << "| N -  Cycle through numberings      |" << endl
      << "| o -  (De)refine elem. (NC shading) |" << endl
      << "| O -  Switch 'o' func. (NC shading) |" << endl
      << "| p/P  Cycle through color palettes  |" << endl
      << "| q -  Quits                         |" << endl
      << "| Q -  Cycle quadrature data mode    |" << endl
      << "| r -  Reset the plot to 3D view     |" << endl
      << "| R -  Reset the plot to 2D view     |" << endl
      << "| s -  Turn on/off unit cube scaling |" << endl
      << "| S -  Take snapshot/Record a movie  |" << endl
      << "| t -  Cycle materials and lights    |" << endl
      << "| u -  Vector sampling; scalar func. |" << endl
      << "| U -  Switch 'u' functionality      |" << endl
      << "| v -  Cycle through vector fields   |" << endl
      << "| V -  Change the arrows scaling     |" << endl
      << "| \\ -  Set light source position     |" << endl
      << "| Alt+a  - Axes number format        |" << endl
      << "| Alt+c  - Colorbar number format    |" << endl
      << "| Ctrl+p - Print to a PDF file       |" << endl
      << "+------------------------------------+" << endl
      << "| Function keys                      |" << endl
      << "+------------------------------------+" << endl
      << "| F1 - X window info and keystrokes  |" << endl
      << "| F2 - Update colors, etc.           |" << endl
      << "| F3/F4 - Shrink/Zoom elements       |" << endl
      << "| F5 - Set level lines               |" << endl
      << "| F6 - Palette options               |" << endl
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
   return os.str();
}

thread_local VisualizationSceneVector  * vsvector;
extern thread_local VisualizationScene * locscene;
extern thread_local VisualizationSceneSolution * vssol;
extern thread_local GeometryRefiner GLVisGeometryRefiner;

thread_local int ianim = 0;
thread_local int ianimmax = 10;

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
         {
            vsvector->RefineFactor--;
         }
         else
         {
            update = 0;
         }
         break;

      case 2:
         vsvector->CycleVec2Scalar(1);
         SendExposeEvent();
         break;
   }

   switch (key_u_func)
   {
      case 0:
      case 1:
         if (update &&
             vsvector->GetShading() == VisualizationSceneSolution::Shading::Noncomforming)
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
   switch (key_u_func)
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

void VisualizationSceneVector::ToggleDrawElems()
{
   const char *modes[] =
   { "none", "vector->scalar function", "det(J0)/det(J)", "det(J)/det(J0)" };

   drawelems = (drawelems+3)%4;

   cout << "Surface elements mode : " << modes[drawelems] << endl;

   if (drawelems != 0 && shading == Shading::Noncomforming)
   {
      DoAutoscaleValue(false);
      PrepareLines();
      PrepareBoundary();
      Prepare();
      PrepareLevelCurves();
      PrepareCP();
   }
}

void VisualizationSceneVector::ToggleVectorField()
{
   drawvector = (drawvector+1)%6;
   PrepareVectorField();
}

const char *Vec2ScalarNames[7] =
{
   "magnitude", "direction", "x-component", "y-component", "divergence",
   "curl", "anisotropy"
};

VisualizationSceneVector::VisualizationSceneVector(Mesh & m,
                                                   Vector & sx, Vector & sy, Mesh *mc)
{
   mesh = &m;
   mesh_coarse = mc;
   solx = &sx;
   soly = &sy;

   sol  = new Vector(mesh -> GetNV());

   VecGridF = NULL;

   Init();
}

VisualizationSceneVector::VisualizationSceneVector(GridFunction &vgf)
{
   FiniteElementSpace *fes = vgf.FESpace();
   if (fes == NULL || vgf.VectorDim() != 2)
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
   {
      vc0.Destroy();
   }
   else
   {
      vc0 += vgf;
   }
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
{
   VecLength, VecDirection, VecDotNx, VecDotNy, VecDivSubst, VecCurlSubst,
   VecAnisotrSubst
};

void VisualizationSceneVector::CycleVec2Scalar(int print)
{
   int i;

   for (i = 0; Vec2Scalar != Vec2ScalarFunctions[i]; i++)
      ;

   if (VecGridF->FESpace()->GetVDim() == 1)
   {
      if (dynamic_cast<const ND_FECollection*>(VecGridF->FESpace()->FEColl()))
      {
         if (i == 5) { i = 0; }
         else if (i == 3) { i = 5; }
         else { i = (i + 1) % 5; }
      }
      else
      {
         i = (i + 1) % 5;
      }
   }
   else
   {
      i = (i + 1) % 7;
   }

   if (print)
   {
      cout << "Vector-to-scalar function: " << Vec2ScalarNames[i] << endl;
   }

   Vec2Scalar = Vec2ScalarFunctions[i];
   extra_caption = Vec2ScalarNames[i];

   for (i = 0; i < mesh->GetNV(); i++)
   {
      (*sol)(i) = Vec2Scalar((*solx)(i), (*soly)(i));
   }

   // update scalar stuff
   DoAutoscaleValue(false);
   PrepareLines();
   PrepareBoundary();
   PrepareLevelCurves();
   PrepareCP();
   Prepare();

   if (i == 0)
   {
      maxlen = maxv;
   }

   // update the colors of the vectors
   if (drawvector > 1)
   {
      PrepareVectorField();
   }
}

void VisualizationSceneVector::NewMeshAndSolution(GridFunction &vgf, Mesh *mc)
{
   delete sol;

   if (VecGridF)
   {
      delete soly;
      delete solx;
   }

   Mesh *old_m = mesh;
   Mesh *new_mesh = vgf.FESpace()->GetMesh();
   mesh = new_mesh;
   mesh_coarse = mc;
   VecGridF = &vgf;

   // If the number of elements changes, recompute the refinement factor
   if (mesh->GetNE() != old_m->GetNE())
   {
      int ref = GetAutoRefineFactor();
      if (TimesToRefine != ref || EdgeRefineFactor != 1)
      {
         TimesToRefine = ref;
         EdgeRefineFactor = 1;
         cout << "Subdivision factors = " << TimesToRefine << ", 1" << endl;
      }
      if (RefineFactor != 1)
      {
         RefineFactor = 1;
         cout << "Vector subdivision factor = 1" << endl;
      }
   }

   solx = new Vector(mesh->GetNV());
   soly = new Vector(mesh->GetNV());

   vgf.GetNodalValues(*solx, 1);
   vgf.GetNodalValues(*soly, 2);

   mesh->GetNodes(vc0);
   if (vc0.Size() != vgf.Size())
   {
      vc0.Destroy();
   }
   else
   {
      vc0 += vgf;
   }

   sol = new Vector(mesh->GetNV());
   for (int i = 0; i < mesh->GetNV(); i++)
   {
      (*sol)(i) = Vec2Scalar((*solx)(i), (*soly)(i));
   }

   VisualizationSceneSolution::NewMeshAndSolution(mesh, mesh_coarse, sol, &vgf);

   if (autoscale)
   {
      if (Vec2Scalar == VecLength)
      {
         maxlen = maxv;
      }
      else
      {
         cout << "VisualizationSceneVector::NewMeshAndSolution() : "
              " maxlen not updated!" << endl;
      }
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
   extra_caption = Vec2ScalarNames[0];

   for (int i = 0; i < mesh->GetNV(); i++)
   {
      (*sol)(i) = Vec2Scalar((*solx)(i), (*soly)(i));
   }

   VisualizationSceneSolution::Init();

   PrepareVectorField();
   // PrepareDisplacedMesh(); // called by PrepareLines()

   vsvector = this;

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      wnd->setOnKeyDown('d', KeyDPressed);
      wnd->setOnKeyDown('D', KeyDPressed);
      wnd->setOnKeyDown('n', KeyNPressed);
      wnd->setOnKeyDown('b', KeyBPressed);
      wnd->setOnKeyDown('v', KeyvPressed);
      wnd->setOnKeyDown('V', KeyVPressed);
      wnd->setOnKeyDown('u', KeyuPressed);
      wnd->setOnKeyDown('U', KeyUPressed);
   }

   // Vec2Scalar is VecLength
   maxlen = maxv;
}

VisualizationSceneVector::~VisualizationSceneVector()
{
   delete sol;

   if (VecGridF)
   {
      delete soly;
      delete solx;
   }
}

void VisualizationSceneVector::GetRefinedValues(const int i,
                                                const IntegrationRule &ir,
                                                Vector &vals, DenseMatrix &tr,
                                                const bool do_shrink)
{
   if (drawelems < 2)
   {
      DenseMatrix vec_vals;
      ElementTransformation *T;
      DenseMatrix gv;
      double ev[2], evec[4];
      int vdim = VecGridF->FESpace()->GetVDim();
      double curl_v[1];
      Vector curl(curl_v, 1);

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
            else if (Vec2Scalar == VecCurlSubst)
            {
               VecGridF->GetCurl(*T, curl);
               vals(j) = curl(0);
            }
            else
            {
               if (vdim == 1)
               {
                  vals(j) = 0.0;
                  continue;
               }
               VecGridF->GetVectorGradient(*T, gv);
               if (Vec2Scalar == VecCurlSubst)
               {
                  vals(j) = gv(1, 0) - gv(0, 1); // curl v in 2D
               }
               else
               {
                  gv.Symmetrize();
                  gv.CalcEigenvalues(ev, evec);
                  vals(j) = ev[1] - ev[0];
                  // vals(j) = ev[1];
                  // vals(j) = (ev[0] <= 0.0) ? ev[0] : ev[1];
                  // vals(j) = (fabs(ev[0]) >= fabs(ev[1])) ? ev[0] : ev[1];
               }
            }
         }
         else
         {
            vals(j) = Vec2Scalar(vec_vals(0, j), vec_vals(1, j));
         }
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

   if (logscale)
      for (int j = 0; j < vals.Size(); j++)
      {
         vals(j) = _LogVal(vals(j));
      }

   if (do_shrink && (shrink != 1.0 || shrinkmat != 1.0))
   {
      ShrinkPoints(tr, i, 0, 0);
   }
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
   int ne = mesh -> GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double zc = 0.5*(bb.z[0]+bb.z[1]);

   // prepare the displaced mesh
   displine_buf.clear();
   gl3::GlBuilder build = displine_buf.createBuilder();
   if (shading != Shading::Noncomforming)
   {
      for (int i = 0; i < ne; i++)
      {
         build.glBegin(GL_LINE_LOOP);
         mesh->GetPointMatrix (i, pointmat);
         mesh->GetElementVertices (i, vertices);

         for (int j = 0; j < pointmat.Size(); j++)
         {
            build.glVertex3d (
               pointmat.Elem(0, j)+
               (*solx)(vertices[j])*(ianim)/ianimmax,
               pointmat.Elem(1, j)+
               (*soly)(vertices[j])*(ianim)/ianimmax,
               zc);
         }
         build.glEnd();
      }
   }
   else if (drawdisp < 2)
   {
      double sc = double(ianim)/ianimmax;
      DenseMatrix vvals, pm;

      for (int i = 0; i < ne; i++)
      {
         RefinedGeometry *RefG =
            GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                        TimesToRefine, EdgeRefineFactor);
         VecGridF->GetVectorValues(i, RefG->RefPts, vvals, pm);

         Array<int> &RE = RefG->RefEdges;
         build.glBegin (GL_LINES);
         for (int k = 0; k+1 < RE.Size(); k++)
         {
            build.glVertex3d (
               pm(0, RE[k]) + sc * vvals(0, RE[k]),
               pm(1, RE[k]) + sc * vvals(1, RE[k]), zc);
            k++;
            build.glVertex3d (
               pm(0, RE[k]) + sc * vvals(0, RE[k]),
               pm(1, RE[k]) + sc * vvals(1, RE[k]), zc);
         }
         build.glEnd();
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

      for (int i = 0; i < ne; i++)
      {
         RefinedGeometry *RefG =
            GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
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

            if (x_min > x) { x_min = x; }
            if (x_max < x) { x_max = x; }
            if (y_min > y) { y_min = y; }
            if (y_max < y) { y_max = y; }
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
         for (int i = 0; i <= nx; i++)
         {
            double a = double(i) / nx;
            levels_x[i] = (1. - a) * x_min + a * x_max;
         }
         double offs = 1e-3 / nx;
         levels_x[0]  = (1. - offs) * x_min + offs * x_max;
         levels_x[nx] = offs * x_min + (1. - offs) * x_max;
         for (int i = 0; i <= ny; i++)
         {
            double a = double(i) / ny;
            levels_y[i] = (1. - a) * y_min + a * y_max;
         }
         offs = 1e-3 / ny;
         levels_y[0]  = (1. - offs) * y_min + offs * y_max;
         levels_y[ny] = offs * y_min + (1. - offs) * y_max;
      }

      for (int i = 0; i < ne; i++)
      {
         RefinedGeometry *RefG =
            GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
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
         {
            vals(j) = vvals(0, j);
         }
         DrawLevelCurves(build, RG, pm, vals, sides, levels_x, 1);
         for (int j = 0; j < vvals.Width(); j++)
         {
            vals(j) = vvals(1, j);
         }
         DrawLevelCurves(build, RG, pm, vals, sides, levels_y, 1);
      }
   }

   updated_bufs.emplace_back(&displine_buf);
}

thread_local double new_maxlen;

void VisualizationSceneVector::DrawVector(double px, double py, double vx,
                                          double vy, double cval)
{
   double zc = (drawvector > 3)?(bb.z[1]):(0.5*(bb.z[0]+bb.z[1]));

   if (drawvector == 1)
   {
      arrow_type = 0;
      arrow_scaling_type = 0;
      Arrow(vector_buf, px, py, zc, vx, vy, 0.0, 1.0, 1./4./3.);
   }
   else if (drawvector > 0)
   {
      double area = (bb.x[1]-bb.x[0])*(bb.y[1]-bb.y[0]);
      double h = sqrt(area/mesh->GetNV()) * ArrowScale;

      arrow_type = 1;
      arrow_scaling_type = 1;

      if (drawvector > 3) { cval = HUGE_VAL; }

      if (drawvector == 2 || drawvector == 4)
      {
         Arrow(vector_buf, px, py, zc, vx, vy, 0.0, h, 0.125, cval);
      }
      else
      {
         double len = VecLength(vx, vy);
         Arrow(vector_buf, px, py, zc, vx, vy, 0.0,
               h*max(0.01, len/maxlen), 0.125, cval);
         if (len > new_maxlen)
         {
            new_maxlen = len;
         }
      }
   }
}

int VisualizationSceneVector::GetFunctionAutoRefineFactor()
{
   if (!VecGridF) { return 1;}

   return VisualizationSceneScalarData::GetFunctionAutoRefineFactor(*VecGridF);
}

void VisualizationSceneVector::PrepareVectorField()
{
   int rerun;
   do
   {
      rerun = 0;

      // glNewList(vectorlist, GL_COMPILE);
      vector_buf.clear();

      if (drawvector > 0)
      {
         int i;

         palette.SetUseLogscale(logscale);
         if (drawvector == 3 || drawvector == 5)
         {
            new_maxlen = 0.0;
         }
         for (i = 0; i < mesh->GetNV(); i++)
         {
            double *v = mesh->GetVertex(i);
            DrawVector(v[0], v[1], (*solx)(i), (*soly)(i), (*sol)(i));
         }

         if (shading == Shading::Noncomforming && RefineFactor > 1)
         {
            DenseMatrix vvals, pm;
            for (i = 0; i < mesh->GetNE(); i++)
            {
               const IntegrationRule *ir =
                  GLVisGeometryRefiner.RefineInterior(
                     mesh->GetElementBaseGeometry(i), RefineFactor);
               if (ir == NULL)
               {
                  continue;
               }
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
                  GLVisGeometryRefiner.RefineInterior(
                     mesh->GetFaceGeometry(i), RefineFactor);
               if (ir == NULL)
               {
                  continue;
               }
               VecGridF->GetFaceVectorValues(i, 0, *ir, vvals, pm);
               for (int j = 0; j < vvals.Width(); j++)
               {
                  DrawVector(pm(0, j), pm(1,j), vvals(0, j), vvals(1, j),
                             Vec2Scalar(vvals(0, j), vvals(1, j)));
               }
            }
         }

         if ((drawvector == 3 || drawvector == 5) && new_maxlen != maxlen)
         {
            maxlen = new_maxlen;
            rerun = 1;
         }
      }

      // glEndList();

   }
   while (rerun);
   updated_bufs.emplace_back(&vector_buf);
}

gl3::SceneInfo VisualizationSceneVector::GetSceneObjs()
{
   if (colorbar)
   {
      // update color bar before we get the base class scene
      PrepareColorBar(minv, maxv, (drawmesh == 2) ? &level : nullptr );
   }
   gl3::SceneInfo scene = VisualizationSceneScalarData::GetSceneObjs();
   gl3::RenderParams params = GetMeshDrawParams();
   params.use_clip_plane = draw_cp;
   double* cp_eqn = CuttingPlane->Equation();
   params.clip_plane_eqn = {cp_eqn[0], cp_eqn[1], cp_eqn[2], cp_eqn[3]};
   params.contains_translucent = false;
   if (drawvector == 2 || drawvector == 3)
   {
      scene.queue.emplace_back(params, &vector_buf);
   }
   params.contains_translucent = matAlpha < 1.0 ||
                                 palette.GetPalette()->IsTranslucent();
   if (drawelems)
   {
      scene.queue.emplace_back(params, &disp_buf);
   }
   params.contains_translucent = false;
   params.mesh_material = BLK_MAT;
   params.static_color = GetLineColor();
   if (draw_cp)
   {
      params.use_clip_plane = false;
      scene.queue.emplace_back(params, &cp_buf);
      params.use_clip_plane = true;
   }
   // disable lighting for objects below
   params.num_pt_lights = 0;

   // draw boundary in 2D
   if (drawbdr)
   {
      scene.queue.emplace_back(params, &bdr_buf);
   }

   // draw lines
   if (drawmesh == 1)
   {
      scene.queue.emplace_back(params, &line_buf);
   }
   else if (drawmesh == 2)
   {
      scene.queue.emplace_back(params, &lcurve_buf);
   }

   // draw numberings
   if (drawnums == Numbering::ELEMENTS)
   {
      scene.queue.emplace_back(params, &e_nums_buf);
   }
   else if (drawnums == Numbering::VERTICES)
   {
      scene.queue.emplace_back(params, &v_nums_buf);
   }
   else if (drawnums == Numbering::EDGES)
   {
      scene.queue.emplace_back(params, &f_nums_buf);
   }
   else if (drawnums == Numbering::DOFS)
   {
      scene.queue.emplace_back(params, &d_nums_buf);
   }
   else { /* Numbering::NONE */ }

   if (drawvector == 1 || drawvector > 3)
   {
      if (drawvector > 3)
      {
         params.static_color = {.3f, .3f, .3f, 1.f};
      }
      scene.queue.emplace_back(params, &vector_buf);
   }

   if (drawdisp > 0)
   {
      if (drawmesh == 1)
      {
         params.static_color = {1.f, 0.f, 0.f, 1.f};
      }
      else
      {
         params.static_color = GetLineColor();
      }
      scene.queue.emplace_back(params, &displine_buf);
   }
   ProcessUpdatedBufs(scene);

   return scene;
}

void VisualizationSceneVector::glTF_Export()
{
   string name = "GLVis_scene_000";

   glTF_Builder bld(name);

   auto palette_mat = AddPaletteMaterial(bld);
   auto pal_lines_mat = AddPaletteLinesMaterial(bld, palette_mat);
   auto black_mat = AddBlackMaterial(bld);
   auto buf = bld.addBuffer("buffer");
   if (drawvector)
   {
      auto vec_node = AddModelNode(bld, "Vectors");
      auto vec_mesh = bld.addMesh("Vectors Mesh");
      bld.addNodeMesh(vec_node, vec_mesh);

      int ntria = AddTriangles(
                     bld,
                     vec_mesh,
                     buf,
                     (drawvector == 1 || drawvector > 3) ? black_mat : palette_mat,
                     vector_buf);
      int nlines = AddLines(
                      bld,
                      vec_mesh,
                      buf,
                      (drawvector == 1 || drawvector > 3) ? black_mat : pal_lines_mat,
                      vector_buf);
      if (ntria == 0 || nlines == 0)
      {
         cout << "glTF export: no vectors found to export!" << endl;
      }
   }
   if (drawelems) { glTF_ExportElements(bld, buf, palette_mat, disp_buf); }
   if (drawmesh)
   {
      glTF_ExportMesh(bld, buf, black_mat,
                      (drawmesh == 1) ? line_buf : lcurve_buf);
   }
   if (drawbdr) { glTF_ExportBoundary(bld, buf, black_mat); }
   if (drawaxes) { glTF_ExportBox(bld, buf, black_mat); }
   bld.writeFile();

   cout << "Exported glTF -> " << name << ".gltf" << endl;
}
