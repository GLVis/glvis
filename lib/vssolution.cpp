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

#include <cstdlib>
#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <vector>

#include "mfem.hpp"
using namespace mfem;
#include "visual.hpp"

using namespace std;


VisualizationSceneSolution *vssol;
extern VisualizationScene  *locscene;
extern GeometryRefiner GLVisGeometryRefiner;

#ifdef GLVIS_ISFINITE
/* This test for INFs or NaNs is the same as the one used in hypre's PCG and
   should work on all IEEE-compliant compilers. Fro detail see "Lecture Notes on
   the Status of IEEE 754" by W. Kahan, http://tinyurl.com/cfz5d88 */
int isfinite(double x)
{
   double ieee_check = 1;

   if (x != 0)
   {
      ieee_check = x/x;   // INF -> NaN conversion
   }

   return (ieee_check == ieee_check);
}
#endif

// Definitions of some more keys

static void SolutionKeyHPressed()
{
   cout << endl
        << "+------------------------------------+" << endl
        << "| Keys                               |" << endl
        << "+------------------------------------+" << endl
        << "| a -  Displays/Hides the axes       |" << endl
        << "| A -  Turns antialiasing on/off     |" << endl
        << "| b/B  Displays/Hides the boundary   |" << endl
        << "| c -  Toggle colorbar and caption   |" << endl
        << "| C -  Change the main plot caption  |" << endl
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
        << "| n/N  Cycle through numberings      |" << endl
        << "| p/P  Cycle through color palettes  |" << endl
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
}

static void KeyF8Pressed()
{
   int attr;
   Array<int> attr_list(&attr, 1);
   const Array<int> &all_attr = vssol->GetMesh()->attributes;

   cout << "El attributes ON: ";
   for (int i = 0; i < all_attr.Size(); i++)
      if (vssol->el_attr_to_show[all_attr[i]-1])
      {
         cout << " " << all_attr[i];
      }
   cout << endl;

   cout << "El attribute to toggle : " << flush;
   cin >> attr;
   vssol->ToggleAttributes(attr_list);
   SendExposeEvent();
}

static void SwitchAttribute(int increment, int &attribute,
                            Array<int> &attribute_marker,
                            bool bdr)
{
   const char *attr_type = bdr ? "bdr" : "element";
   if (attribute_marker.Size() == 0)
   {
      cout << "There are no " << attr_type << " attributes" << endl;
      return;
   }
   if (attribute == -1)
   {
      attribute_marker = 0;
      attribute = (increment >= 0) ? 0 : attribute_marker.Size()-1;
   }
   else
   {
      if (attribute != attribute_marker.Size())
      {
         attribute_marker[attribute] = 0;
      }
      else
      {
         attribute_marker = 0;
      }
      attribute += increment;
   }
   attribute += attribute_marker.Size()+1;
   attribute %= attribute_marker.Size()+1;
   if (attribute != attribute_marker.Size())
   {
      attribute_marker[attribute] = 1;
      cout << "Showing " << attr_type << " attribute " << attribute+1 << endl;
   }
   else
   {
      attribute_marker = 1;
      cout << "Showing all " << attr_type << " attributes" << endl;
   }
   if (bdr)
   {
      vssol->PrepareBoundary();
   }
   else
   {
      vssol->PrepareLines();
      vssol->Prepare();
   }
   SendExposeEvent();
}

static void KeyF9Pressed(GLenum state)
{
   if (!(state & KMOD_SHIFT))
   {
      SwitchAttribute(+1, vssol->attr_to_show, vssol->el_attr_to_show, false);
   }
   else
   {
      SwitchAttribute(+1, vssol->bdr_attr_to_show, vssol->bdr_el_attr_to_show,
                      true);
   }
}

static void KeyF10Pressed(GLenum state)
{
   if (!(state & KMOD_SHIFT))
   {
      SwitchAttribute(-1, vssol->attr_to_show, vssol->el_attr_to_show, false);
   }
   else
   {
      SwitchAttribute(-1, vssol->bdr_attr_to_show, vssol->bdr_el_attr_to_show,
                      true);
   }
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

static void KeyNPressed()
{
   vssol -> ToggleDrawNumberings();
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
         {
            vssol -> TimesToRefine -= vssol -> EdgeRefineFactor;
         }
         else
         {
            update = 0;
         }
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
         {
            update = 0;
         }
         break;
   }
   if (update && vssol -> shading == 2)
   {
      vssol -> DoAutoscale(false);
      vssol -> PrepareLines();
      vssol -> PrepareBoundary();
      vssol -> Prepare();
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
      vssol->PrepareLevelCurves();
      vssol->PrepareNumbering();
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
      vssol->PrepareLevelCurves();
      vssol->PrepareNumbering();
      SendExposeEvent();
   }
}

static void KeyF11Pressed()
{
   if (vssol->shading == 2)
   {
      if (vssol->matc.Width() == 0)
      {
         vssol->ComputeElemAttrCenter();
      }
      vssol->shrinkmat *= 0.9;
      vssol->Prepare();
      vssol->PrepareLines();
      vssol->PrepareBoundary();
      vssol->PrepareLevelCurves();
      vssol->PrepareNumbering();
      SendExposeEvent();
   }
}

static void KeyF12Pressed()
{
   if (vssol->shading == 2)
   {
      if (vssol->matc.Width() == 0)
      {
         vssol->ComputeElemAttrCenter();
      }
      vssol->shrinkmat *= 1.11111111111111111111111;
      vssol->Prepare();
      vssol->PrepareLines();
      vssol->PrepareBoundary();
      vssol->PrepareLevelCurves();
      vssol->PrepareNumbering();
      SendExposeEvent();
   }
}

VisualizationSceneSolution::VisualizationSceneSolution()
{
   v_normals = NULL;
}

VisualizationSceneSolution::VisualizationSceneSolution(
   Mesh &m, Vector &s, Vector *normals)
{
   mesh = &m;
   sol = &s;
   v_normals = normals;

   Init();

   wnd->setOnKeyDown('h', SolutionKeyHPressed);
   wnd->setOnKeyDown('H', SolutionKeyHPressed);
}

void VisualizationSceneSolution::Init()
{
   rsol  = NULL;
   vssol = this;

   drawelems = shading = 1;
   drawmesh  = 0;
   drawnums  = 0;

   shrink = 1.0;
   shrinkmat = 1.0;
   bdrc.SetSize(2,0);
   matc.SetSize(2,0);

   TimesToRefine = EdgeRefineFactor = 1;

   attr_to_show = bdr_attr_to_show = -1;
   el_attr_to_show.SetSize(mesh->attributes.Max());
   el_attr_to_show = 1;
   bdr_el_attr_to_show.SetSize(mesh->bdr_attributes.Size() > 0 ?
                               mesh->bdr_attributes.Max() : 0);
   bdr_el_attr_to_show = 1;

   drawbdr = 0;

   VisualizationSceneScalarData::Init();  // Calls FindNewBox() !!!

   SetUseTexture(1);

   double eps = 1e-6; // move the cutting plane a bit to avoid artifacts
   CuttingPlane = new Plane(-1.0,0.0,0.0,(0.5-eps)*x[0]+(0.5+eps)*x[1]);
   draw_cp = 0;

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      wnd->setOnKeyDown('b', KeyBPressed);
      wnd->setOnKeyDown('B', KeyBPressed);

      wnd->setOnKeyDown('m', KeyMPressed);
      wnd->setOnKeyDown('M', KeyMPressed);

      wnd->setOnKeyDown('n', KeyNPressed);
      wnd->setOnKeyDown('N', KeyNPressed);

      wnd->setOnKeyDown('e', KeyEPressed);
      wnd->setOnKeyDown('E', KeyEPressed);

      wnd->setOnKeyDown('f', KeyFPressed);
      wnd->setOnKeyDown('F', KeyFPressed);

      wnd->setOnKeyDown('i', KeyiPressed);
      wnd->setOnKeyDown('I', KeyIPressed);

      wnd->setOnKeyDown('w', KeywPressed);
      wnd->setOnKeyDown('y', KeyyPressed);
      wnd->setOnKeyDown('Y', KeyYPressed);
      wnd->setOnKeyDown('z', KeyzPressed);
      wnd->setOnKeyDown('Z', KeyZPressed);

      wnd->setOnKeyDown(SDLK_F3, KeyF3Pressed);
      wnd->setOnKeyDown(SDLK_F4, KeyF4Pressed);
      wnd->setOnKeyDown(SDLK_F8, KeyF8Pressed);
      wnd->setOnKeyDown(SDLK_F9,  KeyF9Pressed);
      wnd->setOnKeyDown(SDLK_F10, KeyF10Pressed);
      wnd->setOnKeyDown(SDLK_F11, KeyF11Pressed);
      wnd->setOnKeyDown(SDLK_F12, KeyF12Pressed);
   }

   Prepare();
   PrepareLines();
   PrepareLevelCurves();
   PrepareBoundary();
   PrepareNumbering();
}

VisualizationSceneSolution::~VisualizationSceneSolution()
{
}

void VisualizationSceneSolution::ToggleDrawElems()
{
   const char *modes[] =
   {
      "none", "solution", "kappa + 1/kappa", "kappa", "1/det(J)", "det(J)",
      "attribute"
   };

   drawelems = (drawelems + 6) % 7;

   cout << "Surface elements mode : " << modes[drawelems] << endl;
   if (drawelems < 2)
   {
      extra_caption.clear();
   }
   else
   {
      extra_caption = modes[drawelems];
   }

   if (drawelems != 0 && shading == 2)
   {
      DoAutoscaleValue(false);
      PrepareLines();
      PrepareBoundary();
      Prepare();
      PrepareLevelCurves();
      PrepareCP();
      PrepareNumbering();
   }
}

void VisualizationSceneSolution::NewMeshAndSolution(
   Mesh *new_m, Vector *new_sol, GridFunction *new_u)
{
   // If the number of elements changes, recompute the refinement factor
   if (mesh->GetNE() != new_m->GetNE())
   {
      mesh = new_m;
      int ref = GetAutoRefineFactor();
      if (TimesToRefine != ref || EdgeRefineFactor != 1)
      {
         TimesToRefine = ref;
         EdgeRefineFactor = 1;
         cout << "Subdivision factors = " << TimesToRefine << ", 1" << endl;
      }
   }
   mesh = new_m;
   sol = new_sol;
   rsol = new_u;

   DoAutoscale(false);

   Prepare();
   PrepareLines();
   PrepareLevelCurves();
   PrepareBoundary();
   PrepareCP();
}


void VisualizationSceneSolution::GetRefinedDetJ(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr)
{
   int geom = mesh->GetElementBaseGeometry(i);
   ElementTransformation *T = mesh->GetElementTransformation(i);
   double Jd[4];
   DenseMatrix J(Jd, 2, 2);

   T->Transform(ir, tr);

   vals.SetSize(ir.GetNPoints());
   for (int j = 0; j < ir.GetNPoints(); j++)
   {
      T->SetIntPoint(&ir.IntPoint(j));
      Geometries.JacToPerfJac(geom, T->Jacobian(), J);
      if (drawelems == 6) // attribute
      {
         vals(j) = mesh->GetAttribute(i);
      }
      else if (drawelems >= 4)
      {
         vals(j) = J.Det();
         // if (vals(j) >= 0.0)
         //    vals(j) = sqrt(vals(j));
         // else
         //    vals(j) = -sqrt(-vals(j));
      }
      else
      {
         vals(j) = J.CalcSingularvalue(0)/J.CalcSingularvalue(1);
         if (drawelems == 2)
         {
            vals(j) = vals(j) + 1.0/vals(j);
         }
      }
   }

   if (drawelems == 4)
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

   J.ClearExternalData();
}

void VisualizationSceneSolution::GetRefinedValues(
   int i, const IntegrationRule &ir, Vector &vals, DenseMatrix &tr)
{
   if (drawelems < 2)
   {
      rsol->GetValues(i, ir, vals, tr);
   }
   else
   {
      GetRefinedDetJ(i, ir, vals, tr);
   }

   if (logscale)
      for (int j = 0; j < vals.Size(); j++)
      {
         vals(j) = _LogVal(vals(j));
      }

   if (shrink != 1.0 || shrinkmat != 1.0)
   {
      ShrinkPoints(tr, i, 0, 0);
   }
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
   {
      GetRefinedDetJ(i, ir, vals, tr);
   }

   if (logscale)
   {
      if (have_normals)
         for (int j = 0; j < normals.Width(); j++)
            if (vals(j) >= minv && vals(j) <= maxv)
            {
               normals(0, j) *= log_a/vals(j);
               normals(1, j) *= log_a/vals(j);
            }
      for (int j = 0; j < vals.Size(); j++)
      {
         vals(j) = _LogVal(vals(j));
      }
   }

   if (shrink != 1.0 || shrinkmat != 1.0)
   {
      ShrinkPoints(tr, i, 0, 0);
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

void VisualizationSceneSolution::SetShading(int s, bool print)
{
   if (shading == s || s < 0)
   {
      return;
   }

   if (rsol)
   {
      if (s > 2)
      {
         return;
      }

      if (s == 2 || shading == 2)
      {
         shading = s;
         DoAutoscale(false);
         PrepareLines();
         PrepareBoundary();
         PrepareLevelCurves();
         PrepareCP();
         PrepareNumbering();
      }
      else
      {
         shading = s;
      }
   }
   else
   {
      if (s > 1)
      {
         return;
      }
      shading = s;
   }
   Prepare();

   static const char *shading_type[3] =
   {"flat", "smooth", "non-conforming (with subdivision)"};
   if (print)
   {
      cout << "Shading type : " << shading_type[shading] << endl;
   }
}

void VisualizationSceneSolution::ToggleShading()
{
   if (rsol)
   {
      SetShading((shading + 1) % 3, true);
   }
   else
   {
      SetShading(1 - shading, true);
   }
}

void VisualizationSceneSolution::SetRefineFactors(int tot, int bdr)
{
   if ((tot == TimesToRefine && bdr == EdgeRefineFactor) || tot < 1 || bdr < 1)
   {
      return;
   }

   if (tot % bdr)
   {
      tot += bdr - tot % bdr;
   }

   TimesToRefine = tot;
   EdgeRefineFactor = bdr;

   if (shading == 2)
   {
      DoAutoscale(false);
      PrepareLines();
      PrepareBoundary();
      Prepare();
      PrepareLevelCurves();
      PrepareCP();
   }
}

int VisualizationSceneSolution::GetAutoRefineFactor()
{
   int ne = mesh->GetNE(), ref = 1;

   while (ref < auto_ref_max && ne*(ref+1)*(ref+1) <= auto_ref_max_surf_elem)
   {
      ref++;
   }

   return ref;
}

void VisualizationSceneSolution::AutoRefine()
{
   int ref = GetAutoRefineFactor();

   cout << "Subdivision factors = " << ref << ", 1" << endl;

   SetRefineFactors(ref, 1);
}

void VisualizationSceneSolution::ToggleAttributes(Array<int> &attr_list)
{
   Array<int> &attr_marker = el_attr_to_show;

   for (int i = 0; i < attr_list.Size(); i++)
   {
      int attr = attr_list[i];
      if (attr < 1)
      {
         cout << "Hiding all attributes." << endl;
         attr_marker = 0;
      }
      else if (attr > attr_marker.Size())
      {
         cout << "Showing all attributes." << endl;
         attr_marker = 1;
      }
      else
      {
         attr_marker[attr-1] = !attr_marker[attr-1];
      }
   }
   PrepareLines();
   Prepare();
}

void VisualizationSceneSolution::SetNewScalingFromBox()
{
   if (scaling)
   {
      VisualizationSceneScalarData::SetNewScalingFromBox();
   }
   else
   {
      xscale = x[1]-x[0];
      yscale = y[1]-y[0];
      zscale = z[1]-z[0];
      xscale = (xscale < yscale) ? yscale : xscale;
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
         if ((*sol)(i) < rval[0]) { rval[0] = (*sol)(i); }
         if ((*sol)(i) > rval[1]) { rval[1] = (*sol)(i); }
      }
      rx[0] = rx[1] = coord[0];
      ry[0] = ry[1] = coord[1];

      for (i = 1; i < nv; i++)
      {
         coord = mesh->GetVertex(i);
         if (coord[0] < rx[0]) { rx[0] = coord[0]; }
         if (coord[1] < ry[0]) { ry[0] = coord[1]; }
         if (coord[0] > rx[1]) { rx[1] = coord[0]; }
         if (coord[1] > ry[1]) { ry[1] = coord[1]; }
      }
   }
   else
   {
      int ne = mesh -> GetNE();
      DenseMatrix pointmat;
      Vector values;
      RefinedGeometry *RefG;
      bool log_scale = logscale;

      logscale = false;
      rx[0] = ry[0] = rval[0] = numeric_limits<double>::infinity();
      rx[1] = ry[1] = rval[1] = -rx[0];
      for (i = 0; i < ne; i++)
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine, EdgeRefineFactor);
         GetRefinedValues(i, RefG->RefPts, values, pointmat);
         for (j = 0; j < values.Size(); j++)
         {
            if (isfinite(pointmat(0,j)))
            {
               if (pointmat(0,j) < rx[0]) { rx[0] = pointmat(0,j); }
               if (pointmat(0,j) > rx[1]) { rx[1] = pointmat(0,j); }
            }
            if (isfinite(pointmat(1,j)))
            {
               if (pointmat(1,j) < ry[0]) { ry[0] = pointmat(1,j); }
               if (pointmat(1,j) > ry[1]) { ry[1] = pointmat(1,j); }
            }
            if (isfinite(values(j)))
            {
               if (values(j) < rval[0]) { rval[0] = values(j); }
               if (values(j) > rval[1]) { rval[1] = values(j); }
            }
         }
      }
      logscale = log_scale;
   }
}

void VisualizationSceneSolution::FindNewBox(bool prepare)
{
   FindNewBox(x, y, z);

   minv = z[0];
   maxv = z[1];

   FixValueRange();

   z[0] = minv;
   z[1] = maxv;

   SetNewScalingFromBox(); // UpdateBoundingBox minus PrepareAxes
   UpdateValueRange(prepare);
}

void VisualizationSceneSolution::FindNewValueRange(bool prepare)
{
   double rx[2], ry[2], rv[2];

   FindNewBox(rx, ry, rv);
   minv = rv[0];
   maxv = rv[1];

   FixValueRange();

   UpdateValueRange(prepare);
}

void VisualizationSceneSolution::FindMeshBox(bool prepare)
{
   double rv[2];

   FindNewBox(x, y, rv);

   UpdateBoundingBox(); // SetNewScalingFromBox plus PrepareAxes
}

void VisualizationSceneSolution::ToggleLogscale(bool print)
{
   if (logscale || LogscaleRange())
   {
      // we do not change 'MySetColorLogscale' here. It is set to 0 in
      // Prepare() since we apply logarithmic scaling to the values.
      // In PrepareVectorField() we set 'MySetColorLogscale' to 'logscale'.
      logscale = !logscale;
      SetLogA();
      SetLevelLines(minv, maxv, nl);
      EventUpdateColors(); // Prepare() [+ PrepareVectorField() for vectors]
      PrepareLines();
      PrepareLevelCurves();
      PrepareBoundary();
      PrepareCP();
      if (print)
      {
         PrintLogscale(false);
      }
   }
   else if (print)
   {
      PrintLogscale(true);
   }
}

void DrawNumberedMarker(gl3::GlDrawable& buff, const double x[3], double dx, int n)
{
    gl3::GlBuilder bld = buff.createBuilder();
    bld.glBegin(GL_LINES);
    // glColor4d(0, 0, 0, 0);
    bld.glVertex3d(x[0]-dx, x[1]-dx, x[2]);
    bld.glVertex3d(x[0]+dx, x[1]+dx, x[2]);
    bld.glVertex3d(x[0]+dx, x[1]-dx, x[2]);
    bld.glVertex3d(x[0]-dx, x[1]+dx, x[2]);
    bld.glEnd();

    buff.addText(x[0], x[1], x[2], std::to_string(n));
}

void DrawTriangle(gl3::GlDrawable& buff,
                  const double (&pts)[4][3], const double (&cv)[4],
                  const double minv, const double maxv)
{
   double nor[3];
   if (Compute3DUnitNormal(pts[0], pts[1], pts[2], nor))
   {
      return;
   }
   
   float rgba[3][4];
   std::array<float, 2> texcoord[3];
   std::array<float, 3> fpts[3];
   std::array<float, 3> fnorm = {nor[0], nor[1], nor[2]};

   for (int i = 0; i < 3; i++) {
       texcoord[i][0] = cv[i];
       texcoord[i][1] = MySetColor(cv[i], minv, maxv, rgba[i]);
       fpts[i] = {pts[i][0], pts[i][1], pts[i][2]};
   }
   if (GetUseTexture()) {
       buff.addTriangle(gl3::VertexNormTex{fpts[0], fnorm, texcoord[0]},
                        gl3::VertexNormTex{fpts[1], fnorm, texcoord[1]},
                        gl3::VertexNormTex{fpts[2], fnorm, texcoord[2]});
   } else {
       buff.addTriangle(gl3::VertexNormColor{fpts[0], fnorm, gl3::ColorU8(rgba[0])},
                        gl3::VertexNormColor{fpts[1], fnorm, gl3::ColorU8(rgba[1])},
                        gl3::VertexNormColor{fpts[2], fnorm, gl3::ColorU8(rgba[2])});
   }
}

void DrawQuad(gl3::GlDrawable& buff,
              const double (&pts)[4][3], const double (&cv)[4],
              const double minv, const double maxv)
{
   double nor[3];
   if (Compute3DUnitNormal(pts[0], pts[1], pts[2], nor))
   {
      return;
   }
   
   std::array<float, 2> texcoord[4];
   float rgba[4][4];
   std::array<float, 3> fpts[4];
   std::array<float, 3> fnorm = {nor[0], nor[1], nor[2]};
   
   for (int i = 0; i < 4; i++) { 
       texcoord[i][0] = cv[i];
       texcoord[i][1] = MySetColor(cv[i], minv, maxv, rgba[i]);
       fpts[i] = {pts[i][0], pts[i][1], pts[i][2]};
   }
   if (GetUseTexture()) {
       buff.addQuad(gl3::VertexNormTex{fpts[0], fnorm, texcoord[0]},
                    gl3::VertexNormTex{fpts[1], fnorm, texcoord[1]},
                    gl3::VertexNormTex{fpts[2], fnorm, texcoord[2]},
                    gl3::VertexNormTex{fpts[3], fnorm, texcoord[3]});
   } else {
       buff.addQuad(gl3::VertexNormColor{fpts[0], fnorm, gl3::ColorU8(rgba[0])},
                    gl3::VertexNormColor{fpts[1], fnorm, gl3::ColorU8(rgba[1])},
                    gl3::VertexNormColor{fpts[2], fnorm, gl3::ColorU8(rgba[2])},
                    gl3::VertexNormColor{fpts[3], fnorm, gl3::ColorU8(rgba[3])});
   }
}

void RemoveFPErrors(const DenseMatrix &pts, Vector &vals, DenseMatrix &normals,
                    const int n, const Array<int> &ind, Array<int> &f_ind)
{
   int o = 0;

   f_ind.SetSize(ind.Size());
   for (int i = 0; i < ind.Size(); i += n)
   {
      bool good = true;
      for (int j = 0; j < n; j++)
      {
         f_ind[o+j] = ind[i+j];

         if (!isfinite(pts(0, ind[i+j])) || !isfinite(pts(1, ind[i+j])) ||
             !isfinite(pts(2, ind[i+j])) || !isfinite(vals(ind[i+j])))
            // check normals?
         {
            good = false;
            break;
         }
      }
      if (good)
      {
         o += n;
      }
   }
   f_ind.SetSize(o);
}

void DrawPatch(gl3::GlDrawable& drawable, const DenseMatrix &pts, Vector &vals, DenseMatrix &normals,
               const int n, const Array<int> &ind, const double minv,
               const double maxv, const int normals_opt)
{
   gl3::GlBuilder poly = drawable.createBuilder();
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
               {
                  normals(k, ind[i+j]) += na[k];
               }
      }
   }

   if (n == 3)
   {
      poly.glBegin(GL_TRIANGLES);
   }
   else
   {
      poly.glBegin(GL_QUADS);
   }
   if (normals_opt != 0 && normals_opt != -1)
   {
      if (normals_opt > 0)
      {
         for (int i = 0; i < ind.Size(); i++)
         {
            poly.glNormal3dv(&normals(0, ind[i]));
            MySetColor(poly, vals(ind[i]), minv, maxv);
            poly.glVertex3dv(&pts(0, ind[i]));
         }
      }
      else
      {
         for (int i = ind.Size()-1; i >= 0; i--)
         {
            poly.glNormal3dv(&normals(0, ind[i]));
            MySetColor(poly, vals(ind[i]), minv, maxv);
            poly.glVertex3dv(&pts(0, ind[i]));
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
               poly.glNormal3dv(na);
               for ( ; j < n; j++)
               {
                  MySetColor(poly, vals(ind[i]), minv, maxv);
                  poly.glVertex3dv(&pts(0, ind[i+j]));
               }
            }
            else
            {
               poly.glNormal3d(-na[0], -na[1], -na[2]);
               for (j = n-1; j >= 0; j--)
               {
                  MySetColor(poly, vals(ind[i]), minv, maxv);
                  poly.glVertex3dv(&pts(0, ind[i+j]));
               }
            }
         }
      }
   }
   poly.glEnd();
}



void VisualizationSceneSolution::PrepareWithNormals()
{
   disp_buf.clear();
   gl3::GlBuilder poly = disp_buf.createBuilder();
   Array<int> vertices;
   double *vtx, *nor, val, s;

   for (int i = 0; i < mesh->GetNE(); i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      mesh->GetElementVertices(i, vertices);
      GLenum shape;
      if (vertices.Size() == 3)
      {
         shape = GL_TRIANGLES;
      }
      else
      {
         shape = GL_QUADS;
      }
      poly.glBegin(shape);
      for (int j = 0; j < vertices.Size(); j++)
      {
         vtx = mesh->GetVertex(vertices[j]);
         nor = &(*v_normals)(3*vertices[j]);
         val = (*sol)(vertices[j]);
         if (logscale && val >= minv && val <= maxv)
         {
            s = log_a/val;
            val = _LogVal_(val);
            poly.glNormal3d(s*nor[0], s*nor[1], nor[2]);
         }
         else
         {
            poly.glNormal3dv(nor);
         }
         MySetColor(poly, val, minv, maxv);
         poly.glVertex3d(vtx[0], vtx[1], val);
      }
      poly.glEnd();
   }
   disp_buf.buffer();
}

void VisualizationSceneSolution::PrepareFlat()
{
   int i, j;
   disp_buf.clear();
   int ne = mesh -> GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double pts[4][3], col[4];

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      mesh->GetPointMatrix (i, pointmat);
      mesh->GetElementVertices (i, vertices);

      for (j = 0; j < pointmat.Width(); j++)
      {
         pts[j][0] = pointmat(0, j);
         pts[j][1] = pointmat(1, j);
         pts[j][2] = col[j] = LogVal((*sol)(vertices[j]));
      }
      if (j == 3)
      {
         DrawTriangle(disp_buf, pts, col, minv, maxv);
      }
      else
      {
         DrawQuad(disp_buf, pts, col, minv, maxv);
      }
   }
   disp_buf.buffer();
}

// determines how quads and their level lines are drawn
// when using subdivision:
// 0 - draw a quad
// 1 - draw 2 triangles (split using the '0-2' diagonal)
// 2 - draw 4 triangles (split using both diagonals)
const int split_quads = 1;

void VisualizationSceneSolution::PrepareFlat2()
{
   int i, j, k;
   disp_buf.clear();
   int ne = mesh -> GetNE();
   DenseMatrix pointmat, pts3d, normals;
   Vector values;
   RefinedGeometry *RefG;
   Array<int> fRG;

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
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
      RemoveFPErrors(pts3d, values, normals, sides, RG, fRG);
      DrawPatch(disp_buf, pts3d, values, normals, sides, fRG, minv, maxv, j);
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
            {
               DrawTriangle(pts, col, minv, maxv);
            }
            else
            {
               DrawQuad(pts, col, minv, maxv);
            }
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
   disp_buf.buffer();
}

void VisualizationSceneSolution::Prepare()
{
   MySetColorLogscale = 0;

   switch (shading)
   {
      case 0:
         PrepareFlat();
         return;
      case 2:
         PrepareFlat2();
         return;
      default:
         if (v_normals)
         {
            PrepareWithNormals();
            return;
         }
         break;
   }

   int i, j;

   disp_buf.clear();
   gl3::GlBuilder poly = disp_buf.createBuilder();
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

      if (!el_attr_to_show[mesh -> attributes[d]-1]) { continue; }

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
               p[j][2] = LogVal((*sol)(vertices[j]));
            }

            if (pointmat.Width() == 3)
            {
               j = Compute3DUnitNormal(p[0], p[1], p[2], nor);
            }
            else
            {
               j = Compute3DUnitNormal(p[0], p[1], p[2], p[3], nor);
            }

            if (j == 0)
               for (j = 0; j < pointmat.Size(); j++)
               {
                  nx(vertices[j]) += nor[0];
                  ny(vertices[j]) += nor[1];
                  nz(vertices[j]) += nor[2];
               }
         }

      for (i = 0; i < ne; i++)
      {
         if (mesh -> GetAttribute(i) == mesh -> attributes[d])
         {
            GLenum shape;
            switch (mesh->GetElementType(i))
            {
               case Element::TRIANGLE:
                  shape = GL_TRIANGLES;
                  break;

               case Element::QUADRILATERAL:
                  shape = GL_QUADS;
                  break;
            }
            poly.glBegin(shape);
            mesh->GetPointMatrix (i, pointmat);
            mesh->GetElementVertices (i, vertices);

            for (j = 0; j < pointmat.Size(); j++)
            {
               double z = LogVal((*sol)(vertices[j]));
               MySetColor(poly, z, minv, maxv);
               poly.glNormal3d(nx(vertices[j]), ny(vertices[j]), nz(vertices[j]));
               poly.glVertex3d(pointmat(0, j), pointmat(1, j), z);
            }
            poly.glEnd();
         }
      }
   }
   disp_buf.buffer();
}

void VisualizationSceneSolution::PrepareLevelCurves()
{
   if (shading == 2)
   {
      PrepareLevelCurves2();
      return;
   }

   static int vt[4] = { 0, 1, 2, 3 };
   Array<int> RG(vt, 4), vertices;
   Vector values;
   DenseMatrix pointmat;

   lcurve_buf.clear();
   gl3::GlBuilder build = lcurve_buf.createBuilder();
   for (int i = 0; i < mesh->GetNE(); i++)
   {
      mesh->GetElementVertices(i, vertices);
      mesh->GetPointMatrix(i, pointmat);
      sol->GetSubVector(vertices, values);
      if (logscale)
         for (int j = 0; j < vertices.Size(); j++)
         {
            values(j) = _LogVal(values(j));
         }
      RG.SetSize(vertices.Size());
      DrawLevelCurves(build, RG, pointmat, values, vertices.Size(), level);
   }
   lcurve_buf.buffer();
}

void VisualizationSceneSolution::DrawLevelCurves(
   gl3::GlBuilder& builder, Array<int> &RG, DenseMatrix &pointmat, Vector &values,
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
         DrawPolygonLevelLines(builder, point[0], sides, lvl, logscale);
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
            DrawPolygonLevelLines(builder, point[0], 3, lvl, logscale);
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
         {
            point[2][2] = zc;
         }

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

            DrawPolygonLevelLines(builder, point[0], 3, lvl, logscale);
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

   lcurve_buf.clear();
   gl3::GlBuilder build = lcurve_buf.createBuilder();
   for (i = 0; i < ne; i++)
   {
      RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      GetRefinedValues (i, RefG->RefPts, values, pointmat);
      Array<int> &RG = RefG->RefGeoms;
      int sides = mesh->GetElement(i)->GetNVertices();

      DrawLevelCurves(build, RG, pointmat, values, sides, level);
   }
   lcurve_buf.buffer();
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

   line_buf.clear();
   gl3::GlBuilder lb = line_buf.createBuilder();

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      lb.glBegin(GL_LINE_LOOP);
      mesh->GetPointMatrix (i, pointmat);
      mesh->GetElementVertices (i, vertices);

      for (j = 0; j < pointmat.Size(); j++)
         lb.glVertex3d(pointmat(0, j), pointmat(1, j),
                    LogVal((*sol)(vertices[j])));
      lb.glEnd();
   }

   line_buf.buffer();
}

double VisualizationSceneSolution::GetElementLengthScale(int k)
{
   DenseMatrix pointmat;
   Array<int> vertices;

   mesh->GetPointMatrix(k, pointmat);
   mesh->GetElementVertices(k, vertices);

   // Get length scale for x mark
   double xmax = -numeric_limits<double>::infinity();
   double ymax = -numeric_limits<double>::infinity();
   double xmin = numeric_limits<double>::infinity();
   double ymin = numeric_limits<double>::infinity();

   int nv = vertices.Size();
   for (int j = 0; j < nv; j++)
   {
      double x = pointmat(0,j);
      double y = pointmat(1,j);
      if (x > xmax) { xmax = x; }
      if (x < xmin) { xmin = x; }
      if (y > ymax) { ymax = y; }
      if (y < ymin) { ymin = y; }
   }
   double dx = xmax-xmin;
   double dy = ymax-ymin;
   double ds = std::min<double>(dx,dy);

   return ds;
}

void VisualizationSceneSolution::PrepareElementNumbering()
{
   int ne = mesh -> GetNE();

   if (ne > MAX_RENDER_NUMBERING)
   {
      cout << "Element numbering disabled when #elements > "
           << MAX_RENDER_NUMBERING << endl;
      cout << "Rendering the text would be very slow." << endl;
      return;
   }

   if (2 == shading)
   {
      PrepareElementNumbering2();
   }
   else
   {
      PrepareElementNumbering1();
   }
}

void VisualizationSceneSolution::PrepareElementNumbering1()
{
   e_nums_buf.clear();

   DenseMatrix pointmat;
   Array<int> vertices;

   int ne = mesh->GetNE();
   for (int k = 0; k < ne; k++)
   {
      mesh->GetPointMatrix (k, pointmat);
      mesh->GetElementVertices (k, vertices);
      int nv = vertices.Size();

      ShrinkPoints(pointmat, k, 0, 0);

      double xs = 0.0;
      double ys = 0.0;
      double us = 0.0;
      for (int j = 0; j < nv; j++)
      {
         xs += pointmat(0,j);
         ys += pointmat(1,j);
         us += LogVal((*sol)(vertices[j]));
      }
      xs /= nv;
      ys /= nv;
      us /= nv;

      double ds = GetElementLengthScale(k);
      double dx = 0.05*ds;

      double xx[3] = {xs,ys,us};
      DrawNumberedMarker(e_nums_buf,xx,dx,k);
   }

   e_nums_buf.buffer();
}

void VisualizationSceneSolution::PrepareElementNumbering2()
{
   IntegrationRule center_ir(1);
   DenseMatrix pointmat;
   Vector values;

   e_nums_buf.clear();

   int ne = mesh->GetNE();
   for (int i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      center_ir.IntPoint(0) =
         Geometries.GetCenter(mesh->GetElementBaseGeometry(i));
      GetRefinedValues (i, center_ir, values, pointmat);

      double xc = pointmat(0,0);
      double yc = pointmat(1,0);
      double uc = values(0);

      double ds = GetElementLengthScale(i);
      double dx = 0.05*ds;

      double xx[3] = {xc,yc,uc};
      DrawNumberedMarker(e_nums_buf,xx,dx,i);
   }

   e_nums_buf.buffer();
}

void VisualizationSceneSolution::PrepareVertexNumbering()
{
   int nv = mesh->GetNV();

   if (nv > MAX_RENDER_NUMBERING)
   {
      cout << "Vertex numbering disabled when #vertices > "
           << MAX_RENDER_NUMBERING << endl;
      cout << "Rendering the text would be very slow." << endl;
      return;
   }

   if (2 == shading)
   {
      PrepareVertexNumbering2();
   }
   else
   {
      PrepareVertexNumbering1();
   }
}

void VisualizationSceneSolution::PrepareVertexNumbering1()
{
   v_nums_buf.clear();

   DenseMatrix pointmat;
   Array<int> vertices;

   // Draw the vertices for each element.  This is redundant, except
   // when the elements or domains are shrunk.

   const int ne = mesh->GetNE();
   for (int k = 0; k < ne; k++)
   {
      mesh->GetPointMatrix (k, pointmat);
      mesh->GetElementVertices (k, vertices);
      int nv = vertices.Size();

      ShrinkPoints(pointmat, k, 0, 0);

      double ds = GetElementLengthScale(k);
      double xs = 0.05*ds;

      for (int j = 0; j < nv; j++)
      {
         double x = pointmat(0,j);
         double y = pointmat(1,j);
         double u = LogVal((*sol)(vertices[j]));

         double xx[3] = {x,y,u};
         DrawNumberedMarker(v_nums_buf,xx,xs,vertices[j]);
      }
   }

   v_nums_buf.buffer();
}

void VisualizationSceneSolution::PrepareVertexNumbering2()
{
   DenseMatrix pointmat;
   Vector values;
   Array<int> vertices;

   v_nums_buf.clear();

   const int ne = mesh->GetNE();
   for (int i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      mesh->GetElementVertices (i, vertices);

      const IntegrationRule &vert_ir =
         *Geometries.GetVertices(mesh->GetElementBaseGeometry(i));

      GetRefinedValues (i, vert_ir, values, pointmat);

      double ds = GetElementLengthScale(i);
      double xs = 0.05*ds;

      for (int j = 0; j < values.Size(); j++)
      {
         double xv = pointmat(0, j);
         double yv = pointmat(1, j);

         double u = values[j];

         double xx[3] = {xv,yv,u};
         DrawNumberedMarker(v_nums_buf,xx,xs,vertices[j]);
      }
   }

   v_nums_buf.buffer();
}

void VisualizationSceneSolution::PrepareNumbering()
{
   PrepareElementNumbering();
   PrepareVertexNumbering();
}

void VisualizationSceneSolution::PrepareLines2()
{
   int i, j, k, ne = mesh -> GetNE();
   Vector values;
   DenseMatrix pointmat;
   RefinedGeometry *RefG;

   line_buf.clear();
   gl3::GlBuilder lb = line_buf.createBuilder();

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

      RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      GetRefinedValues (i, RefG->RefPts, values, pointmat);
      Array<int> &RG = RefG->RefGeoms;
      int sides = mesh->GetElement(i)->GetNVertices();

      for (k = 0; k < RG.Size()/sides; k++)
      {
         lb.glBegin(GL_LINE_LOOP);

         for (j = 0; j < sides; j++)
            lb.glVertex3d(pointmat(0, RG[sides*k+j]),
                       pointmat(1, RG[sides*k+j]),
                       values(RG[sides*k+j]));
         lb.glEnd();
      }
   }

   line_buf.buffer();
}

void VisualizationSceneSolution::PrepareLines3()
{
   int i, k, ne = mesh -> GetNE();
   Vector values;
   DenseMatrix pointmat;
   RefinedGeometry *RefG;

   line_buf.clear();
   gl3::GlBuilder lb = line_buf.createBuilder();

   for (i = 0; i < ne; i++)
   {
      if (!el_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }
      RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                         TimesToRefine, EdgeRefineFactor);
      GetRefinedValues (i, RefG->RefPts, values, pointmat);
      Array<int> &RE = RefG->RefEdges;

      lb.glBegin (GL_LINES);
      for (k = 0; k < RE.Size()/2; k++)
      {
         lb.glVertex3d (pointmat(0, RE[2*k]),
                     pointmat(1, RE[2*k]),
                     values(RE[2*k]));
         lb.glVertex3d (pointmat(0, RE[2*k+1]),
                     pointmat(1, RE[2*k+1]),
                     values(RE[2*k+1]));
      }
      lb.glEnd();
   }

   line_buf.buffer();
}

void VisualizationSceneSolution::UpdateValueRange(bool prepare)
{
   bool had_logscale = logscale;
   logscale = logscale && LogscaleRange();
   SetLogA();
   SetLevelLines(minv, maxv, nl);
   // preserve the current box z-size
   zscale *= (z[1]-z[0])/(maxv-minv);
   z[0] = minv;
   z[1] = maxv;
   PrepareAxes();
   if (prepare)
   {
      UpdateLevelLines();
      EventUpdateColors();
      if (had_logscale)
      {
         PrepareLines();
         PrepareBoundary();
         PrepareCP();
      }
   }
}

void VisualizationSceneSolution::PrepareBoundary()
{
   int i, j, ne = mesh->GetNBE();
   Array<int> vertices;
   DenseMatrix pointmat;

   bdr_buf.clear();
   gl3::GlBuilder bl = bdr_buf.createBuilder();
   if (shading != 2)
   {
      bl.glBegin(GL_LINES);
      for (i = 0; i < ne; i++)
      {
         if (!bdr_el_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }
         mesh->GetBdrElementVertices(i, vertices);
         mesh->GetBdrPointMatrix(i, pointmat);
         for (j = 0; j < pointmat.Size(); j++)
            bl.glVertex3d(pointmat(0, j), pointmat(1, j),
                       LogVal((*sol)(vertices[j])));
      }
      bl.glEnd();
   }
   else // shading == 2
   {
      int en;
      FaceElementTransformations *T;
      RefinedGeometry *RefG =
         GLVisGeometryRefiner.Refine(Geometry::SEGMENT, TimesToRefine,
                                     EdgeRefineFactor);
      IntegrationRule &ir = RefG->RefPts;
      IntegrationRule eir(ir.GetNPoints());
      Vector vals;
      double shr = shrink;
      shrink = 1.0;

      for (i = 0; i < ne; i++)
      {
         if (!bdr_el_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }
         en = mesh->GetBdrElementEdgeIndex(i);
         T = mesh->GetFaceElementTransformations(en, 4);
         T->Loc1.Transform(ir, eir);
         GetRefinedValues(T->Elem1No, eir, vals, pointmat);
         bl.glBegin(GL_LINE_STRIP);
         for (j = 0; j < vals.Size(); j++)
         {
            bl.glVertex3d(pointmat(0, j), pointmat(1, j), vals(j));
         }
         bl.glEnd();

         if (T->Elem2No >= 0)
         {
            T = mesh->GetFaceElementTransformations(en, 8);
            T->Loc2.Transform(ir, eir);
            GetRefinedValues(T->Elem2No, eir, vals, pointmat);
            bl.glBegin(GL_LINE_STRIP);
            for (j = 0; j < vals.Size(); j++)
            {
               bl.glVertex3d(pointmat(0, j), pointmat(1, j), vals(j));
            }
            bl.glEnd();
         }
      }
      shrink = shr;
   }

   bdr_buf.buffer();
}

void VisualizationSceneSolution::PrepareCP()
{
   Vector values;
   DenseMatrix pointmat;
   Array<int> ind;

   if (draw_cp == 0)
   {
      return;
   }

   cp_buf.clear();
   gl3::GlBuilder bld = cp_buf.createBuilder();
   bld.glBegin(GL_LINES);

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
            {
               n++;
            }
         }
         if (n == 0 || n == pointmat.Width())
         {
            continue;
         }

         mesh->GetElementVertices(i, vertices);
         values.SetSize(vertices.Size());
         ind.SetSize(vertices.Size());
         for (int j = 0; j < values.Size(); j++)
         {
            values(j) = LogVal((*sol)(vertices[j]));
            ind[j] = j;
         }

         DrawCPLine(bld, pointmat, values, ind);
      }
   }
   else
   {
      RefinedGeometry *RefG;

      for (int i = 0; i < mesh->GetNE(); i++)
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine, EdgeRefineFactor);
         GetRefinedValues (i, RefG->RefPts, values, pointmat);
         Array<int> &RG = RefG->RefGeoms;
         int sides = mesh->GetElement(i)->GetNVertices();

         ind.SetSize(sides);
         for (int k = 0; k < RG.Size()/sides; k++)
         {
            for (int j = 0; j < sides; j++)
            {
               ind[j] = RG[k*sides+j];
            }

            int n = 0;
            for (int j = 0; j < sides; j++)
            {
               const double s =
                  CuttingPlane->Transform(pointmat(0, ind[j]),
                                          pointmat(1, ind[j]), 0.0);
               if (s >= 0.0)
               {
                  n++;
               }
            }
            if (n == 0 || n == sides)
            {
               continue;
            }

            DrawCPLine(bld, pointmat, values, ind);
         }
      }
   }

   bld.glEnd();
   cp_buf.buffer();
}

void VisualizationSceneSolution::DrawCPLine(
   gl3::GlBuilder& bld, DenseMatrix &pointmat, Vector &values, Array<int> &ind)
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

         bld.glVertex3d((1.-a) * xs + a * xt,
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
   gl->enableDepthTest();

   Set_Background();
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // model transformation
   ModelView();

   glPolygonOffset (1, 1);
   glEnable (GL_POLYGON_OFFSET_FILL);
   //glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   gl->disableClipPlane();
   gl->disableLight();

#if 0
   // Testing: moved the drawing of the colorbar at the end. If there are no
   // undesired effects we can delete this disabled code. See also below.

   // draw colorbar
   if (colorbar)
   {
      if (drawmesh == 2)
      {
         DrawColorBar(minv,maxv,&level);
      }
      else
      {
         DrawColorBar(minv,maxv);
      }
   }
#endif

   if (draw_cp)
   {
      gl->setClipPlane(CuttingPlane->Equation());
      gl->enableClipPlane();
   }

   Set_Material();
   if (light)
   {
      gl->enableLight();
   }

   if (MatAlpha < 1.0)
   {
      Set_Transparency();
   }

   // draw elements
   if (drawelems)
   {
      disp_buf.draw();
   }

   if (MatAlpha < 1.0)
   {
      Remove_Transparency();
   }

   if (light)
   {
      gl->disableLight();
   }
   Set_Black_Material();

   // ruler may have mixture of polygons and lines
   if (draw_cp)
   {
      gl->disableClipPlane();
      DrawRuler(logscale);
      cp_buf.draw();
      gl->enableClipPlane();
   }
   else
   {
      DrawRuler(logscale);
   }
   if (drawbdr)
   {
      bdr_buf.draw();
   }

   // draw lines
   if (drawmesh == 1)
   {
      line_buf.draw();
   }
   else if (drawmesh == 2)
   {
      lcurve_buf.draw();
   }

   // draw numberings
   if (drawnums)
   {
      if (1 == drawnums)
      {
         e_nums_buf.draw();
      }
      else if (2 == drawnums)
      {
         v_nums_buf.draw();
      }
   }

   if (draw_cp)
   {
      gl->disableClipPlane();
   }

   // draw axes
   if (drawaxes)
   {
      axes_buf.draw();
      DrawCoordinateCross();
   }

#if 1
   // Testing: moved the drawing of the colorbar from the beginning. If there
   // are no undesired effects we can remove this comment and "#if 1". We can
   // also do the same in vector, 3D, and vector 3D modes.

   // draw colorbar
   if (colorbar)
   {
      if (drawmesh == 2)
      {
         DrawColorBar(minv,maxv,&level);
      }
      else
      {
         DrawColorBar(minv,maxv);
      }
   }
#endif
}
