// Copyright (c) 2010-2022, Lawrence Livermore National Security, LLC. Produced
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
#include <string>
#include <limits>
#include <cmath>
#include <vector>

#include "mfem.hpp"
#include "visual.hpp"
#include "palettes.hpp"
#include "gltf.hpp"

using namespace mfem;
using namespace std;


thread_local VisualizationSceneSolution *vssol;
extern thread_local VisualizationScene  *locscene;
extern thread_local GeometryRefiner GLVisGeometryRefiner;

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

std::string VisualizationSceneSolution::GetHelpString() const
{
   std::stringstream os;
   os << endl
      << "+------------------------------------+" << endl
      << "| Keys                               |" << endl
      << "+------------------------------------+" << endl
      << "| a -  Displays/Hides the axes       |" << endl
      << "| A -  Turns antialiasing on/off     |" << endl
      << "| b/B  Toggle 2D boundary            |" << endl
      << "| c -  Toggle colorbar and caption   |" << endl
      << "| C -  Change the main plot caption  |" << endl
      << "| e -  Displays/Hides the elements   |" << endl
      << "| f -  Smooth/Nonconf/Flat shading   |" << endl
      << "| g -  Toggle background             |" << endl
      << "| h -  Displays help menu            |" << endl
      << "| i -  Toggle the cutting plane      |" << endl
      << "| j -  Turn on/off perspective       |" << endl
      << "| k/K  Adjust the transparency level |" << endl
      << "| ,/<  Adjust color transparency     |" << endl
      << "| l -  Turns on/off the light        |" << endl
      << "| L -  Toggle logarithmic scale      |" << endl
      << "| m -  Displays/Hides the mesh       |" << endl
      << "| n/N  Cycle through numberings      |" << endl
      << "| o -  (De)refine elem. (NC shading) |" << endl
      << "| O -  Switch 'o' func. (NC shading) |" << endl
      << "| p/P  Cycle through color palettes  |" << endl
      << "| q -  Quits                         |" << endl
      << "| r -  Reset the plot to 3D view     |" << endl
      << "| R -  Reset the plot to 2D view     |" << endl
      << "| s -  Turn on/off unit cube scaling |" << endl
      << "| S -  Take snapshot/Record a movie  |" << endl
      << "| t -  Cycle materials and lights    |" << endl
      << "| y/Y  Rotate the cutting plane      |" << endl
      << "| z/Z  Move the cutting plane        |" << endl
      << "| \\ -  Set light source position     |" << endl
      << "| Ctrl+o - Element ordering curve    |" << endl
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

void VisualizationSceneSolution::ToggleDrawBdr()
{
   drawbdr = (drawbdr+1)%3;

   // Switch the colorbar range to correspond to the boundary attributes when
   // showing them in color (drawbdr == 2) and restore it back when not.
   if (drawbdr == 2)
   {
      if (drawelems == 1 || drawelems == 0)
      {
         minv_sol = minv;
         maxv_sol = maxv;
         have_sol_range = true;
      }
      drawelems = 0;
      minv = 1;
      maxv = mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 1;
      FixValueRange();
      UpdateValueRange(true);
      PrepareBoundary();
   }
   else if (drawbdr == 1)
   {
      PrepareBoundary();
   }
   else if (drawbdr == 0)
   {
      drawelems = 1;
      if (have_sol_range)
      {
         minv = minv_sol;
         maxv = maxv_sol;
         bb.z[0] = minv;
         bb.z[1] = maxv;
         SetNewScalingFromBox(); // UpdateBoundingBox minus PrepareAxes
         UpdateValueRange(true);
      }
      else
      {
         FindNewValueRange(true);
      }
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

static void KeyoPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      vssol -> ToggleDrawOrdering();
      vssol -> PrepareOrderingCurve();
      SendExposeEvent();
   }
   else
   {
      vssol->ToggleRefinements();
   }
}

static void KeyOPressed(GLenum state)
{
   (void)state;
   vssol->ToggleRefinementFunction();
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

void KeyiPressed()
{
   vssol->ToggleDrawCP();
   SendExposeEvent();
}

void KeyIPressed()
{
   // no-op, available
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
      vssol->PrepareOrderingCurve();
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
}

void VisualizationSceneSolution::Init()
{
   rsol  = NULL;
   vssol = this;

   drawelems = shading = 1;
   drawmesh  = 0;
   draworder = 0;
   drawnums  = 0;

   refine_func = 0;
   have_sol_range = false;

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

   palette.SetIndex(2); // use the 'jet-like' palette in 2D

   double eps = 1e-6; // move the cutting plane a bit to avoid artifacts
   CuttingPlane = new Plane(-1.0,0.0,0.0,(0.5-eps)*bb.x[0]+(0.5+eps)*bb.x[1]);
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

      wnd->setOnKeyDown('o', KeyoPressed);
      wnd->setOnKeyDown('O', KeyOPressed);

      wnd->setOnKeyDown('e', KeyEPressed);
      wnd->setOnKeyDown('E', KeyEPressed);

      wnd->setOnKeyDown('f', KeyFPressed);
      wnd->setOnKeyDown('F', KeyFPressed);

      wnd->setOnKeyDown('i', KeyiPressed);
      wnd->setOnKeyDown('I', KeyIPressed);

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
   PrepareOrderingCurve();
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

   if (drawbdr == 2) { return; }

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

   if (drawelems == 0)
   {
      minv_sol = minv;
      maxv_sol = maxv;
      have_sol_range = true;
   }
   else if (shading == 2)
   {
      if (drawelems == 1 && have_sol_range)
      {
         minv = minv_sol;
         maxv = maxv_sol;
         bb.z[0] = minv;
         bb.z[1] = maxv;
         SetNewScalingFromBox(); // UpdateBoundingBox minus PrepareAxes
         UpdateValueRange(false);
      }
      else
      {
         DoAutoscaleValue(false);
      }
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

   have_sol_range = false;
   DoAutoscale(false);

   Prepare();
   PrepareLines();
   PrepareLevelCurves();
   PrepareBoundary();
   PrepareCP();
   PrepareNumbering();
   PrepareOrderingCurve();
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
      {
         for (int j = 0; j < normals.Width(); j++)
         {
            if (vals(j) >= minv && vals(j) <= maxv)
            {
               normals(0, j) *= log_a/vals(j);
               normals(1, j) *= log_a/vals(j);
            }
         }
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
         have_sol_range = false;
         DoAutoscale(false);
         PrepareLines();
         PrepareBoundary();
         PrepareLevelCurves();
         PrepareCP();
         PrepareNumbering();
         PrepareOrderingCurve();
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

void VisualizationSceneSolution::ToggleRefinements()
{
   int update = 1;
   switch (refine_func)
   {
      case 0:
         TimesToRefine += EdgeRefineFactor;
         break;
      case 1:
         if (TimesToRefine > EdgeRefineFactor)
         {
            TimesToRefine -= EdgeRefineFactor;
         }
         else
         {
            update = 0;
         }
         break;
      case 2:
         TimesToRefine /= EdgeRefineFactor;
         EdgeRefineFactor++;
         TimesToRefine *= EdgeRefineFactor;
         break;
      case 3:
         if (EdgeRefineFactor > 1)
         {
            TimesToRefine /= EdgeRefineFactor;
            EdgeRefineFactor--;
            TimesToRefine *= EdgeRefineFactor;
         }
         else
         {
            update = 0;
         }
         break;
   }
   if (update && shading == 2)
   {
      have_sol_range = false;
      DoAutoscale(false);
      PrepareLines();
      PrepareBoundary();
      Prepare();
      PrepareLevelCurves();
      PrepareCP();
      SendExposeEvent();
   }
   cout << "Subdivision factors = " << TimesToRefine << ", " << EdgeRefineFactor
        << endl;
}

void VisualizationSceneSolution::ToggleRefinementFunction()
{
   refine_func = (refine_func+1)%4;
   cout << "Key 'o' will: ";
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
      have_sol_range = false;
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
      xscale = bb.x[1]-bb.x[0];
      yscale = bb.y[1]-bb.y[0];
      zscale = bb.z[1]-bb.z[0];
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
   FindNewBox(bb.x, bb.y, bb.z);

   minv = bb.z[0];
   maxv = bb.z[1];

   FixValueRange();

   bb.z[0] = minv;
   bb.z[1] = maxv;

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

   FindNewBox(bb.x, bb.y, rv);

   UpdateBoundingBox(); // SetNewScalingFromBox plus PrepareAxes
}

void VisualizationSceneSolution::ToggleLogscale(bool print)
{
   if (logscale || LogscaleRange())
   {
      // we do not change the palette logscale setting here. It is set to 0 in
      // Prepare() since we apply logarithmic scaling to the values.
      // In PrepareVectorField() we call 'palette.SetUseLogscale(logscale)'.
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

void VisualizationSceneSolution::EventUpdateColors()
{
   Prepare();
   PrepareOrderingCurve();
}

void VisualizationSceneSolution::EventUpdateBackground()
{
   PrepareNumbering();
}

void DrawNumberedMarker(gl3::GlDrawable& buff, const double x[3], double dx,
                        int n)
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
   updated_bufs.emplace_back(&disp_buf);
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
   updated_bufs.emplace_back(&disp_buf);
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
   updated_bufs.emplace_back(&disp_buf);
}

void VisualizationSceneSolution::Prepare()
{
   palette.SetUseLogscale(0);

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
            GLenum shape = GL_NONE;
            switch (mesh->GetElementType(i))
            {
               case Element::TRIANGLE:
                  shape = GL_TRIANGLES;
                  break;
               case Element::QUADRILATERAL:
                  shape = GL_QUADS;
                  break;
               default:
                  MFEM_ABORT("Invalid 2D element type");
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
   updated_bufs.emplace_back(&disp_buf);
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
   updated_bufs.emplace_back(&lcurve_buf);
}

void VisualizationSceneSolution::DrawLevelCurves(
   gl3::GlBuilder& builder, Array<int> &RG, DenseMatrix &pointmat, Vector &values,
   int sides, Array<double> &lvl, int flat)
{
   double point[4][4];
   // double zc = 0.5*(z[0]+z[1]);
   double zc = bb.z[1];

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
   updated_bufs.emplace_back(&lcurve_buf);
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

   updated_bufs.emplace_back(&line_buf);
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

   updated_bufs.emplace_back(&e_nums_buf);
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

   updated_bufs.emplace_back(&e_nums_buf);
}

void VisualizationSceneSolution::PrepareVertexNumbering()
{
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

   updated_bufs.emplace_back(&v_nums_buf);
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

   updated_bufs.emplace_back(&v_nums_buf);
}

void VisualizationSceneSolution::PrepareEdgeNumbering()
{
   f_nums_buf.clear();

   DenseMatrix p;
   Array<int> vertices;
   Array<int> edges;
   Array<int> edges_ori;

   const int ne = mesh->GetNE();
   for (int k = 0; k < ne; k++)
   {
      mesh->GetElementEdges(k, edges, edges_ori);

      double ds = GetElementLengthScale(k);
      double xs = 0.05 * ds;

      for (int i = 0; i < edges.Size(); i++)
      {
         mesh->GetEdgeVertices(edges[i], vertices);

         p.SetSize(mesh->Dimension(), vertices.Size());
         p.SetCol(0, mesh->GetVertex(vertices[0]));
         p.SetCol(1, mesh->GetVertex(vertices[1]));

         ShrinkPoints(p, k, 0, 0);

         const double m[2] = {0.5 * (p(0,0) + p(0,1)), 0.5 * (p(1,0) + p(1,1))};
         // TODO: figure out something better...
         double u = LogVal(0.5 * ((*sol)(vertices[0]) + (*sol)(vertices[1])));

         double xx[3] = {m[0], m[1], u};
         DrawNumberedMarker(f_nums_buf, xx, xs, edges[i]);
      }
   }

   updated_bufs.emplace_back(&f_nums_buf);
}

void VisualizationSceneSolution::PrepareOrderingCurve()
{
   bool color = draworder < 3;
   order_buf.clear();
   order_noarrow_buf.clear();
   PrepareOrderingCurve1(order_buf, true, color);
   PrepareOrderingCurve1(order_noarrow_buf, false, color);
   updated_bufs.emplace_back(&order_buf);
   updated_bufs.emplace_back(&order_noarrow_buf);
}

void VisualizationSceneSolution::PrepareOrderingCurve1(gl3::GlDrawable& buf,
                                                       bool arrows,
                                                       bool color)
{
   gl3::GlBuilder builder = buf.createBuilder();
   DenseMatrix pointmat;
   Array<int> vertices;

   DenseMatrix pointmat1;
   Array<int> vertices1;

   int ne = mesh->GetNE();
   for (int k = 0; k < ne-1; k++)
   {
      mesh->GetPointMatrix (k, pointmat);
      mesh->GetElementVertices (k, vertices);
      mesh->GetPointMatrix (k+1, pointmat1);
      mesh->GetElementVertices (k+1, vertices1);
      int nv = vertices.Size();
      int nv1 = vertices1.Size();

      ShrinkPoints(pointmat, k, 0, 0);
      ShrinkPoints(pointmat1, k+1, 0, 0);

      double xs = 0.0;
      double ys = 0.0;
      double us = 0.0;
      for (int j = 0; j < nv; j++)
      {
         xs += pointmat(0,j);
         ys += pointmat(1,j);
         us += maxv + double(k)/ne*(maxv-minv);
      }
      xs /= nv;
      ys /= nv;
      us /= nv;

      double xs1 = 0.0;
      double ys1 = 0.0;
      double us1 = 0.0;
      for (int j = 0; j < nv1; j++)
      {
         xs1 += pointmat1(0,j);
         ys1 += pointmat1(1,j);
         us1 += maxv + double(k+1)/ne*(maxv-minv);
      }
      xs1 /= nv1;
      ys1 /= nv1;
      us1 /= nv1;

      double dx = xs1-xs;
      double dy = ys1-ys;
      double du = us1-us;
      double ds = sqrt(dx*dx+dy*dy+du*du);

      if (color)
      {
         double cval = minv+double(k)/ne*(maxv-minv);
         MySetColor(builder, cval, minv, maxv);
      }

      if (arrows)
      {
         Arrow3(builder,
                xs,ys,us,
                dx,dy,du,
                ds,0.05);
      }
      else
      {
         Arrow3(builder,
                xs,ys,us,
                dx,dy,du,
                ds,0.0);
      }
   }
}

void VisualizationSceneSolution::PrepareNumbering()
{
   PrepareElementNumbering();
   PrepareEdgeNumbering();
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

   updated_bufs.emplace_back(&line_buf);
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

   updated_bufs.emplace_back(&line_buf);
}

void VisualizationSceneSolution::UpdateValueRange(bool prepare)
{
   bool had_logscale = logscale;
   logscale = logscale && LogscaleRange();
   SetLogA();
   SetLevelLines(minv, maxv, nl);
   // preserve the current box z-size
   zscale *= (bb.z[1]-bb.z[0])/(maxv-minv);
   bb.z[0] = minv;
   bb.z[1] = maxv;
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
         if (drawbdr == 2)
         {
            const double val = mesh->GetBdrAttribute(i);
            MySetColor(bl, val, minv, maxv);
            for (j = 0; j < vals.Size(); j++)
            {
               bl.glVertex3d(pointmat(0, j), pointmat(1, j), val);
            }
         }
         else
         {
            for (j = 0; j < vals.Size(); j++)
            {
               bl.glVertex3d(pointmat(0, j), pointmat(1, j), vals(j));
            }
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

   updated_bufs.emplace_back(&bdr_buf);
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
   updated_bufs.emplace_back(&cp_buf);
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

gl3::SceneInfo VisualizationSceneSolution::GetSceneObjs()
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
   params.contains_translucent = matAlpha < 1.0;
   if (drawelems)
   {
      // draw elements
      scene.queue.emplace_back(params, &disp_buf);
   }
   // draw orderings -- color modes
   params.contains_translucent = false;
   if (draworder == 1)
   {
      scene.queue.emplace_back(params, &order_noarrow_buf);
   }
   else if (draworder == 2)
   {
      scene.queue.emplace_back(params, &order_buf);
   }
   params.contains_translucent = matAlpha < 1.0;
   params.mesh_material = VisualizationScene::BLK_MAT;
   // everything below will be drawn in "black"
   params.static_color = GetLineColor();
   if (draw_cp)
   {
      // draw cutting plane
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
   if (drawnums == 1)
   {
      scene.queue.emplace_back(params, &e_nums_buf);
   }
   else if (drawnums == 2)
   {
      scene.queue.emplace_back(params, &f_nums_buf);
   }
   else if (drawnums == 3)
   {
      scene.queue.emplace_back(params, &v_nums_buf);
   }

   // draw orderings -- "black" modes
   if (draworder == 3)
   {
      scene.queue.emplace_back(params, &order_noarrow_buf);
   }
   else if (draworder == 4)
   {
      scene.queue.emplace_back(params, &order_buf);
   }

   return scene;
}

void VisualizationSceneSolution::glTF_ExportBoundary(
   glTF_Builder &bld,
   glTF_Builder::buffer_id buffer,
   glTF_Builder::material_id black_mat)
{
   auto bdr_node = AddModelNode(bld, "Boundary");
   auto bdr_mesh = bld.addMesh("Boundary Mesh");
   bld.addNodeMesh(bdr_node, bdr_mesh);

   int nlines = AddLines(
                   bld,
                   bdr_mesh,
                   buffer,
                   black_mat,
                   bdr_buf);
   if (nlines == 0)
   {
      cout << "glTF export: no boundary found to export!" << endl;
   }
}

void VisualizationSceneSolution::glTF_Export()
{
   string name = "GLVis_scene_000";

   glTF_Builder bld(name);

   auto palette_mat = AddPaletteMaterial(bld);
   auto black_mat = AddBlackMaterial(bld);
   auto buf = bld.addBuffer("buffer");
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
