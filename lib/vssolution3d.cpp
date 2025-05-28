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
#include <cmath>
#include <limits>

#include <mfem.hpp>

#include "palettes.hpp"
#include "vssolution3d.hpp"

using namespace std;
using namespace mfem;

thread_local VisualizationSceneSolution3d *vssol3d;
extern thread_local GeometryRefiner GLVisGeometryRefiner;

// Reference geometries with a cut in the middle, based on subdivision of
// GLVisGeometryRefiner in 3-4 quads. Updated when cut_lambda is updated, see
// keys Ctrl+F3/F4. We need these variables because the GLVisGeometryRefiner
// caches its RefinedGeometry objects.
thread_local IntegrationRule cut_QuadPts;
thread_local Array<int> cut_QuadGeoms;
thread_local IntegrationRule cut_TriPts;
thread_local Array<int> cut_TriGeoms;

// Definitions of some more keys

std::string VisualizationSceneSolution3d::GetHelpString() const
{
   std::stringstream os;
   os << endl
      << "+------------------------------------+" << endl
      << "| Keys                               |" << endl
      << "+------------------------------------+" << endl
      << "| a -  Displays/Hides the axes       |" << endl
      << "| A -  Turns antialiasing on/off     |" << endl
      << "| c -  Toggle colorbar and caption   |" << endl
      << "| C -  Change the main plot caption  |" << endl
      << "| e -  Displays/Hides the elements   |" << endl
      << "| E -  Toggle the elements in the CP |" << endl
      << "| f -  Smooth/Flat/discont. shading  |" << endl
      << "| g -  Toggle background             |" << endl
      << "| G -  Export to glTF format         |" << endl
      << "| h -  Displays help menu            |" << endl
      << "| i -  Toggle cutting plane          |" << endl
      << "| I -  Toggle cutting plane algorithm|" << endl
      << "| j -  Turn on/off perspective       |" << endl
      << "| k/K  Adjust the transparency level |" << endl
      << "| ,/<  Adjust color transparency     |" << endl
      << "| l -  Turns on/off the light        |" << endl
      << "| L -  Toggle logarithmic scale      |" << endl
      << "| m -  Displays/Hides the mesh       |" << endl
      << "| M -  Toggle the mesh in the CP     |" << endl
      << "| o/O  (De)refine elem, disc shading |" << endl
      << "| p/P  Cycle through color palettes  |" << endl
      << "| q -  Quits                         |" << endl
      << "| Q -  Cycle quadrature data mode    |" << endl
      << "| r -  Reset the plot to 3D view     |" << endl
      << "| R -  Reset the plot to 2D view     |" << endl
      << "| s -  Turn on/off unit cube scaling |" << endl
      << "| S -  Take snapshot/Record a movie  |" << endl
      << "| t -  Cycle materials and lights    |" << endl
      << "| u/U  Move the level surface        |" << endl
      << "| v/V  Add/Delete a level surface    |" << endl
      << "| w/W  Move bdr elements up/down     |" << endl
      << "| x/X  Rotate cutting plane (phi)    |" << endl
      << "| y/Y  Rotate cutting plane (theta)  |" << endl
      << "| z/Z  Translate cutting plane       |" << endl
      << "| \\ -  Set light source position     |" << endl
      << "| Alt+a  - Axes number format        |" << endl
      << "| Alt+c  - Colorbar number format    |" << endl
      << "| Ctrl+o - Element ordering curve    |" << endl
      << "| Ctrl+p - Print to a PDF file       |" << endl
      << "+------------------------------------+" << endl
      << "| Function keys                      |" << endl
      << "+------------------------------------+" << endl
      << "| F1 - X window info and keystrokes  |" << endl
      << "| F2 - Update colors, etc.           |" << endl
      << "| F3/F4 - Shrink/Zoom bdr elements   |" << endl
      << "| Ctrl+F3/F4 - Cut face bdr elements |" << endl
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

static void KeyiPressed()
{
   vssol3d -> ToggleCuttingPlane();
   SendExposeEvent();
}

static void KeyIPressed()
{
   vssol3d -> ToggleCPAlgorithm();
   SendExposeEvent();
}

void VisualizationSceneSolution3d::PrepareOrderingCurve()
{
   bool color = draworder < 3;
   order_buf.clear();
   order_noarrow_buf.clear();
   PrepareOrderingCurve1(order_buf, true, color);
   PrepareOrderingCurve1(order_noarrow_buf, false, color);
   updated_bufs.emplace_back(&order_buf);
   updated_bufs.emplace_back(&order_noarrow_buf);
}


void VisualizationSceneSolution3d::PrepareOrderingCurve1(gl3::GlDrawable& buf,
                                                         bool arrows,
                                                         bool color)
{
   // make the lines of the ordering curve thicker
   double ThicknessFactor = 2.0;
   double MS_Thickness = 2.0;
   double LineWidth;
   if (GetMultisample() > 0)
   {
      LineWidth = GetLineWidthMS();
      SetLineWidthMS(MS_Thickness);
   }
   else
   {
      LineWidth = GetLineWidth();
      SetLineWidth(ThicknessFactor*LineWidth);
   }

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
      double zs = 0.0;
      for (int j = 0; j < nv; j++)
      {
         xs += pointmat(0,j);
         ys += pointmat(1,j);
         zs += pointmat(2,j);
      }
      xs /= nv;
      ys /= nv;
      zs /= nv;

      double xs1 = 0.0;
      double ys1 = 0.0;
      double zs1 = 0.0;
      for (int j = 0; j < nv1; j++)
      {
         xs1 += pointmat1(0,j);
         ys1 += pointmat1(1,j);
         zs1 += pointmat1(2,j);
      }
      xs1 /= nv1;
      ys1 /= nv1;
      zs1 /= nv1;

      double dx = xs1-xs;
      double dy = ys1-ys;
      double dz = zs1-zs;
      double ds = sqrt(dx*dx+dy*dy+dz*dz);

      double cval = HUGE_VAL;
      if (color)
      {
         cval = minv+double(k)/ne*(maxv-minv);
      }

      if (arrows)
      {
         Arrow3(buf,
                xs,ys,zs,
                dx,dy,dz,
                ds,0.05,cval);
      }
      else
      {
         Arrow3(buf,
                xs,ys,zs,
                dx,dy,dz,
                ds,0.0,cval);
      }
   }

   if (GetMultisample() > 0)
   {
      SetLineWidthMS(LineWidth);
   }
   else
   {
      SetLineWidth(LineWidth);
   }
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

static void KeyoPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      vssol3d -> ToggleDrawOrdering();
      vssol3d -> PrepareOrderingCurve();
      SendExposeEvent();
   }
   else
   {
      if (vssol3d -> TimesToRefine < 32)
      {
         cout << "Subdivision factor = " << ++vssol3d->TimesToRefine << endl;
         if (vssol3d -> GetShading() ==
             VisualizationSceneScalarData::Shading::Noncomforming)
         {
            vssol3d->DoAutoscale(false);
            vssol3d -> Prepare();
            vssol3d -> PrepareLines();
            vssol3d -> CPPrepare();
            vssol3d -> PrepareLevelSurf();
            SendExposeEvent();
         }
      }
   }
}

static void KeyOPressed()
{
   if (vssol3d -> TimesToRefine > 1)
   {
      cout << "Subdivision factor = " << --vssol3d->TimesToRefine << endl;
      if (vssol3d -> GetShading() ==
          VisualizationSceneScalarData::Shading::Noncomforming)
      {
         vssol3d->DoAutoscale(false);
         vssol3d -> Prepare();
         vssol3d -> PrepareLines();
         vssol3d -> CPPrepare();
         vssol3d -> PrepareLevelSurf();
         SendExposeEvent();
      }
   }
}

static void KeywPressed()
{
   if (vssol3d -> GetShading() ==
       VisualizationSceneScalarData::Shading::Noncomforming)
   {
      if ( fabs(vssol3d -> FaceShiftScale += 0.01) < 0.001 )
      {
         vssol3d -> FaceShiftScale = 0.0;
      }
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
   if (vssol3d -> GetShading() ==
       VisualizationSceneScalarData::Shading::Noncomforming)
   {
      if ( fabs(vssol3d -> FaceShiftScale -= 0.01) < 0.001 )
      {
         vssol3d -> FaceShiftScale = 0.0;
      }
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

static void KeyF3Pressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      if (vssol3d->cut_lambda <= 0.95)
      {
         vssol3d->cut_lambda += 0.05;
         vssol3d->cut_updated = false;
      }
      if (fabs(vssol3d->cut_lambda-1.0) < 1e-12) // snap to 1
      {
         vssol3d->cut_lambda = 1.0;
      }
      vssol3d->Prepare();
      SendExposeEvent();
   }
   else
   {
      if (vssol3d->GetShading() ==
          VisualizationSceneScalarData::Shading::Noncomforming)
      {
         if (vssol3d->GetMesh()->Dimension() == 3 && vssol3d->bdrc.Width() == 0)
         {
            vssol3d->ComputeBdrAttrCenter();
         }
         if (vssol3d->GetMesh()->Dimension() == 2 && vssol3d->matc.Width() == 0)
         {
            vssol3d->ComputeElemAttrCenter();
         }
         vssol3d->shrink *= 0.9;
         if (magic_key_pressed)
         {
            vssol3d -> Scale(1.11111111111111111111111);
         }
         SendExposeEvent();
         vssol3d->Prepare();
         vssol3d->PrepareLines();
         SendExposeEvent();
      }
   }
}

static void KeyF4Pressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      if (vssol3d->cut_lambda >= 0.05)
      {
         vssol3d->cut_lambda -= 0.05;
         vssol3d->cut_updated = false;
      }
      if (fabs(vssol3d->cut_lambda-0.0) < 1e-12) // snap to 0
      {
         vssol3d->cut_lambda = 0.0;
      }
      vssol3d->Prepare();
      SendExposeEvent();
   }
   else
   {
      if (vssol3d->GetShading() ==
          VisualizationSceneScalarData::Shading::Noncomforming)
      {
         if (vssol3d->GetMesh()->Dimension() == 3 && vssol3d->bdrc.Width() == 0)
         {
            vssol3d->ComputeBdrAttrCenter();
         }
         if (vssol3d->GetMesh()->Dimension() == 2 && vssol3d->matc.Width() == 0)
         {
            vssol3d->ComputeElemAttrCenter();
         }
         vssol3d->shrink *= 1.11111111111111111111111;
         if (magic_key_pressed)
         {
            vssol3d -> Scale(0.9);
         }
         vssol3d->Prepare();
         vssol3d->PrepareLines();
         SendExposeEvent();
      }
   }
}

static void KeyF11Pressed()
{
   if (vssol3d->GetShading() ==
       VisualizationSceneScalarData::Shading::Noncomforming)
   {
      if (vssol3d->matc.Width() == 0)
      {
         vssol3d->ComputeElemAttrCenter();
      }
      vssol3d->shrinkmat *= 0.9;
      if (magic_key_pressed)
      {
         vssol3d -> Scale(1.11111111111111111111111);
      }
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

static void KeyF12Pressed()
{
   if (vssol3d->GetShading() ==
       VisualizationSceneScalarData::Shading::Noncomforming)
   {
      if (vssol3d->matc.Width() == 0)
      {
         vssol3d->ComputeElemAttrCenter();
      }
      vssol3d->shrinkmat *= 1.11111111111111111111111;
      if (magic_key_pressed)
      {
         vssol3d -> Scale(0.9);
      }
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

static void KeyF8Pressed()
{
   Mesh &mesh = *vssol3d->GetMesh();
   int dim = mesh.Dimension();
   const Array<int> &all_attr = ((dim == 3) ? mesh.bdr_attributes :
                                 mesh.attributes);
   Array<int> &attr_marker = vssol3d->bdr_attr_to_show;
   int attr;
   Array<int> attr_list(&attr, 1);

   cout << ((dim == 3) ? "Bdr a" : "A") << "ttributes ON: ";
   for (int i = 0; i < all_attr.Size(); i++)
      if (attr_marker[all_attr[i]-1])
      {
         cout << " " << all_attr[i];
      }
   cout << endl;

   cout << ((dim == 3) ? "Bdr a" : "A") << "ttribute to toggle : " << flush;
   cin >> attr;
   vssol3d->ToggleAttributes(attr_list);
   SendExposeEvent();
}

static void KeyF9Pressed()
{
   Mesh &mesh = *vssol3d->GetMesh();
   int dim = mesh.Dimension();
   const Array<int> &attr_list = ((dim == 3) ? mesh.bdr_attributes :
                                  mesh.attributes);
   Array<int> &attr_marker = vssol3d->bdr_attr_to_show;
   int attr, n, j, i;

   if (attr_list.Size() == 0)
   {
      return;
   }
   for (i = j = n = 0; i < attr_list.Size(); i++)
      if (attr_marker[attr_list[i]-1])
      {
         j = i, n++;
      }
   if (n == 1)
   {
      j = (j + 1) % (attr_list.Size() + 1);
   }
   else
   {
      j = 0;
   }
   if (j == attr_list.Size())
   {
      attr_marker = 1;
      cout << "Showing all " << ((dim == 3) ? "bdr " : "") << "attributes "
           << endl;
   }
   else
   {
      attr = attr_list[j];
      attr_marker = 0;
      attr_marker[attr-1] = 1;
      cout << "Showing " << ((dim == 3) ? "bdr " : "") << "attribute "
           << attr << endl;
   }
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

static void KeyF10Pressed()
{
   Mesh &mesh = *vssol3d->GetMesh();
   int dim = mesh.Dimension();
   const Array<int> &attr_list = ((dim == 3) ? mesh.bdr_attributes :
                                  mesh.attributes);
   Array<int> &attr_marker = vssol3d->bdr_attr_to_show;
   int attr, n, j, i;

   if (attr_list.Size() == 0)
   {
      return;
   }
   n = 0;
   for (i = j = n = 0; i < attr_list.Size(); i++)
      if (attr_marker[attr_list[i]-1])
      {
         j = i, n++;
      }
   if (n == 1)
   {
      j = (j + attr_list.Size()) % (attr_list.Size() + 1);
   }
   else
   {
      j = attr_list.Size() - 1;
   }
   if (j == attr_list.Size())
   {
      attr_marker = 1;
      cout << "Showing all " << ((dim == 3) ? "bdr " : "") << "attributes "
           << endl;
   }
   else
   {
      attr = attr_list[j];
      attr_marker = 0;
      attr_marker[attr-1] = 1;
      cout << "Showing " << ((dim == 3) ? "bdr " : "") << "attribute "
           << attr << endl;
   }
   vssol3d -> PrepareLines();
   vssol3d -> Prepare();
   SendExposeEvent();
}

VisualizationSceneSolution3d::VisualizationSceneSolution3d()
{}

VisualizationSceneSolution3d::VisualizationSceneSolution3d(Mesh &m, Vector &s,
                                                           Mesh *mc)
{
   mesh = &m;
   mesh_coarse = mc;
   sol = &s;
   GridF = NULL;

   Init();
}


void VisualizationSceneSolution3d::Init()
{
   vssol3d = this;

   cplane = 0;
   cp_drawmesh = 0; cp_drawelems = 1;
   drawlsurf = 0;
   cp_algo = 0;

   drawelems = 1;
   shading = Shading::Smooth;
   drawmesh = 0;
   draworder = 0;
   scaling = 0;

   shrink = 1.0;
   shrinkmat = 1.0;
   bdrc.SetSize(3,0);
   matc.SetSize(3,0);

   TimesToRefine = 1;
   FaceShiftScale = 0.0;

   if (mesh->Dimension() == 3)
   {
      if (!mesh->bdr_attributes.Size())
      {
         MFEM_WARNING("Missing boundary attributes!");
      }
      bdr_attr_to_show.SetSize(mesh->bdr_attributes.Size() > 0 ?
                               mesh->bdr_attributes.Max() : 0);
   }
   else
   {
      if (!mesh->attributes.Size())
      {
         MFEM_WARNING("Missing element attributes!");
      }
      bdr_attr_to_show.SetSize(mesh->attributes.Size() > 0 ?
                               mesh->attributes.Max() : 0);
   }
   bdr_attr_to_show = 1;

   VisualizationSceneScalarData::Init(); // calls FindNewBox

   FindNewValueRange(false);

   node_pos = new double[mesh->GetNV()];

   palette.SetIndex(12); // use the 'vivid' palette in 3D

   double eps = 1e-6; // move the cutting plane a bit to avoid artifacts
   CuttingPlane = new Plane(-1.0,0.0,0.0,(0.5-eps)*bb.x[0]+(0.5+eps)*bb.x[1]);

   nlevels = 1;

   FindNodePos();

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      wnd->setOnKeyDown('m', KeymPressed);
      wnd->setOnKeyDown('M', KeyMPressed);

      wnd->setOnKeyDown('e', KeyePressed);
      wnd->setOnKeyDown('E', KeyEPressed);

      wnd->setOnKeyDown('f', KeyFPressed);
      // wnd->setOnKeyDown('F', KeyFPressed);

      wnd->setOnKeyDown('i', KeyiPressed);
      wnd->setOnKeyDown('I', KeyIPressed);

      wnd->setOnKeyDown('o', KeyoPressed);
      wnd->setOnKeyDown('O', KeyOPressed);

      wnd->setOnKeyDown('w', KeywPressed);
      wnd->setOnKeyDown('W', KeyWPressed);

      wnd->setOnKeyDown('x', KeyxPressed);
      wnd->setOnKeyDown('X', KeyXPressed);

      wnd->setOnKeyDown('y', KeyyPressed);
      wnd->setOnKeyDown('Y', KeyYPressed);

      wnd->setOnKeyDown('z', KeyzPressed);
      wnd->setOnKeyDown('Z', KeyZPressed);

      wnd->setOnKeyDown('u', KeyuPressed);
      wnd->setOnKeyDown('U', KeyUPressed);

      wnd->setOnKeyDown('v', KeyvPressed);
      wnd->setOnKeyDown('V', KeyVPressed);

      wnd->setOnKeyDown(SDLK_F3, KeyF3Pressed);
      wnd->setOnKeyDown(SDLK_F4, KeyF4Pressed);
      wnd->setOnKeyDown(SDLK_NUMLOCKCLEAR, ToggleMagicKey);

      wnd->setOnKeyDown(SDLK_F8, KeyF8Pressed);
      wnd->setOnKeyDown(SDLK_F9, KeyF9Pressed);
      wnd->setOnKeyDown(SDLK_F10, KeyF10Pressed);

      wnd->setOnKeyDown(SDLK_F11, KeyF11Pressed);
      wnd->setOnKeyDown(SDLK_F12, KeyF12Pressed);
   }
   Prepare();
   PrepareLines();
   CPPrepare();
   PrepareLevelSurf();
   PrepareOrderingCurve();
}

VisualizationSceneSolution3d::~VisualizationSceneSolution3d()
{
   delete [] node_pos;
}

void VisualizationSceneSolution3d::NewMeshAndSolution(
   Mesh *new_m, Mesh *new_mc, Vector *new_sol, GridFunction *new_u)
{
   if (mesh->GetNV() != new_m->GetNV())
   {
      delete [] node_pos;
      node_pos = new double[new_m->GetNV()];
   }

   Mesh *old_m = mesh;
   mesh = new_m;
   mesh_coarse = new_mc;
   sol = new_sol;
   GridF = new_u;

   // If the number of surface elements changes, recompute the refinement factor
   if (mesh->Dimension() != old_m->Dimension() ||
       (mesh->Dimension() == 2 && mesh->GetNE() != old_m->GetNE()) ||
       (mesh->Dimension() == 3 && mesh->GetNBE() != old_m->GetNBE()))
   {
      int ref = GetAutoRefineFactor();
      if (TimesToRefine != ref)
      {
         TimesToRefine = ref;
         cout << "Subdivision factor = " << TimesToRefine << endl;
      }
   }
   FindNodePos();

   DoAutoscale(false);

   Prepare();
   PrepareLines();
   CPPrepare();
   PrepareLevelSurf();
   PrepareOrderingCurve();
}

void VisualizationSceneSolution3d::SetShading(Shading s, bool print)
{
   if (shading == s || s <= Shading::Min)
   {
      return;
   }

   if (s >= Shading::Max || (GridF == NULL && s > Shading::Smooth))
   {
      return;
   }
   Shading os = shading;
   shading = s;
   if (GridF != NULL && (s == Shading::Noncomforming ||
                         os == Shading::Noncomforming))
   {
      DoAutoscale(false);
      PrepareLines();
      CPPrepare();
   }
   Prepare();
   PrepareLevelSurf();

   static const char *shading_type[3] =
   {"flat", "smooth", "non-conforming (with subdivision)"};
   if (print)
   {
      cout << "Shading type : " << shading_type[(int)shading] << endl;
   }
}

void VisualizationSceneSolution3d::ToggleShading()
{
   if (GridF)
   {
      VisualizationSceneScalarData::ToggleShading();
   }
   else
   {
      SetShading((Shading)(1 - (int)shading), true);
   }
}

void VisualizationSceneSolution3d::SetRefineFactors(int f, int ignored)
{
   if (TimesToRefine == f || f < 1)
   {
      return;
   }

   TimesToRefine = f;

   if (shading == Shading::Noncomforming)
   {
      DoAutoscale(false);
      Prepare();
      PrepareLines();
      CPPrepare();
      PrepareOrderingCurve();
   }
}

int VisualizationSceneSolution3d::GetFunctionAutoRefineFactor()
{
   if (!GridF) { return 1; }

   return VisualizationSceneScalarData::GetFunctionAutoRefineFactor(*GridF);
}

void VisualizationSceneSolution3d::AutoRefine()
{
   int ref = GetAutoRefineFactor();

   cout << "Subdivision factor = " << ref << endl;

   SetRefineFactors(ref, 1);
}

void VisualizationSceneSolution3d::ToggleAttributes(Array<int> &attr_list)
{
   int dim = mesh->Dimension();
   Array<int> &attr_marker = bdr_attr_to_show;

   for (int i = 0; i < attr_list.Size(); i++)
   {
      int attr = attr_list[i];
      if (attr < 1)
      {
         cout << "Hiding all" << ((dim == 3) ? " bdr" : "") << " attributes."
              << endl;
         attr_marker = 0;
      }
      else if (attr > attr_marker.Size())
      {
         cout << "Showing all" << ((dim == 3) ? " bdr" : "") << " attributes."
              << endl;
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

void VisualizationSceneSolution3d::FindNewBox(bool prepare)
{
   int nv = mesh -> GetNV();

   double *coord = mesh->GetVertex(0);

   bb.x[0] = bb.x[1] = coord[0];
   bb.y[0] = bb.y[1] = coord[1];
   bb.z[0] = bb.z[1] = coord[2];

   for (int i = 1; i < nv; i++)
   {
      coord = mesh->GetVertex(i);
      if (coord[0] < bb.x[0]) { bb.x[0] = coord[0]; }
      if (coord[1] < bb.y[0]) { bb.y[0] = coord[1]; }
      if (coord[2] < bb.z[0]) { bb.z[0] = coord[2]; }
      if (coord[0] > bb.x[1]) { bb.x[1] = coord[0]; }
      if (coord[1] > bb.y[1]) { bb.y[1] = coord[1]; }
      if (coord[2] > bb.z[1]) { bb.z[1] = coord[2]; }
   }

   if (shading == Shading::Noncomforming)
   {
      int dim = mesh->Dimension();
      int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
      int fn, fo;
      DenseMatrix pointmat;
      RefinedGeometry *RefG;
      IntegrationRule eir;
      FaceElementTransformations *Tr;
      ElementTransformation *T;

      for (int i = 0; i < ne; i++)
      {
         if (dim == 3)
         {
            mesh->GetBdrElementFace(i, &fn, &fo);
            RefG = GLVisGeometryRefiner.Refine(mesh->GetFaceGeometry(fn),
                                               TimesToRefine);
            Tr = mesh->GetFaceElementTransformations(fn, 5);
            eir.SetSize(RefG->RefPts.GetNPoints());
            Tr->Loc1.Transform(RefG->RefPts, eir);
            Tr->Elem1->Transform(eir, pointmat);
         }
         else
         {
            T = mesh->GetElementTransformation(i);
            RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                               TimesToRefine);
            T->Transform(RefG->RefPts, pointmat);
         }
         for (int j = 0; j < pointmat.Width(); j++)
         {
            if (pointmat(0,j) < bb.x[0]) { bb.x[0] = pointmat(0,j); }
            if (pointmat(1,j) < bb.y[0]) { bb.y[0] = pointmat(1,j); }
            if (pointmat(2,j) < bb.z[0]) { bb.z[0] = pointmat(2,j); }
            if (pointmat(0,j) > bb.x[1]) { bb.x[1] = pointmat(0,j); }
            if (pointmat(1,j) > bb.y[1]) { bb.y[1] = pointmat(1,j); }
            if (pointmat(2,j) > bb.z[1]) { bb.z[1] = pointmat(2,j); }
         }
      }
   }

   UpdateBoundingBox();
}

void VisualizationSceneSolution3d::FindNewValueRange(bool prepare)
{
   int map_type = (GridF) ?
                  GridF->FESpace()->FEColl()->GetMapType(mesh->Dimension()) :
                  FiniteElement::VALUE;

   if (shading < Shading::Noncomforming || map_type != (int)FiniteElement::VALUE)
   {
      minv = sol->Min();
      maxv = sol->Max();
   }
   else
   {
      minv = GridF->Min();
      maxv = GridF->Max();
   }
   FixValueRange();
   UpdateValueRange(prepare);
}

void VisualizationSceneSolution3d::EventUpdateColors()
{
   Prepare();
   PrepareCuttingPlane();
   PrepareLevelSurf();
   PrepareOrderingCurve();
   if (shading == Shading::Noncomforming && drawmesh != 0 && FaceShiftScale != 0.0)
   {
      PrepareLines();
   }
}

void VisualizationSceneSolution3d::UpdateValueRange(bool prepare)
{
   logscale = logscale && LogscaleRange();
   palette.SetUseLogscale(logscale);
   SetLogA();
   SetLevelLines(minv, maxv, nl);
   if (prepare)
   {
      UpdateLevelLines();
      EventUpdateColors();
   }
}

void VisualizationSceneSolution3d::FindNodePos()
{
   int i, nnodes = mesh -> GetNV();

   for (i = 0; i < nnodes; i++)
   {
      node_pos[i] = CuttingPlane -> Transform (mesh -> GetVertex (i));
   }
}

void VisualizationSceneSolution3d::ToggleDrawMesh()
{
   drawmesh = (drawmesh+1)%3;
   PrepareLines();
}

void VisualizationSceneSolution3d::ToggleCuttingPlane()
{
   if (cplane == 2 && cp_drawmesh == 3)
   {
      cp_drawmesh = 2;
   }

   cplane = (cplane+1)%3;
#ifdef GLVIS_DEBUG
   cout << "cplane = " << cplane << endl;
#endif
   CPPrepare();
   if (cplane == 0 || cplane == 2)
   {
      Prepare();
      PrepareLines();
      PrepareOrderingCurve();
   }
}

void VisualizationSceneSolution3d::ToggleCPDrawElems()
{
   cp_drawelems = 1-cp_drawelems;
   PrepareCuttingPlane();
}

void VisualizationSceneSolution3d::ToggleCPDrawMesh()
{
   if (cplane == 1)
   {
      cp_drawmesh = (cp_drawmesh+1)%3;
   }
   else if (cplane == 2)
   {
      cp_drawmesh = (cp_drawmesh+1)%4;
   }
   PrepareCuttingPlaneLines();
}

void VisualizationSceneSolution3d::ToggleCPAlgorithm()
{
   cp_algo = (cp_algo+1)%2;
   if (shading == Shading::Noncomforming && cplane == 1)
   {
      CPPrepare();
   }
}

void VisualizationSceneSolution3d::MoveLevelSurf(int move)
{
   drawlsurf += move;
   if (drawlsurf < 0)
   {
      drawlsurf = 0;
   }
   if (drawlsurf > 49)
   {
      drawlsurf = 49;
   }
   PrepareLevelSurf();
}

void VisualizationSceneSolution3d::NumberOfLevelSurf(int c)
{
   nlevels += c;
   if (nlevels < 1)
   {
      nlevels = 1;
   }
   PrepareLevelSurf();
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
   FaceElementTransformations *Tr;
   normals.SetSize(3, ir.GetNPoints());
   ElementTransformation *ETr;
   IntegrationPointTransformation *LTr;
   if (side == 0)
   {
      Tr = mesh->GetFaceElementTransformations(FaceNo, 5);
      ETr = Tr->Elem1;
      LTr = &Tr->Loc1;
   }
   else
   {
      Tr = mesh->GetFaceElementTransformations(FaceNo, 10);
      ETr = Tr->Elem2;
      LTr = &Tr->Loc2;
   }
   LTr->Transform(ir, eir);
   for (int i = 0; i < normals.Width(); i++)
   {
      LTr->Transf.SetIntPoint(&ir.IntPoint(i));
      CalcOrtho(LTr->Transf.Jacobian(), lnor);
      ETr->SetIntPoint(&eir.IntPoint(i));
      const DenseMatrix &Jac = ETr->Jacobian();
      CalcInverse(Jac, JInv);
      normals.GetColumnReference(i, nr);
      JInv.MultTranspose(lnor, nr);
   }
   if (side)
   {
      normals *= -1.;
   }
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
         RefG = GLVisGeometryRefiner.Refine(Geometry::TRIANGLE, TimesToRefine);
         ip_transf.Transf.SetFE (&TriangleFE);
         break;
      case 4:
         RefG = GLVisGeometryRefiner.Refine(Geometry::SQUARE, TimesToRefine);
         ip_transf.Transf.SetFE (&QuadrilateralFE);
         break;
      case 5:
         DrawRefinedSurf (3, points, elem, func, 0); // draw (0,1,2)
         for (i = 0; i < 3; i++)
         {
            points[4+i] = points[i];   // move point 0 to point 1
         }
         DrawRefinedSurf (4, points+4, elem, func, 1); // draw (0,2,3,4)
         return;
      case 6:
         DrawRefinedSurf (4, points, elem, func, 0); // draw (0,1,2,3)
         for (i = 0; i < 3; i++)
         {
            points[8+i] = points[i];   // move point 0 to point 2
         }
         DrawRefinedSurf (4, points+8, elem, func, 1); // draw (0,3,4,5)
         return;
      default:
         return;
   }
   DenseMatrix &pm = ip_transf.Transf.GetPointMat();
   pm.SetSize (3, n);
   for (i = 0; i < n; i++)
      for (j = 0; j < 3; j++)
      {
         pm(j,i) = points[4*i+j];
      }
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
   {
      return;
   }

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
         {
            pts[i][j] = pointmat(j,RG[i]);
         }
      if (n > 3)
      {
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], pts[3], norm);
      }
      else
      {
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], norm);
      }
      if (j)
      {
         // could not compute normal
         cerr << "WARNING: VisualizationSceneSolution3d::LiftRefinedSurf"
              << endl;
         return;
      }
   }

   double bbox_diam = sqrt ( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                             (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                             (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   for (i = 0; i < pointmat.Width(); i++)
   {
      double val = sc * (values(i) - minv) / (maxv - minv);
      for (j = 0; j < 3; j++)
      {
         pointmat(j, i) +=  val * norm[j];
      }
   }
}

void VisualizationSceneSolution3d::DrawRefinedSurf(
   int n, DenseMatrix &pointmat, Vector &values, Array<int> &RefGeoms)
{
   double norm[3], pts[4][3];
   gl3::GlBuilder draw = cplane_buf.createBuilder();

   for (int i = 0; i < RefGeoms.Size()/n; i++)
   {
      int *RG = &(RefGeoms[i*n]);
      int j;
      for (j = 0; j < n; j++)
         for (int l = 0; l < 3; l++)
         {
            pts[j][l] = pointmat(l, RG[j]);
         }
      if (n > 3)
      {
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], pts[3], norm);
      }
      else
      {
         j = Compute3DUnitNormal (pts[0], pts[1], pts[2], norm);
      }
      if (!j)
      {
         draw.glBegin (GL_POLYGON);
         draw.glNormal3dv (norm);
         for (j = 0; j < n; j++)
         {
            MySetColor (draw, values(RG[j]), minv, maxv);
            draw.glVertex3dv (pts[j]);
         }
         draw.glEnd();
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
   gl3::GlBuilder line = cplines_buf.createBuilder();
   if (part == 0)
   {
      k_end = (k_end/n) * (n-1);
   }
   if (part == 1)
   {
      k_start = (k_end/n);
   }

   line.glBegin(GL_LINES);
   for (k = k_start; k < k_end; k++)
   {
      int RE = RefEdges[k];

      line.glVertex3d (pointmat(0, RE), pointmat(1, RE),
                       pointmat(2, RE));
   }
   line.glEnd();
}

void VisualizationSceneSolution3d::DrawBdrElCoarseSurfEdges(
   gl3::GlBuilder &line, int be, DenseMatrix &pointmat, const IntegrationRule *ir,
   Array<int> *idxs)
{
   const int dim = mesh->Dimension();
   int f, e1;
   if (dim == 3)
   {
      int o;
      mesh->GetBdrElementFace(be, &f, &o);
      e1 = -1;
   }
   else
   {
      f = -1;
      e1 = be;
   }
   DrawCoarseSurfEdges(line, f, e1, -1, pointmat, ir, idxs);
}

void VisualizationSceneSolution3d::DrawFaceCoarseSurfEdges(
   gl3::GlBuilder &line, int f, DenseMatrix &pointmat, const IntegrationRule *ir,
   Array<int> *idxs)
{
   DrawCoarseSurfEdges(line, f, -1, -1, pointmat, ir, idxs);
}

void VisualizationSceneSolution3d::DrawCoarseSurfEdges(
   gl3::GlBuilder &line, int f, int e1, int e2, DenseMatrix &pointmat,
   const IntegrationRule *ir,
   Array<int> *idxs)
{
   MFEM_ASSERT(mesh_coarse, "Cannot be used without the coarse mesh");
   MFEM_ASSERT(mesh->GetLastOperation() == Mesh::Operation::REFINE,
               "Not a refined mesh");

   const int dim = mesh->Dimension();
   FaceElementTransformations *ftr;
   if (f >= 0)
   {
      ftr = mesh->GetFaceElementTransformations(f);
      e1 = ftr->Elem1No;
      e2 = ftr->Elem2No;
   }
   auto &ref = mesh->GetRefinementTransforms();
   IsoparametricTransformation trans;
   static const BiLinear2DFiniteElement fe_face;
   static const TriLinear3DFiniteElement fe;
   trans.SetFE(&fe);
   DenseMatrix emb_pointmat;

   // we assume that mesh_course is used only for tensor finite elements,
   // like for representation of quadratures, so in 2D it is square
   const int geom = (dim == 3)?(Geometry::Type::CUBE):(Geometry::Type::SQUARE);
   const int mat = ref.embeddings[e1].matrix;
   const DenseMatrix &emb_mat = ref.point_matrices[geom](mat);
   trans.SetPointMat(emb_mat);
   if (!ir)
   {
      if (dim == 3)
      {
         IntegrationRule nodes3d(4);
         ftr->Loc1.Transform(fe_face.GetNodes(), nodes3d);
         trans.Transform(nodes3d, emb_pointmat);
      }
      else
      {
         trans.Transform(fe_face.GetNodes(), emb_pointmat);
      }
   }

   line.glBegin(GL_LINES);

   const int k_max = (idxs)?(idxs->Size()/2):(4);

   for (int k = 0; k < k_max; k++)
   {
      int j, jp1;
      Vector emb_ip1, emb_ip2;
      if (ir && idxs)
      {
         j = (*idxs)[2*k];
         jp1 = (*idxs)[2*k+1];
         if (dim == 3)
         {
            IntegrationPoint ip1_3d, ip2_3d;
            ftr->Loc1.Transform((*ir)[j], ip1_3d);
            ftr->Loc1.Transform((*ir)[jp1], ip2_3d);
            trans.Transform(ip1_3d, emb_ip1);
            trans.Transform(ip2_3d, emb_ip2);
         }
         else
         {
            trans.Transform((*ir)[j], emb_ip1);
            trans.Transform((*ir)[jp1], emb_ip2);
         }
      }
      else
      {
         j = k;
         jp1 = (k+1) % 4;
         emb_pointmat.GetColumnReference(j, emb_ip1);
         emb_pointmat.GetColumnReference(jp1, emb_ip2);
      }

      // check if we are on the outer edge
      int inter = 0;
      for (int d = 0; d < 3; d++)
         if ((emb_ip1(d) != 0. && emb_ip1(d) != 1.)
             || (emb_ip2(d) != 0. && emb_ip2(d) != 1.))
         { inter++; }
      if (e2 >= 0 && ref.embeddings[e1].parent == ref.embeddings[e2].parent)
      { inter--; }
      if (inter != 1) { continue; }

      line.glVertex3dv(&pointmat(0, j));
      line.glVertex3dv(&pointmat(0, jp1));
   }

   line.glEnd();
}

void VisualizationSceneSolution3d::DrawRefinedSurfLevelLines(
   int n, DenseMatrix &pointmat, Vector &values, Array<int> &RefGeoms)
{
   int j, k;
   int *RG;
   double point[4][4];

   gl3::GlBuilder line = cplines_buf.createBuilder();

   for (k = 0; k < RefGeoms.Size()/n; k++)
   {
      RG = &(RefGeoms[k*n]);

      for (j = 0; j < n; j++)
      {
         for (int i = 0; i < 3; i++)
         {
            point[j][i] = pointmat(i, RG[j]);
         }
         point[j][3] = values(RG[j]);
      }
      DrawPolygonLevelLines(line, point[0], n, level, false);
   }
}

void VisualizationSceneSolution3d::PrepareFlat()
{
   int i, j;
   disp_buf.clear();
   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double p[4][3], c[4];

   for (i = 0; i < ne; i++)
   {
      if (dim == 3)
      {
         if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }

         if (cplane == 2)
         {
            // for cplane == 2, get vertices of the volume element, not bdr
            int f, o, e1, e2;
            mesh->GetBdrElementFace(i, &f, &o);
            mesh->GetFaceElements(f, &e1, &e2);
            mesh->GetElementVertices(e1, vertices);
         }
         else
         {
            mesh->GetBdrElementVertices(i, vertices);
         }
      }
      else
      {
         if (!bdr_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

         mesh->GetElementVertices(i, vertices);
      }

      if (cplane == 2 && CheckPositions(vertices)) { continue; }

      if (dim == 3)
      {
         mesh->GetBdrPointMatrix(i, pointmat);
      }
      else
      {
         mesh->GetPointMatrix(i, pointmat);
      }

      for (j = 0; j < pointmat.Width(); j++)
      {
         p[j][0] = pointmat(0, j);
         p[j][1] = pointmat(1, j);
         p[j][2] = pointmat(2, j);
         c[j] = (*sol)(vertices[j]);
      }
      if (j == 3)
      {
         if (cut_lambda > 0)
         {
            DrawCutTriangle(disp_buf, p, c, minv, maxv);
         }
         else
         {
            DrawTriangle(disp_buf, p, c, minv, maxv);
         }
      }
      else if (j == 4)
      {
         if (cut_lambda > 0)
         {
            DrawCutQuad(disp_buf, p, c, minv, maxv);
         }
         else
         {
            DrawQuad(disp_buf, p, c, minv, maxv);
         }
      }
      else if (j == 2)
      {
         DrawLine(disp_buf, p, c, minv, maxv);
      }
      else
      {
         mfem_error("VisualizationSceneSolution3d::PrepareFlat() :Unknown geometry.");
      }
   }

   updated_bufs.emplace_back(&disp_buf);
}


// Cut the reference square by subdividing it into 4 trapezoids with a central
// square removed: (fl,fl)-(fr,fl)-(fr,fr)-(fl,fr). The input RefG corresponds
// to the reference square. The value of lambda controls the cut: 0 = no cut, 1
// = full cut. See keys Ctrl+F3/F4.
static void CutReferenceSquare(RefinedGeometry *RefG, double lambda,
                               IntegrationRule &RefPts, Array<int> &RefGeoms)
{
   // lambda * vertex + (1-lambda) * center
   double fl = (1.0-lambda)/2.0; // left corner of the cut frame
   double fr = (1.0+lambda)/2.0; // right corner of the cut frame

   int np = RefG->RefPts.Size();
   RefPts.SetSize(4*np);
   for (int i = 0; i < np; i++)
   {
      double X = RefG->RefPts[i].x;
      double Y = RefG->RefPts[i].y;

      // First order unit square basis functions
      double phi3 = (1.0-X)*Y,       phi2 = X*Y;       // (0,1)-(1,1)
      double phi0 = (1.0-X)*(1.0-Y), phi1 = X*(1.0-Y); // (0,0)-(1,0)

      // bottom trapezoid: (0,0)-(1,1)-(fr,fl)-(fl,fl)
      RefPts[i].x      = phi1 + fr * phi2 + fl * phi3;
      RefPts[i].y      =        fl * phi2 + fl * phi3;

      // right trapezoid: (fr,fl)-(1,0)-(1,1)-(fr,fr)
      RefPts[i+np].x   = fr * phi0 + phi1 + phi2 + fr * phi3;
      RefPts[i+np].y   = fl * phi0 +        phi2 + fr * phi3;

      // top trapezoid: (fl,fr)-(fr,fr)-(1,1)-(0,1)
      RefPts[i+2*np].x = fl * phi0 + fr * phi1 + phi2;
      RefPts[i+2*np].y = fr * phi0 + fr * phi1 + phi2 + phi3;

      // left trapezoid: (0,0)-(fl,fl)-(fl,fr)-(0,1)
      RefPts[i+3*np].x = fl * phi1 + fl * phi2;
      RefPts[i+3*np].y = fl * phi1 + fr * phi2 + phi3;

      RefPts[i].z      = RefG->RefPts[i].z;
      RefPts[i+np].z   = RefG->RefPts[i].z;
      RefPts[i+2*np].z = RefG->RefPts[i].z;
      RefPts[i+3*np].z = RefG->RefPts[i].z;
   }

   int ne = RefG->RefGeoms.Size();
   RefGeoms.SetSize(4*ne);
   for (int i = 0; i < ne; i++)
   {
      RefGeoms[i]      = RefG->RefGeoms[i];
      RefGeoms[i+ne]   = RefG->RefGeoms[i] + np;
      RefGeoms[i+2*ne] = RefG->RefGeoms[i] + 2*np;
      RefGeoms[i+3*ne] = RefG->RefGeoms[i] + 3*np;
   }
}

// Cut the reference triangle by subdividing it into 3 trapezoids with a central
// triangle removed: (fl,fl)-(fr,fl)-(fl,fr). Note that the input RefG
// corresponds to a reference square, not reference triangle. The value of
// lambda controls the cut: 0 = no cut, 1 = full cut. See keys Ctrl+F3/F4.
static void CutReferenceTriangle(RefinedGeometry *RefG, double lambda,
                                 IntegrationRule &RefPts, Array<int> &RefGeoms)
{
   // lambda * vertex + (1-lambda) * center
   double fl = (1.0-lambda)/3.0;     // left corner of the cut frame
   double fr = (1.0+2.0*lambda)/3.0; // right corner of the cut frame

   int np = RefG->RefPts.Size();
   RefPts.SetSize(3*np);
   for (int i = 0; i < np; i++)
   {
      double X = RefG->RefPts[i].x;
      double Y = RefG->RefPts[i].y;

      // First order unit square basis functions
      double phi3 = (1.0-X)*Y,       phi2 = X*Y;       // (0,1)-(1,1)
      double phi0 = (1.0-X)*(1.0-Y), phi1 = X*(1.0-Y); // (0,0)-(1,0)

      // bottom trapezoid: (0,0)-(1,0)-(fr,fl)-(fl,fl)
      RefPts[i].x      = phi1 + fr * phi2 + fl * phi3;
      RefPts[i].y      =        fl * phi2 + fl * phi3;

      // diagonal trapezoid: (fr,fl)-(1,0)-(0,1)-(fl,fr)
      RefPts[i+np].x   = fr * phi0 + phi1        + fl * phi3;
      RefPts[i+np].y   = fl * phi0        + phi2 + fr * phi3;

      // left trapezoid: (0,0)-(fl,fl)-(fl,fr)-(0,1)
      RefPts[i+2*np].x = fl * phi1 + fl * phi2;
      RefPts[i+2*np].y = fl * phi1 + fr * phi2 + phi3;

      RefPts[i].z      = RefG->RefPts[i].z;
      RefPts[i+np].z   = RefG->RefPts[i].z;
      RefPts[i+2*np].z = RefG->RefPts[i].z;
   }

   int ne = RefG->RefGeoms.Size();
   RefGeoms.SetSize(3*ne);
   for (int i = 0; i < ne; i++)
   {
      RefGeoms[i]      = RefG->RefGeoms[i];
      RefGeoms[i+ne]   = RefG->RefGeoms[i] + np;
      RefGeoms[i+2*ne] = RefG->RefGeoms[i] + 2*np;
   }
}

// Call CutReferenceTriangle and CutReferenceSquare to update the global
// variables cut_TriPts, cut_TriGeoms, cut_QuadPts, cut_QuadGeoms.
void CutReferenceElements(int TimesToRefine, double lambda)
{
   RefinedGeometry *RefG =
      GLVisGeometryRefiner.Refine(Geometry::SQUARE, TimesToRefine);
   CutReferenceTriangle(RefG, lambda, cut_TriPts, cut_TriGeoms);
   CutReferenceSquare(RefG, lambda, cut_QuadPts, cut_QuadGeoms);
}

void VisualizationSceneSolution3d::PrepareFlat2()
{
   const int dim = mesh->Dimension();

   int fn, fo, di, have_normals = 0;
   double bbox_diam, vmin, vmax;
   disp_buf.clear();

   int nbe = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat, normals;
   Vector values, normal;
   RefinedGeometry * RefG;
   Array<int> vertices;
   double norm[3];
   IsoparametricTransformation T;

   bbox_diam = sqrt ( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                      (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                      (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   vmin = numeric_limits<double>::infinity();
   vmax = -vmin;
   for (int i = 0; i < nbe; i++)
   {
      int sides;
      switch ((dim == 3) ? mesh->GetBdrElementType(i) : mesh->GetElementType(i))
      {
         case Element::TRIANGLE:
            sides = 3;
            break;
         // TODO: can we improve this? It is quite a hack because sides
         case Element::SEGMENT:
            sides = 0;
            break;

         case Element::QUADRILATERAL:
         default:
            sides = 4;
            break;
      }

      if (dim == 3)
      {
         if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }

         if (cplane == 2)
         {
            // for cplane == 2, get vertices of the volume element, not bdr
            int f, o, e1, e2;
            mesh->GetBdrElementFace(i, &f, &o);
            mesh->GetFaceElements(f, &e1, &e2);
            mesh->GetElementVertices(e1, vertices);
         }
         else
         {
            mesh->GetBdrElementVertices(i, vertices);
         }
      }
      else
      {
         if (!bdr_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }
         mesh->GetElementVertices(i, vertices);
      }

      if (cplane == 2 && CheckPositions(vertices)) { continue; }

      if (dim == 3)
      {
         mesh -> GetBdrElementFace (i, &fn, &fo);
         RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceGeometry (fn),
                                            TimesToRefine);
         if (!cut_updated)
         {
            // Update the cut version of the reference geometries
            CutReferenceElements(TimesToRefine, cut_lambda);
            cut_updated = true;
         }

         // di = GridF -> GetFaceValues (fn, 2, RefG->RefPts, values, pointmat);
         // this assumes the interior boundary faces are properly oriented ...
         di = fo % 2;
         if (di == 1 && !mesh->FaceIsInterior(fn))
         {
            di = 0;
         }

         IntegrationRule &RefPts = (cut_lambda > 0) ?
                                   ((sides == 3) ? cut_TriPts : cut_QuadPts) :
                                   RefG->RefPts;
         GridF -> GetFaceValues (fn, di, RefPts, values, pointmat);
         GetFaceNormals(fn, di, RefPts,normals);
         have_normals = 1;
         ShrinkPoints(pointmat, i, fn, di);
      }
      else
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine);
         if (!cut_updated)
         {
            // Update the cut version of the reference geometries
            CutReferenceElements(TimesToRefine, cut_lambda);
            cut_updated = true;
         }
         const IntegrationRule &ir = (cut_lambda > 0 && dim > 1) ?
                                     ((sides == 3) ? cut_TriPts : cut_QuadPts) :
                                     RefG->RefPts;
         GridF->GetValues(i, ir, values, pointmat);
         mesh->GetElementTransformation(i, &T);
         di = 0;
         ShrinkPoints(pointmat, i, 0, 0);

         // Compute normals. Skip in 1D.
         if (dim > 1)
         {
            normals.SetSize(3, values.Size());
            for (int j = 0; j < values.Size(); j++)
            {
               T.SetIntPoint(&ir.IntPoint(j));
               const DenseMatrix &J = T.Jacobian();
               normals.GetColumnReference(j, normal);
               CalcOrtho(J, normal);
               normal /= normal.Norml2();
            }
            have_normals = 1;
         }
      }

      vmin = fmin(vmin, values.Min());
      vmax = fmax(vmax, values.Max());

      // compute an average normal direction for the current face
      if (sc != 0.0 && have_normals)
      {
         for (int j = 0; j < 3; j++)
         {
            norm[j] = 0.0;
         }
         Normalize(normals);
         for (int k = 0; k < normals.Width(); k++)
            for (int j = 0; j < 3; j++)
            {
               norm[j] += normals(j, k);
            }
         Normalize(norm);
         for (int k = 0; k < pointmat.Width(); k++)
         {
            double val = sc * (values(k) - minv) / (maxv - minv);
            for (int j = 0; j < 3; j++)
            {
               pointmat(j, k) += val * norm[j];
            }
         }
         have_normals = 0;
      }

      have_normals = have_normals && (dim > 1) ? 2 : 0;
      if (di)
      {
         have_normals = -1 - have_normals;
      }
      // Comment the above lines and use the below version in order to remove
      // the 3D dark artifacts (indicating wrong boundary element orientation)
      // have_normals = have_normals ? 1 : 0;

      Array<int> &RefGeoms = (cut_lambda > 0) ?
                             ((sides == 3) ? cut_TriGeoms : cut_QuadGeoms) :
                             RefG->RefGeoms;
      int psides = (cut_lambda > 0) ? 4 : sides;
      if (dim == 1) { psides = 2; } // Hack to trigger line rendering.
      DrawPatch(disp_buf, pointmat, values, normals, psides, RefGeoms,
                minv, maxv, have_normals);
   }
   updated_bufs.emplace_back(&disp_buf);

   cout << "VisualizationSceneSolution3d::PrepareFlat2() : [min,max] = ["
        << vmin << "," << vmax << "]" << endl;
}

void VisualizationSceneSolution3d::Prepare()
{
   int j;

   if (!drawelems)
   {
      return;
   }

   switch (shading)
   {
      case Shading::Flat:
         PrepareFlat();
         return;
      case Shading::Noncomforming:
         PrepareFlat2();
         return;
      default:
         break;
   }

   disp_buf.clear();
   gl3::GlBuilder poly = disp_buf.createBuilder();

   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
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
      {
         be_to_ba.AddAColumnInRow(i);
      }
      be_to_ba.MakeJ();
      if (dim == 3)
      {
         for (int i = 0; i < ne; i++)
         {
            be_to_ba.AddConnection(i, mesh->GetBdrAttribute(i)-1);
         }
      }
      else
      {
         for (int i = 0; i < ne; i++)
         {
            be_to_ba.AddConnection(i, mesh->GetAttribute(i)-1);
         }
      }
      be_to_ba.ShiftUpI();

      Transpose(be_to_ba, ba_to_be);
   }

   const Array<int> &attributes =
      ((dim == 3) ? mesh->bdr_attributes : mesh->attributes);
   for (int d = 0; d < attributes.Size(); d++)
   {
      const int attr = attributes[d]-1;

      if (!bdr_attr_to_show[attr]) { continue; }

      const int nelem = ba_to_be.RowSize(attr);
      const int *elem = ba_to_be.GetRow(attr);

      for (int i = 0; i < nelem; i++)
      {
         if (dim == 3)
         {
            mesh->GetBdrElementVertices(elem[i], vertices);
         }
         else
         {
            mesh->GetElementVertices(elem[i], vertices);
         }
         for (j = 0; j < vertices.Size() && (dim > 1); j++)
         {
            nx(vertices[j]) = ny(vertices[j]) = nz(vertices[j]) = 0.;
         }
      }

      // Compute normals. Skip in 1D.
      for (int i = 0; i < nelem && (dim > 1); i++)
      {
         if (dim == 3)
         {
            mesh->GetBdrPointMatrix(elem[i], pointmat);
            mesh->GetBdrElementVertices(elem[i], vertices);
         }
         else
         {
            mesh->GetPointMatrix(elem[i], pointmat);
            mesh->GetElementVertices(elem[i], vertices);
         }

         if (pointmat.Width() == 3)
         {
            j = Compute3DUnitNormal(&pointmat(0,0), &pointmat(0,1),
                                    &pointmat(0,2), nor);
         }
         else
         {
            j = Compute3DUnitNormal(&pointmat(0,0), &pointmat(0,1),
                                    &pointmat(0,2), &pointmat(0,3), nor);
         }
         if (j == 0)
         {
            for (j = 0; j < pointmat.Size(); j++)
            {
               nx(vertices[j]) += nor[0];
               ny(vertices[j]) += nor[1];
               nz(vertices[j]) += nor[2];
            }
         }
      }

      for (int i = 0; i < nelem; i++)
      {
         if (dim == 3)
         {
            if (cplane == 2)
            {
               // for cplane == 2, get vertices of the volume element, not bdr
               int f, o, e1, e2;
               mesh->GetBdrElementFace(i, &f, &o);
               mesh->GetFaceElements(f, &e1, &e2);
               mesh->GetElementVertices(e1, vertices);

               if (CheckPositions(vertices)) { continue; }
            }
            else
            {
               mesh->GetBdrElementVertices(elem[i], vertices);
            }
         }
         else
         {
            mesh->GetElementVertices(elem[i], vertices);
         }

         GLenum elemType = GL_NONE;
         switch ((dim == 3) ? mesh->GetBdrElementType(elem[i]) :
                 mesh->GetElementType(elem[i]))
         {
            case Element::TRIANGLE:
               elemType = GL_TRIANGLES;
               break;
            case Element::QUADRILATERAL:
               elemType = GL_QUADS;
               break;
            case Element::SEGMENT:
               elemType = GL_LINES;
               break;
            default:
               MFEM_ABORT("Invalid boundary element type");
               break;
         }

         poly.glBegin(elemType);
         if (dim == 3)
         {
            mesh->GetBdrPointMatrix(elem[i], pointmat);
         }
         else
         {
            mesh->GetPointMatrix(elem[i], pointmat);
         }

         for (j = 0; j < pointmat.Size(); j++)
         {
            MySetColor(poly, (*sol)(vertices[j]), minv, maxv);
            if (dim > 1) { poly.glNormal3d(nx(vertices[j]), ny(vertices[j]), nz(vertices[j])); }
            poly.glVertex3dv(&pointmat(0, j));
         }
         poly.glEnd();
      }
   }

   updated_bufs.emplace_back(&disp_buf);
}

void VisualizationSceneSolution3d::PrepareLines()
{
   if (!drawmesh)
   {
      return;
   }

   if (shading == Shading::Noncomforming)
   {
      PrepareLines2();
      return;
   }

   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat;

   line_buf.clear();

   Array<int> vertices;

   for (int i = 0; i < ne; i++)
   {
      if (dim == 3)
      {
         if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }

         if (cplane == 2)
         {
            // for cplane == 2, get vertices of the volume element, not bdr
            int f, o, e1, e2;
            mesh->GetBdrElementFace(i, &f, &o);
            mesh->GetFaceElements(f, &e1, &e2);
            mesh->GetElementVertices(e1, vertices);
         }
         else
         {
            mesh->GetBdrElementVertices(i, vertices);
         }
      }
      else
      {
         if (!bdr_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

         mesh->GetElementVertices(i, vertices);
      }

      if (cplane == 2 && CheckPositions(vertices)) { continue; }

      double point[4][4];
      if (dim == 3)
      {
         mesh->GetBdrPointMatrix(i, pointmat);
      }
      else
      {
         mesh->GetPointMatrix(i, pointmat);
      }

      gl3::GlBuilder line = line_buf.createBuilder();
      switch (drawmesh)
      {
         case 1:
         {
            if (mesh_coarse)
            {
               DrawBdrElCoarseSurfEdges(line, i, pointmat);
            }
            else
            {
               line.glBegin(GL_LINE_LOOP);

               for (int j = 0; j < pointmat.Size(); j++)
               {
                  line.glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j));
               }
               line.glEnd();
            }
            break;
         }
         case 2:
            for (int j = 0; j < pointmat.Size(); j++)
            {
               for (int k = 0; k < 3; k++)
               {
                  point[j][k] = pointmat(k,j);
               }
               point[j][3] = (*sol)(vertices[j]);
            }
            DrawPolygonLevelLines(line, point[0], pointmat.Size(), level, false);
            break;
      }
   }
   updated_bufs.emplace_back(&line_buf);
}

void VisualizationSceneSolution3d::PrepareLines2()
{
   int fn, fo, di = 0;
   double bbox_diam;

   line_buf.clear();

   int dim = mesh->Dimension();
   int nbe = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat, normals;
   Vector values, normal;
   RefinedGeometry * RefG;
   Array<int> vertices;
   IsoparametricTransformation T;

   bbox_diam = sqrt ( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                      (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                      (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   for (int i = 0; i < nbe; i++)
   {
      if (dim == 3)
      {
         if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }

         if (cplane == 2)
         {
            // for cplane == 2, get vertices of the volume element, not bdr
            int f, o, e1, e2;
            mesh->GetBdrElementFace(i, &f, &o);
            mesh->GetFaceElements(f, &e1, &e2);
            mesh->GetElementVertices(e1, vertices);
         }
         else
         {
            mesh->GetBdrElementVertices(i, vertices);
         }
      }
      else
      {
         if (!bdr_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

         mesh->GetElementVertices(i, vertices);
      }

      if (cplane == 2 && CheckPositions(vertices)) { continue; }

      if (dim == 3)
      {
         mesh -> GetBdrElementFace (i, &fn, &fo);
         RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceGeometry (fn),
                                            TimesToRefine);
         // di = GridF -> GetFaceValues (fn, 2, RefG->RefPts, values, pointmat);
         di = fo % 2;
         if (di == 1 && !mesh->FaceIsInterior(fn))
         {
            di = 0;
         }
         GridF -> GetFaceValues (fn, di, RefG->RefPts, values, pointmat);
         ShrinkPoints(pointmat, i, fn, di);
      }
      else
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine);
         GridF->GetValues(i, RefG->RefPts, values, pointmat);
         ShrinkPoints(pointmat, i, 0, 0);
      }

      if (sc != 0.0)
      {
         if (dim == 3)
         {
            GetFaceNormals(fn, di, RefG->RefPts, normals);
         }
         else
         {
            normals.SetSize(3, values.Size());
            mesh->GetElementTransformation(i, &T);
            for (int j = 0; j < values.Size(); j++)
            {
               T.SetIntPoint(&RefG->RefPts.IntPoint(j));
               const DenseMatrix &J = T.Jacobian();
               normals.GetColumnReference(j, normal);
               CalcOrtho(J, normal);
               normal /= normal.Norml2();
            }
         }
         double norm[3];
         for (int j = 0; j < 3; j++)
         {
            norm[j] = 0.0;
         }
         for (int k = 0; k < normals.Width(); k++)
         {
            for (int j = 0; j < 3; j++)
            {
               norm[j] += normals(j, k);
            }
         }
         double len = sqrt(InnerProd(norm, norm));
         if (len > 0.0)
         {
            len = 1.0 / len;
         }
         for (int j = 0; j < 3; j++)
         {
            norm[j] *= len;
         }
         for (int k = 0; k < pointmat.Width(); k++)
         {
            double val = sc * (values(k) - minv) / (maxv - minv);
            for (int j = 0; j < 3; j++)
            {
               pointmat(j, k) += val * norm[j];
            }
         }
      }
      gl3::GlBuilder line = line_buf.createBuilder();
      if (drawmesh == 1)
      {
         Array<int> &REdges = RefG->RefEdges;

         line.glBegin(GL_LINES);

         if (mesh_coarse)
         {
            DrawBdrElCoarseSurfEdges(line, i, pointmat, &RefG->RefPts, &REdges);
         }
         else
         {
            for (int k = 0; k < REdges.Size(); k++)
            {
               line.glVertex3dv(&pointmat(0, REdges[k]));
            }
         }
         line.glEnd();
      }
      else if (drawmesh == 2)
      {
         double point[4][4];
         int sides;
         switch ((dim == 3) ? mesh->GetBdrElementType(i) :
                 mesh->GetElementType(i))
         {
            case Element::TRIANGLE:
               sides = 3;
               break;

            case Element::QUADRILATERAL:
            default:
               sides = 4;
               break;
         }
         for (int k = 0; k < RefG->RefGeoms.Size()/sides; k++)
         {
            int *RG = &(RefG->RefGeoms[k*sides]);

            for (int j = 0; j < sides; j++)
            {
               for (int ii = 0; ii < 3; ii++)
               {
                  point[j][ii] = pointmat(ii, RG[j]);
               }
               point[j][3] = values(RG[j]);
            }
            DrawPolygonLevelLines(line, point[0], sides, level, false);
         }
      }
   }
   updated_bufs.emplace_back(&line_buf);
}

static void CutElement(const Geometry::Type geom, const int *vert_flags,
                       const int **edge_vert_ptr, int *cut_edges,
                       int *num_cut_edges, int *n2_cut_edges)
{
   static const int tet_edges[12]= {0,3, 0,2, 0,1, 1,2, 1,3, 2,3};
   static const int pyr_edges[16] =
   {0,1, 1,2, 3,2, 0,3, 0,4, 1,4, 2,4, 3,4};
   static const int pyr_cutting[16][3] =
   {
      { 3,  4,  6}, { 9, -1, 10}, { 4,  6,  1}, {11, -1, 12},
      {13, -1, 14}, { 6,  1,  3}, {15, -1,  8}, { 1,  3,  4},
      {-1, 10,  0}, { 7, 15, -1}, {-1, 12,  2}, { 0,  9, -1},
      {-1, 14,  5}, { 2, 11, -1}, {-1,  8,  7}, { 5, 13, -1}
   };
   static const int pri_edges[18] =
   {0,1, 1,2, 2,0, 3,4, 4,5, 5,3, 0,3, 1,4, 2,5};
   static const int pri_cutting[18][3] =
   {
      { 3, -1,  5}, {13,  7, 14}, { 5, -1,  1}, {15,  9, 16},
      { 1, -1,  3}, {17, 11, 12},
      {14,  0, 13}, {10, -1,  8}, {16,  2, 15}, { 6, -1, 10},
      {12,  4, 17}, { 8, -1,  6},
      { 7, 14,  0}, { 4, 17, 11}, { 9, 16,  2}, { 0, 13,  7},
      {11, 12,  4}, { 2, 15,  9}
   };
   static const int hex_edges[24] =
   { 0,1, 1,2, 3,2, 0,3, 4,5, 5,6, 7,6, 4,7, 0,4, 1,5, 2,6, 3,7 };
   static const int hex_cutting[24][3] =
   {
      { 3,  4,  6}, {17,  9, 18}, { 4,  6,  1}, {19, 11, 20},
      {21, 12, 22}, { 6,  1,  3}, {23, 14, 16}, { 1,  3,  4},
      {18,  0, 17}, {15, 13, 10}, {20,  2, 19}, { 8, 15, 13},
      {10,  8, 15}, {22,  5, 21}, {13, 10,  8}, {16,  7, 23},
      { 9, 18,  0}, { 7, 23, 14}, {11, 20,  2}, { 0, 17,  9},
      {12, 22,  5}, { 2, 19, 11}, {14, 16,  7}, { 5, 21, 12}
   };
   const int *ev;
   int n, n2;

   n = n2 = 0;
   switch (geom)
   {
      case Geometry::TETRAHEDRON:
      {
         ev = tet_edges;
         for (int j = 0; j < 6; j++, ev += 2)
            if (vert_flags[ev[0]] != vert_flags[ev[1]])
            {
               cut_edges[n++] = j;
            }
         ev = tet_edges;
      }
      break;

      case Geometry::PYRAMID:
      {
         int emark[8];
         ev = pyr_edges;
         for (int j = 0; j < 8; j++, ev += 2)
         {
            emark[j] = vert_flags[ev[1]] - vert_flags[ev[0]];
         }
         do
         {
            int j;
            for (j = 0; j < 8; j++)
            {
               if (emark[j]) { break; }
            }
            if (j == 8)
            {
               break;
            }
            int k = 2 * j;
            if (emark[j] > 0)
            {
               k++;
            }
            do
            {
               int m;
               for (j = 0; j < 3; j++)
               {
                  m = pyr_cutting[k][j];
                  if (m >= 0)
                  {
                     ev = pyr_edges + 2 * (m / 2);
                     if ((m%2 == 0 && vert_flags[ev[0]] > vert_flags[ev[1]]) ||
                         (m%2 == 1 && vert_flags[ev[1]] > vert_flags[ev[0]]))
                     {
                        break;
                     }
                  }
               }
               cut_edges[n2++] = k/2;
               emark[k/2] = 0;
               k = m;
            }
            while (k/2 != cut_edges[n]);
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
         ev = pyr_edges;
      }
      break;

      case Geometry::PRISM:
      {
         int emark[9];
         ev = pri_edges;
         for (int j = 0; j < 9; j++, ev += 2)
         {
            emark[j] = vert_flags[ev[1]] - vert_flags[ev[0]];
         }
         do
         {
            int j;
            for (j = 0; j < 9; j++)
            {
               if (emark[j]) { break; }
            }
            if (j == 9)
            {
               break;
            }
            int k = 2 * j;
            if (emark[j] > 0)
            {
               k++;
            }
            do
            {
               int m;
               for (j = 0; j < 3; j++)
               {
                  m = pri_cutting[k][j];
                  if (m >= 0)
                  {
                     ev = pri_edges + 2 * (m / 2);
                     if ((m%2 == 0 && vert_flags[ev[0]] > vert_flags[ev[1]]) ||
                         (m%2 == 1 && vert_flags[ev[1]] > vert_flags[ev[0]]))
                     {
                        break;
                     }
                  }
               }
               cut_edges[n2++] = k/2;
               emark[k/2] = 0;
               k = m;
            }
            while (k/2 != cut_edges[n]);
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
         ev = pri_edges;
      }
      break;

      case Geometry::CUBE:
      {
         int emark[12];
         ev = hex_edges;
         for (int j = 0; j < 12; j++, ev += 2)
         {
            emark[j] = vert_flags[ev[1]] - vert_flags[ev[0]];
         }
         do
         {
            int j;
            for (j = 0; j < 12; j++)
            {
               if (emark[j]) { break; }
            }
            if (j == 12)
            {
               break;
            }
            int k = 2 * j;
            if (emark[j] > 0)
            {
               k++;
            }
            do
            {
               int m;
               for (j = 0; j < 3; j++)
               {
                  m = hex_cutting[k][j];
                  ev = hex_edges + 2 * (m / 2);
                  if ((m%2 == 0 && vert_flags[ev[0]] > vert_flags[ev[1]]) ||
                      (m%2 == 1 && vert_flags[ev[1]] > vert_flags[ev[0]]))
                  {
                     break;
                  }
               }
               cut_edges[n2++] = k/2;
               emark[k/2] = 0;
               k = m;
            }
            while (k/2 != cut_edges[n]);
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
      }
      break;

      default:
         ev = NULL; // suppress a warning
         break;
   }

   *edge_vert_ptr = ev;
   *num_cut_edges = n;
   *n2_cut_edges = n2;
}

void VisualizationSceneSolution3d::CuttingPlaneFunc(int func)
{
   int m, n, n2;
   int flag[8], cut_edges[6];
   const int *ev;
   double t, point[6][4], norm[3];

   DenseMatrix pointmat;

   Array<int> nodes;
   for (int i = 0; i < mesh -> GetNE(); i++)
   {
      n = n2 = 0; // n will be the number of intersection points
      mesh -> GetElementVertices(i, nodes);
      for (int j = 0; j < nodes.Size(); j++)
      {
         if (node_pos[nodes[j]] >= 0.0)
         {
            flag[j] = 1;
         }
         else
         {
            flag[j] = -1;
         }
      }

      CutElement(mesh->GetElementBaseGeometry(i), flag,
                 &ev, cut_edges, &n, &n2);

      while (n > 2)
      {
         if (shading != Shading::Noncomforming)
         {
            mesh -> GetPointMatrix (i, pointmat);
         }
         else
         {
            const IntegrationRule *ir;
            ir = Geometries.GetVertices (mesh -> GetElementBaseGeometry(i));
            pointmat.SetSize (3, ir -> GetNPoints());
            for (int j = 0; j < ir -> GetNPoints(); j++)
            {
               const IntegrationPoint &ip = ir -> IntPoint (j);
               pointmat(0,j) = ip.x;
               pointmat(1,j) = ip.y;
               pointmat(2,j) = ip.z;
            }
         }
         for (int j = 0; j < n; j++)
         {
            const int *en = ev + 2*cut_edges[j];
            t = node_pos[ nodes[en[1]] ];
            t = t / ( t - node_pos[ nodes[en[0]] ] );
            for (int k = 0; k < 3; k++)
            {
               point[j][k] = t*pointmat(k,en[0]) + (1-t)*pointmat(k,en[1]);
            }
            point[j][3] = t*(*sol)(nodes[en[0]]) + (1-t)*(*sol)(nodes[en[1]]);
         }

         switch (func)
         {
            case 1:  // PrepareCuttingPlane()
            {
               if (shading == Shading::Noncomforming)
               {
                  // changes point for n > 4
                  DrawRefinedSurf(n, point[0], i, 1);
               }
               else
               {
                  m = n;
                  int no_norm;
                  while (1)
                  {
                     if (m > 3)
                     {
                        no_norm = Compute3DUnitNormal(point[0], point[1], point[2],
                                                      point[3], norm);
                        if (no_norm && m > 4)
                        {
                           for (int j = 3; j < m; j++)
                              for (int k = 0; k < 4; k++)
                              {
                                 point[j-2][k] = point[j][k];
                              }
                           m -= 2;
                           continue;
                        }
                     }
                     else
                        no_norm = Compute3DUnitNormal(point[0], point[1], point[2],
                                                      norm);
                     break;
                  }

                  gl3::GlBuilder draw = cplane_buf.createBuilder();
                  if (!no_norm)
                  {
                     draw.glBegin(GL_POLYGON);
                     draw.glNormal3dv(norm);
                     for (int j = 0; j < m; j++)
                     {
                        MySetColor(draw, point[j][3], minv, maxv);
                        draw.glVertex3dv(point[j]);
                     }
                     draw.glEnd();
                  }
               }
            }
            break;

            case 2:  // PrepareCuttingPlaneLines() with mesh
            {
               if (shading == Shading::Noncomforming)
               {
                  // changes point for n > 4
                  DrawRefinedSurf(n, point[0], i, 2);
               }
               else
               {
                  // glBegin (GL_POLYGON);
                  gl3::GlBuilder line = cplines_buf.createBuilder();
                  line.glBegin(GL_LINE_LOOP);
                  for (int j = 0; j < n; j++)
                  {
                     line.glVertex3dv(point[j]);
                  }
                  line.glEnd();
               }
            }
            break;

            case 3:  // PrepareCuttingPlaneLines() with level lines
            {
               if (shading == Shading::Noncomforming)
               {
                  // changes point for n > 4
                  DrawRefinedSurf(n, point[0], i, 3);
               }
               else
               {
                  gl3::GlBuilder line = cplines_buf.createBuilder();
                  DrawPolygonLevelLines(line, point[0], n, level, false);
               }
            }
            break;
         }

         for (int j = 0; j < n2; j++)
         {
            cut_edges[j] = cut_edges[j+n];
         }
         n = n2;
         n2 = 0;
      }
   }
}

void VisualizationSceneSolution3d::CutRefinedElement(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vert_dist, const Vector &vals,
   const Geometry::Type geom, const int *elems, int num_elems, int func)
{
   double sc = 0.0;
   if (FaceShiftScale != 0.0)
   {
      double bbox_diam = sqrt ( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                                (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                                (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
      sc = FaceShiftScale * bbox_diam;
   }
   const int nv = Geometry::NumVerts[geom];

   gl3::GlBuilder bld = target.createBuilder();

   for (int i = 0; i < num_elems; i++)
   {
      Geometry::Type egeom = geom;
      int vert_flag[8], cut_edges[8];
      const int *elem = elems + i*nv;
      const int *edge_vert;
      int n = 0, n2, nev = 0;
      for (int j = 0; j < nv; j++)
      {
         if (elem[j] < 0)
         {
            // This appears to be a tetrahedral subelement of a refined pyramid
            egeom = Geometry::TETRAHEDRON;
            break;
         }
         nev++;
         vert_flag[j] = (vert_dist(elem[j]) >= 0.0) ? n++, 1 : 0;
      }
      if (n == 0 || n == nev) { continue; }

      CutElement(egeom, vert_flag, &edge_vert, cut_edges, &n, &n2);
      // n  = number of intersected edges
      // n2 = number of intersected edges, second polygon

      while (n > 2)
      {
         // 'pts' describe the intersecting polygon: triangle, quad, etc
         double pts[6][4]; // up to 6 points x (3 coordinates + 1 value)

         for (int j = 0; j < n; j++)
         {
            const int *ev = edge_vert + 2*cut_edges[j];
            double t = vert_dist(elem[ev[0]]);
            t = t / (t - vert_dist(elem[ev[1]]));
            for (int d = 0; d < 3; d++)
            {
               pts[j][d] = (1-t)*verts(d,elem[ev[0]]) + t*verts(d,elem[ev[1]]);
            }
            pts[j][3] = (1-t)*vals(elem[ev[0]]) + t*vals(elem[ev[1]]);
         }
         if (sc != 0.0)
         {
            const double *dir = CuttingPlane->Equation();
            for (int j = 0; j < n; j++)
            {
               // dir points into the visible side, so we add a minus to val:
               const double val = -sc * (pts[j][3] - minv) / (maxv - minv);
               for (int d = 0; d < 3; d++)
               {
                  pts[j][d] += val*dir[d];
               }
            }
         }

         if (func == 0) // draw surface
         {
            double norm[3];
            int no_norm, m = n;
            while (1)
            {
               if (m > 3)
               {
                  no_norm = Compute3DUnitNormal(pts[0], pts[1], pts[2], pts[3],
                                                norm);
                  if (no_norm && m > 4)
                  {
                     for (int j = 3; j < m; j++)
                     {
                        for (int k = 0; k < 4; k++)
                        {
                           pts[j-2][k] = pts[j][k];
                        }
                     }
                     m -= 2;
                     continue;
                  }
               }
               else
               {
                  no_norm = Compute3DUnitNormal(pts[0], pts[1], pts[2], norm);
               }
               break;
            }
            if (!no_norm)
            {
               bld.glBegin(GL_POLYGON);
               bld.glNormal3dv(norm);
               for (int j = 0; j < m; j++)
               {
                  MySetColor(bld, pts[j][3], minv, maxv);
                  bld.glVertex3dv(pts[j]);
               }
               bld.glEnd();
            }
         }
         else // draw level lines
         {
            DrawPolygonLevelLines(bld, pts[0], n, level, false);
         }

         for (int j = 0; j < n2; j++)
         {
            cut_edges[j] = cut_edges[j+n];
         }
         n = n2;
         n2 = 0;
      }
   }
}

void VisualizationSceneSolution3d::CutRefinedFace(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vert_dist, const Vector &vals,
   const Geometry::Type geom, const int *faces, int num_faces)
{
   double sc = 0.0;
   if (FaceShiftScale != 0.0)
   {
      double bbox_diam = sqrt ( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                                (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                                (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
      sc = FaceShiftScale * bbox_diam;
   }
   const int nv = Geometry::NumVerts[geom];

   for (int i = 0; i < num_faces; i++)
   {
      int vert_flag[4], cut_edges[4];
      const int *face = faces + i*nv;
      int n = 0;
      for (int j = 0; j < nv; j++)
      {
         vert_flag[j] = (vert_dist(face[j]) >= 0.0) ? n++, 1 : 0;
      }
      if (n == 0 || n == nv) { continue; }
      n = 0;
      for (int j = 0; j < nv; j++)
      {
         const int j1 = (j+1)%nv;
         if (vert_flag[j] != vert_flag[j1])
         {
            cut_edges[n++] = j;
         }
      }
      // n = number of intersected edges (2, or 4)
      double pts[4][4]; // up to 4 points x (3 coordinates + 1 value)
      for (int j = 0; j < n; j++)
      {
         const int v0 = cut_edges[j];
         const int v1 = (v0+1)%nv;
         double t = vert_dist(face[v0]);
         t = t / (t - vert_dist(face[v1]));
         for (int d = 0; d < 3; d++)
         {
            pts[j][d] = (1-t)*verts(d,face[v0]) + t*verts(d,face[v1]);
         }
         if (sc != 0.0)
         {
            pts[j][3] = (1-t)*vals(face[v0]) + t*vals(face[v1]);
         }
      }
      if (sc != 0.0)
      {
         const double *dir = CuttingPlane->Equation();
         for (int j = 0; j < n; j++)
         {
            // dir points into the visible side, so we add a minus to val:
            const double val = -sc * (pts[j][3] - minv) / (maxv - minv);
            for (int d = 0; d < 3; d++)
            {
               pts[j][d] += val*dir[d];
            }
         }
      }

      target.addLine(gl3::Vertex::create(pts[0]), gl3::Vertex::create(pts[1]));
      if (n == 4)
      {
         target.addLine(gl3::Vertex::create(pts[2]), gl3::Vertex::create(pts[3]));
      }
   }
}

void VisualizationSceneSolution3d::PrepareCuttingPlane()
{
   cplane_buf.clear();
   if (cp_drawelems && cplane && mesh->Dimension() == 3)
   {
      if (cplane == 2)
      {
         PrepareCuttingPlane2();
      }
      else if (cp_algo == 1)
      {
         CuttingPlaneFunc(1);
      }
      else
      {
         Vector vals, vert_dist;
         DenseMatrix pointmat;
         for (int i = 0; i < mesh->GetNE(); i++)
         {
            const Geometry::Type geom = mesh->GetElementBaseGeometry(i);
            RefinedGeometry *RefG =
               GLVisGeometryRefiner.Refine(geom, TimesToRefine);
            GridF->GetValues(i, RefG->RefPts, vals, pointmat);
            vert_dist.SetSize(pointmat.Width());
            for (int j = 0; j < pointmat.Width(); j++)
            {
               vert_dist(j) = CuttingPlane->Transform(&pointmat(0,j));
            }
            Array<int> &RG = RefG->RefGeoms;

            const int func = 0; // draw surface
            CutRefinedElement(cplane_buf, pointmat, vert_dist, vals, geom,
                              RG, RG.Size()/Geometry::NumVerts[geom], func);
         }
      }
   }
   updated_bufs.emplace_back(&cplane_buf);
}

void VisualizationSceneSolution3d::PrepareCuttingPlane2()
{
   int i, j, n = 0;
   double p[4][3], c[4], *coord;
   DenseMatrix pointmat, normals;
   Vector values;
   RefinedGeometry *RefG;

   Array<int> nodes;
   Array<int> partition (mesh -> GetNE());

   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = 0; // n will be the number of nodes behind the cutting plane
      mesh -> GetElementVertices(i, nodes);
      for (j = 0; j < nodes.Size(); j++)
      {
         if (node_pos[nodes[j]] >= 0.0) { n++; }
      }
      partition[i] = (n == nodes.Size()) ? 0 : 1;
   }

   for (i = 0; i < mesh -> GetNFaces(); i++)
   {
      int e1, e2;
      mesh -> GetFaceElements (i, &e1, &e2);
      if (e2 >= 0 && partition[e1] != partition[e2])
      {
         if (shading != Shading::Noncomforming)
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
            {
               DrawTriangle(cplane_buf, p, c, minv, maxv);
            }
            else
            {
               DrawQuad(cplane_buf, p, c, minv, maxv);
            }
         }
         else // shading == 2
         {
            RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceGeometry (i),
                                               TimesToRefine);
            // partition[e1] is 0 if e1 is behind the cutting plane
            // and 1 otherwise
            int dir = partition[e1];
            GridF -> GetFaceValues (i, dir, RefG->RefPts, values, pointmat);
            GetFaceNormals(i, dir, RefG->RefPts, normals);
            switch (mesh -> GetFaceGeometry (i))
            {
               case Geometry::TRIANGLE:  n = 3; break;
               case Geometry::SQUARE:    n = 4; break;
               default:
                  MFEM_ABORT("Invalid element type");
                  break;
            }
            // DrawRefinedSurf (n, pointmat, values, RefG->RefGeoms);
            DrawPatch(cplane_buf, pointmat, values, normals, n, RefG->RefGeoms,
                      minv, maxv, dir ? -3 : 2);
         } // end shading == 2
      }
   }
}

void VisualizationSceneSolution3d::PrepareCuttingPlaneLines()
{
   cplines_buf.clear();

   if (cp_drawmesh && cplane && mesh->Dimension() == 3)
   {
      if (cplane == 2 && cp_drawmesh != 3)
      {
         PrepareCuttingPlaneLines2();
      }
      else
      {
         if (cp_drawmesh == 1 && cp_algo == 1)
         {
            CuttingPlaneFunc(2);
         }
         else if (cp_drawmesh == 1)
         {
            Vector vert_dist, vals;
            DenseMatrix pointmat;
            int num_faces = mesh->GetNFaces();
            if (mesh->NURBSext)
            {
               // Note: for NURBS meshes, the methods
               //       Mesh::GetFaceTransformation() and
               //       GridFunction::GetFaceValues() are not supported.
               cout << _MFEM_FUNC_NAME
                    << ": NURBS mesh: cut faces will not be drawn!" << endl;
               num_faces = 0;
            }
            for (int i = 0; i < num_faces; i++)
            {
               const Geometry::Type geom = mesh->GetFaceGeometry(i);
               RefinedGeometry *RefG =
                  GLVisGeometryRefiner.Refine(geom, TimesToRefine);
               if (FaceShiftScale == 0.0)
               {
                  ElementTransformation *T = mesh->GetFaceTransformation(i);
                  T->Transform(RefG->RefPts, pointmat);
               }
               else
               {
                  const int side = 2;
                  GridF->GetFaceValues(i, side, RefG->RefPts, vals, pointmat);
                  // For discontinuous grid function, we should draw two edges.
               }
               vert_dist.SetSize(pointmat.Width());
               for (int j = 0; j < pointmat.Width(); j++)
               {
                  vert_dist(j) = CuttingPlane->Transform(&pointmat(0,j));
               }
               Array<int> &RG = RefG->RefGeoms;

               CutRefinedFace(cplines_buf, pointmat, vert_dist, vals, geom,
                              RG, RG.Size()/Geometry::NumVerts[geom]);
            }
         }
         else if (cp_algo == 1)
         {
            CuttingPlaneFunc(3);
         }
         else
         {
            Vector vals, vert_dist;
            DenseMatrix pointmat;
            for (int i = 0; i < mesh->GetNE(); i++)
            {
               const Geometry::Type geom = mesh->GetElementBaseGeometry(i);
               RefinedGeometry *RefG =
                  GLVisGeometryRefiner.Refine(geom, TimesToRefine);
               GridF->GetValues(i, RefG->RefPts, vals, pointmat);
               vert_dist.SetSize(pointmat.Width());
               for (int j = 0; j < pointmat.Width(); j++)
               {
                  vert_dist(j) = CuttingPlane->Transform(&pointmat(0,j));
               }
               Array<int> &RG = RefG->RefGeoms;

               const int func = 1; // draw level lines
               CutRefinedElement(cplines_buf, pointmat, vert_dist, vals, geom,
                                 RG, RG.Size()/Geometry::NumVerts[geom], func);
            }
         }
      }
   }

   updated_bufs.emplace_back(&cplines_buf);
}

void VisualizationSceneSolution3d::PrepareCuttingPlaneLines2()
{
   int i, j, n = 0;
   double *coord;
   DenseMatrix pointmat;
   Vector values;
   RefinedGeometry *RefG;

   Array<int> nodes;
   Array<int> partition (mesh -> GetNE());
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = 0;  // n will be the number of nodes behind the cutting plane
      mesh -> GetElementVertices(i,nodes);
      for (j=0; j<nodes.Size(); j++)
      {
         if (node_pos[nodes[j]] >= 0.0) { n++; }
      }
      partition[i] = (n == nodes.Size()) ? 0 : 1;
   }

   for (i = 0; i < mesh -> GetNFaces(); i++)
   {
      int e1, e2;
      mesh -> GetFaceElements (i, &e1, &e2);
      if (e2 >= 0 && partition[e1] != partition[e2])
      {
         if (shading != Shading::Noncomforming)
         {
            mesh -> GetFaceVertices (i, nodes);
            pointmat.SetSize(4, nodes.Size());
            for (j = 0; j < nodes.Size(); j++)
            {
               coord = mesh -> GetVertex(nodes[j]);
               pointmat(0, j) = coord[0];
               pointmat(1, j) = coord[1];
               pointmat(2, j) = coord[2];
               pointmat(3, j) = (*sol)(nodes[j]);
            }
            gl3::GlBuilder line = cplines_buf.createBuilder();
            switch (cp_drawmesh)
            {
               case 1:
               {
                  if (mesh_coarse)
                  {
                     DrawFaceCoarseSurfEdges(line, i, pointmat);
                  }
                  else
                  {
                     // glBegin(GL_POLYGON);
                     line.glBegin(GL_LINE_LOOP);
                     for (j = 0; j < nodes.Size(); j++)
                     {
                        line.glVertex3dv(&pointmat(0,j));
                     }
                     line.glEnd();
                  }
                  break;
               }
               case 2:
                  DrawPolygonLevelLines(line, pointmat.GetData(), nodes.Size(), level, false);
                  break;
            }
         }
         else // shading == 2
         {
            RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceGeometry (i),
                                               TimesToRefine);
            // partition[e1] is 0 if e1 is behind the cutting plane
            // and 1 otherwise
            int di = partition[e1];
            GridF -> GetFaceValues (i, di, RefG->RefPts, values, pointmat);
            switch (mesh -> GetFaceGeometry (i))
            {
               case Geometry::TRIANGLE:  n = 3; break;
               case Geometry::SQUARE:    n = 4; break;
               default:
                  MFEM_ABORT("Invalid element type");
                  break;
            }
            switch (cp_drawmesh)
            {
               case 1:
               {
                  if (mesh_coarse)
                  {
                     gl3::GlBuilder line = cplines_buf.createBuilder();
                     DrawFaceCoarseSurfEdges(line, i, pointmat, &RefG->RefPts, &RefG->RefEdges);
                  }
                  else
                  {
                     DrawRefinedSurfEdges (n, pointmat, values, RefG->RefEdges);
                  }
                  break;
               }
               case 2:
                  DrawRefinedSurfLevelLines (n, pointmat, values,
                                             RefG->RefGeoms);
                  break;
            }
         } // end shading == 2
      }
   }
}

thread_local int triangle_counter;
thread_local int quad_counter;

void VisualizationSceneSolution3d::DrawTetLevelSurf(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vals, const int *ind,
   const Array<double> &surf_levels, const DenseMatrix *grad)
{
   double t, lvl, normal[3], vert[4][3], norm[4][3];
   int i, j, l, pos[4];
   bool flipped;

   gl3::GlBuilder draw = target.createBuilder();

   for (l = 0; l < surf_levels.Size(); l++)
   {
      lvl = surf_levels[l];

      for (i = 0; i < 4; i++)
      {
         pos[i] = ind[i];
      }
      i = 0;
      j = 4;
      flipped = false;
      do
      {
         // (i < j) is true
         while (vals(pos[i]) < lvl)
         {
            i++;
            if (i == j)
            {
               goto step_one;
            }
         }
         // (i < j) && (vals[pos[i]] >= lvl) is true
         do
         {
            j--;
            if (i == j)
            {
               goto step_one;
            }
         }
         while (vals(pos[j]) >= lvl);
         // (i < j) && (vals[pos[i]] >= lvl) && (vals[pos[j]] < lvl) is true
         Swap<int>(pos[i], pos[j]);
         flipped = !flipped;
         i++;
      }
      while (i < j);
   step_one:
      if (flipped)
      {
         if (i >= 2)
         {
            Swap<int>(pos[0], pos[1]);
         }
         else
         {
            Swap<int>(pos[2], pos[3]);
         }
      }

      if (j == 3)
      {
         Swap<int>(pos[0], pos[3]);
         j = 1;
      }
      if (j == 1)
      {
         for (int k = 0; k < 3; k++)
         {
            int p0 = pos[0];
            int p1 = pos[k+1];
            t = (lvl - vals(p0)) / (vals(p1) - vals(p0));
            for (int d = 0; d < 3; d++)
            {
               vert[k][d] = (1.0 - t) * verts(d, p0) + t * verts(d, p1);
            }
            if (grad)
               for (int d = 0; d < 3; d++)
               {
                  norm[k][d] = (1.0 - t) * (*grad)(d, p0) + t * (*grad)(d, p1);
               }
         }

         if (grad == NULL)
         {
            if (!Compute3DUnitNormal(vert[0], vert[1], vert[2], normal))
            {
               MySetColor(draw, lvl, minv, maxv);
               draw.glNormal3dv(normal);
               draw.glBegin(GL_TRIANGLES);
               for (int k = 0; k < 3; k++)
               {
                  draw.glVertex3dv(vert[k]);
               }
               draw.glEnd();
               triangle_counter++;
            }
         }
         else
         {
            MySetColor(draw, lvl, minv, maxv);
            draw.glBegin(GL_TRIANGLES);
            for (int k = 0; k < 3; k++)
            {
               Normalize(norm[k]);
               draw.glNormal3dv(norm[k]);
               draw.glVertex3dv(vert[k]);
            }
            draw.glEnd();
            triangle_counter++;
         }
      }
      else if (j == 2)
      {
         static const int idx[4][2] = { {0, 2}, {0, 3}, {1, 3}, {1, 2} };

         for (int k = 0; k < 4; k++)
         {
            int p0 = pos[idx[k][0]];
            int p1 = pos[idx[k][1]];
            t = (lvl - vals(p0)) / (vals(p1) - vals(p0));
            for (int d = 0; d < 3; d++)
            {
               vert[k][d] = (1.0 - t) * verts(d, p0) + t * verts(d, p1);
            }
            if (grad)
               for (int d = 0; d < 3; d++)
               {
                  norm[k][d] = (1.0 - t) * (*grad)(d, p0) + t * (*grad)(d, p1);
               }
         }

         if (grad == NULL)
         {
            if (!Compute3DUnitNormal(vert[0], vert[1], vert[2], vert[3],
                                     normal))
            {
               MySetColor(draw, lvl, minv, maxv);
               draw.glNormal3dv(normal);
               draw.glBegin(GL_QUADS);
               for (int k = 0; k < 4; k++)
               {
                  draw.glVertex3dv(vert[k]);
               }
               draw.glEnd();
               quad_counter++;
            }
         }
         else
         {
            MySetColor(draw, lvl, minv, maxv);
            draw.glBegin(GL_QUADS);
            for (int k = 0; k < 4; k++)
            {
               Normalize(norm[k]);
               draw.glNormal3dv(norm[k]);
               draw.glVertex3dv(vert[k]);
            }
            draw.glEnd();
            quad_counter++;
         }
      }
   }
}

// static method
int VisualizationSceneSolution3d::GetPyramidFaceSplits(
   const Array<bool> &quad_diag, const Array<int> &faces,
   const Array<int> &ofaces)
{
   int fs = 0;
   bool diag = quad_diag[faces[0]];
   if ((ofaces[0]/2)%2) // orientations 2,3,6,7
   {
      diag = !diag;
   }
   fs = 2*fs + diag;
   return fs;
}

void VisualizationSceneSolution3d::DrawRefinedPyramidLevelSurf(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vals, const int *RG, const int np,
   const int face_splits, const DenseMatrix *grad)
{
#if 1
   static const int pyr_tets[2][4] =
   {
      { 0, 1, 2, 4 }, { 0, 2, 3, 4 }
   };
   for (int k = 0; k < np; k++)
   {
      const int *hv = &RG[5*k];
      if (hv[4] > 0)
      {
         for (int j = 0; j < 2; j++)
         {
            int m_ind[4];
            for (int i = 0; i < 4; i++)
            {
               m_ind[i] = hv[pyr_tets[j][i]];
            }
            DrawTetLevelSurf(target, verts, vals, m_ind, levels, grad);
         }
      }
      else
      {
         DrawTetLevelSurf(target, verts, vals, hv, levels, grad);
      }
   }
#else
   // MLS: Not sure how to adapt this to Pyramids so skip it for now
   static const int pri_tets[8-2][3][4] =
   {
      // 0 = 000 (see below; prism is split into 6 tets)
      { { 0, 1, 2, 5 }, { 0, 1, 5, 4 }, { 0, 3, 4, 5 } }, // 1 = 001
      { { 0, 1, 2, 4 }, { 0, 2, 3, 4 }, { 2, 3, 4, 5 } }, // 2 = 010
      { { 0, 1, 2, 4 }, { 0, 2, 5, 4 }, { 0, 3, 4, 5 } }, // 3 = 011
      { { 0, 1, 2, 3 }, { 1, 2, 3, 5 }, { 1, 3, 4, 5 } }, // 4 = 100
      { { 0, 1, 2, 5 }, { 0, 1, 5, 3 }, { 1, 3, 4, 5 } }, // 5 = 101
      { { 0, 1, 2, 3 }, { 1, 2, 3, 4 }, { 2, 3, 4, 5 } }  // 6 = 110
      // 7 = 111 (see below; prism is split into 6 tets)
   };
   static const int pri_tets_0[6][4] =
   { {0,1,2,6}, {0,1,6,4}, {0,2,3,6}, {0,6,3,4}, {2,3,6,5}, {3,4,6,5} };
   static const int pri_tets_7[6][4] =
   { {0,1,2,6}, {0,2,5,6}, {0,1,6,3}, {0,3,6,5}, {1,3,4,6}, {3,4,6,5} };
   const int fs2 = (~face_splits&7)/2 + 4*((~face_splits&7)%2);
   const int n = (np == 1) ? 1 : TimesToRefine;

   double vs_data[7], pm_data[3*7], gd_data[3*7];
   Vector vs(vs_data, 7);
   DenseMatrix pm(pm_data, 3, 7), gd(gd_data, 3, 7);

   for (int k = 0, l0 = -1, l1 = -1; k < np; k++)
   {
      const int pk = k % (n*n);
      if (pk == 0)
      {
         l0 = 0; l1 = 2*n-1;
      }
      else if (pk == l1)
      {
         const int s = l1-l0;
         l0 = l1;
         l1 += (s-2);
      }
      const int fsl = ((pk-l0)%2 == 0) ? face_splits : fs2;
      // The algorithm for choosing 'fsl' used above assumes the refined prisms
      // are listed in certain order -- see the prism case in
      // mfem::GeometryRefiner::Refine.
      const int *pv = &RG[6*k];
      if (fsl == 0)
      {
         for (int j = 0; j < 6; j++)
         {
            const int idx = pv[j];
            for (int d = 0; d < 3; d++)
            {
               pm(d,j) = verts(d,idx);
               if (grad) { gd(d,j) = (*grad)(d,idx); }
            }
            vs(j) = vals(idx);
         }
         for (int d = 0; d < 3; d++)
         {
            pm(d,6) = 0.5*(pm(d,1) + pm(d,5));
            if (grad) { gd(d,6) = 0.5*(gd(d,1) + gd(d,5)); }
         }
         vs(6) = 0.5*(vs(1) + vs(5));
         const DenseMatrix *gd_ = grad ? &gd : NULL;
         for (int k = 0; k < 6; k++)
         {
            DrawTetLevelSurf(pm, vs, pri_tets_0[k], levels, gd_);
         }
      }
      else if (fsl == 7)
      {
         for (int j = 0; j < 6; j++)
         {
            const int idx = pv[j];
            for (int d = 0; d < 3; d++)
            {
               pm(d,j) = verts(d,idx);
               if (grad) { gd(d,j) = (*grad)(d,idx); }
            }
            vs(j) = vals(idx);
         }
         for (int d = 0; d < 3; d++)
         {
            pm(d,6) = 0.5*(pm(d,2) + pm(d,4));
            if (grad) { gd(d,6) = 0.5*(gd(d,2) + gd(d,4)); }
         }
         vs(6) = 0.5*(vs(2) + vs(4));
         const DenseMatrix *gd_ = grad ? &gd : NULL;
         for (int k = 0; k < 6; k++)
         {
            DrawTetLevelSurf(pm, vs, pri_tets_7[k], levels, gd_);
         }
      }
      else
      {
         int m_ind[4];
         for (int j = 0; j < 3; j++)
         {
            for (int i = 0; i < 4; i++)
            {
               m_ind[i] = pv[pri_tets[fsl-1][j][i]];
            }
            DrawTetLevelSurf(verts, vals, m_ind, levels, grad);
         }
      }
   }
#endif
}

// static method
int VisualizationSceneSolution3d::GetWedgeFaceSplits(
   const Array<bool> &quad_diag, const Array<int> &faces,
   const Array<int> &ofaces)
{
   int fs = 0;
   for (int lf = 2; lf < 5; lf++)
   {
      bool diag = quad_diag[faces[lf]];
      if ((ofaces[lf]/2)%2) // orientations 2,3,6,7
      {
         diag = !diag;
      }
      fs = 2*fs + diag;
   }
   return fs;
}

void VisualizationSceneSolution3d::DrawRefinedWedgeLevelSurf(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vals, const int *RG, const int np,
   const int face_splits, const DenseMatrix *grad)
{
#if 0
   static const int pri_tets[3][4] =
   {
      { 0, 1, 2, 5 }, { 0, 1, 5, 3 }, { 1, 3, 4, 5 }
   };
   for (int k = 0; k < np; k++)
   {
      const int *hv = &RG[6*k];
      for (int j = 0; j < 3; j++)
      {
         int m_ind[4];
         for (int i = 0; i < 4; i++)
         {
            m_ind[i] = hv[pri_tets[j][i]];
         }
         DrawTetLevelSurf(verts, vals, m_ind, levels, grad);
      }
   }
#else
   static const int pri_tets[8-2][3][4] =
   {
      // 0 = 000 (see below; prism is split into 6 tets)
      { { 0, 1, 2, 5 }, { 0, 1, 5, 4 }, { 0, 3, 4, 5 } }, // 1 = 001
      { { 0, 1, 2, 4 }, { 0, 2, 3, 4 }, { 2, 3, 4, 5 } }, // 2 = 010
      { { 0, 1, 2, 4 }, { 0, 2, 5, 4 }, { 0, 3, 4, 5 } }, // 3 = 011
      { { 0, 1, 2, 3 }, { 1, 2, 3, 5 }, { 1, 3, 4, 5 } }, // 4 = 100
      { { 0, 1, 2, 5 }, { 0, 1, 5, 3 }, { 1, 3, 4, 5 } }, // 5 = 101
      { { 0, 1, 2, 3 }, { 1, 2, 3, 4 }, { 2, 3, 4, 5 } }  // 6 = 110
      // 7 = 111 (see below; prism is split into 6 tets)
   };
   static const int pri_tets_0[6][4] =
   { {0,1,2,6}, {0,1,6,4}, {0,2,3,6}, {0,6,3,4}, {2,3,6,5}, {3,4,6,5} };
   static const int pri_tets_7[6][4] =
   { {0,1,2,6}, {0,2,5,6}, {0,1,6,3}, {0,3,6,5}, {1,3,4,6}, {3,4,6,5} };
   const int fs2 = (~face_splits&7)/2 + 4*((~face_splits&7)%2);
   const int n = (np == 1) ? 1 : TimesToRefine;

   double vs_data[7], pm_data[3*7], gd_data[3*7];
   Vector vs(vs_data, 7);
   DenseMatrix pm(pm_data, 3, 7), gd(gd_data, 3, 7);

   for (int k = 0, l0 = -1, l1 = -1; k < np; k++)
   {
      const int pk = k % (n*n);
      if (pk == 0)
      {
         l0 = 0; l1 = 2*n-1;
      }
      else if (pk == l1)
      {
         const int s = l1-l0;
         l0 = l1;
         l1 += (s-2);
      }
      const int fsl = ((pk-l0)%2 == 0) ? face_splits : fs2;
      // The algorithm for choosing 'fsl' used above assumes the refined prisms
      // are listed in certain order -- see the prism case in
      // mfem::GeometryRefiner::Refine.
      const int *pv = &RG[6*k];
      if (fsl == 0)
      {
         for (int j = 0; j < 6; j++)
         {
            const int idx = pv[j];
            for (int d = 0; d < 3; d++)
            {
               pm(d,j) = verts(d,idx);
               if (grad) { gd(d,j) = (*grad)(d,idx); }
            }
            vs(j) = vals(idx);
         }
         for (int d = 0; d < 3; d++)
         {
            pm(d,6) = 0.5*(pm(d,1) + pm(d,5));
            if (grad) { gd(d,6) = 0.5*(gd(d,1) + gd(d,5)); }
         }
         vs(6) = 0.5*(vs(1) + vs(5));
         const DenseMatrix *gd_ = grad ? &gd : NULL;
         for (int j = 0; j < 6; j++)
         {
            DrawTetLevelSurf(target, pm, vs, pri_tets_0[j], levels, gd_);
         }
      }
      else if (fsl == 7)
      {
         for (int j = 0; j < 6; j++)
         {
            const int idx = pv[j];
            for (int d = 0; d < 3; d++)
            {
               pm(d,j) = verts(d,idx);
               if (grad) { gd(d,j) = (*grad)(d,idx); }
            }
            vs(j) = vals(idx);
         }
         for (int d = 0; d < 3; d++)
         {
            pm(d,6) = 0.5*(pm(d,2) + pm(d,4));
            if (grad) { gd(d,6) = 0.5*(gd(d,2) + gd(d,4)); }
         }
         vs(6) = 0.5*(vs(2) + vs(4));
         const DenseMatrix *gd_ = grad ? &gd : NULL;
         for (int j = 0; j < 6; j++)
         {
            DrawTetLevelSurf(target, pm, vs, pri_tets_7[j], levels, gd_);
         }
      }
      else
      {
         int m_ind[4];
         for (int j = 0; j < 3; j++)
         {
            for (int i = 0; i < 4; i++)
            {
               m_ind[i] = pv[pri_tets[fsl-1][j][i]];
            }
            DrawTetLevelSurf(target, verts, vals, m_ind, levels, grad);
         }
      }
   }
#endif
}

// static method
int VisualizationSceneSolution3d::GetHexFaceSplits(
   const Array<bool> &quad_diag, const Array<int> &faces,
   const Array<int> &ofaces)
{
   int fs = 0;
   for (int lf = 0; lf < 6; lf++)
   {
      bool diag = quad_diag[faces[lf]];
      if ((ofaces[lf]/2)%2) // orientations 2,3,6,7
      {
         diag = !diag;
      }
      fs = 2*fs + diag;
   }
   return fs;
}

void VisualizationSceneSolution3d::DrawRefinedHexLevelSurf(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vals, const int *RG, const int nh,
   const int face_splits, const DenseMatrix *grad)
{
#if 0
   static const int hex_tets[6][4] =
   {
      { 0, 1, 2, 6 }, { 0, 5, 1, 6 }, { 0, 4, 5, 6 },
      { 0, 2, 3, 6 }, { 0, 3, 7, 6 }, { 0, 7, 4, 6 }
   };
   for (int k = 0; k < nh; k++)
   {
      const int *hv = &RG[8*k];
      for (int j = 0; j < 6; j++)
      {
         int m_ind[4];
         for (int i = 0; i < 4; i++)
         {
            m_ind[i] = hv[hex_tets[j][i]];
         }
         DrawTetLevelSurf(verts, vals, m_ind, levels, grad);
      }
   }
#else
   const int n = (nh == 1) ? 1 : TimesToRefine;

   double vs_data[9], pm_data[3*9], gd_data[3*9];
   Vector vs(vs_data, 9);
   DenseMatrix pm(pm_data, 3, 9), gd(gd_data, 3, 9);

   for (int k = 0; k < nh; k++)
   {
      const int ix = k%n, iy = (k/n)%n, iz = k/(n*n);
      int fsl = face_splits;
      if (ix != 0)
      {
         // Copy the right face bit (face 2, bit mask 8) to the left face bit
         // (face 4, bit mask 2) and flip it.
         fsl = (fsl&(63-2)) + (2-(fsl&8)/4);
      }
      if (iy != 0)
      {
         // Copy the back face bit (face 3, bit mask 4) to the front face bit
         // (face 1, bit mask 16) and flip it.
         fsl = (fsl&(63-16)) + (16-(fsl&4)*4);
      }
      if (iz != 0)
      {
         // Copy the top face bit (face 5, bit mask 1) to the bottom face bit
         // (face 0, bit mask 32) and flip it.
         fsl = (fsl&(63-32)) + (32-(fsl&1)*32);
      }

      const int *hv_orig = &RG[8*k], *hv;
      int hv1[8], hv2[8];
#if 1
      // Find a pair of opposite faces that are split into triangles using
      // two parallel diagonals (if such pair exists).
      if (1-(fsl&1) == (fsl&32)/32)
      {
         // The bottom and top faces are split in the same direction.
         hv = hv_orig;
      }
      else if (4-(fsl&4) == (fsl&16)/4)
      {
         // The front and back faces are split in the same direction.
         // Rotate the hex around the x-axis to bring front-back to bottom-top.
         static const int rot_x[8] = { 4, 5, 1, 0, 7, 6, 2, 3 };
         for (int j = 0; j < 8; j++)
         {
            hv1[j] = hv_orig[rot_x[j]];
         }
         hv = hv1;
         // fsl bits change: a|b|c|d|e|f -> b|f|1-c|a|1-e|d
         fsl = (fsl&16)*2 + (fsl&1)*16 + (8-(fsl&8)) + (fsl&32)/8 +
               (2-(fsl&2)) + (fsl&4)/4;
      }
      else if (2-(fsl&2) == (fsl&8)/4)
      {
         // The left and right faces are split in the same direction.
         // Rotate the hex around the y-axis to bring left-right to top-bottom.
         static const int rot_y[8] = { 1, 5, 6, 2, 0, 4, 7, 3 };
         for (int j = 0; j < 8; j++)
         {
            hv1[j] = hv_orig[rot_y[j]];
         }
         hv = hv1;
         // fsl bit change: a|b|c|d|e|f -> 1-c|1-b|1-f|1-d|1-a|1-e
         fsl = ~fsl&63;
         fsl = (fsl&8)*4 + (fsl&16) + (fsl&1)*8 + (fsl&4) + (fsl&32)/16 +
               (fsl&2)/2;
      }
      else
#endif
      {
         // All opposite faces are split in opposite directions.
         // Split the hex into 12 tets using the hex center as a vertex in all
         // 12 tets.
         for (int j = 0; j < 8; j++)
         {
            const int idx = hv_orig[j];
            for (int d = 0; d < 3; d++)
            {
               pm(d,j) = verts(d,idx);
               if (grad) { gd(d,j) = (*grad)(d,idx); }
            }
            vs(j) = vals(idx);
         }
         for (int d = 0; d < 3; d++)
         {
            pm(d,8) = 0.0;
            if (grad) { gd(d,8) = 0.0; }
         }
         vs(8) = 0.0;
         for (int j = 0; j < 8; j++)
         {
            for (int d = 0; d < 3; d++)
            {
               pm(d,8) += pm(d,j);
               if (grad) { gd(d,8) += gd(d,j); }
            }
            vs(8) += vs(j);
         }
         for (int d = 0; d < 3; d++)
         {
            pm(d,8) *= 0.125;
            if (grad) { gd(d,8) *= 0.125; }
         }
         vs(8) *= 0.125;

         typedef Mesh::hex_t hex_t;
         for (int j = 0; j < hex_t::NumFaces; j++)
         {
            int tv[8];
            const int *fv = hex_t::FaceVert[j];
            const bool diag = fsl&(32>>j);
            if (diag == 0)
            {
               tv[0] = fv[2];
               tv[1] = fv[1];
               tv[2] = fv[0];
               tv[3] = 8;
               tv[4] = fv[0];
               tv[5] = fv[3];
               tv[6] = fv[2];
               tv[7] = 8;
            }
            else
            {
               tv[0] = fv[1];
               tv[1] = fv[0];
               tv[2] = fv[3];
               tv[3] = 8;
               tv[4] = fv[3];
               tv[5] = fv[2];
               tv[6] = fv[1];
               tv[7] = 8;
            }
            const DenseMatrix *gp = grad ? &gd : NULL;
            DrawTetLevelSurf(target, pm, vs, &tv[0], levels, gp);
            DrawTetLevelSurf(target, pm, vs, &tv[4], levels, gp);
         }

         continue;
      }

      // Rotate the hex so that the diagonal edge splitting the bottom face is
      // the edge 0-2.
      if ((fsl&32) == 0)
      {
         // Rotate the hex around the z-axis -- left-right becomes front-back.
         static const int rot_z[8] = { 3, 0, 1, 2, 7, 4, 5, 6 };
         for (int j = 0; j < 8; j++)
         {
            hv2[j] = hv[rot_z[j]];
         }
         hv = hv2;
         // fsl bit change: a|b|c|d|e|f -> 1-a|e|b|c|d|1-f
         fsl = (32-(fsl&32)) + (fsl&2)*8 + (fsl&(16+8+4))/2 + (1-(fsl&1));
      }

      // Split the hex into two prisms using the diagonal face 0-2-6-4.
      const int pv[2][6] =
      {
         { hv[0], hv[1], hv[2], hv[4], hv[5], hv[6] },
         { hv[2], hv[3], hv[0], hv[6], hv[7], hv[4] }
      };
      // Choose the shorter diagonal on the face 0-2-6-4.
      const double l06 = Distance(&verts(0,hv[0]), &verts(0,hv[6]), 3);
      const double l24 = Distance(&verts(0,hv[2]), &verts(0,hv[4]), 3);
      const bool diag = (l06 > 1.01*l24);
      const int fs1 = (fsl&(16+8))/4 + !diag; // a|b|c|d|e|f -> b|c|1-diag
      const int fs2 = (fsl&(4+2)) + diag; // a|b|c|d|e|f -> d|e|diag
      DrawRefinedWedgeLevelSurf(target, verts, vals, pv[0], 1, fs1, grad);
      DrawRefinedWedgeLevelSurf(target, verts, vals, pv[1], 1, fs2, grad);
   }
#endif
}

void VisualizationSceneSolution3d::PrepareLevelSurf()
{
   static const int ident[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

   Vector vals;
   DenseMatrix pointmat, grad;
   Array<int> vertices;

   lsurf_buf.clear();
   if (drawlsurf == 0 || mesh->Dimension() != 3)
   {
      //  Create empty list
      return;
   }

   triangle_counter = quad_counter = 0;

   levels.SetSize(nlevels);
   for (int l = 0; l < nlevels; l++)
   {
      double lvl = ((double)(50*l+drawlsurf) / (nlevels*50));
      levels[l] = ULogVal(lvl);
   }

   // For every quad face, choose the shorter diagonal to split the quad into
   // two triangles. Elements adjacent to that quad face (wedge or hex) will use
   // the same diagonal when subdividing the element.
   Array<bool> quad_diag;
   if (mesh->HasGeometry(Geometry::SQUARE))
   {
      quad_diag.SetSize(mesh->GetNFaces());
      for (int fi = 0; fi < mesh->GetNFaces(); fi++)
      {
         const Element *face = mesh->GetFace(fi);
         if (face->GetType() != Element::QUADRILATERAL) { continue; }
         ElementTransformation *T = mesh->GetFaceTransformation(fi);
         T->Transform(*Geometries.GetVertices(Geometry::SQUARE), pointmat);
         const double l02 = Distance(&pointmat(0,0), &pointmat(0,2), 3);
         const double l13 = Distance(&pointmat(0,1), &pointmat(0,3), 3);
         quad_diag[fi] = (l02 > 1.01*l13);
      }
   }

   Array<int> faces, ofaces;

   if (shading != Shading::Noncomforming)
   {
      for (int ie = 0; ie < mesh->GetNE(); ie++)
      {
         mesh->GetPointMatrix(ie, pointmat);
         mesh->GetElementVertices(ie, vertices);
         vals.SetSize(vertices.Size());
         for (int j = 0; j < vertices.Size(); j++)
         {
            vals(j) = (*sol)(vertices[j]);
         }

         switch (mesh->GetElementType(ie))
         {
            case Element::TETRAHEDRON:
               DrawTetLevelSurf(lsurf_buf, pointmat, vals, ident, levels);
               break;
            case Element::PYRAMID:
            {
               mesh->GetElementFaces(ie, faces, ofaces);
               const int fs = GetPyramidFaceSplits(quad_diag, faces, ofaces);
               DrawRefinedPyramidLevelSurf(lsurf_buf, pointmat, vals, ident, 1, fs);
            }
            case Element::WEDGE:
            {
               mesh->GetElementFaces(ie, faces, ofaces);
               const int fs = GetWedgeFaceSplits(quad_diag, faces, ofaces);
               DrawRefinedWedgeLevelSurf(lsurf_buf, pointmat, vals, ident, 1, fs);
            }
            break;
            case Element::HEXAHEDRON:
            {
               mesh->GetElementFaces(ie, faces, ofaces);
               const int fs = GetHexFaceSplits(quad_diag, faces, ofaces);
               DrawRefinedHexLevelSurf(lsurf_buf, pointmat, vals, ident, 1, fs);
            }
            break;
            default:
               MFEM_ABORT("Unrecognized 3D element type \""
                          << mesh->GetElementType(ie) << "\"");
         }
      }
   }
   else // shading == 2
   {
      RefinedGeometry *RefG;
#define GLVIS_SMOOTH_LEVELSURF_NORMALS
      const DenseMatrix *gp = NULL;

      for (int ie = 0; ie < mesh->GetNE(); ie++)
      {
         const Geometry::Type geom = mesh->GetElementBaseGeometry(ie);

         RefG = GLVisGeometryRefiner.Refine(geom, TimesToRefine);
         GridF->GetValues(ie, RefG->RefPts, vals, pointmat);
#ifdef GLVIS_SMOOTH_LEVELSURF_NORMALS
         const int map_type = GridF->FESpace()->GetFE(ie)->GetMapType();
         if (map_type == FiniteElement::MapType::VALUE)
         {
            GridF->GetGradients(ie, RefG->RefPts, grad);
            gp = &grad;
         }
         else if (map_type == FiniteElement::MapType::INTEGRAL)
         {
            FiniteElementSpace *fes = GridF->FESpace();
            const FiniteElement *fe = fes->GetFE(ie);
            const int ndof = fe->GetDof();
            const int ndim = fe->GetDim();
            ElementTransformation *Trans = fes->GetElementTransformation(ie);
            const IntegrationRule &ir = RefG->RefPts;
            DenseMatrix dshape(ndof, ndim);
            Vector lval, gh(ndim), gcol;

            GridF->GetElementDofValues(ie, lval);

            // Local projection to value-based FE
            const IntegrationRule &nodes = fe->GetNodes();
            for (int n = 0; n < nodes.GetNPoints(); n++)
            {
               const IntegrationPoint &ip = nodes.IntPoint(n);
               Trans->SetIntPoint(&ip);
               lval(n) /= Trans->Weight(); // value = dof / |J|
            }

            // Gradient calculation
            grad.SetSize(fe->GetDim(), ir.GetNPoints());
            for (int q = 0; q < ir.GetNPoints(); q++)
            {
               const IntegrationPoint &ip = ir.IntPoint(q);
               fe->CalcDShape(ip, dshape);
               dshape.MultTranspose(lval, gh);
               Trans->SetIntPoint(&ip);
               grad.GetColumnReference(q, gcol);
               const DenseMatrix &Jinv = Trans->InverseJacobian();
               Jinv.MultTranspose(gh, gcol);
            }
            gp = &grad;
         }
         else
         {
            MFEM_ABORT("Unknown mapping type");
         }
#endif

         Array<int> &RG = RefG->RefGeoms;
         const int nv = mesh->GetElement(ie)->GetNVertices();
         const int nre = RG.Size()/nv;

         if (geom == Geometry::TETRAHEDRON)
         {
            for (int k = 0; k < nre; k++)
            {
               DrawTetLevelSurf(lsurf_buf, pointmat, vals, &RG[nv*k], levels, gp);
            }
         }
         else if (geom == Geometry::PYRAMID)
         {
            mesh->GetElementFaces(ie, faces, ofaces);
            const int fs = GetPyramidFaceSplits(quad_diag, faces, ofaces);
            DrawRefinedPyramidLevelSurf(lsurf_buf, pointmat, vals, RG, nre, fs, gp);
         }
         else if (geom == Geometry::PRISM)
         {
            mesh->GetElementFaces(ie, faces, ofaces);
            const int fs = GetWedgeFaceSplits(quad_diag, faces, ofaces);
            DrawRefinedWedgeLevelSurf(lsurf_buf, pointmat, vals, RG, nre, fs, gp);
         }
         else if (geom == Geometry::CUBE)
         {
            mesh->GetElementFaces(ie, faces, ofaces);
            const int fs = GetHexFaceSplits(quad_diag, faces, ofaces);
            DrawRefinedHexLevelSurf(lsurf_buf, pointmat, vals, RG, nre, fs, gp);
         }
      }
   }

   updated_bufs.emplace_back(&lsurf_buf);

#ifdef GLVIS_DEBUG
   cout << "VisualizationSceneSolution3d::PrepareLevelSurf() : "
        << triangle_counter << " triangles + " << quad_counter
        << " quads used" << endl;
#endif
}

gl3::SceneInfo VisualizationSceneSolution3d::GetSceneObjs()
{
   if (colorbar)
   {
      Array<double>* cb_level = nullptr,
                     * cb_levels = nullptr;
      if (drawlsurf) { cb_levels = &levels; }
      if (drawmesh == 2 || cp_drawmesh >= 2) { cb_level = &level; }
      PrepareColorBar(minv, maxv, cb_level, cb_levels);
   }
   gl3::SceneInfo scene = VisualizationSceneScalarData::GetSceneObjs();
   gl3::RenderParams params = GetMeshDrawParams();
   params.use_clip_plane = cplane;
   double* cp_eqn = CuttingPlane->Equation();
   params.clip_plane_eqn = {cp_eqn[0], cp_eqn[1], cp_eqn[2], cp_eqn[3]};
   params.contains_translucent = matAlpha < 1.0 ||
                                 palette.GetPalette()->IsTranslucent();

   if (drawlsurf)
   {
      scene.queue.emplace_back(params, &lsurf_buf);
   }
   if (drawelems)
   {
      scene.queue.emplace_back(params, &disp_buf);
   }
   if (cplane && cp_drawelems)
   {
      params.use_clip_plane = false;
      scene.queue.emplace_back(params, &cplane_buf);
      params.use_clip_plane = true;
   }
   params.contains_translucent = false;
   // draw orderings -- color modes
   if (draworder == 1)
   {
      scene.queue.emplace_back(params, &order_noarrow_buf);
   }
   else if (draworder == 2)
   {
      scene.queue.emplace_back(params, &order_buf);
   }
   params.mesh_material = VisualizationScene::BLK_MAT;
   // everything below will be drawn in "black"
   params.static_color = GetLineColor();
   params.num_pt_lights = 0;
   if (drawmesh)
   {
      scene.queue.emplace_back(params, &line_buf);
   }
   if (cp_drawmesh)
   {
      params.use_clip_plane = false;
      scene.queue.emplace_back(params, &cplines_buf);
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
   ProcessUpdatedBufs(scene);
   return scene;
}

void VisualizationSceneSolution3d::glTF_Export()
{
   string name = "GLVis_scene_000";

   glTF_Builder bld(name);

   auto palette_mat = AddPaletteMaterial(bld);
   auto black_mat = AddBlackMaterial(bld);
   auto buf = bld.addBuffer("buffer");
   if (drawelems) { glTF_ExportElements(bld, buf, palette_mat, disp_buf); }
   if (drawmesh) { glTF_ExportMesh(bld, buf, black_mat, line_buf); }
   if (cplane && cp_drawelems)
   {
      auto cp_elems_node = AddModelNode(bld, "CP Elements");
      auto cp_elems_mesh = bld.addMesh("CP Elements Mesh");
      bld.addNodeMesh(cp_elems_node, cp_elems_mesh);

      int ntria = AddTriangles(
                     bld,
                     cp_elems_mesh,
                     buf,
                     palette_mat,
                     cplane_buf);
      if (ntria == 0)
      {
         cout << "glTF export: no cp elements found to export!" << endl;
      }
   }
   if (cp_drawmesh)
   {
      auto cp_lines_node = AddModelNode(bld, "CP Lines");
      auto cp_lines_mesh = bld.addMesh("CP Lines Mesh");
      bld.addNodeMesh(cp_lines_node, cp_lines_mesh);

      int nlines = AddLines(
                      bld,
                      cp_lines_mesh,
                      buf,
                      black_mat,
                      cplines_buf);
      if (nlines == 0)
      {
         cout << "glTF export: no cp mesh/level lines found to export!" << endl;
      }
   }
   if (drawlsurf)
   {
      auto lsurf_node = AddModelNode(bld, "Level Surface");
      auto lsurf_mesh = bld.addMesh("Level Surface Mesh");
      bld.addNodeMesh(lsurf_node, lsurf_mesh);

      int ntria = AddTriangles(
                     bld,
                     lsurf_mesh,
                     buf,
                     palette_mat,
                     lsurf_buf);
      if (ntria == 0)
      {
         cout << "glTF export: no level surface elements found to export!"
              << endl;
      }
   }
   if (drawaxes) { glTF_ExportBox(bld, buf, black_mat); }
   bld.writeFile();

   cout << "Exported glTF -> " << name << ".gltf" << endl;
}
