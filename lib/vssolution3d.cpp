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
#include <cmath>
#include <limits>

#include "mfem.hpp"
using namespace mfem;
#include "visual.hpp"
#include "palettes.hpp"
using namespace std;


VisualizationSceneSolution3d *vssol3d;
extern GeometryRefiner GLVisGeometryRefiner;

// Definitions of some more keys

static void Solution3dKeyHPressed()
{
   cout << endl
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
        << "| p/P  Cycle through color palettes  |" << endl
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
        << "| Ctrl+p - Print to a PDF file       |" << endl
        << "+------------------------------------+" << endl
        << "| Function keys                      |" << endl
        << "+------------------------------------+" << endl
        << "| F1 - X window info and keystrokes  |" << endl
        << "| F2 - Update colors, etc.           |" << endl
        << "| F3/F4 - Shrink/Zoom bdr elements   |" << endl
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
      cout << "Subdivision factor = " << ++vssol3d->TimesToRefine << endl;
      if (vssol3d -> GetShading() == 2)
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

static void KeyOPressed()
{
   if (vssol3d -> TimesToRefine > 1)
   {
      cout << "Subdivision factor = " << --vssol3d->TimesToRefine << endl;
      if (vssol3d -> GetShading() == 2)
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
   if (vssol3d -> GetShading() == 2)
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
   if (vssol3d -> GetShading() == 2)
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

static void KeyF3Pressed()
{
   if (vssol3d->GetShading() == 2)
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
      vssol3d->Prepare();
      vssol3d->PrepareLines();
      SendExposeEvent();
   }
}

static void KeyF4Pressed()
{
   if (vssol3d->GetShading() == 2)
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

static void KeyF11Pressed()
{
   if (vssol3d->GetShading() == 2)
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
   if (vssol3d->GetShading() == 2)
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

// 'v' must define three vectors that add up to the zero vector
// that is v[0], v[1], and v[2] are the sides of a triangle
int UnitCrossProd(double v[][3], double nor[])
{
   // normalize the three vectors
   for (int i = 0; i < 3; i++)
      if (Normalize(v[i]))
      {
         return 1;
      }

   // find the pair that forms an angle closest to pi/2, i.e. having
   // the longest cross product:
   double cp[3][3], max_a = 0.;
   int k = 0;
   for (int i = 0; i < 3; i++)
   {
      CrossProd(v[(i+1)%3], v[(i+2)%3], cp[i]);
      double a = sqrt(InnerProd(cp[i], cp[i]));
      if (max_a < a)
      {
         max_a = a, k = i;
      }
   }
   if (max_a == 0.)
   {
      return 1;
   }
   for (int i = 0; i < 3; i++)
   {
      nor[i] = cp[k][i] / max_a;
   }

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

   wnd->setOnKeyDown('h', Solution3dKeyHPressed);
   wnd->setOnKeyDown('H', Solution3dKeyHPressed);
}


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

   paletteSet(12); // use the 'vivid' palette in 3D
   SetUseTexture(1);

   double eps = 1e-6; // move the cutting plane a bit to avoid artifacts
   CuttingPlane = new Plane(-1.0,0.0,0.0,(0.5-eps)*x[0]+(0.5+eps)*x[1]);

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

      wnd->setOnKeyDown('i', KeyIPressed);
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
}

VisualizationSceneSolution3d::~VisualizationSceneSolution3d()
{
   delete [] node_pos;
}

void VisualizationSceneSolution3d::NewMeshAndSolution(
   Mesh *new_m, Vector *new_sol, GridFunction *new_u)
{
   if (mesh->GetNV() != new_m->GetNV())
   {
      delete [] node_pos;
      node_pos = new double[new_m->GetNV()];
   }
   // If the number of surface elements changes, recompute the refinement factor
   if (mesh->Dimension() != new_m->Dimension() ||
       (mesh->Dimension() == 2 && mesh->GetNE() != new_m->GetNE()) ||
       (mesh->Dimension() == 3 && mesh->GetNBE() != new_m->GetNBE()))
   {
      mesh = new_m;
      int ref = GetAutoRefineFactor();
      if (TimesToRefine != ref)
      {
         TimesToRefine = ref;
         cout << "Subdivision factor = " << TimesToRefine << endl;
      }
   }
   mesh = new_m;
   sol = new_sol;
   GridF = new_u;
   FindNodePos();

   DoAutoscale(false);

   Prepare();
   PrepareLines();
   CPPrepare();
   PrepareLevelSurf();
}

void VisualizationSceneSolution3d::SetShading(int s, bool print)
{
   if (shading == s || s < 0)
   {
      return;
   }

   if (s > 2 || (GridF == NULL && s > 1))
   {
      return;
   }
   int os = shading;
   shading = s;
   if (GridF != NULL && (s == 2 || os == 2))
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
      cout << "Shading type : " << shading_type[shading] << endl;
   }
}

void VisualizationSceneSolution3d::ToggleShading()
{
   if (GridF)
   {
      SetShading((shading+1)%3, true);
   }
   else
   {
      SetShading(1-shading, true);
   }
}

void VisualizationSceneSolution3d::SetRefineFactors(int f, int ignored)
{
   if (TimesToRefine == f || f < 1)
   {
      return;
   }

   TimesToRefine = f;

   if (shading == 2)
   {
      DoAutoscale(false);
      Prepare();
      PrepareLines();
      CPPrepare();
   }
}

int VisualizationSceneSolution3d::GetAutoRefineFactor()
{
   int ne = mesh->GetNBE(), ref = 1;
   if (mesh->Dimension() == 2)
   {
      ne = mesh->GetNE();
   }

   while (ref < auto_ref_max && ne*(ref+1)*(ref+1) <= auto_ref_max_surf_elem)
   {
      ref++;
   }

   return ref;
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

   x[0] = x[1] = coord[0];
   y[0] = y[1] = coord[1];
   z[0] = z[1] = coord[2];

   for (int i = 1; i < nv; i++)
   {
      coord = mesh->GetVertex(i);
      if (coord[0] < x[0]) { x[0] = coord[0]; }
      if (coord[1] < y[0]) { y[0] = coord[1]; }
      if (coord[2] < z[0]) { z[0] = coord[2]; }
      if (coord[0] > x[1]) { x[1] = coord[0]; }
      if (coord[1] > y[1]) { y[1] = coord[1]; }
      if (coord[2] > z[1]) { z[1] = coord[2]; }
   }

   if (shading == 2)
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
            RefG = GLVisGeometryRefiner.Refine(mesh->GetFaceBaseGeometry(fn),
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
            if (pointmat(0,j) < x[0]) { x[0] = pointmat(0,j); }
            if (pointmat(1,j) < y[0]) { y[0] = pointmat(1,j); }
            if (pointmat(2,j) < z[0]) { z[0] = pointmat(2,j); }
            if (pointmat(0,j) > x[1]) { x[1] = pointmat(0,j); }
            if (pointmat(1,j) > y[1]) { y[1] = pointmat(1,j); }
            if (pointmat(2,j) > z[1]) { z[1] = pointmat(2,j); }
         }
      }
   }

   UpdateBoundingBox();
}

void VisualizationSceneSolution3d::FindNewValueRange(bool prepare)
{
   if (shading < 2)
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
   if (shading == 2 && drawmesh != 0 && FaceShiftScale != 0.0)
   {
      PrepareLines();
   }
}

void VisualizationSceneSolution3d::UpdateValueRange(bool prepare)
{
   logscale = logscale && LogscaleRange();
   MySetColorLogscale = logscale;
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

int Normalize(DenseMatrix &normals)
{
   int err = 0;
   for (int i = 0; i < normals.Width(); i++)
   {
      err += Normalize(&normals(0, i));
   }
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

   double bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                             (y[1]-y[0])*(y[1]-y[0]) +
                             (z[1]-z[0])*(z[1]-z[0]) );
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
          DrawTriangle(disp_buf, p, c, minv, maxv);
      }
      else
      {
          DrawQuad(disp_buf, p, c, minv, maxv);
      }
   }
   disp_buf.buffer();
}

void VisualizationSceneSolution3d::PrepareFlat2()
{
   int i, k, fn, fo, di, have_normals;
   double bbox_diam, vmin, vmax;
   disp_buf.clear();

   int dim = mesh->Dimension();
   int nbe = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat, normals;
   Vector values, normal;
   RefinedGeometry * RefG;
   Array<int> vertices;
   double norm[3];
   IsoparametricTransformation T;

   bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                      (y[1]-y[0])*(y[1]-y[0]) +
                      (z[1]-z[0])*(z[1]-z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   vmin = numeric_limits<double>::infinity();
   vmax = -vmin;
   for (i = 0; i < nbe; i++) {
      if (dim == 3) {
         if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) { continue; }

         if (cplane == 2) {
            // for cplane == 2, get vertices of the volume element, not bdr
            int f, o, e1, e2;
            mesh->GetBdrElementFace(i, &f, &o);
            mesh->GetFaceElements(f, &e1, &e2);
            mesh->GetElementVertices(e1, vertices);
         } else {
            mesh->GetBdrElementVertices(i, vertices);
         }
      } else {
         if (!bdr_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }
         mesh->GetElementVertices(i, vertices);
      }
      if (cplane == 2 && CheckPositions(vertices)) { continue; }
      if (dim == 3) {
         mesh -> GetBdrElementFace (i, &fn, &fo);
         RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceBaseGeometry (fn),
                                            TimesToRefine);
         // di = GridF -> GetFaceValues (fn, 2, RefG->RefPts, values, pointmat);

         // this assumes the interior boundary faces are properly oriented ...
         di = fo % 2;
         if (di == 1 && !mesh->FaceIsInterior(fn)) {
            di = 0;
         }
         GridF -> GetFaceValues (fn, di, RefG->RefPts, values, pointmat);
         GetFaceNormals(fn, di, RefG->RefPts, normals);
         have_normals = 1;
         ShrinkPoints(pointmat, i, fn, di);
      } else {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine);
         const IntegrationRule &ir = RefG->RefPts;
         GridF->GetValues(i, ir, values, pointmat);
         normals.SetSize(3, values.Size());
         mesh->GetElementTransformation(i, &T);
         for (int j = 0; j < values.Size(); j++) {
            T.SetIntPoint(&ir.IntPoint(j));
            const DenseMatrix &J = T.Jacobian();
            normals.GetColumnReference(j, normal);
            CalcOrtho(J, normal);
            normal /= normal.Norml2();
         }
         have_normals = 1;
         di = 0;
         ShrinkPoints(pointmat, i, 0, 0);
      }

      vmin = fmin(vmin, values.Min());
      vmax = fmax(vmax, values.Max());

      int sides;
      switch ((dim == 3) ? mesh->GetBdrElementType(i) : mesh->GetElementType(i))
      {
         case Element::TRIANGLE:
            sides = 3;
            break;

         case Element::QUADRILATERAL:
         default:
            sides = 4;
            break;
      }

      // compute an average normal direction for the current face
      if (sc != 0.0)
      {
         for (int i = 0; i < 3; i++)
         {
            norm[i] = 0.0;
         }
         Normalize(normals);
         for (k = 0; k < normals.Width(); k++)
            for (int j = 0; j < 3; j++)
            {
               norm[j] += normals(j, k);
            }
         Normalize(norm);
         for (int i = 0; i < pointmat.Width(); i++)
         {
            double val = sc * (values(i) - minv) / (maxv - minv);
            for (int j = 0; j < 3; j++)
            {
               pointmat(j, i) += val * norm[j];
            }
         }
         have_normals = 0;
      }

      have_normals = have_normals ? 2 : 0;
      if (di)
      {
         have_normals = -1 - have_normals;
      }
      // Comment the above lines and use the below version in order to remove
      // the 3D dark artifacts (indicating wrong boundary element orientation)
      // have_normals = have_normals ? 1 : 0;
      DrawPatch(disp_buf, pointmat, values, normals, sides, RefG->RefGeoms,
                minv, maxv, have_normals);
   }
   disp_buf.buffer();

   cout << "VisualizationSceneSolution3d::PrepareFlat2() : [min,max] = ["
        << vmin << "," << vmax << "]" << endl;
}

void VisualizationSceneSolution3d::Prepare()
{
   int i,j;

   if (!drawelems)
   {
      return;
   }

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

      for (i = 0; i < nelem; i++)
      {
         if (dim == 3)
         {
            mesh->GetBdrElementVertices(elem[i], vertices);
         }
         else
         {
            mesh->GetElementVertices(elem[i], vertices);
         }
         for (j = 0; j < vertices.Size(); j++)
         {
            nx(vertices[j]) = ny(vertices[j]) = nz(vertices[j]) = 0.;
         }
      }
      for (i = 0; i < nelem; i++)
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
         
         GLenum elemType;
         switch ((dim == 3) ? mesh->GetBdrElementType(elem[i]) :
                 mesh->GetElementType(elem[i]))
         {
            case Element::TRIANGLE:
               elemType = GL_TRIANGLES;
               break;

            case Element::QUADRILATERAL:
               elemType = GL_QUADS;
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
            poly.glNormal3d(nx(vertices[j]), ny(vertices[j]), nz(vertices[j]));
            poly.glVertex3dv(&pointmat(0, j));
         }
         poly.glEnd();
      }
   }
   disp_buf.buffer();
}

void VisualizationSceneSolution3d::PrepareLines()
{
   if (!drawmesh)
   {
      return;
   }

   if (shading == 2)
   {
      PrepareLines2();
      return;
   }

   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   int i, j, k;
   DenseMatrix pointmat;

   line_buf.clear();

   Array<int> vertices;

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
            line.glBegin(GL_LINE_LOOP);

            for (j = 0; j < pointmat.Size(); j++)
            {
               line.glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j));
            }
            line.glEnd();
            break;

         case 2:
            for (j = 0; j < pointmat.Size(); j++)
            {
               for (k = 0; k < 3; k++)
               {
                  point[j][k] = pointmat(k,j);
               }
               point[j][3] = (*sol)(vertices[j]);
            }
            DrawPolygonLevelLines(line, point[0], pointmat.Size(), level, false);
            break;
      }
   }
   line_buf.buffer();
}

void VisualizationSceneSolution3d::PrepareLines2()
{
   int i, j, k, fn, fo, di = 0;
   double bbox_diam;

   line_buf.clear();

   int dim = mesh->Dimension();
   int nbe = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat, normals;
   Vector values, normal;
   RefinedGeometry * RefG;
   Array<int> vertices;
   IsoparametricTransformation T;

   bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                      (y[1]-y[0])*(y[1]-y[0]) +
                      (z[1]-z[0])*(z[1]-z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   for (i = 0; i < nbe; i++)
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
         RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceBaseGeometry (fn),
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
         for (int i = 0; i < 3; i++)
         {
            norm[i] = 0.0;
         }
         for (k = 0; k < normals.Width(); k++)
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
         for (int i = 0; i < 3; i++)
         {
            norm[i] *= len;
         }
         for (int i = 0; i < pointmat.Width(); i++)
         {
            double val = sc * (values(i) - minv) / (maxv - minv);
            for (int j = 0; j < 3; j++)
            {
               pointmat(j, i) += val * norm[j];
            }
         }
      }
      gl3::GlBuilder line = line_buf.createBuilder();
      if (drawmesh == 1)
      {
         Array<int> &REdges = RefG->RefEdges;

         line.glBegin(GL_LINES);
         for (k = 0; k < REdges.Size(); k++)
         {
            line.glVertex3dv(&pointmat(0, REdges[k]));
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
         for (k = 0; k < RefG->RefGeoms.Size()/sides; k++)
         {
            int *RG = &(RefG->RefGeoms[k*sides]);

            for (j = 0; j < sides; j++)
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
   line_buf.buffer();
}

void VisualizationSceneSolution3d::CuttingPlaneFunc(int func)
{
   int i, j, k, m, n, n2;
   int flag[8], oedges[6];
   static const int tet_edges[12]= {0,3, 0,2, 0,1, 1,2, 1,3, 2,3};
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
   double t, point[6][4], norm[3];

   DenseMatrix pointmat;

   Array<int> nodes;
   for (i = 0; i < mesh -> GetNE(); i++)
   {
      n = n2 = 0; // n will be the number of intersection points
      mesh -> GetElementVertices(i, nodes);
      for (j = 0; j < nodes.Size(); j++)
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
      switch (mesh -> GetElementType(i))
      {
         case Element::TETRAHEDRON:
            ev = tet_edges;
            for (j = 0; j < 6; j++, ev += 2)
               if (flag[ev[0]] != flag[ev[1]])
               {
                  oedges[n++] = j;
               }
            ev = tet_edges;
            break;
         case Element::HEXAHEDRON:
            int emark[12];
            ev = hex_edges;
            for (j = 0; j < 12; j++, ev += 2)
            {
               emark[j] = flag[ev[1]] - flag[ev[0]];
            }
            do
            {
               for (j = 0; j < 12; j++)
               {
                  if (emark[j]) { break; }
               }
               if (j == 12)
               {
                  break;
               }
               k = 2 * j;
               if (emark[j] > 0)
               {
                  k++;
               }
               do
               {
                  for (j = 0; j < 3; j++)
                  {
                     m = hex_cutting[k][j];
                     ev = hex_edges + 2 * (m / 2);
                     if ((m % 2 == 0 && flag[ev[0]] > flag[ev[1]]) ||
                         (m % 2 == 1 && flag[ev[1]] > flag[ev[0]]))
                     {
                        break;
                     }
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
         {
            mesh -> GetPointMatrix (i, pointmat);
         }
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
            t = node_pos[ nodes[en[1]] ];
            t = t / ( t - node_pos[ nodes[en[0]] ] );
            for (k = 0; k < 3; k++)
            {
               point[j][k] = t*pointmat(k,en[0]) + (1-t)*pointmat(k,en[1]);
            }
            point[j][3] = t*(*sol)(nodes[en[0]]) + (1-t)*(*sol)(nodes[en[1]]);
         }

         switch (func)
         {
            case 1:  // PrepareCuttingPlane()
            {
               if (shading == 2)
               {
                  // changes point for n > 4
                  DrawRefinedSurf(n, point[0], i, 1);
               }
               else
               {
                  m = n;
                  while (1)
                  {
                     if (m > 3)
                     {
                        j = Compute3DUnitNormal(point[0], point[1], point[2],
                                                point[3], norm);
                        if (j && m > 4)
                        {
                           for (int j = 3; j < m; j++)
                              for (int i = 0; i < 4; i++)
                              {
                                 point[j-2][i] = point[j][i];
                              }
                           m -= 2;
                           continue;
                        }
                     }
                     else
                        j = Compute3DUnitNormal(point[0], point[1], point[2],
                                                norm);
                     break;
                  }

                  gl3::GlBuilder draw = cplane_buf.createBuilder();
                  if (!j)
                  {
                     draw.glBegin(GL_POLYGON);
                     draw.glNormal3dv(norm);
                     for (j = 0; j < m; j++)
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
               if (shading == 2)
               {
                  // changes point for n > 4
                  DrawRefinedSurf(n, point[0], i, 2);
               }
               else
               {
                  // glBegin (GL_POLYGON);
                  gl3::GlBuilder line = cplines_buf.createBuilder();
                  line.glBegin(GL_LINE_LOOP);
                  for (j = 0; j < n; j++)
                  {
                     line.glVertex3dv(point[j]);
                  }
                  line.glEnd();
               }
            }
            break;

            case 3:  // PrepareCuttingPlaneLines() with level lines
            {
               if (shading == 2)
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

         for (j = 0; j < n2; j++)
         {
            oedges[j] = oedges[j+n];
         }
         n = n2;
         n2 = 0;
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
      else
      {
         CuttingPlaneFunc(1);
      }
   }
   cplane_buf.buffer();
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
            RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceBaseGeometry (i),
                                               TimesToRefine);
            // partition[e1] is 0 if e1 is behind the cutting plane
            // and 1 otherwise
            int dir = partition[e1];
            GridF -> GetFaceValues (i, dir, RefG->RefPts, values, pointmat);
            GetFaceNormals(i, dir, RefG->RefPts, normals);
            switch (mesh -> GetFaceBaseGeometry (i))
            {
               case Geometry::TRIANGLE:  n = 3; break;
               case Geometry::SQUARE:    n = 4; break;
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
         if (cp_drawmesh == 1)
         {
            CuttingPlaneFunc(2);
         }
         else
         {
            CuttingPlaneFunc(3);
         }
      }
   }

   cplines_buf.buffer();
}

void VisualizationSceneSolution3d::PrepareCuttingPlaneLines2()
{
   int i, j, n = 0;
   double point[4][4], *coord;
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
         if (shading != 2)
         {
            mesh -> GetFaceVertices (i, nodes);
            for (j = 0; j < nodes.Size(); j++)
            {
               coord = mesh -> GetVertex(nodes[j]);
               point[j][0] = coord[0];
               point[j][1] = coord[1];
               point[j][2] = coord[2];
               point[j][3] = (*sol)(nodes[j]);
            }
            gl3::GlBuilder line = cplines_buf.createBuilder();
            switch (cp_drawmesh)
            {
               case 1:
                  // glBegin(GL_POLYGON);
                  line.glBegin(GL_LINE_LOOP);
                  for (j = 0; j < nodes.Size(); j++)
                  {
                     line.glVertex3dv (point[j]);
                  }
                  line.glEnd();
                  break;
               case 2:
                  DrawPolygonLevelLines(line, point[0], nodes.Size(), level, false);
                  break;
            }
         }
         else // shading == 2
         {
            RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceBaseGeometry (i),
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
}

int triangle_counter;
int quad_counter;

void VisualizationSceneSolution3d::DrawTetLevelSurf(
   gl3::GlDrawable& target,
   const DenseMatrix &verts, const Vector &vals, const int *ind,
   const Array<double> &levels, const DenseMatrix *grad)
{
   double t, lvl, normal[3], vert[4][3], norm[4][3];
   int i, j, l, pos[4];
   bool flipped;

   gl3::GlBuilder draw = target.createBuilder();

   for (l = 0; l < levels.Size(); l++)
   {
      lvl = levels[l];

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
               draw.glNormal3dv(norm[k]);
               draw.glVertex3dv(vert[k]);
            }
            draw.glEnd();
            quad_counter++;
         }
      }
   }
}

void VisualizationSceneSolution3d::PrepareLevelSurf()
{
   static const int tet_id[4] = { 0, 1, 2, 3 };
   static const int hex_tets[6][4] =
   {
      { 0, 1, 2, 6 }, { 0, 5, 1, 6 }, { 0, 4, 5, 6 },
      { 0, 2, 3, 6 }, { 0, 3, 7, 6 }, { 0, 7, 4, 6 }
   };

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

   if (shading != 2)
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

         if (mesh->GetElementType(ie) == Element::TETRAHEDRON)
         {
            DrawTetLevelSurf(lsurf_buf, pointmat, vals, tet_id, levels);
         }
         else if (mesh->GetElementType(ie) == Element::HEXAHEDRON)
         {
            for (int k = 0; k < 6; k++)
            {
               DrawTetLevelSurf(lsurf_buf, pointmat, vals, hex_tets[k], levels);
            }
         }
      }
   }
   else // shading == 2
   {
      RefinedGeometry *RefG;

      for (int ie = 0; ie < mesh->GetNE(); ie++)
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(ie),
                                            TimesToRefine);
         GridF->GetValues(ie, RefG->RefPts, vals, pointmat);
#define GLVIS_SMOOTH_LEVELSURF_NORMALS
#ifdef GLVIS_SMOOTH_LEVELSURF_NORMALS
         GridF->GetGradients(ie, RefG->RefPts, grad);
#endif

         Array<int> &RG = RefG->RefGeoms;
         int nv = mesh->GetElement(ie)->GetNVertices();

         for (int k = 0; k < RG.Size()/nv; k++)
         {
            if (nv == 4)
            {
#ifndef GLVIS_SMOOTH_LEVELSURF_NORMALS
               DrawTetLevelSurf(lsurf_buf, pointmat, vals, &RG[nv*k], levels);
#else
               DrawTetLevelSurf(lsurf_buf, pointmat, vals, &RG[nv*k], levels, &grad);
#endif
            }
            else if (nv == 8)
            {
               int m_ind[4];
               for (int j = 0; j < 6; j++)
               {
                  for (int i = 0; i < 4; i++)
                  {
                     m_ind[i] = RG[nv*k+hex_tets[j][i]];
                  }
#ifndef GLVIS_SMOOTH_LEVELSURF_NORMALS
                  DrawTetLevelSurf(lsurf_buf, pointmat, vals, m_ind, levels);
#else
                  DrawTetLevelSurf(lsurf_buf, pointmat, vals, m_ind, levels, &grad);
#endif
               }
            }
         }
      }
   }

   lsurf_buf.buffer();

#ifdef GLVIS_DEBUG
   cout << "VisualizationSceneSolution3d::PrepareLevelSurf() : "
        << triangle_counter << " triangles + " << quad_counter
        << " quads used" << endl;
#endif
}

void VisualizationSceneSolution3d::Draw()
{
   gl->enableDepthTest();

   Set_Background();
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   // model transformation
   ModelView();

   // draw colored faces
   glPolygonOffset (1, 1);
   glEnable (GL_POLYGON_OFFSET_FILL);
   //glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   gl->disableClipPlane();
   // draw colorbar
   gl->disableLight();
   if (colorbar)
   {
      if (drawmesh == 2 || cp_drawmesh >= 2)
      {
         if (drawlsurf)
         {
            DrawColorBar(minv,maxv,&level,&levels);
         }
         else
         {
            DrawColorBar(minv,maxv,&level);
         }
      }
      else
      {
         if (drawlsurf)
         {
            DrawColorBar(minv,maxv,NULL,&levels);
         }
         else
         {
            DrawColorBar(minv,maxv);
         }
      }
   }

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

   if (drawlsurf)
   {
      lsurf_buf.draw();
      // Set_Black_Material();
      // glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
   }

   // draw elements
   if (drawelems)
   {
       disp_buf.draw();
   }

   if (cplane && cp_drawelems)
   {
      gl->disableClipPlane();
      cplane_buf.draw();
      gl->enableClipPlane();
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
   if (cplane)
   {
      gl->disableClipPlane();
      DrawRuler();
      if (cp_drawmesh)
      {
         cplines_buf.draw();
      }
      gl->enableClipPlane();
   }
   else
   {
      DrawRuler();
   }
   // draw lines
   if (drawmesh)
   {
      line_buf.draw();
   }

   if (cplane)
   {
      gl->disableClipPlane();
   }

   // draw axes
   if (drawaxes)
   {
      axes_buf.draw();
      DrawCoordinateCross();
   }
}
