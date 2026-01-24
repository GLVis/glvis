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

#include "vsvector3d.hpp"

using namespace mfem;
using namespace std;

// Reference geometry with a cut in the middle, which subdivides GeometryRefiner
// when cut_lambda is updated, see keys Ctrl+F3/F4. These variables are defined
// in lib/vssolution3d.cpp.
extern thread_local IntegrationRule cut_QuadPts;
extern thread_local Array<int> cut_QuadGeoms;
extern thread_local IntegrationRule cut_TriPts;
extern thread_local Array<int> cut_TriGeoms;
extern void CutReferenceElements(int TimesToRefine, double lambda);


std::string VisualizationSceneVector3d::GetHelpString() const
{
   std::stringstream os;
   os << endl
      << "+------------------------------------+" << endl
      << "| Keys                               |" << endl
      << "+------------------------------------+" << endl
      << "| a -  Displays/Hides the axes       |" << endl
      << "| A -  Turns antialiasing on/off     |" << endl
      << "| b -  Displacements step back       |" << endl
      << "| c -  Toggle colorbar and caption   |" << endl
      << "| C -  Change the main plot caption  |" << endl
      << "| d -  Displays/Hides displacements  |" << endl
      << "| e -  Displays/Hides the elements   |" << endl
      << "| E -  Toggle the elements in the CP |" << endl
      << "| f -  Smooth/Flat shading           |" << endl
      << "| F -  Display mag./x/y/z component  |" << endl
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
      << "| n -  Displacements step forward    |" << endl
      << "| o/O  (De)refine elem, disc shading |" << endl
      << "| p/P  Cycle through color palettes  |" << endl
      << "| q -  Quits                         |" << endl
      << "| Q -  Cycle quadrature data mode    |" << endl
      << "| r -  Reset the plot to 3D view     |" << endl
      << "| R -  Reset the plot to 2D view     |" << endl
      << "| s -  Turn on/off unit cube scaling |" << endl
      << "| S -  Take snapshot/Record a movie  |" << endl
      << "| t -  Cycle materials and lights    |" << endl
      << "| u/U  Move the level field vectors  |" << endl
      << "| v/V  Vector field                  |" << endl
      << "| w/W  Add/Delete level field vector |" << endl
      << "| x/X  Rotate cutting plane (phi)    |" << endl
      << "| y/Y  Rotate cutting plane (theta)  |" << endl
      << "| z/Z  Translate cutting plane       |" << endl
      << "| \\ -  Set light source position     |" << endl
      << "| Alt+a  - Axes number format        |" << endl
      << "| Alt+c  - Colorbar number format    |" << endl
      << "| Ctrl+p - Print to a PDF file       |" << endl
      << "+------------------------------------+" << endl
      << "| Function keys                      |" << endl
      << "+------------------------------------+" << endl
      << "| F1 - X window info and keystrokes  |" << endl
      << "| F2 - Update colors, etc.           |" << endl
      << "| F3/F4 - Shrink/Zoom bdr elements   |" << endl
      << "| Ctrl+F3/F4 - Cut face bdr elements |" << endl
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

thread_local VisualizationSceneVector3d  *vsvector3d;
extern thread_local VisualizationScene *locscene;
extern thread_local GeometryRefiner GLVisGeometryRefiner;

static void KeyDPressed()
{
   vsvector3d -> ToggleDisplacements();
   SendExposeEvent();
}

static void KeyNPressed()
{
   if (vsvector3d -> drawdisp)
      vsvector3d -> ianimd =( (vsvector3d ->ianimd + 1) %
                              (vsvector3d -> ianimmax + 1) );
   else
      vsvector3d -> ianim =( (vsvector3d ->ianim + 1) %
                             (vsvector3d -> ianimmax + 1) );
   vsvector3d -> NPressed();
}

static void KeyBPressed()
{
   if (vsvector3d -> drawdisp)
      vsvector3d ->ianimd = ((vsvector3d ->ianimd +
                              vsvector3d ->ianimmax) %
                             (vsvector3d ->ianimmax + 1));
   else
      vsvector3d ->ianim = ((vsvector3d ->ianim +
                             vsvector3d ->ianimmax) %
                            (vsvector3d ->ianimmax + 1));
   vsvector3d -> NPressed();
}

static void KeyrPressed()
{
   locscene -> spinning = 0;
   RemoveIdleFunc(MainLoop);
   vsvector3d -> CenterObject();
   locscene -> ViewAngle = 45.0;
   locscene -> ViewScale = 1.0;
   locscene -> ViewCenterX = 0.0;
   locscene -> ViewCenterY = 0.0;
   vsvector3d -> ianim = vsvector3d -> ianimd = 0;
   vsvector3d -> Prepare();
   vsvector3d -> PrepareLines();
   vsvector3d -> PrepareDisplacedMesh();
   vsvector3d -> key_r_state = 0;
   SendExposeEvent();
}

static void KeyRPressed()
{
   locscene->spinning = 0;
   RemoveIdleFunc(MainLoop);
   vsvector3d -> ianim = vsvector3d -> ianimd = 0;
   vsvector3d -> Prepare();
   vsvector3d -> PrepareLines();
   vsvector3d -> PrepareDisplacedMesh();
   vsvector3d->Toggle2DView();
   SendExposeEvent();
}

void VisualizationSceneVector3d::NPressed()
{
   if (drawdisp)
   {
      PrepareDisplacedMesh();
   }
   else
   {
      Prepare();
      PrepareLines();
   }

   SendExposeEvent();
}

static void KeyuPressed()
{
   vsvector3d -> ToggleVectorFieldLevel(+1);
   SendExposeEvent();
}

static void KeyUPressed()
{
   vsvector3d -> ToggleVectorFieldLevel(-1);
   SendExposeEvent();
}

void VisualizationSceneVector3d::ToggleVectorFieldLevel(int v)
{
   int i;
   for (i = 0; i < vflevel.Size(); i++)
      if (vflevel[i] == 0 && v == -1)
      {
         vflevel[i] = nl;
      }
      else
      {
         vflevel[i] = (vflevel[i] + v) % (nl+1);
      }
   for (i = 0; i < vflevel.Size(); i++)
   {
      dvflevel[i] = level[vflevel[i]];
   }
   vsvector3d -> PrepareVectorField();
}

static void KeywPressed()
{
   vsvector3d -> AddVectorFieldLevel();
   SendExposeEvent();
}

static void KeyWPressed()
{
   vsvector3d -> RemoveVectorFieldLevel();
   SendExposeEvent();
}

void VisualizationSceneVector3d::AddVectorFieldLevel()
{
   int next = vflevel[vflevel.Size()-1];
   next = (next + 1) % (nl+1);
   vflevel.Append(next);
   dvflevel.Append(level[next]);
   vsvector3d -> PrepareVectorField();
}

void VisualizationSceneVector3d::RemoveVectorFieldLevel()
{
   vflevel.DeleteLast();
   dvflevel.DeleteLast();
   vsvector3d -> PrepareVectorField();
}

static void KeyvPressed()
{
   vsvector3d -> ToggleVectorField(1);
   SendExposeEvent();
}

static void KeyVPressed()
{
   vsvector3d -> ToggleVectorField(-1);
   SendExposeEvent();
}

static void VectorKeyFPressed()
{
   vsvector3d->ToggleScalarFunction();
   SendExposeEvent();
}

void VisualizationSceneVector3d::ToggleVectorField(int i)
{
   drawvector = (drawvector+i+6)%6;
   PrepareVectorField();
}

static const char *scal_func_name[] =
{"magnitude", "x-component", "y-component", "z-component"};

void VisualizationSceneVector3d::SetScalarFunction()
{
   FiniteElementSpace *fes = (VecGridF) ? VecGridF->FESpace() : NULL;

   switch (scal_func)
   {
      case 0: // magnitude
         for (int i = 0; i < sol->Size(); i++)
            (*sol)(i) = sqrt((*solx)(i) * (*solx)(i) +
                             (*soly)(i) * (*soly)(i) +
                             (*solz)(i) * (*solz)(i) );
         if (GridF)
         {
            Array<int> dofs(3);
            for (int i = 0; i < GridF->Size(); i++)
            {
               dofs.SetSize(1);
               dofs[0] = i;
               fes->DofsToVDofs(dofs);
               double x = (*VecGridF)(dofs[0]);
               double y = (*VecGridF)(dofs[1]);
               double z = (*VecGridF)(dofs[2]);

               (*GridF)(i) = sqrt(x*x+y*y+z*z);
            }
         }
         break;
      case 1: // x-component
         *sol = *solx;
         if (GridF)
            for (int i = 0; i < GridF->Size(); i++)
            {
               (*GridF)(i) = (*VecGridF)(fes->DofToVDof(i, 0));
            }
         break;
      case 2: // y-component
         *sol = *soly;
         if (GridF)
            for (int i = 0; i < GridF->Size(); i++)
            {
               (*GridF)(i) = (*VecGridF)(fes->DofToVDof(i, 1));
            }
         break;
      case 3: // z-component
         *sol = *solz;
         if (GridF)
            for (int i = 0; i < GridF->Size(); i++)
            {
               (*GridF)(i) = (*VecGridF)(fes->DofToVDof(i, 2));
            }
         break;
   }
   extra_caption = scal_func_name[scal_func];
}

void VisualizationSceneVector3d::ToggleScalarFunction()
{
   scal_func = (scal_func + 1) % 4;
   cout << "Displaying " << scal_func_name[scal_func] << endl;
   SetScalarFunction();
   FindNewValueRange(true);
}

VisualizationSceneVector3d::VisualizationSceneVector3d(Mesh &m, Vector &sx,
                                                       Vector &sy, Vector &sz, Mesh *mc)
{
   mesh = &m;
   mesh_coarse = mc;
   solx = &sx;
   soly = &sy;
   solz = &sz;

   sol = new Vector(mesh->GetNV());

   sfes = NULL;
   VecGridF = NULL;

   Init();
}

VisualizationSceneVector3d::VisualizationSceneVector3d(GridFunction &vgf,
                                                       Mesh *mc)
{
   FiniteElementSpace *fes = vgf.FESpace();
   if (fes == NULL || fes->GetVDim() != 3)
   {
      cout << "VisualizationSceneVector3d::VisualizationSceneVector3d" << endl;
      exit(1);
   }

   VecGridF = &vgf;

   mesh = fes->GetMesh();
   mesh_coarse = mc;

   sfes = new FiniteElementSpace(mesh, fes->FEColl(), 1, fes->GetOrdering());
   GridF = new GridFunction(sfes);

   solx = new Vector(mesh->GetNV());
   soly = new Vector(mesh->GetNV());
   solz = new Vector(mesh->GetNV());

   vgf.GetNodalValues(*solx, 1);
   vgf.GetNodalValues(*soly, 2);
   vgf.GetNodalValues(*solz, 3);

   sol = new Vector(mesh->GetNV());

   Init();
}

void VisualizationSceneVector3d::Init()
{
   key_r_state = 0;

   drawdisp = 0;
   drawvector = 0;
   scal_func = 0;

   ianim = ianimd = 0;
   ianimmax = 10;

   SetScalarFunction();

   VisualizationSceneSolution3d::Init();

   mesh_volume = 0.0;
   if (mesh)
   {
      for (int i=0; i<mesh->GetNE(); i++)
      {
         mesh_volume += mesh->GetElementVolume(i);
      }
   }

   PrepareVectorField();
   PrepareDisplacedMesh();

   vflevel.Append(0);
   dvflevel.Append(level[0]);

   vsvector3d = this;

   // static int init = 0;
   // if (!init)
   {
      // init = 1;

      wnd->setOnKeyDown('d', KeyDPressed);
      wnd->setOnKeyDown('D', KeyDPressed);

      wnd->setOnKeyDown('n', KeyNPressed);
      wnd->setOnKeyDown('N', KeyNPressed);

      wnd->setOnKeyDown('b', KeyBPressed);
      wnd->setOnKeyDown('B', KeyBPressed);

      wnd->setOnKeyDown('r', KeyrPressed); // adds another function to 'r' and 'R'
      wnd->setOnKeyDown('R', KeyRPressed); // the other function is in vsdata.cpp

      wnd->setOnKeyDown('u', KeyuPressed); // Keys u, U are also used in
      wnd->setOnKeyDown('U', KeyUPressed); // VisualizationSceneSolution3d

      wnd->setOnKeyDown('w', KeywPressed); // Keys w, W are also used in
      wnd->setOnKeyDown('W', KeyWPressed); // VisualizationSceneSolution3d

      wnd->setOnKeyDown('v', KeyvPressed); // Keys v, V are also used in
      wnd->setOnKeyDown('V', KeyVPressed); // VisualizationSceneSolution3d

      wnd->setOnKeyDown('F', VectorKeyFPressed);
   }
}

int VisualizationSceneVector3d::GetFunctionAutoRefineFactor()
{
   if (!VecGridF) { return 1; }

   return VisualizationSceneScalarData::GetFunctionAutoRefineFactor(*VecGridF);
}

VisualizationSceneVector3d::~VisualizationSceneVector3d()
{
   delete sol;

   if (VecGridF)
   {
      delete solz;
      delete soly;
      delete solx;
      delete GridF;
      delete sfes;
   }
}

void VisualizationSceneVector3d::NewMeshAndSolution(
   Mesh *new_m, Mesh *new_mc, GridFunction *new_v)
{
   delete sol;
   if (VecGridF)
   {
      delete solz;
      delete soly;
      delete solx;
      delete GridF;
      delete sfes;
   }
   if (mesh->GetNV() != new_m->GetNV())
   {
      delete [] node_pos;
      node_pos = new double[new_m->GetNV()];
   }

   Mesh *old_m = mesh;
   VecGridF = new_v;
   mesh = new_m;
   mesh_coarse = new_mc;

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

   FiniteElementSpace *new_fes = new_v->FESpace();

   FindNodePos();

   sfes = new FiniteElementSpace(mesh, new_fes->FEColl(), 1,
                                 new_fes->GetOrdering());
   GridF = new GridFunction(sfes);

   solx = new Vector(mesh->GetNV());
   soly = new Vector(mesh->GetNV());
   solz = new Vector(mesh->GetNV());

   VecGridF->GetNodalValues(*solx, 1);
   VecGridF->GetNodalValues(*soly, 2);
   VecGridF->GetNodalValues(*solz, 3);

   sol = new Vector(mesh->GetNV());

   SetScalarFunction();

   DoAutoscale(false);

   Prepare();
   PrepareLines();
   CPPrepare();
   PrepareLevelSurf();

   PrepareVectorField();
   PrepareDisplacedMesh();
}

void VisualizationSceneVector3d::PrepareFlat()
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

         mesh->GetBdrPointMatrix(i, pointmat);
         mesh->GetBdrElementVertices(i, vertices);
      }
      else
      {
         if (!bdr_attr_to_show[mesh->GetAttribute(i)-1]) { continue; }

         mesh->GetPointMatrix(i, pointmat);
         mesh->GetElementVertices(i, vertices);
      }

      for (j = 0; j < pointmat.Width(); j++)
      {
         pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
         pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
         pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
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
         mfem_error("VisualizationSceneVector3d::PrepareFlat() :Unknown geometry.");
      }
   }
   updated_bufs.emplace_back(&disp_buf);
}

void VisualizationSceneVector3d::PrepareFlat2()
{
   int fn, fo, di = 0, have_normals = 0;
   double bbox_diam, vmin, vmax;
   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();

   DenseMatrix pointmat, normals;
   DenseMatrix vec_vals;
   Vector values, normal;
   RefinedGeometry * RefG;
   Array<int> vertices;
   double norm[3];
   IsoparametricTransformation T;

   bbox_diam = sqrt( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                     (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                     (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   disp_buf.clear();

   vmin = numeric_limits<double>::infinity();
   vmax = -vmin;
   for (int i = 0; i < ne; i++)
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

         // di = GridF->GetFaceValues(fn, 2, RefG->RefPts, values, pointmat);
         di = fo % 2;
         if (di == 1 && !mesh->FaceIsInterior(fn))
         {
            di = 0;
         }

         IntegrationRule &RefPts = (cut_lambda > 0) ?
                                   ((sides == 3) ? cut_TriPts : cut_QuadPts) :
                                   RefG->RefPts;
         GridF->GetFaceValues(fn, di, RefPts, values, pointmat);
         if (ianim > 0)
         {
            VecGridF->GetFaceVectorValues(fn, di, RefPts, vec_vals,
                                          pointmat);
            pointmat.Add(double(ianim)/ianimmax, vec_vals);
            have_normals = 0;
         }
         else
         {
            GetFaceNormals(fn, di, RefPts, normals);
            have_normals = 1;
         }
         ShrinkPoints(pointmat, i, fn, di);
      }
      else // dim < 3
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine);
         if (!cut_updated)
         {
            // Update the cut version of the reference geometries
            CutReferenceElements(TimesToRefine, cut_lambda);
            cut_updated = true;
         }
         IntegrationRule &RefPts = (cut_lambda > 0  && dim > 1) ?
                                   ((sides == 3) ? cut_TriPts : cut_QuadPts) :
                                   RefG->RefPts;
         GridF->GetValues(i, RefPts, values, pointmat);
         if (ianim > 0)
         {
            VecGridF->GetVectorValues(i, RefPts, vec_vals, pointmat);
            pointmat.Add(double(ianim)/ianimmax, vec_vals);
            have_normals = 0;
         }
         else
         {
            // Compute normals. Skip in 1D.
            if (dim > 1)
            {
               const IntegrationRule &ir = (cut_lambda > 0 && dim > 1) ?
                                           ((sides == 3) ? cut_TriPts : cut_QuadPts) :
                                           RefG->RefPts;
               normals.SetSize(3, values.Size());
               mesh->GetElementTransformation(i, &T);
               for (int j = 0; j < values.Size(); j++)
               {
                  T.SetIntPoint(&ir.IntPoint(j));
                  normals.GetColumnReference(j, normal);
                  CalcOrtho(T.Jacobian(), normal);
                  normal /= normal.Norml2();
               }
               have_normals = 1;
               di = 0;
            }
         }
         ShrinkPoints(pointmat, i, 0, 0);
      }

      vmin = fmin(vmin, values.Min());
      vmax = fmax(vmax, values.Max());

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

      Array<int> &RefGeoms = (cut_lambda > 0 && dim > 1) ?
                             ((sides == 3) ? cut_TriGeoms : cut_QuadGeoms) :
                             RefG->RefGeoms;
      int psides = (cut_lambda > 0) ? 4 : sides;
      if (dim == 1) { psides = 2; } // Hack to trigger line rendering.
      DrawPatch(disp_buf, pointmat, values, normals, psides, RefGeoms,
                minv, maxv, have_normals);
   }
   updated_bufs.emplace_back(&disp_buf);
   cout << "VisualizationSceneVector3d::PrepareFlat2() : [min,max] = ["
        << vmin << "," << vmax << "]" << endl;
}

void VisualizationSceneVector3d::Prepare()
{
   int i,j;

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
   gl3::GlBuilder draw = disp_buf.createBuilder();

   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   int nv = mesh -> GetNV();
   DenseMatrix pointmat;
   Array<int> vertices;
   double nor[3];

   Vector nx(nv);
   Vector ny(nv);
   Vector nz(nv);

   const Array<int> &attributes =
      ((dim == 3) ? mesh->bdr_attributes : mesh->attributes);
   for (int d = 0; d < attributes.Size(); d++)
   {
      if (!bdr_attr_to_show[attributes[d]-1]) { continue; }

      nx=0.;
      ny=0.;
      nz=0.;

      for (i = 0; i < ne; i++)
      {
         int attr =
            (dim == 3) ? mesh->GetBdrAttribute(i) : mesh->GetAttribute(i);
         if (attr != attributes[d]) { continue; }

         if (dim == 3)
         {
            mesh->GetBdrPointMatrix (i, pointmat);
            mesh->GetBdrElementVertices (i, vertices);
         }
         else
         {
            mesh->GetPointMatrix(i, pointmat);
            mesh->GetElementVertices(i, vertices);
         }

         for (j = 0; j < pointmat.Size(); j++)
         {
            pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
            pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
            pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
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

      for (i = 0; i < ne; i++)
      {
         int attr =
            (dim == 3) ? mesh->GetBdrAttribute(i) : mesh->GetAttribute(i);
         if (attr != attributes[d]) { continue; }

         int el_type =
            (dim == 3) ? mesh->GetBdrElementType(i) : mesh->GetElementType(i);
         switch (el_type)
         {
            case Element::TRIANGLE:
               draw.glBegin (GL_TRIANGLES);
               break;
            case Element::QUADRILATERAL:
               draw.glBegin (GL_QUADS);
               break;
            case Element::SEGMENT:
               draw.glBegin(GL_LINES);
               break;
            default:
               MFEM_ABORT("Invalid boundary element type");
               break;
         }
         if (dim == 3)
         {
            mesh->GetBdrPointMatrix (i, pointmat);
            mesh->GetBdrElementVertices (i, vertices);
         }
         else
         {
            mesh->GetPointMatrix(i, pointmat);
            mesh->GetElementVertices(i, vertices);
         }
         for (j = 0; j < pointmat.Size(); j++)
         {
            pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
            pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
            pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
         }
         for (j = 0; j < pointmat.Size(); j++)
         {
            MySetColor(draw, (*sol)(vertices[j]), minv, maxv);
            draw.glNormal3d(nx(vertices[j]), ny(vertices[j]), nz(vertices[j]));
            draw.glVertex3dv(&pointmat(0, j));
         }
         draw.glEnd();
      }
   }
   updated_bufs.emplace_back(&disp_buf);
}

void VisualizationSceneVector3d::PrepareLines()
{
   if (!drawmesh) { return; }

   if (shading == Shading::Noncomforming)
   {
      PrepareLines2();
      return;
   }

   int dim = mesh->Dimension();
   int i, j, ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;
   double point[4][4];

   line_buf.clear();

   for (i = 0; i < ne; i++)
   {
      int attr = (dim == 3) ? mesh->GetBdrAttribute(i) : mesh->GetAttribute(i);
      if (!bdr_attr_to_show[attr-1]) { continue; }

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
         mesh->GetBdrElementVertices(i, vertices);
         mesh->GetBdrPointMatrix(i, pointmat);
      }
      else
      {
         mesh->GetElementVertices(i, vertices);
         mesh->GetPointMatrix(i, pointmat);
      }

      if (cplane == 2 && CheckPositions(vertices)) { continue; }

      for (j = 0; j < pointmat.Size(); j++)
      {
         pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
         pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
         pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
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

               for (j = 0; j < pointmat.Size(); j++)
               {
                  line.glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j));
               }
               line.glEnd();
            }
            break;
         }
         case 2:
            for (j = 0; j < pointmat.Size(); j++)
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

void VisualizationSceneVector3d::PrepareLines2()
{
   int fn, fo, di = 0;
   double bbox_diam;

   line_buf.clear();
   gl3::GlBuilder line = line_buf.createBuilder();

   int dim = mesh->Dimension();
   int ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat;
   DenseMatrix vec_vals;
   Vector values;
   RefinedGeometry * RefG;
   Array<int> vertices;

   bbox_diam = sqrt ( (bb.x[1]-bb.x[0])*(bb.x[1]-bb.x[0]) +
                      (bb.y[1]-bb.y[0])*(bb.y[1]-bb.y[0]) +
                      (bb.z[1]-bb.z[0])*(bb.z[1]-bb.z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   for (int i = 0; i < ne; i++)
   {
      int attr = (dim == 3) ? mesh->GetBdrAttribute(i) : mesh->GetAttribute(i);
      if (!bdr_attr_to_show[attr-1]) { continue; }

      if (dim == 3)
      {
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
            mesh->GetBdrElementVertices (i, vertices);
         }
      }
      else
      {
         mesh->GetElementVertices(i, vertices);
      }

      if (cplane == 2 && CheckPositions(vertices)) { continue; }

      if (dim == 3)
      {
         mesh -> GetBdrElementFace (i, &fn, &fo);
         RefG = GLVisGeometryRefiner.Refine(mesh -> GetFaceGeometry (fn),
                                            TimesToRefine);
         // di = GridF->GetFaceValues(fn, 2, RefG->RefPts, values, pointmat);
         di = fo % 2;
         if (di == 1 && !mesh->FaceIsInterior(fn))
         {
            di = 0;
         }
         GridF -> GetFaceValues (fn, di, RefG->RefPts, values, pointmat);
         VecGridF -> GetFaceVectorValues (fn, di, RefG->RefPts,
                                          vec_vals, pointmat);
         if (ianim > 0)
         {
            pointmat.Add (double(ianim)/ianimmax, vec_vals);
         }
         ShrinkPoints(pointmat, i, fn, di);
      }
      else
      {
         RefG = GLVisGeometryRefiner.Refine(mesh->GetElementBaseGeometry(i),
                                            TimesToRefine);
         GridF->GetValues(i, RefG->RefPts, values, pointmat);
         VecGridF->GetVectorValues(i, RefG->RefPts, vec_vals, pointmat);
         if (ianim > 0)
         {
            pointmat.Add(double(ianim)/ianimmax, vec_vals);
         }
         ShrinkPoints(pointmat, i, 0, 0);
      }

      int *RG = &(RefG->RefGeoms[0]);
      double pts[][3] =
      {
         { pointmat(0,RG[0]), pointmat(1,RG[0]), pointmat(2,RG[0]) },
         { pointmat(0,RG[1]), pointmat(1,RG[1]), pointmat(2,RG[1]) },
         { pointmat(0,RG[2]), pointmat(1,RG[2]), pointmat(2,RG[2]) }
      };
      double norm[3];
      Compute3DUnitNormal (pts[0], pts[1], pts[2], norm);
      if (di != 0 && sc != 0.0)
      {
         norm[0] = -norm[0];
         norm[1] = -norm[1];
         norm[2] = -norm[2];
      }

      if (drawmesh == 1)
      {
         Array<int> &REdges = RefG->RefEdges;
         if (mesh_coarse)
         {
            if (sc == 0.0)
            {
               DrawBdrElCoarseSurfEdges(line, i, pointmat, &RefG->RefPts, &REdges);
            }
            else
            {
               DenseMatrix pointmat_shift = pointmat;
               for (int j = 0; j < REdges.Size(); j++)
               {
                  double val = sc * (values(REdges[j]) - minv) / (maxv - minv);
                  for (int d = 0; d < 3; d++)
                  {
                     pointmat_shift(d, REdges[j]) += val*norm[d];
                  }
               }
               DrawBdrElCoarseSurfEdges(line, i, pointmat_shift, &RefG->RefPts, &REdges);
            }
         }
         else
         {
            line.glBegin (GL_LINES);
            for (int k = 0; k < REdges.Size()/2; k++)
            {
               int *RE = &(REdges[2*k]);

               if (sc == 0.0)
               {
                  for (int j = 0; j < 2; j++)
                     line.glVertex3d (pointmat(0, RE[j]),
                                      pointmat(1, RE[j]),
                                      pointmat(2, RE[j]));
               }
               else
               {
                  for (int j = 0; j < 2; j++)
                  {
                     double val = sc * (values(RE[j]) - minv) / (maxv - minv);
                     line.glVertex3d (pointmat(0, RE[j])+val*norm[0],
                                      pointmat(1, RE[j])+val*norm[1],
                                      pointmat(2, RE[j])+val*norm[2]);
                  }
               }
            }
            line.glEnd();
         }
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
            RG = &(RefG->RefGeoms[k*sides]);

            if (sc == 0.0)
            {
               for (int j = 0; j < sides; j++)
               {
                  for (int ii = 0; ii < 3; ii++)
                  {
                     point[j][ii] = pointmat(ii, RG[j]);
                  }
                  point[j][3] = values(RG[j]);
               }
            }
            else
            {
               for (int j = 0; j < sides; j++)
               {
                  double val = (values(RG[j]) - minv) / (maxv - minv);
                  val *= sc;
                  for (int ii = 0; ii < 3; ii++)
                  {
                     point[j][ii] = pointmat(ii, RG[j])+val*norm[ii];
                  }
                  point[j][3] = values(RG[j]);
               }
            }
            DrawPolygonLevelLines(line, point[0], sides, level, false);
         }
      }
   }
   updated_bufs.emplace_back(&line_buf);
}

void VisualizationSceneVector3d::PrepareDisplacedMesh()
{
   int dim = mesh->Dimension();
   int i, j, ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
   DenseMatrix pointmat;
   Array<int> vertices;

   // prepare the displaced mesh
   displine_buf.clear();
   gl3::GlBuilder builder = displine_buf.createBuilder();

   for (i = 0; i < ne; i++)
   {
      builder.glBegin(GL_LINE_LOOP);
      if (dim == 3)
      {
         mesh->GetBdrPointMatrix (i, pointmat);
         mesh->GetBdrElementVertices (i, vertices);
      }
      else
      {
         mesh->GetPointMatrix(i, pointmat);
         mesh->GetElementVertices(i, vertices);
      }

      for (j = 0; j < pointmat.Size(); j++)
      {
         pointmat(0, j) += (*solx)(vertices[j])*(ianimd)/ianimmax;
         pointmat(1, j) += (*soly)(vertices[j])*(ianimd)/ianimmax;
         pointmat(2, j) += (*solz)(vertices[j])*(ianimd)/ianimmax;
      }

      for (j = 0; j < pointmat.Size(); j++)
      {
         builder.glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j) );
      }
      builder.glEnd();
   }
   updated_bufs.emplace_back(&displine_buf);
}

void ArrowsDrawOrNot (Array<int> l[], int nv, Vector & sol,
                      int nl, Array<double> & level)
{
   static int first_time = 1;
   static int nll = nl;

   if (!first_time && nll == nl)
   {
      return;
   }
   else
   {
      first_time = 1;
      nll = nl;
   }

   int i,j;
   double v;
   double eps = (level[1] - level[0])/2;

   for (i = 0; i <= nl; i++)
   {
      l[i].SetSize(0);
   }

   for (j = 0; j < nv; j++)
   {
      v = sol(j);
      for (i = 0; i <= nl; i++)
      {
         if (fabs(v-level[i]) < eps)
         {
            l[i].Append(j);
            break;
         }
         if (v < level[i] - eps)
         {
            break;
         }
      }
   }
}

int ArrowDrawOrNot (double v, int nl, Array<double> & level)
{
   double eps = (level[nl] - level[0])/10;
   for (int i = 0; i <= nl; i++)
   {
      if (fabs(v-level[i]) < eps)
      {
         return 1;
      }
      if (v < level[i] - eps)
      {
         return 0;
      }
   }
   return 0;
}

void VisualizationSceneVector3d::DrawVector(gl3::GlDrawable& buf,
                                            int type, double v0, double v1,
                                            double v2, double sx, double sy,
                                            double sz, double s)
{
   static int nv = mesh -> GetNV();
   static double bb_vol = (bb.x[1]-bb.x[0])*(bb.y[1]-bb.y[0])*(bb.z[1]-bb.z[0]);
   static double volume = std::max(bb_vol, mesh_volume);
   static double h      = pow(volume/nv, 0.333);
   static double hh     = pow(volume, 0.333) / 10;

   switch (type)
   {
      case 1:
      {
         arrow_type = 0;
         arrow_scaling_type = 0;
         // glColor3f(0, 0, 0); // color is set in Draw()
         Arrow(buf,v0,v1,v2,sx,sy,sz,s);
      }
      break;

      case 2:
      {
         arrow_type = 1;
         arrow_scaling_type = 1;
         Arrow(buf,v0,v1,v2,sx,sy,sz,h,0.125,s);
      }
      break;

      case 3:
      {
         arrow_type = 1;
         arrow_scaling_type = 1;
         Arrow(buf,v0,v1,v2,sx,sy,sz,h*s/maxv,0.125,s);
      }
      break;

      case 4:
      case 5:
      {
         arrow_type = 1;
         arrow_scaling_type = 1;
         // glColor3f(0.3, 0.3, 0.3); // color is set in Draw
         Arrow(buf,v0,v1,v2,sx,sy,sz,hh*s/maxv,0.125);
      }
      break;
   }
}

void VisualizationSceneVector3d::PrepareVectorField()
{
   int i, nv = mesh -> GetNV();
   double *vertex;

   vector_buf.clear();

   switch (drawvector)
   {
      case 0:
         break;

      case 1:
         for (i = 0; i < nv; i++)
            if (drawmesh != 2 || ArrowDrawOrNot((*sol)(i), nl, level))
            {
               vertex = mesh->GetVertex(i);
               DrawVector(vector_buf, drawvector, vertex[0], vertex[1], vertex[2],
                          (*solx)(i), (*soly)(i), (*solz)(i), (*sol)(i));
            }
         break;

      case 2:
      {
         arrow_type = 1;
         arrow_scaling_type = 1;
         for (i = 0; i < nv; i++)
            if (drawmesh != 2 || ArrowDrawOrNot((*sol)(i), nl, level))
            {
               vertex = mesh->GetVertex(i);
               DrawVector(vector_buf, drawvector, vertex[0], vertex[1], vertex[2],
                          (*solx)(i), (*soly)(i), (*solz)(i), (*sol)(i));
            }
      }
      break;

      case 3:
      {
         arrow_type = 1;
         arrow_scaling_type = 1;

         for (i = 0; i < nv; i++)
            if (drawmesh != 2 || ArrowDrawOrNot((*sol)(i), nl, level))
            {
               vertex = mesh->GetVertex(i);
               DrawVector(vector_buf, drawvector, vertex[0], vertex[1], vertex[2],
                          (*solx)(i), (*soly)(i), (*solz)(i), (*sol)(i));
            }
      }
      break;

      case 4:
      {
         Array<int> *l = new Array<int>[nl+1];
         ArrowsDrawOrNot(l, nv, *sol, nl, level);

         int j,k;

         for (k = 0; k < vflevel.Size(); k++)
         {
            i = vflevel[k];
            for (j = 0; j < l[i].Size(); j++)
            {
               vertex = mesh->GetVertex( l[i][j] );
               DrawVector(vector_buf, drawvector, vertex[0], vertex[1], vertex[2],
                          (*solx)(l[i][j]), (*soly)(l[i][j]), (*solz)(l[i][j]),
                          (*sol)(l[i][j]));
            }
         }

         delete [] l;
      }
      break;

      case 5:
      {
         int dim = mesh->Dimension();
         int j, k, ne = (dim == 3) ? mesh->GetNBE() : mesh->GetNE();
         Array<int> vertices;
         Array<bool> vert_marker(nv);

         vert_marker = false;
         for (k = 0; k < ne; k++)
         {
            if (dim == 3)
            {
               mesh->GetBdrElementVertices(k, vertices);
            }
            else
            {
               mesh->GetElementVertices(k, vertices);
            }
            for (j = 0; j < vertices.Size(); j++)
            {
               i = vertices[j];
               if (vert_marker[i]) { continue; }
               vertex = mesh->GetVertex(i);
               DrawVector(vector_buf, drawvector, vertex[0], vertex[1], vertex[2],
                          (*solx)(i), (*soly)(i), (*solz)(i), (*sol)(i));
               vert_marker[i] = true;
            }
         }
      }
      break;
   }
   updated_bufs.emplace_back(&vector_buf);
}

void VisualizationSceneVector3d::PrepareCuttingPlane()
{
   if (cp_drawelems == 0 || cplane != 1 || drawvector == 0 ||
       mesh->Dimension() != 3)
   {
      VisualizationSceneSolution3d::PrepareCuttingPlane();
      return;
   }

   int i, j, k, n = 0;
   int flag[4], ind[6][2]= {{0,3},{0,2},{0,1},{1,2},{1,3},{2,3}};
   double t, point[4][4], val[4][3];

   cplane_buf.clear();
   gl3::GlBuilder builder = cplane_buf.createBuilder();

   DenseMatrix pointmat(3,4);
   double * coord;

   Array<int> nodes;
   for (i = 0; i < mesh->GetNE(); i++)
   {
      if (mesh->GetElementType(i) != Element::TETRAHEDRON)
      {
         continue;
      }
      // TODO: support for hex elements as in
      // VisualizationSceneSolution3d::CuttingPlaneFunc

      n = 0; // n will be the number of intersection points
      mesh -> GetElementVertices(i,nodes);
      mesh -> GetPointMatrix (i,pointmat);
      for (j=0; j<4; j++)
         if (node_pos[nodes[j]] == 0.0)
         {
            flag[j] = 1;
            coord = mesh -> GetVertex(nodes[j]);
            for (k=0; k<3; k++)
            {
               point[n][k] = coord[k];
            }

            point[n][3] = (*sol)(nodes[j]);
            val[n][0] = (*solx)(nodes[j]);
            val[n][1] = (*soly)(nodes[j]);
            val[n][2] = (*solz)(nodes[j]);
            n++;
         }
         else if (node_pos[nodes[j]] < 0.)
         {
            flag[j] = -1;
         }
         else
         {
            flag[j] = 0;
         }

      for (j=0; j<6; j++)
         if (flag[ind[j][0]] != 1 && flag[ind[j][1]] != 1 &&
             flag[ind[j][0]] != flag[ind[j][1]])
         {
            t = node_pos[ nodes[ind[j][1]] ] / (node_pos[ nodes[ind[j][1]] ] -
                                                node_pos[ nodes[ind[j][0]] ] );
            for (k=0; k<3; k++)
               point[n][k] = t*(pointmat(k,ind[j][0]) -
                                pointmat(k,ind[j][1])) +
                             pointmat(k,ind[j][1]);

            point[n][3] = t * ((*sol)(nodes[ind[j][0]]) -
                               (*sol)(nodes[ind[j][1]])) +
                          (*sol)(nodes[ind[j][1]]);
            val[n][0] = t * ((*solx)(nodes[ind[j][0]]) -
                             (*solx)(nodes[ind[j][1]])) +
                        (*solx)(nodes[ind[j][1]]);
            val[n][1] = t * ((*soly)(nodes[ind[j][0]]) -
                             (*soly)(nodes[ind[j][1]])) +
                        (*soly)(nodes[ind[j][1]]);
            val[n][2] = t * ((*solz)(nodes[ind[j][0]]) -
                             (*solz)(nodes[ind[j][1]])) +
                        (*solz)(nodes[ind[j][1]]);
            n++;
         }

      if (n > 2)
      {
         double v10[] = { point[1][0]-point[0][0],
                          point[1][1]-point[0][1],
                          point[1][2]-point[0][2]
                        };
         double v21[] = { point[2][0]-point[1][0],
                          point[2][1]-point[1][1],
                          point[2][2]-point[1][2]
                        };

         double norm[] = { v10[1]*v21[2]-v10[2]*v21[1],
                           v10[2]*v21[0]-v10[0]*v21[2],
                           v10[0]*v21[1]-v10[1]*v21[0]
                         };

         double * eq = CuttingPlane -> Equation();

         if ( eq[0]*norm[0]+eq[1]*norm[1]+eq[2]*norm[2] > 0.0 )
         {
            if (drawvector != 5)
            {
               builder.glBegin (GL_POLYGON);
               for (j=0; j<n; j++)
               {
                  MySetColor ( builder, point[j][3], maxv, minv);
                  builder.glNormal3dv (CuttingPlane -> Equation());
                  builder.glVertex3d (point[j][0],point[j][1],point[j][2]);
               }
            }
            builder.glEnd();
            if (drawvector)
               for (j=0; j<n; j++)
                  DrawVector(cplane_buf, drawvector, point[n][0], point[n][1], point[n][2],
                             val[n][0],val[n][1], val[n][2], point[n][3]);
         }
         else
         {
            if (drawvector != 5)
            {
               builder.glBegin (GL_POLYGON);
               for (j=n-1; j>=0; j--)
               {
                  MySetColor ( builder, point[j][3], minv, maxv);
                  builder.glNormal3dv (CuttingPlane -> Equation());
                  builder.glVertex3d (point[j][0],point[j][1],point[j][2]);
               }
            }
            builder.glEnd();
            if (drawvector)
               for (j=n-1; j>=0; j--)
                  DrawVector(cplane_buf, drawvector, point[n][0], point[n][1], point[n][2],
                             val[n][0],val[n][1], val[n][2], point[n][3]);
         }
      }
   }
   updated_bufs.emplace_back(&cplane_buf);
}

gl3::SceneInfo VisualizationSceneVector3d::GetSceneObjs()
{
   if (colorbar)
   {
      Array<double>* cb_level = nullptr;
      if (drawvector == 4)
      {
         cb_level = &dvflevel;
      }
      else if (drawmesh == 2 || cp_drawmesh >= 2)
      {
         cb_level = &level;
      }
      PrepareColorBar(minv, maxv, cb_level);
   }
   gl3::SceneInfo scene = VisualizationSceneScalarData::GetSceneObjs();
   gl3::RenderParams params = GetMeshDrawParams();
   params.use_clip_plane = cplane;
   double* cp_eqn = CuttingPlane->Equation();
   params.clip_plane_eqn = {cp_eqn[0], cp_eqn[1], cp_eqn[2], cp_eqn[3]};
   params.contains_translucent = matAlpha < 1.0 ||
                                 palette.GetPalette()->IsTranslucent();
   // draw vector field
   if (drawvector == 2 || drawvector == 3)
   {
      scene.queue.emplace_back(params, &vector_buf);
   }
   // draw elements
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
   params.num_pt_lights = 0;

   params.mesh_material = VisualizationScene::BLK_MAT;

   if (drawvector == 1)
   {
      params.static_color = GetLineColor();
      scene.queue.emplace_back(params, &vector_buf);
   }
   else if (drawvector > 3)
   {
      params.static_color = {0.3,0.3,0.3,1.0};
      scene.queue.emplace_back(params, &vector_buf);
      params.static_color = GetLineColor();
   }

   // draw lines
   if (drawmesh)
   {
      scene.queue.emplace_back(params, &line_buf);
   }

   // draw displacement
   if (drawdisp)
   {
      params.static_color = {1.f, 0.f, 0.f, 1.f};
      scene.queue.emplace_back(params, &displine_buf);
      params.static_color = GetLineColor();
   }
   // ruler may have mixture of polygons and lines
   if (cp_drawmesh)
   {
      params.use_clip_plane = false;
      scene.queue.emplace_back(params, &cplines_buf);
   }
   ProcessUpdatedBufs(scene);
   return scene;
}
