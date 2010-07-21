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

#include "mfem.hpp"
#include "visual.hpp"

static void VectorKeyHPressed (){
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
        << "| E -  Toggle the elements in the CP |" << endl
        << "| f -  Smooth/Flat shading           |" << endl
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
        << "| p -  Cycle through color palettes  |" << endl
        << "| P -  Print to PostScript file      |" << endl
        << "| q -  Quits                         |" << endl
        << "| r -  Reset the plot to 3D view     |" << endl
        << "| R -  Reset the plot to 2D view     |" << endl
        << "| s -  Turn on/off unit cube scaling |" << endl
        << "| S -  Take snapshot/Record a movie  |" << endl
        << "| t -  Cycle materials and lights    |" << endl
        << "| \\ -  Set light source position     |" << endl
        << "| u/U  Move the level field vectors  |" << endl
        << "| v -  Vector field                  |" << endl
        << "| w/W  Add/Delete level field vector |" << endl
        << "| x/X  Rotate clipping plane (phi)   |" << endl
        << "| y/Y  Rotate clipping plane (theta) |" << endl
        << "| z/Z  Translate clipping plane      |" << endl
        << "+------------------------------------+" << endl
        << "| Function keys                      |" << endl
        << "+------------------------------------+" << endl
        << "| F1 - X window info and keystrokes  |" << endl
        << "| F2 - Update colors, etc.           |" << endl
        << "| F3/F4 - Shrink/Zoom subdomains     |" << endl
        << "| F6 - Palete options                |" << endl
        << "| F7 - Manually set min/max value    |" << endl
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

VisualizationSceneVector3d  * vsvector3d;
extern VisualizationScene * locscene;
//extern VisualizationSceneVector3d * vssol3d;

static void KeyDPressed(){
   vsvector3d -> ToggleDisplacements();
   SendExposeEvent();
}

static void KeyNPressed ()
{
   if (vsvector3d -> drawdisp)
      vsvector3d -> ianimd =( (vsvector3d ->ianimd + 1) %
                              (vsvector3d -> ianimmax + 1) );
   else
      vsvector3d -> ianim =( (vsvector3d ->ianim + 1) %
                             (vsvector3d -> ianimmax + 1) );
   vsvector3d -> NPressed();
}

static void KeyBPressed ()
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

static void KeyrPressed (){
   locscene -> spinning = 0;
   auxIdleFunc((void (*)())0);
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

static void KeyRPressed (){
   vsvector3d -> ianim = vsvector3d -> ianimd = 0;
   vsvector3d -> Prepare();
   vsvector3d -> PrepareLines();
   vsvector3d -> PrepareDisplacedMesh();
}

void VisualizationSceneVector3d::NPressed ()
{
   if (drawdisp)
      PrepareDisplacedMesh();
   else
   {
      Prepare();
      PrepareLines();
   }

   SendExposeEvent();
}

static void KeyuPressed(){
   vsvector3d -> ToggleVectorFieldLevel(+1);
   SendExposeEvent();
}

static void KeyUPressed(){
   vsvector3d -> ToggleVectorFieldLevel(-1);
   SendExposeEvent();
}

void VisualizationSceneVector3d::ToggleVectorFieldLevel (int v)
{
   int i;
   for (i = 0; i < vflevel.Size(); i++)
      if (vflevel[i] == 0 && v == -1)
         vflevel[i] = nl;
      else
         vflevel[i] = (vflevel[i] + v) % (nl+1);
   for (i = 0; i < vflevel.Size(); i++)
      dvflevel[i] = level[vflevel[i]];
   vsvector3d -> PrepareVectorField();
}

static void KeywPressed(){
   vsvector3d -> AddVectorFieldLevel();
   SendExposeEvent();
}

static void KeyWPressed(){
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

static void KeyvPressed(){
   vsvector3d -> ToggleVectorField(1);
   SendExposeEvent();
}

static void KeyVPressed(){
   vsvector3d -> ToggleVectorField(-1);
   SendExposeEvent();
}

void VisualizationSceneVector3d::ToggleVectorField(int i){
   if (drawvector == 0 && i == -1)
      drawvector = 5;
   else
      drawvector = (drawvector+i)%6;
   PrepareVectorField();
}

// VisualizationSceneVector3d::VisualizationSceneVector3d(istream & in)
// {
//   mesh = new DMesh<2>();
//   solx = new Vector();
//   soly = new Vector();
//   solz = new Vector();
//
//   mesh -> Load (in);
//   solx -> Load (in, mesh -> GetNV());
//   soly -> Load (in, mesh -> GetNV());
//   solz -> Load (in, mesh -> GetNV());
//
//   sol  = new Vector(mesh -> GetNV());
//
//   Init();
// }

VisualizationSceneVector3d
::VisualizationSceneVector3d(Mesh & m, Vector & sx, Vector & sy, Vector & sz)
{
   mesh = &m;
   solx = &sx;
   soly = &sy;
   solz = &sz;

   sol  = new Vector(mesh -> GetNV());

   sfes = NULL;
   VecGridF = NULL;

   Init();
}

VisualizationSceneVector3d::VisualizationSceneVector3d (GridFunction &vgf)
{
   FiniteElementSpace *fes = vgf.FESpace();
   if (fes == NULL || fes->GetVDim() != 3)
   {
      cout << "VisualizationSceneVector3d::VisualizationSceneVector3d" << endl;
      exit(1);
   }

   VecGridF = &vgf;

   mesh = fes->GetMesh();

   sfes = new FiniteElementSpace (mesh, fes->FEColl(), 1, fes->GetOrdering());
   GridF = new GridFunction (sfes);

   int ndofs = GridF->Size();

   Array<int> dofs(3);
   for (int i = 0; i < ndofs; i++)
   {
      dofs.SetSize(1);
      dofs[0] = i;
      fes->DofsToVDofs (dofs);
      double x = (*VecGridF)(dofs[0]);
      double y = (*VecGridF)(dofs[1]);
      double z = (*VecGridF)(dofs[2]);

      (*GridF)(i) = sqrt(x*x+y*y+z*z);
   }

   solx = new Vector(mesh -> GetNV());
   soly = new Vector(mesh -> GetNV());
   solz = new Vector(mesh -> GetNV());

   vgf.GetNodalValues (*solx, 1);
   vgf.GetNodalValues (*soly, 2);
   vgf.GetNodalValues (*solz, 3);

   sol  = new Vector(mesh -> GetNV());

   Init();
}

void VisualizationSceneVector3d::Init()
{
   int i, nv = mesh -> GetNV();

   key_r_state = 0;

   drawdisp = 0;
   drawvector = 0;

   ianim = ianimd = 0;
   ianimmax = 10;

   for(i=0; i<nv; i++)
      (*sol)(i) = sqrt((*solx)(i) * (*solx)(i) +
                       (*soly)(i) * (*soly)(i) +
                       (*solz)(i) * (*solz)(i) );

   vectorlist = glGenLists(1);
   displinelist = glGenLists(1);

   VisualizationSceneSolution3d :: Init();

   PrepareVectorField();
   PrepareDisplacedMesh();

   vflevel.Append(0);
   dvflevel.Append(level[0]);

   vsvector3d = this;

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

      auxKeyFunc (AUX_r, KeyrPressed);
      auxKeyFunc (AUX_R, KeyRPressed);

      auxKeyFunc (AUX_u, KeyuPressed);
      auxKeyFunc (AUX_U, KeyUPressed);

      auxKeyFuncReplace (AUX_w, KeywPressed); // Keys w, W are also used in
      auxKeyFuncReplace (AUX_W, KeyWPressed); // VisualizationSceneSolution3d

      auxKeyFuncReplace (AUX_v, KeyvPressed); // Keys v, V are also used in
      auxKeyFuncReplace (AUX_V, KeyVPressed); // VisualizationSceneSolution3d
   }
}

VisualizationSceneVector3d::~VisualizationSceneVector3d ()
{
   glDeleteLists (vectorlist, 1);
   glDeleteLists (displinelist, 1);

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


void VisualizationSceneVector3d :: PrepareFlat (){
   int i, j;

   glNewList (displlist, GL_COMPILE);

   Set_Material ();
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   int nbe = mesh -> GetNBE();
   DenseMatrix pointmat;
   Array<int> vertices;

   for (i = 0; i < nbe; i++)
   {
      if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) continue;

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
      mesh->GetBdrElementVertices (i, vertices);

      for (j = 0; j < pointmat.Size(); j++) {
         pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
         pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
         pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
      }

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
         MySetColor ( (*sol)(vertices[j]) , minv, maxv);
         glVertex3d (pointmat(0, j),
                     pointmat(1, j),
                     pointmat(2, j));
      }
      glEnd ();
   }
   glEndList ();
}

void VisualizationSceneVector3d::PrepareFlat2()
{
   int i, j, k, fn, fo, ft, di;
   double bbox_diam, vmin, vmax, mm;
   int nbe = mesh -> GetNBE();

   DenseMatrix pointmat;
   DenseMatrix vec_vals;
   Vector values;
   RefinedGeometry * RefG;
   Array<int> vertices;

   bbox_diam = sqrt ( (x[1]-x[0])*(x[1]-x[0]) +
                      (y[1]-y[0])*(y[1]-y[0]) +
                      (z[1]-z[0])*(z[1]-z[0]) );
   double sc = FaceShiftScale * bbox_diam;

   glNewList (displlist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

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
      di = fo % 2;
      GridF -> GetFaceValues (fn, di, RefG->RefPts, values, pointmat);
      VecGridF -> GetFaceVectorValues (fn, di, RefG->RefPts,
                                       vec_vals, pointmat);
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

      if (ianim > 0)
         pointmat.Add (double(ianim)/ianimmax, vec_vals);

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
      Compute3DUnitNormal (pts[0], pts[1], pts[2], norm);
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
               MySetColor ( values(RG[j]) , minv , maxv);
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
            Compute3DUnitNormal (pts[0], pts[1], pts[2], nnorm);
            glNormal3dv (nnorm);
            for (j = 0; j < 3; j++)
            {
               MySetColor ( values(RG[j]) , minv , maxv);
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
   cout << "VisualizationSceneVector3d::PrepareFlat2() : [min,max] = ["
        << vmin << "," << vmax << "]" << endl;
}

void VisualizationSceneVector3d :: Prepare (){
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
   Set_Material ();
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

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

            for (j = 0; j < pointmat.Size(); j++) {
               pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
               pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
               pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
            }

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
            mesh->GetBdrElementVertices (i, vertices);

            for (j = 0; j < pointmat.Size(); j++) {
               pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
               pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
               pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
            }

            for (j = 0; j < pointmat.Size(); j++)
            {
               MySetColor ( (*sol)(vertices[j]) , minv , maxv);
               glNormal3d (nx(vertices[j]),ny(vertices[j]),nz(vertices[j]));
               glVertex3d (pointmat(0, j),
                           pointmat(1, j),
                           pointmat(2, j));
            }
            glEnd ();
         }
   }
   glEndList ();
}

void VisualizationSceneVector3d::PrepareLines()
{
   if (!drawmesh) return;

   if (shading == 2)
   {
      PrepareLines2();
      return;
   }

   int i, j, ne = mesh -> GetNBE();
   DenseMatrix pointmat;
   Array<int> vertices;

   glNewList(linelist, GL_COMPILE);

   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   for (i = 0; i < ne; i++)
   {
      if (!bdr_attr_to_show[mesh->GetBdrAttribute(i)-1]) continue;

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
      mesh->GetBdrElementVertices (i, vertices);

      for (j = 0; j < pointmat.Size(); j++) {
         pointmat(0, j) += (*solx)(vertices[j])*(ianim)/ianimmax;
         pointmat(1, j) += (*soly)(vertices[j])*(ianim)/ianimmax;
         pointmat(2, j) += (*solz)(vertices[j])*(ianim)/ianimmax;
      }

      for (j = 0; j < pointmat.Size(); j++)
         glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j) );
      glEnd ();
   }
   glEndList ();
}

void VisualizationSceneVector3d::PrepareLines2()
{
   int i, j, k, fn, fo, di;
   double bbox_diam;

   glNewList (linelist, GL_COMPILE);

   int nbe = mesh -> GetNBE();
   DenseMatrix pointmat;
   DenseMatrix vec_vals;
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
      VecGridF -> GetFaceVectorValues (fn, di, RefG->RefPts,
                                       vec_vals, pointmat);

      if (ianim > 0)
         pointmat.Add (double(ianim)/ianimmax, vec_vals);

      int *RG = &(RefG->RefGeoms[0]);
      double pts[][3] = { { pointmat(0,RG[0]), pointmat(1,RG[0]),
                            pointmat(2,RG[0]) },
                          { pointmat(0,RG[1]), pointmat(1,RG[1]),
                            pointmat(2,RG[1]) },
                          { pointmat(0,RG[2]), pointmat(1,RG[2]),
                            pointmat(2,RG[2]) } };
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
         glBegin (GL_LINES);
         for (k = 0; k < REdges.Size()/2; k++)
         {
            int *RE = &(REdges[2*k]);

            if (sc == 0.0)
            {
               for (j = 0; j < 2; j++)
                  glVertex3d (pointmat(0, RE[j]), pointmat(1, RE[j]),
                              pointmat(2, RE[j]));
            }
            else
            {
               for (j = 0; j < 2; j++)
               {
                  double val = sc * (values(RE[j]) - minv) / (maxv - minv);
                  glVertex3d (pointmat(0, RE[j])+val*norm[0],
                              pointmat(1, RE[j])+val*norm[1],
                              pointmat(2, RE[j])+val*norm[2]);
               }
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

void VisualizationSceneVector3d::PrepareDisplacedMesh(){
   int i, j, ne = mesh -> GetNBE();
   DenseMatrix pointmat;
   Array<int> vertices;

   // prepare the displaced mesh
   glNewList(displinelist, GL_COMPILE);

   Set_Material ();
   glColor3f (1, 0, 0);
   glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

   for (i = 0; i < ne; i++){
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
      mesh->GetBdrElementVertices (i, vertices);

      for (j = 0; j < pointmat.Size(); j++) {
         pointmat(0, j) += (*solx)(vertices[j])*(ianimd)/ianimmax;
         pointmat(1, j) += (*soly)(vertices[j])*(ianimd)/ianimmax;
         pointmat(2, j) += (*solz)(vertices[j])*(ianimd)/ianimmax;
      }

      for (j = 0; j < pointmat.Size(); j++)
         glVertex3d (pointmat(0, j), pointmat(1, j), pointmat(2, j) );
      glEnd ();
   }
   glEndList ();
}

void ArrowsDrawOrNot (Array<int> l[], int nv, Vector & sol, int nl, Array<double> & level)
{
   static int first_time = 1;
   static int nll = nl;

   if (!first_time && nll == nl)
      return;
   else {
      first_time = 1;
      nll = nl;
   }

   int i,j;
   double v;
   double eps = (level[1] - level[0])/2;

   for (i = 0; i <= nl; i++)
      l[i].SetSize(0);

   for (j = 0; j < nv; j++) {
      v = sol(j);
      for (i = 0; i <= nl; i++) {
         if (fabs(v-level[i]) < eps) {
            l[i].Append(j);
            break;
         }
         if (v < level[i] - eps)
            break;
      }
   }
}

int ArrowDrawOrNot (double v, int nl, Array<double> & level)
{
   double eps = (level[nl] - level[0])/10;
   for (int i = 0; i <= nl; i++) {
      if (fabs(v-level[i]) < eps)
         return 1;
      if (v < level[i] - eps)
         return 0;
   }
   return 0;
}

void VisualizationSceneVector3d::DrawVector (int type, double v0, double v1, double v2,
                                             double sx, double sy, double sz, double s)
{
   static int nv = mesh -> GetNV();
   static double volume = (x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0]);
   static double h      = pow(volume/nv, 0.333);
   static double hh     = pow(volume, 0.333) / 10;

   switch (type) {
   case 1:
   {
      arrow_type = 0;
      arrow_scaling_type = 0;
      glColor3f(0, 0, 0);
      Arrow(v0,v1,v2,sx,sy,sz,s);
   }
   break;

   case 2:
   {
      arrow_type = 1;
      arrow_scaling_type = 1;
      MySetColor(s, minv, maxv);
      Arrow(v0,v1,v2,sx,sy,sz,h,0.125);
   }
   break;

   case 3:
   {
      arrow_type = 1;
      arrow_scaling_type = 1;
      MySetColor(s,maxv,minv);
      Arrow(v0,v1,v2,sx,sy,sz,h*s/maxv,0.125);
   }
   break;

   case 4:
   case 5:
   {
      arrow_type = 1;
      arrow_scaling_type = 1;
      glColor3f(0.3, 0.3, 0.3);
      Arrow(v0,v1,v2,sx,sy,sz,hh*s/maxv,0.125);
   }
   break;
   }
}

void VisualizationSceneVector3d::PrepareVectorField(){
   int i, nv = mesh -> GetNV();
   double *vertex;

   glNewList(vectorlist, GL_COMPILE);

   Set_Material ();

   switch (drawvector)
   {
   case 0:
      break;

   case 1:
      for(i=0; i<nv; i++)
         if (drawmesh != 2 || ArrowDrawOrNot ((*sol)(i),nl,level)) {
            vertex = mesh -> GetVertex( i );
            DrawVector(drawvector,vertex[0],  vertex[1], vertex[2],
                       (*solx)(i), (*soly)(i), (*solz)(i),(*sol)(i));
         }
      break;

   case 2:
   {
      arrow_type = 1;
      arrow_scaling_type = 1;
      for(i=0; i<nv; i++)
         if (drawmesh != 2 || ArrowDrawOrNot ((*sol)(i),nl,level)) {
            vertex = mesh -> GetVertex( i );
            DrawVector(drawvector,vertex[0],  vertex[1], vertex[2],
                       (*solx)(i), (*soly)(i), (*solz)(i),(*sol)(i));
         }
   }
   break;

   case 3:
   {
      arrow_type = 1;
      arrow_scaling_type = 1;

      for(i=0; i<nv; i++)
         if (drawmesh != 2 || ArrowDrawOrNot ((*sol)(i),nl,level)) {
            vertex = mesh -> GetVertex( i );
            DrawVector(drawvector,vertex[0],  vertex[1], vertex[2],
                       (*solx)(i), (*soly)(i), (*solz)(i),(*sol)(i));
         }
   }
   break;

   case 4: {
      Array<int> *l = new Array<int>[nl+1];
      ArrowsDrawOrNot (l, nv, *sol, nl, level);

      int j,k;
      /*
        double volume = (x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])
        * xscale * yscale * zscale;
        double h = pow(volume, 0.333) / 10;
      */

      for (k = 0; k < vflevel.Size(); k++) {
         i = vflevel[k];
         for (j = 0; j < l[i].Size(); j++) {
            vertex = mesh -> GetVertex( l[i][j] );
            DrawVector(drawvector,vertex[0],  vertex[1], vertex[2],
                       (*solx)(l[i][j]), (*soly)(l[i][j]), (*solz)(l[i][j]),(*sol)(l[i][j]));
         }
      }

      delete [] l;
   }
      break;

   case 5: {
      int j, k, ne = mesh -> GetNBE();
      Array<int> vertices;

      for (k = 0; k < ne; k++) {
         mesh->GetBdrElementVertices (k, vertices);
         for (j = 0; j < vertices.Size(); j++) {
            i = vertices[j];
            vertex = mesh -> GetVertex( i );
            DrawVector(drawvector,vertex[0],  vertex[1], vertex[2],
                       (*solx)(i), (*soly)(i), (*solz)(i),(*sol)(i));
         }
      }
   }
      break;
   }
   glEndList ();

}

void VisualizationSceneVector3d::PrepareCuttingPlane(){

   if (cp_drawelems == 0) return;
   if (cplane == 0) return;
   if (cplane == 2)
   {
      PrepareCuttingPlane2();
      return;
   }

   /*
     double volume = (x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0])
     * xscale * yscale * zscale;
     double h = pow(volume/nv, 0.333);
   */

   int i, j, k, n = 0;
   int flag[4], ind[6][2]={{0,3},{0,2},{0,1},{1,2},{1,3},{2,3}};
   double t, point[4][4], val[4][3];

   glNewList(cplanelist, GL_COMPILE);
   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

   Set_Material ();

   DenseMatrix pointmat(3,4);
   double * coord;

   Array<int> nodes;
   for(i=0; i < mesh -> GetNE(); i++)
   {
      n = 0;             // n will be the number of intersection points
      mesh -> GetElementVertices(i,nodes);
      mesh -> GetPointMatrix (i,pointmat);
      for(j=0; j<4; j++)
         if (node_pos[nodes[j]] == 0.0)
         {
            flag[j] = 1;
            coord = mesh -> GetVertex(nodes[j]);
            for(k=0; k<3; k++)
               point[n][k] = coord[k];

            point[n][3] = (*sol)(nodes[j]);
            val[n][0] = (*solx)(nodes[j]);
            val[n][1] = (*soly)(nodes[j]);
            val[n][2] = (*solz)(nodes[j]);
            n++;
         }
         else if (node_pos[nodes[j]] < 0.)
            flag[j] = -1;
         else
            flag[j] = 0;

      for(j=0; j<6; j++)
         if (flag[ind[j][0]] != 1 && flag[ind[j][1]] != 1 &&
             flag[ind[j][0]] != flag[ind[j][1]])
         {
            t = node_pos[ nodes[ind[j][1]] ] / (node_pos[ nodes[ind[j][1]] ] -
                                                node_pos[ nodes[ind[j][0]] ] );
            for(k=0; k<3; k++)
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

      if (n > 2){
         double v10[] = { point[1][0]-point[0][0],
                          point[1][1]-point[0][1],
                          point[1][2]-point[0][2] };
         double v21[] = { point[2][0]-point[1][0],
                          point[2][1]-point[1][1],
                          point[2][2]-point[1][2] };

         double norm[] = { v10[1]*v21[2]-v10[2]*v21[1],
                           v10[2]*v21[0]-v10[0]*v21[2],
                           v10[0]*v21[1]-v10[1]*v21[0] };

         double * eq = CuttingPlane -> Equation();

         if ( eq[0]*norm[0]+eq[1]*norm[1]+eq[2]*norm[2] > 0.0 )
         {
            if (drawvector != 5) {
               glBegin (GL_POLYGON);
               for(j=0; j<n; j++)
               {
                  MySetColor ( point[j][3] , maxv , minv);
                  glNormal3dv (CuttingPlane -> Equation());
                  glVertex3d (point[j][0],point[j][1],point[j][2]);
               }
            }
            glEnd();
            if (drawvector)
               for(j=0; j<n; j++)
                  DrawVector(drawvector, point[n][0], point[n][1], point[n][2],
                             val[n][0],val[n][1], val[n][2], point[n][3]);
         }
         else
         {
            if (drawvector != 5) {
               glBegin (GL_POLYGON);
               for(j=n-1; j>=0; j--)
               {
                  MySetColor ( point[j][3] , minv , maxv);
                  glNormal3dv (CuttingPlane -> Equation());
                  glVertex3d (point[j][0],point[j][1],point[j][2]);
               }
            }
            glEnd();
            if (drawvector)
               for(j=n-1; j>=0; j--)
                  DrawVector(drawvector, point[n][0], point[n][1], point[n][2],
                             val[n][0],val[n][1], val[n][2], point[n][3]);
         }
      }
   }
   glEndList();
}

void VisualizationSceneVector3d :: Draw (){
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
      if (drawvector == 4)
         DrawColorBar(minv,maxv,&dvflevel);
      else if (drawmesh == 2 || cp_drawmesh >= 2)
         DrawColorBar(minv,maxv,&level);
      else
         DrawColorBar(minv,maxv);
   if (light)
      glEnable(GL_LIGHTING);

   Set_Black_Material ();
   // draw axes
   if (drawaxes)
   {
      glCallList(axeslist);
      DrawCoordinateCross();
   }
   DrawRuler();

   // define and enable the clipping plane
   if (cplane)
   {
      glClipPlane(GL_CLIP_PLANE0,CuttingPlane->Equation());
      glEnable(GL_CLIP_PLANE0);
      Set_Black_Material ();
      if ( cp_drawmesh )
         glCallList(cplanelineslist);
   }

   Set_Black_Material ();

   // draw lines
   if (drawmesh)
      glCallList(linelist);

   if (drawdisp)
      glCallList(displinelist);

   if (MatAlpha < 1.0)
      Set_Transparency ();

   // draw vector field
   glCallList(vectorlist);

   // draw elements
   if (drawelems)
      glCallList(displlist);

   if (cplane && cp_drawelems)
      glCallList(cplanelist);

   if (MatAlpha < 1.0)
      Remove_Transparency ();

   glFlush();
   glXSwapBuffers (auxXDisplay(), auxXWindow());
}
