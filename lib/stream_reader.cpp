// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "stream_reader.hpp"

using namespace std;
using namespace mfem;

StreamState stream_state;

void Extrude1DMeshAndSolution(Mesh **mesh_p, GridFunction **grid_f_p,
                              Vector *sol)
{
   Mesh *mesh = *mesh_p;

   if (mesh->Dimension() != 1 || mesh->SpaceDimension() != 1)
   {
      return;
   }

   // find xmin and xmax over the vertices of the 1D mesh
   double xmin = numeric_limits<double>::infinity();
   double xmax = -xmin;
   for (int i = 0; i < mesh->GetNV(); i++)
   {
      const double x = mesh->GetVertex(i)[0];
      if (x < xmin)
      {
         xmin = x;
      }
      if (x > xmax)
      {
         xmax = x;
      }
   }

   Mesh *mesh2d = Extrude1D(mesh, 1, 0.1*(xmax - xmin));

   if (grid_f_p && *grid_f_p)
   {
      GridFunction *grid_f_2d =
         Extrude1DGridFunction(mesh, mesh2d, *grid_f_p, 1);

      delete *grid_f_p;
      *grid_f_p = grid_f_2d;
   }
   if (sol && sol->Size() == mesh->GetNV())
   {
      Vector sol2d(mesh2d->GetNV());
      for (int i = 0; i < mesh->GetNV(); i++)
      {
         sol2d(2*i+0) = sol2d(2*i+1) = (*sol)(i);
      }
      *sol = sol2d;
   }

   delete mesh;
   *mesh_p = mesh2d;
}


void SetMeshSolution(Mesh *mesh, GridFunction *&grid_f, bool save_coloring)
{
   auto & state = stream_state;

   if (1) // checkerboard solution
   {
      FiniteElementCollection *cfec;
      if (mesh->Dimension() == 1)
      {
         cfec = new L2_FECollection(0, 1);
      }
      else if (mesh->Dimension() == 2)
      {
         cfec = new Const2DFECollection;
      }
      else
      {
         cfec = new Const3DFECollection;
      }
      FiniteElementSpace *cfes = new FiniteElementSpace(mesh, cfec);
      grid_f = new GridFunction(cfes);
      grid_f->MakeOwner(cfec);
      {
         Array<int> coloring;
         srandom(time(0));
         double a = double(random()) / (double(RAND_MAX) + 1.);
         int el0 = (int)floor(a * mesh->GetNE());
         cout << "Generating coloring starting with element " << el0+1
              << " / " << mesh->GetNE() << endl;
         mesh->GetElementColoring(coloring, el0);
         for (int i = 0; i < coloring.Size(); i++)
         {
            (*grid_f)(i) = coloring[i];
         }
         cout << "Number of colors: " << grid_f->Max() + 1 << endl;
      }
      grid_f->GetNodalValues(state.sol);
      state.is_gf = 1;
      if (save_coloring)
      {
         const char col_fname[] = "GLVis_coloring.gf";
         ofstream fgrid(col_fname);
         cout << "Saving the coloring function -> " << flush;
         grid_f->Save(fgrid);
         cout << col_fname << endl;
      }
   }
   else // zero solution
   {
      state.sol.SetSize (mesh -> GetNV());
      state.sol = 0.0;
   }
}


// Read the content of an input stream (e.g. from socket/file)
int ReadStream(istream &is, const string &data_type)
{
   auto & state = stream_state;

   // 0 - scalar data, 1 - vector data, 2 - mesh only, (-1) - unknown
   int field_type = 0;

   delete state.mesh; state.mesh = NULL;
   delete state.grid_f; state.grid_f = NULL;
   state.keys.clear();
   if (data_type == "fem2d_data")
   {
      state.mesh = new Mesh(is, 0, 0, state.fix_elem_orient);
      state.sol.Load(is, state.mesh->GetNV());
   }
   else if (data_type == "vfem2d_data" || data_type == "vfem2d_data_keys")
   {
      field_type = 1;
      state.mesh = new Mesh(is, 0, 0, state.fix_elem_orient);
      state.solu.Load(is, state.mesh->GetNV());
      state.solv.Load(is, state.mesh->GetNV());
      if (data_type == "vfem2d_data_keys")
      {
         is >> state.keys;
      }
   }
   else if (data_type == "fem3d_data")
   {
      state.mesh = new Mesh(is, 0, 0, state.fix_elem_orient);
      state.sol.Load(is, state.mesh->GetNV());
   }
   else if (data_type == "vfem3d_data" || data_type == "vfem3d_data_keys")
   {
      field_type = 1;
      state.mesh = new Mesh(is, 0, 0, state.fix_elem_orient);
      state.solu.Load(is, state.mesh->GetNV());
      state.solv.Load(is, state.mesh->GetNV());
      state.solw.Load(is, state.mesh->GetNV());
      if (data_type == "vfem3d_data_keys")
      {
         is >> state.keys;
      }
   }
   else if (data_type == "fem2d_gf_data" || data_type == "fem2d_gf_data_keys")
   {
      state.mesh = new Mesh(is, 1, 0, state.fix_elem_orient);
      state.grid_f = new GridFunction(state.mesh, is);
      if (data_type == "fem2d_gf_data_keys")
      {
         is >> state.keys;
      }
   }
   else if (data_type == "vfem2d_gf_data" || data_type == "vfem2d_gf_data_keys")
   {
      field_type = 1;
      state.mesh = new Mesh(is, 1, 0, state.fix_elem_orient);
      state.grid_f = new GridFunction(state.mesh, is);
      if (data_type == "vfem2d_gf_data_keys")
      {
         is >> state.keys;
      }
   }
   else if (data_type == "fem3d_gf_data" || data_type == "fem3d_gf_data_keys")
   {
      state.mesh = new Mesh(is, 1, 0, state.fix_elem_orient);
      state.grid_f = new GridFunction(state.mesh, is);
      if (data_type == "fem3d_gf_data_keys")
      {
         is >> state.keys;
      }
   }
   else if (data_type == "vfem3d_gf_data" || data_type == "vfem3d_gf_data_keys")
   {
      field_type = 1;
      state.mesh = new Mesh(is, 1, 0, state.fix_elem_orient);
      state.grid_f = new GridFunction(state.mesh, is);
      if (data_type == "vfem3d_gf_data_keys")
      {
         is >> state.keys;
      }
   }
   else if (data_type == "solution")
   {
      state.mesh = new Mesh(is, 1, 0, state.fix_elem_orient);
      state.grid_f = new GridFunction(state.mesh, is);
      field_type = (state.grid_f->VectorDim() == 1) ? 0 : 1;
   }
   else if (data_type == "mesh")
   {
      state.mesh = new Mesh(is, 1, 0, state.fix_elem_orient);
      SetMeshSolution(state.mesh, state.grid_f, state.save_coloring);
      field_type = 2;
   }
   else if (data_type == "raw_scalar_2d")
   {
      Array<Array<double> *> vertices;
      Array<Array<int> *> elements;
      Array<int> elem_types;
      string ident;
      int num_patches, num_vert, num_elem, n;
      is >> ws >> ident; // 'patches'
      is >> num_patches;
      // cout << ident << ' ' << num_patches << endl;
      vertices.SetSize(num_patches);
      vertices = NULL;
      elements.SetSize(num_patches);
      elements = NULL;
      elem_types.SetSize(num_patches);
      elem_types = 0;
      int tot_num_vert = 0;
      int tot_num_elem = 0;
      int mesh_type = 0;
      for (int i = 0; i < num_patches; i++)
      {
         is >> ws >> ident; // 'vertices'
         is >> num_vert;
         // cout << '\n' << ident << ' ' << num_vert << endl;
         // read vertices in the format: x y z nx ny nz
         vertices[i] = new Array<double>(6*num_vert);
         Array<double> &verts = *vertices[i];
         for (int j = 0; j < verts.Size(); j++)
         {
            is >> verts[j];
         }

         is >> ws >> ident; // 'triangles' or 'quads'
         if (ident == "triangles")
         {
            n = 3, mesh_type |= 1;
         }
         else
         {
            n = 4, mesh_type |= 2;
         }
         elem_types[i] = n;
         is >> num_elem;
         // cout << ident << ' ' << num_elem << endl;
         elements[i] = new Array<int>(n*num_elem);
         Array<int> &elems = *elements[i];
         for (int j = 0; j < elems.Size(); j++)
         {
            is >> elems[j];
            elems[j] += tot_num_vert;
         }
         tot_num_vert += num_vert;
         tot_num_elem += num_elem;
      }

      state.mesh = new Mesh(2, tot_num_vert, tot_num_elem, 0);
      state.sol.SetSize(tot_num_vert);
      state.normals.SetSize(3*tot_num_vert);

      int v_off = 0;
      for (int i = 0; i < num_patches; i++)
      {
         Array<double> &verts = *vertices[i];
         num_vert = verts.Size()/6;
         for (int j = 0; j < num_vert; j++)
         {
            state.mesh->AddVertex(&verts[6*j]);
            state.sol(v_off) = verts[6*j+2];
            state.normals(3*v_off+0) = verts[6*j+3];
            state.normals(3*v_off+1) = verts[6*j+4];
            state.normals(3*v_off+2) = verts[6*j+5];
            v_off++;
         }

         n = elem_types[i];
         Array<int> &elems = *elements[i];
         num_elem = elems.Size()/n;
         // int attr = 1;
         int attr = i + 1;
         if (n == 3)
            for (int j = 0; j < num_elem; j++)
            {
               state.mesh->AddTriangle(&elems[3*j], attr);
            }
         else
            for (int j = 0; j < num_elem; j++)
            {
               state.mesh->AddQuad(&elems[4*j], attr);
            }
      }

      if (mesh_type == 1)
      {
         state.mesh->FinalizeTriMesh(1, 0, state.fix_elem_orient);
      }
      else if (mesh_type == 2)
      {
         state.mesh->FinalizeQuadMesh(1, 0, state.fix_elem_orient);
      }
      else
      {
         mfem_error("Input data contains mixture of triangles and quads!");
      }

      state.mesh->GenerateBoundaryElements();

      for (int i = num_patches; i > 0; )
      {
         i--;
         delete elements[i];
         delete vertices[i];
      }

      field_type = 0;
   }
   else
   {
      field_type = -1;
      cerr << "Unknown data format" << endl;
      cerr << data_type << endl;
   }

   if (field_type >= 0 && field_type <= 2)
   {
      if (state.grid_f)
      {
         Extrude1DMeshAndSolution(&state.mesh, &state.grid_f, NULL);
      }
      else
      {
         Extrude1DMeshAndSolution(&state.mesh, NULL, &state.sol);
      }
   }

   return field_type;
}
