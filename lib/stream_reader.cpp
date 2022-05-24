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

#include "stream_reader.hpp"
#include "visual.hpp"

#include <cstdlib>

using namespace std;
using namespace mfem;

void StreamState::Extrude1DMeshAndSolution()
{
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

   Mesh *mesh2d = Extrude1D(mesh.get(), 1, 0.1*(xmax - xmin));

   if (grid_f)
   {
      GridFunction *grid_f_2d =
         Extrude1DGridFunction(mesh.get(), mesh2d, grid_f.get(), 1);

      grid_f.reset(grid_f_2d);
   }
   else if (sol.Size() == mesh->GetNV())
   {
      Vector sol2d(mesh2d->GetNV());
      for (int i = 0; i < mesh->GetNV(); i++)
      {
         sol2d(2*i+0) = sol2d(2*i+1) = sol(i);
      }
      sol = sol2d;
   }

   mesh.reset(mesh2d);
}


void StreamState::SetMeshSolution()
{
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
      FiniteElementSpace *cfes = new FiniteElementSpace(mesh.get(), cfec);
      grid_f.reset(new GridFunction(cfes));
      grid_f->MakeOwner(cfec);
      {
         Array<int> coloring;
         srand(time(0));
         double a = double(rand()) / (double(RAND_MAX) + 1.);
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
      grid_f->GetNodalValues(sol);
      is_gf = 1;
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
      sol.SetSize (mesh -> GetNV());
      sol = 0.0;
   }
}


// Read the content of an input stream (e.g. from socket/file)
int StreamState::ReadStream(istream &is, const string &data_type)
{
   // 0 - scalar data, 1 - vector data, 2 - mesh only, (-1) - unknown
   int field_type = 0;
   keys.clear();
   if (data_type == "fem2d_data")
   {
      mesh.reset(new Mesh(is, 0, 0, fix_elem_orient));
      sol.Load(is, mesh->GetNV());
   }
   else if (data_type == "vfem2d_data" || data_type == "vfem2d_data_keys")
   {
      field_type = 1;
      mesh.reset(new Mesh(is, 0, 0, fix_elem_orient));
      solu.Load(is, mesh->GetNV());
      solv.Load(is, mesh->GetNV());
      if (data_type == "vfem2d_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem3d_data")
   {
      mesh.reset(new Mesh(is, 0, 0, fix_elem_orient));
      sol.Load(is, mesh->GetNV());
   }
   else if (data_type == "vfem3d_data" || data_type == "vfem3d_data_keys")
   {
      field_type = 1;
      mesh.reset(new Mesh(is, 0, 0, fix_elem_orient));
      solu.Load(is, mesh->GetNV());
      solv.Load(is, mesh->GetNV());
      solw.Load(is, mesh->GetNV());
      if (data_type == "vfem3d_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem2d_gf_data" || data_type == "fem2d_gf_data_keys")
   {
      mesh.reset(new Mesh(is, 1, 0, fix_elem_orient));
      grid_f.reset(new GridFunction(mesh.get(), is));
      if (data_type == "fem2d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "vfem2d_gf_data" || data_type == "vfem2d_gf_data_keys")
   {
      field_type = 1;
      mesh.reset(new Mesh(is, 1, 0, fix_elem_orient));
      grid_f.reset(new GridFunction(mesh.get(), is));
      if (data_type == "vfem2d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem3d_gf_data" || data_type == "fem3d_gf_data_keys")
   {
      mesh.reset(new Mesh(is, 1, 0, fix_elem_orient));
      grid_f.reset(new GridFunction(mesh.get(), is));
      if (data_type == "fem3d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "vfem3d_gf_data" || data_type == "vfem3d_gf_data_keys")
   {
      field_type = 1;
      mesh.reset(new Mesh(is, 1, 0, fix_elem_orient));
      grid_f.reset(new GridFunction(mesh.get(), is));
      if (data_type == "vfem3d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "solution")
   {
      mesh.reset(new Mesh(is, 1, 0, fix_elem_orient));
      grid_f.reset(new GridFunction(mesh.get(), is));
      field_type = (grid_f->VectorDim() == 1) ? 0 : 1;
   }
   else if (data_type == "mesh")
   {
      mesh.reset(new Mesh(is, 1, 0, fix_elem_orient));
      SetMeshSolution();
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

      mesh.reset(new Mesh(2, tot_num_vert, tot_num_elem, 0));
      sol.SetSize(tot_num_vert);
      normals.SetSize(3*tot_num_vert);

      int v_off = 0;
      for (int i = 0; i < num_patches; i++)
      {
         Array<double> &verts = *vertices[i];
         num_vert = verts.Size()/6;
         for (int j = 0; j < num_vert; j++)
         {
            mesh->AddVertex(&verts[6*j]);
            sol(v_off) = verts[6*j+2];
            normals(3*v_off+0) = verts[6*j+3];
            normals(3*v_off+1) = verts[6*j+4];
            normals(3*v_off+2) = verts[6*j+5];
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
               mesh->AddTriangle(&elems[3*j], attr);
            }
         else
            for (int j = 0; j < num_elem; j++)
            {
               mesh->AddQuad(&elems[4*j], attr);
            }
      }

      if (mesh_type == 1)
      {
         mesh->FinalizeTriMesh(1, 0, fix_elem_orient);
      }
      else if (mesh_type == 2)
      {
         mesh->FinalizeQuadMesh(1, 0, fix_elem_orient);
      }
      else
      {
         mfem_error("Input data contains mixture of triangles and quads!");
      }

      mesh->GenerateBoundaryElements();

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
      Extrude1DMeshAndSolution();
   }

   return field_type;
}

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
std::unique_ptr<GridFunction>
ProjectVectorFEGridFunction(std::unique_ptr<GridFunction> gf)
{
   if ((gf->VectorDim() == 3) && (gf->FESpace()->GetVDim() == 1))
   {
      int p = gf->FESpace()->GetOrder(0);
      cout << "Switching to order " << p
           << " discontinuous vector grid function..." << endl;
      int dim = gf->FESpace()->GetMesh()->Dimension();
      FiniteElementCollection *d_fec =
         (p == 1 && dim == 3) ?
         (FiniteElementCollection*)new LinearDiscont3DFECollection :
         (FiniteElementCollection*)new L2_FECollection(p, dim, 1);
      FiniteElementSpace *d_fespace =
         new FiniteElementSpace(gf->FESpace()->GetMesh(), d_fec, 3);
      GridFunction *d_gf = new GridFunction(d_fespace);
      d_gf->MakeOwner(d_fec);
      gf->ProjectVectorFieldOn(*d_gf);
      gf.reset(d_gf);
   }
   return gf;
}

bool StreamState::SetNewMeshAndSolution(StreamState new_state,
                                        VisualizationScene* vs)
{
   if (new_state.mesh->SpaceDimension() == mesh->SpaceDimension() &&
       new_state.grid_f->VectorDim() == grid_f->VectorDim())
   {
      std::unique_ptr<mfem::Mesh> new_m = std::move(new_state.mesh);
      std::unique_ptr<mfem::GridFunction> new_g = std::move(new_state.grid_f);
      if (new_m->SpaceDimension() == 2)
      {
         if (new_g->VectorDim() == 1)
         {
            VisualizationSceneSolution *vss =
               dynamic_cast<VisualizationSceneSolution *>(vs);
            new_g->GetNodalValues(sol);
            vss->NewMeshAndSolution(new_m.get(), &sol, new_g.get());
         }
         else
         {
            VisualizationSceneVector *vsv =
               dynamic_cast<VisualizationSceneVector *>(vs);
            vsv->NewMeshAndSolution(*new_g);
         }
      }
      else
      {
         if (new_g->VectorDim() == 1)
         {
            VisualizationSceneSolution3d *vss =
               dynamic_cast<VisualizationSceneSolution3d *>(vs);
            new_g->GetNodalValues(sol);
            vss->NewMeshAndSolution(new_m.get(), &sol, new_g.get());
         }
         else
         {
            new_g = ProjectVectorFEGridFunction(std::move(new_g));

            VisualizationSceneVector3d *vss =
               dynamic_cast<VisualizationSceneVector3d *>(vs);
            vss->NewMeshAndSolution(new_m.get(), new_g.get());
         }
      }
      grid_f = std::move(new_g);
      mesh = std::move(new_m);
      return true;
   }
   else
   {
      return false;
   }
}
