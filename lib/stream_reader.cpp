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

#include "stream_reader.hpp"

using namespace std;
using namespace mfem;

int StreamReader::ReadStream(
   istream &is, const string &data_type)
{
   data.SetMesh(NULL);
   data.keys.clear();

   if (data_type == "fem2d_data")
   {
      data.type = DataState::FieldType::SCALAR;
      data.SetMesh(new Mesh(is, 0, 0, data.save_coloring));
      data.sol.Load(is, data.mesh->GetNV());
   }
   else if (data_type == "vfem2d_data" || data_type == "vfem2d_data_keys")
   {
      data.type = DataState::FieldType::VECTOR;
      data.SetMesh(new Mesh(is, 0, 0, data.save_coloring));
      data.solu.Load(is, data.mesh->GetNV());
      data.solv.Load(is, data.mesh->GetNV());
      if (data_type == "vfem2d_data_keys")
      {
         is >> data.keys;
      }
   }
   else if (data_type == "fem3d_data")
   {
      data.type = DataState::FieldType::SCALAR;
      data.SetMesh(new Mesh(is, 0, 0, data.save_coloring));
      data.sol.Load(is, data.mesh->GetNV());
   }
   else if (data_type == "vfem3d_data" || data_type == "vfem3d_data_keys")
   {
      data.type = DataState::FieldType::VECTOR;
      data.SetMesh(new Mesh(is, 0, 0, data.save_coloring));
      data.solu.Load(is, data.mesh->GetNV());
      data.solv.Load(is, data.mesh->GetNV());
      data.solw.Load(is, data.mesh->GetNV());
      if (data_type == "vfem3d_data_keys")
      {
         is >> data.keys;
      }
   }
   else if (data_type == "fem2d_gf_data" || data_type == "fem2d_gf_data_keys")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetGridFunction(new GridFunction(data.mesh.get(), is));
      if (data_type == "fem2d_gf_data_keys")
      {
         is >> data.keys;
      }
   }
   else if (data_type == "vfem2d_gf_data" || data_type == "vfem2d_gf_data_keys")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetGridFunction(new GridFunction(data.mesh.get(), is));
      if (data_type == "vfem2d_gf_data_keys")
      {
         is >> data.keys;
      }
   }
   else if (data_type == "fem3d_gf_data" || data_type == "fem3d_gf_data_keys")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetGridFunction(new GridFunction(data.mesh.get(), is));
      if (data_type == "fem3d_gf_data_keys")
      {
         is >> data.keys;
      }
   }
   else if (data_type == "vfem3d_gf_data" || data_type == "vfem3d_gf_data_keys")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetGridFunction(new GridFunction(data.mesh.get(), is));
      if (data_type == "vfem3d_gf_data_keys")
      {
         is >> data.keys;
      }
   }
   else if (data_type == "solution")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetGridFunction(new GridFunction(data.mesh.get(), is));
   }
   else if (data_type == "quadrature")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetQuadFunction(new QuadratureFunction(data.mesh.get(), is));
      data.SetQuadSolution();
   }
   else if (data_type == "mesh")
   {
      data.SetMesh(new Mesh(is, 1, 0, data.save_coloring));
      data.SetMeshSolution();
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

      data.SetMesh(new Mesh(2, tot_num_vert, tot_num_elem, 0));
      data.sol.SetSize(tot_num_vert);
      data.normals.SetSize(3*tot_num_vert);

      int v_off = 0;
      for (int i = 0; i < num_patches; i++)
      {
         Array<double> &verts = *vertices[i];
         num_vert = verts.Size()/6;
         for (int j = 0; j < num_vert; j++)
         {
            data.mesh->AddVertex(&verts[6*j]);
            data.sol(v_off) = verts[6*j+2];
            data.normals(3*v_off+0) = verts[6*j+3];
            data.normals(3*v_off+1) = verts[6*j+4];
            data.normals(3*v_off+2) = verts[6*j+5];
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
               data.mesh->AddTriangle(&elems[3*j], attr);
            }
         else
            for (int j = 0; j < num_elem; j++)
            {
               data.mesh->AddQuad(&elems[4*j], attr);
            }
      }

      if (mesh_type == 1)
      {
         data.mesh->FinalizeTriMesh(1, 0, data.save_coloring);
      }
      else if (mesh_type == 2)
      {
         data.mesh->FinalizeQuadMesh(1, 0, data.save_coloring);
      }
      else
      {
         mfem_error("Input data contains mixture of triangles and quads!");
      }

      data.mesh->GenerateBoundaryElements();

      for (int i = num_patches; i > 0; )
      {
         i--;
         delete elements[i];
         delete vertices[i];
      }

      data.type = DataState::FieldType::SCALAR;
   }
   else
   {
      cerr << "Unknown data format" << endl;
      cerr << data_type << endl;
      return 1;
   }

   data.ExtrudeMeshAndSolution();
   return 0;
}

int StreamReader::ReadStreams(const StreamCollection &input_streams)
{
   const int nproc = input_streams.size();
   std::vector<Mesh*> mesh_array(nproc);
   std::vector<GridFunction*> gf_array(nproc);
   std::vector<QuadratureFunction*> qf_array(nproc);

   std::string data_type;

   int gf_count = 0;
   int qf_count = 0;

   for (int p = 0; p < nproc; p++)
   {
#ifdef GLVIS_DEBUG
      cout << "connection[" << p << "]: reading initial data ... " << flush;
#endif
      istream &isock = *input_streams[p];
      // assuming the "parallel nproc p" part of the stream has been read
      isock >> ws >> data_type >> ws; // "*_data" / "mesh" / "solution"
#ifdef GLVIS_DEBUG
      cout << " type " << data_type << " ... " << flush;
#endif
      mesh_array[p] = new Mesh(isock, 1, 0, data.save_coloring);
      if (!data.keep_attr)
      {
         // set element and boundary attributes to proc+1
         for (int i = 0; i < mesh_array[p]->GetNE(); i++)
         {
            mesh_array[p]->GetElement(i)->SetAttribute(p+1);
         }
         for (int i = 0; i < mesh_array[p]->GetNBE(); i++)
         {
            mesh_array[p]->GetBdrElement(i)->SetAttribute(p+1);
         }
      }
      gf_array[p] = NULL;
      if (data_type == "quadrature")
      {
         qf_array[p] = new QuadratureFunction(mesh_array[p], isock);
         qf_count++;
      }
      else if (data_type != "mesh")
      {
         gf_array[p] = new GridFunction(mesh_array[p], isock);
         gf_count++;
      }
#ifdef GLVIS_DEBUG
      cout << "done." << endl;
#endif
   }

   if ((gf_count > 0 && gf_count != nproc)
       || (qf_count > 0 && qf_count != nproc))
   {
      mfem_error("Input streams contain a mixture of data types!");
   }

   data.SetMesh(new Mesh(mesh_array.data(), nproc));
   if (gf_count > 0)
   {
      data.SetGridFunction(new GridFunction(data.mesh.get(), gf_array.data(), nproc));
      if (!data.keep_attr) { data.ComputeDofsOffsets(gf_array); }
   }
   else if (qf_count > 0)
   {
      data.SetQuadFunction(qf_array);
   }
   else
   {
      data.SetMeshSolution();
   }

   for (int p = 0; p < nproc; p++)
   {
      delete mesh_array[nproc-1-p];
      delete gf_array[nproc-1-p];
   }

   data.ExtrudeMeshAndSolution();

   return 0;
}

void StreamReader::WriteStream(std::ostream &os)
{
   os.precision(8);
   if (data.quad_f)
   {
      os << "quadrature\n";
      if (data.mesh_quad.get())
      {
         data.mesh_quad->Print(os);
      }
      else
      {
         data.mesh->Print(os);
      }
      data.quad_f->Save(os);
   }
   else if (data.grid_f)
   {
      os << "solution\n";
      data.mesh->Print(os);
      data.grid_f->Save(os);
   }
}
