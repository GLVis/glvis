// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
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

#include <array>
#include <algorithm>

using namespace std;
using namespace mfem;

enum class Command
{
   Mesh,
   Solution,
   Quadrature,
   Fem2D,
   VFem2D,
   VFem2D_keys,
   Fem3D,
   VFem3D,
   VFem3D_keys,
   Fem2D_GF,
   Fem2D_GF_keys,
   VFem2D_GF,
   VFem2D_GF_keys,
   Fem3D_GF,
   Fem3D_GF_keys,
   VFem3D_GF,
   VFem3D_GF_keys,
   RawScalar2D,
   //----------
   Max
};

class StreamCommands
{
   struct CmdItem
   {
      const char *keyword;
      bool keys;
      const char *params;
      const char *desc;

      bool operator==(const string &key) const { return key == keyword; }
   };
   array<CmdItem,(size_t)Command::Max> commands;

public:
   StreamCommands();

   decltype(commands)::const_iterator begin() const { return commands.begin(); }
   decltype(commands)::const_iterator end() const { return commands.end(); }
   CmdItem& operator[](Command cmd) { return commands[(size_t)cmd]; }
   const CmdItem& operator[](Command cmd) const { return commands[(size_t)cmd]; }
};
static const StreamCommands commands;

StreamCommands::StreamCommands()
{
   (*this)[Command::Mesh]                 = {"mesh", false, "<mesh>", "Visualize the mesh."};
   (*this)[Command::Solution]             = {"solution", false, "<mesh> <solution>", "Visualize the solution."};
   (*this)[Command::Quadrature]           = {"quadrature", false, "<mesh> <quadrature>", "Visualize the quadrature."};
   (*this)[Command::Fem2D]                = {"fem2d_data", false, "<mesh> <data>", "Visualize the 2D scalar data."};
   (*this)[Command::VFem2D]               = {"vfem2d_data", false, "<mesh> <data_x> <data_y>", "Visualize the 2D vector data."};
   (*this)[Command::VFem2D_keys]          = {"vfem2d_data_keys", true, "<mesh> <data_x> <data_y> <keys>", "Visualize the 2D vector data and apply control keys."};
   (*this)[Command::Fem3D]                = {"fem3d_data", false, "<mesh> <data>", "Visualize the 3D scalar data."};
   (*this)[Command::VFem3D]               = {"vfem3d_data", false, "<mesh> <data_x> <data_y> <data_z>", "Visualize the 3D vector data."};
   (*this)[Command::VFem3D_keys]          = {"vfem3d_data_keys", true, "<mesh> <data_x> <data_y> <data_z> <keys>", "Visualize the 3D vector data and apply control keys."};
   (*this)[Command::Fem2D_GF]             = {"fem2d_gf_data", false, "<mesh> <solution>", "Visualize the 2D scalar grid function."};
   (*this)[Command::Fem2D_GF_keys]        = {"fem2d_gf_data_keys", true, "<mesh> <solution> <keys>", "Visualize the 2D scalar grid function and apply control keys."};
   (*this)[Command::VFem2D_GF]            = {"vfem2d_gf_data", false, "<mesh> <solution>", "Visualize the 2D vector grid function."};
   (*this)[Command::VFem2D_GF_keys]       = {"vfem2d_gf_data_keys", true, "<mesh> <solution> <keys>", "Visualize the 2D vector grid function and apply control keys."};
   (*this)[Command::Fem3D_GF]             = {"fem3d_gf_data", false, "<mesh> <solution>", "Visualize the 3D scalar grid function."};
   (*this)[Command::Fem3D_GF_keys]        = {"fem3d_gf_data_keys", true, "<mesh> <solution> <keys>", "Visualize the 3D scalar grid function and apply control keys."};
   (*this)[Command::VFem3D_GF]            = {"vfem3d_gf_data", false, "<mesh> <solution>", "Visualize the 3D vector grid function."};
   (*this)[Command::VFem3D_GF_keys]       = {"vfem3d_gf_data_keys", true, "<mesh> <solution> <keys>", "Visualize the 3D vector grid function and apply control keys."};
   (*this)[Command::RawScalar2D]          = {"raw_scalar_2d", false, "<data>", "Visualize the 2D scalar data (see stream_reader.cpp)."};
}

void StreamReader::PrintCommands()
{
   cout << "Available commands are:" << endl;

   for (const auto &ci : commands)
   {
      cout << "\t" << ci.keyword << " " << ci.params << " - " << ci.desc << endl;
   }
}

bool StreamReader::SupportsDataType(const string &data_type)
{
   auto it = find(commands.begin(), commands.end(), data_type);
   return it != commands.end();
}

bool StreamReader::CheckStreamIsComplex(std::istream &solin,
                                        bool parallel)
{
   const char *header = (parallel)?("ParComplexGridFunction"):
                        ("ComplexGridFunction");
   solin >> ws;
   // Returning of the characters to stream(buffer) does not work reliably,
   // so compare only the initial character, which fortunately does not
   // coincide with the header of FiniteElementSpace.
   return (solin.peek() == header[0]);
}

int StreamReader::ReadStream(
   istream &is, const string &data_type)
{
   data.SetMesh(NULL);
   data.keys.clear();

   auto it = find(commands.begin(), commands.end(), data_type);
   if (it == commands.end())
   {
      cerr << "Unknown data format " << data_type << endl;
      PrintCommands();
      return 1;
   }

   Command cmd = (Command)(it - commands.begin());
   switch (cmd)
   {
      case Command::Fem2D:
      {
         Vector sol;
         data.SetMesh(new Mesh(is, 0, 0, data.fix_elem_orient));
         sol.Load(is, data.mesh->GetNV());
         data.SetScalarData(std::move(sol));
      }
      break;
      case Command::VFem2D:
      case Command::VFem2D_keys:
      {
         Vector solx, soly;
         data.SetMesh(new Mesh(is, 0, 0, data.fix_elem_orient));
         solx.Load(is, data.mesh->GetNV());
         soly.Load(is, data.mesh->GetNV());
         data.SetVectorData(std::move(solx), std::move(soly));
      }
      break;
      case Command::Fem3D:
      {
         Vector sol;
         data.SetMesh(new Mesh(is, 0, 0, data.fix_elem_orient));
         sol.Load(is, data.mesh->GetNV());
         data.SetScalarData(std::move(sol));
      }
      break;
      case Command::VFem3D:
      case Command::VFem3D_keys:
      {
         Vector solx, soly, solz;
         data.SetMesh(new Mesh(is, 0, 0, data.fix_elem_orient));
         solx.Load(is, data.mesh->GetNV());
         soly.Load(is, data.mesh->GetNV());
         solz.Load(is, data.mesh->GetNV());
         data.SetVectorData(std::move(solx), std::move(soly), std::move(solz));
      }
      break;
      case Command::Fem2D_GF:
      case Command::Fem2D_GF_keys:
      case Command::VFem2D_GF:
      case Command::VFem2D_GF_keys:
      case Command::Fem3D_GF:
      case Command::Fem3D_GF_keys:
      case Command::VFem3D_GF:
      case Command::VFem3D_GF_keys:
      case Command::Solution:
         data.SetMesh(new Mesh(is, 1, 0, data.fix_elem_orient));
         if (CheckStreamIsComplex(is))
         {
            data.SetCmplxGridFunction(new ComplexGridFunction(data.mesh.get(), is));
         }
         else
         {
            data.SetGridFunction(new GridFunction(data.mesh.get(), is));
         }
         break;
      case Command::Quadrature:
         data.SetMesh(new Mesh(is, 1, 0, data.fix_elem_orient));
         data.SetQuadFunction(new QuadratureFunction(data.mesh.get(), is));
         break;
      case Command::Mesh:
         data.SetMesh(new Mesh(is, 1, 0, data.fix_elem_orient));
         data.SetMeshSolution();
         break;
      case Command::RawScalar2D:
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
         Vector sol(tot_num_vert);
         Vector normals(3*tot_num_vert);

         int v_off = 0;
         for (int i = 0; i < num_patches; i++)
         {
            Array<double> &verts = *vertices[i];
            num_vert = verts.Size()/6;
            for (int j = 0; j < num_vert; j++)
            {
               data.mesh->AddVertex(&verts[6*j]);
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
            data.mesh->FinalizeTriMesh(1, 0, data.fix_elem_orient);
         }
         else if (mesh_type == 2)
         {
            data.mesh->FinalizeQuadMesh(1, 0, data.fix_elem_orient);
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

         data.SetScalarData(std::move(sol));
         data.SetNormals(std::move(normals));
      }
      break;
      case Command::Max: //dummy
         break;
   }

   if (commands[cmd].keys)
   {
      is >> data.keys;
   }

   data.ExtrudeMeshAndSolution();
   return 0;
}

int StreamReader::ReadStreams(const StreamCollection &input_streams)
{
   const int nproc = input_streams.size();
   std::vector<Mesh*> mesh_array(nproc);
   std::vector<GridFunction*> gf_array(nproc);
   std::vector<ComplexGridFunction*> cgf_array(nproc);
   std::vector<QuadratureFunction*> qf_array(nproc);

   std::string data_type;

   int gf_count = 0;
   int cgf_count = 0;
   int qf_count = 0;

   for (int p = 0; p < nproc; p++)
   {
#ifdef GLVIS_DEBUG
      cout << "connection[" << p << "]: reading data ... " << flush;
#endif
      istream &isock = *input_streams[p];
      // assuming the "parallel nproc p" part of the stream has been read
      isock >> ws >> data_type >> ws; // "*_data" / "mesh" / "solution"
#ifdef GLVIS_DEBUG
      cout << " type " << data_type << " ... " << flush;
#endif
      mesh_array[p] = new Mesh(isock, 1, 0, data.fix_elem_orient);
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
         if (CheckStreamIsComplex(isock, true))
         {
            isock.ignore(3);// ignore 'Par' prefix to load as serial
            cgf_array[p] = new ComplexGridFunction(mesh_array[p], isock);
            cgf_count++;
         }
         else
         {
            gf_array[p] = new GridFunction(mesh_array[p], isock);
            gf_count++;
         }
      }
#ifdef GLVIS_DEBUG
      cout << "done." << endl;
#endif
   }

   if ((gf_count > 0 && gf_count != nproc)
       || (cgf_count > 0 && cgf_count != nproc)
       || (qf_count > 0 && qf_count != nproc))
   {
      mfem_error("Input streams contain a mixture of data types!");
   }

   data.SetMesh(new Mesh(mesh_array.data(), nproc));
   if (gf_count > 0)
   {
      data.SetGridFunction(gf_array, nproc);
   }
   else if (cgf_count > 0)
   {
      data.SetCmplxGridFunction(cgf_array);
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
   else if (data.cgrid_f)
   {
      os << "solution\n";
      data.mesh->Print(os);
      data.cgrid_f->Save(os);
   }
   else if (data.grid_f)
   {
      os << "solution\n";
      data.mesh->Print(os);
      data.grid_f->Save(os);
   }
}
