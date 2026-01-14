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

#include "file_reader.hpp"

#include <vector>
#include <general/text.hpp>

using namespace std;
using namespace mfem;

int FileReader::ReadSerial(FileReader::FileType ft, const char *mesh_file,
                           const char *sol_file, int component)
{
   // get the mesh from a file
   named_ifgzstream meshin(mesh_file);
   if (!meshin)
   {
      cerr << "Can not open mesh file " << mesh_file << ". Exit.\n";
      return 1;
   }

   data.SetMesh(new Mesh(meshin, 1, 0, data.fix_elem_orient));

   if (ft == FileType::MESH)
   {
      data.SetMeshSolution();
   }
   else
   {
      // get the solution from file
      shared_ptr<ifgzstream> solin;
      if (!strcmp(mesh_file, sol_file))
      {
         solin.reset(&meshin, [](ifgzstream*) {}); // no op deleter
      }
      else
      {
         solin.reset(new ifgzstream(sol_file));
         if (!(*solin))
         {
            cerr << "Can not open solution file " << sol_file << ". Exit.\n";
            return 2;
         }
      }

      switch (ft)
      {
         case FileType::GRID_FUNC:
         {
            if (CheckStreamIsComplex(*solin))
            {
               data.SetCmplxGridFunction(new ComplexGridFunction(data.mesh.get(), *solin),
                                         component);
            }
            else
            {
               data.SetGridFunction(new GridFunction(data.mesh.get(), *solin), component);
            }
         }
         break;
         case FileType::QUAD_FUNC:
            data.SetQuadFunction(new QuadratureFunction(data.mesh.get(), *solin),
                                 component);
            break;
         case FileType::SCALAR_SOL:
         {
            // get rid of NetGen's info line
            char buff[128];
            solin->getline(buff,128);
            data.sol.Load(*solin, data.mesh->GetNV());
            data.type = DataState::FieldType::SCALAR;
         }
         break;
         case FileType::VECTOR_SOL:
            data.solu.Load(*solin, data.mesh->GetNV());
            data.solv.Load(*solin, data.mesh->GetNV());
            if (data.mesh->SpaceDimension() == 3)
            {
               data.solw.Load(*solin, data.mesh->GetNV());
            }
            data.type = DataState::FieldType::VECTOR;
            break;
         default:
            cerr << "Unknown file type. Exit" << endl;
            exit(1);
      }
   }

   data.ExtrudeMeshAndSolution();
   return 0;
}

int FileReader::ReadParallel(int np, FileType ft, const char *mesh_file,
                             const char *sol_file, int component)
{
   int read_err;

   switch (ft)
   {
      case FileType::MESH:
         read_err = ReadParMeshAndGridFunction(np, mesh_file, NULL);
         break;
      case FileType::GRID_FUNC:
         read_err = ReadParMeshAndGridFunction(np, mesh_file, sol_file, component);
         break;
      case FileType::QUAD_FUNC:
         read_err = ReadParMeshAndQuadFunction(np, mesh_file, sol_file, component);
         break;
      default:
         cerr << "Unknown file type. Exit" << endl;
         exit(1);
   }

   if (!read_err)
   {
      data.ExtrudeMeshAndSolution();
   }

   return read_err;
}

bool FileReader::CheckStreamIsComplex(std::istream &solin, bool parallel) const
{
   string buff;
   auto pos = solin.tellg();
   solin >> ws;
   getline(solin, buff);
   solin.seekg(pos);
   mfem::filter_dos(buff);
   const char *header = (parallel)?("ParComplexGridFunction"):
                        ("ComplexGridFunction");
   return (buff == header);
}

int FileReader::ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                                           const char *sol_prefix, int component)
{
   data.SetMesh(NULL);

   // are the solutions bundled together with the mesh files?
   bool same_file = false;
   if (sol_prefix)
   {
      same_file = !strcmp(sol_prefix, mesh_prefix);
      data.SetGridFunction(NULL);
   }

   std::vector<Mesh *> mesh_array(np);
   std::vector<GridFunction *> gf_array(np);
   std::vector<ComplexGridFunction *> cgf_array(np);

   int read_err = 0;
   for (int p = 0; p < np; p++)
   {
      ostringstream fname;
      fname << mesh_prefix << '.' << setfill('0') << setw(pad_digits) << p;
      named_ifgzstream meshfile(fname.str().c_str());
      if (!meshfile)
      {
         cerr << "Could not open mesh file: " << fname.str() << '!' << endl;
         read_err = 1;
         break;
      }

      mesh_array[p] = new Mesh(meshfile, 1, 0, data.fix_elem_orient);

      if (!data.keep_attr)
      {
         // set element and boundary attributes to be the processor number + 1
         for (int i = 0; i < mesh_array[p]->GetNE(); i++)
         {
            mesh_array[p]->GetElement(i)->SetAttribute(p+1);
         }
         for (int i = 0; i < mesh_array[p]->GetNBE(); i++)
         {
            mesh_array[p]->GetBdrElement(i)->SetAttribute(p+1);
         }
      }

      // read the solution
      if (sol_prefix)
      {
         if (!same_file)
         {
            ostringstream sol_fname;
            sol_fname << sol_prefix << '.' << setfill('0') << setw(pad_digits) << p;
            ifgzstream solfile(sol_fname.str().c_str());
            if (!solfile)
            {
               cerr << "Could not open solution file "
                    << sol_fname.str() << '!' << endl;
               read_err = 2;
               break;
            }

            if (CheckStreamIsComplex(solfile, true))
            {
               solfile >> ws;
               solfile.ignore(3);// ignore 'Par' prefix to load as serial
               cgf_array[p] = new ComplexGridFunction(mesh_array[p], solfile);
            }
            else
            {
               gf_array[p] = new GridFunction(mesh_array[p], solfile);
            }
         }
         else  // mesh and solution in the same file
         {
            if (CheckStreamIsComplex(meshfile))
            {
               cgf_array[p] = new ComplexGridFunction(mesh_array[p], meshfile);
            }
            else
            {
               gf_array[p] = new GridFunction(mesh_array[p], meshfile);
            }
         }
      }
   }

   if (!read_err)
   {
      bool bcmplx{}, breal{};
      for (int p = 0; p < np; p++)
      {
         if (cgf_array[p]) { bcmplx = true; }
         if (gf_array[p]) { breal = true; }
      }
      if (bcmplx && breal)
      {
         cerr << "Inconsistent input files" << endl;
         read_err = 3;
      }
   }

   if (!read_err)
   {
      // create the combined mesh and gf
      data.SetMesh(new Mesh(mesh_array.data(), np));
      if (sol_prefix)
      {
         if (cgf_array[0])
         {
            std::vector<GridFunction *> r_array(np), i_array(np);
            for (int p = 0; p < np; p++)
            {
               r_array[p] = &(cgf_array[p]->real());
               i_array[p] = &(cgf_array[p]->imag());
            }
            GridFunction *rgf = new GridFunction(data.mesh.get(), r_array.data(), np);
            GridFunction *igf = new GridFunction(data.mesh.get(), i_array.data(), np);
            ComplexGridFunction *cgf = new ComplexGridFunction(rgf->FESpace());
            cgf->MakeOwner(rgf->OwnFEC());
            rgf->MakeOwner(NULL);
            cgf->real() = *rgf;
            cgf->imag() = *igf;
            delete rgf;
            delete igf;
            data.SetCmplxGridFunction(cgf, component);
         }
         else
            data.SetGridFunction(new GridFunction(data.mesh.get(), gf_array.data(), np),
                                 component);
      }
      else
      {
         data.SetMeshSolution();
      }
   }

   for (int p = 0; p < np; p++)
   {
      delete gf_array[np-1-p];
      delete cgf_array[np-1-p];
      delete mesh_array[np-1-p];
   }

   return read_err;
}

int FileReader::ReadParMeshAndQuadFunction(int np, const char *mesh_prefix,
                                           const char *sol_prefix, int component)
{
   data.SetMesh(NULL);

   // are the solutions bundled together with the mesh files?
   bool same_file = false;
   if (sol_prefix)
   {
      same_file = !strcmp(sol_prefix, mesh_prefix);
      data.SetGridFunction(NULL);
   }

   std::vector<Mesh*> mesh_array(np);
   std::vector<QuadratureFunction*> qf_array(np);

   int read_err = 0;
   for (int p = 0; p < np; p++)
   {
      ostringstream fname;
      fname << mesh_prefix << '.' << setfill('0') << setw(pad_digits) << p;
      named_ifgzstream meshfile(fname.str().c_str());
      if (!meshfile)
      {
         cerr << "Could not open mesh file: " << fname.str() << '!' << endl;
         read_err = 1;
         break;
      }

      mesh_array[p] = new Mesh(meshfile, 1, 0, data.fix_elem_orient);

      if (!data.keep_attr)
      {
         // set element and boundary attributes to be the processor number + 1
         for (int i = 0; i < mesh_array[p]->GetNE(); i++)
         {
            mesh_array[p]->GetElement(i)->SetAttribute(p+1);
         }
         for (int i = 0; i < mesh_array[p]->GetNBE(); i++)
         {
            mesh_array[p]->GetBdrElement(i)->SetAttribute(p+1);
         }
      }

      // read the solution
      if (sol_prefix)
      {
         if (!same_file)
         {
            ostringstream sol_fname;
            sol_fname << sol_prefix << '.' << setfill('0') << setw(pad_digits) << p;
            ifgzstream solfile(sol_fname.str().c_str());
            if (!solfile)
            {
               cerr << "Could not open solution file "
                    << sol_fname.str() << '!' << endl;
               read_err = 2;
               break;
            }

            qf_array[p] = new QuadratureFunction(mesh_array[p], solfile);
         }
         else  // mesh and solution in the same file
         {
            qf_array[p] = new QuadratureFunction(mesh_array[p], meshfile);
         }
      }
   }

   if (!read_err)
   {
      // create the combined mesh and gf
      data.SetMesh(new Mesh(mesh_array.data(), np));
      if (sol_prefix)
      {
         data.SetQuadFunction(qf_array, component);
      }
      else
      {
         data.SetMeshSolution();
      }
   }

   for (int p = 0; p < np; p++)
   {
      delete qf_array[np-1-p];
      delete mesh_array[np-1-p];
   }

   return read_err;
}
