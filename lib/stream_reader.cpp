// Copyright (c) 2010-2024, Lawrence Livermore National Security, LLC. Produced
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

/// Helper function for extrusion of 1D quadrature functions to 2D
QuadratureFunction* Extrude1DQuadFunction(Mesh *mesh, Mesh *mesh2d,
                                          QuadratureFunction *qf, int ny);

StreamState &StreamState::operator=(StreamState &&ss)
{
   internal = std::move(ss.internal);

   sol = std::move(ss.sol);
   solu = std::move(ss.solu);
   solv = std::move(ss.solv);
   solw = std::move(ss.solw);
   normals = std::move(ss.normals);
   keys = std::move(ss.keys);

   is_gf = ss.is_gf;
   is_qf = ss.is_qf;
   fix_elem_orient = ss.fix_elem_orient;
   save_coloring = ss.save_coloring;
   keep_attr = ss.keep_attr;

   quad_sol = ss.quad_sol;

   return *this;
}

void StreamState::SetMesh(mfem::Mesh *mesh_)
{
   internal.mesh.reset(mesh_);
   SetMesh(std::move(internal.mesh));
}

void StreamState::SetMesh(std::unique_ptr<mfem::Mesh> &&pmesh)
{
   internal.mesh = std::move(pmesh);
   internal.mesh_quad.reset();
   if (grid_f && grid_f->FESpace()->GetMesh() != mesh.get()) { SetGridFunction(NULL); }
   if (quad_f && quad_f->GetSpace()->GetMesh() != mesh.get()) { SetQuadFunction(NULL); }
}

void StreamState::SetGridFunction(mfem::GridFunction *gf)
{
   internal.grid_f.reset(gf);
   SetGridFunction(std::move(internal.grid_f));
}

void StreamState::SetGridFunction(std::unique_ptr<mfem::GridFunction> &&pgf)
{
   internal.grid_f = std::move(pgf);
   internal.quad_f.reset();
   quad_sol = QuadSolution::NONE;
}

void StreamState::SetQuadFunction(mfem::QuadratureFunction *qf)
{
   if (quad_f.get() != qf)
   {
      internal.grid_f.reset();
      quad_sol = QuadSolution::NONE;
   }
   internal.quad_f.reset(qf);
}

void StreamState::SetQuadFunction(std::unique_ptr<mfem::QuadratureFunction>
                                  &&pqf)
{
   if (quad_f.get() != pqf.get())
   {
      internal.grid_f.reset();
      quad_sol = QuadSolution::NONE;
   }
   internal.quad_f = std::move(pqf);
}

void StreamState::ExtrudeMeshAndSolution()
{
   Extrude1DMeshAndSolution();
   Extrude2D3VMeshAndSolution();
}

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

      internal.grid_f.reset(grid_f_2d);
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

   if (!mesh_quad) { internal.mesh.swap(internal.mesh_quad); }
   internal.mesh.reset(mesh2d);
}

void StreamState::Extrude2D3VMeshAndSolution()
{
   if (mesh->SpaceDimension() == 3 || !grid_f || grid_f->VectorDim() < 3) { return; }

   // not all vector elements can be embedded in 3D, so conversion to L2 elements
   // is performed already here
   ProjectVectorFEGridFunction();

   Mesh *mesh3d = new Mesh(*mesh);
   const FiniteElementSpace *fesx = mesh->GetNodalFESpace();
   const FiniteElementCollection *fecx = (fesx)?(fesx->FEColl()):(NULL);
   if (fecx)
   {
      mesh3d->SetCurvature(fecx->GetOrder(),
                           fecx->GetContType() == FiniteElementCollection::DISCONTINUOUS,
                           3);
   }
   else
   {
      mesh3d->SetCurvature(1, false, 3);
   }

   FiniteElementSpace *fes2d = grid_f->FESpace();
   FiniteElementSpace *fes3d = new FiniteElementSpace(*fes2d, mesh3d);
   GridFunction *gf3d = new GridFunction(fes3d);
   *gf3d = *grid_f;
   grid_f->MakeOwner(NULL);
   gf3d->MakeOwner(const_cast<FiniteElementCollection*>(fes2d->FEColl()));

   SetGridFunction(gf3d);
   delete fes2d;
   SetMesh(mesh3d);
}

void StreamState::CollectQuadratures(QuadratureFunction *qf_array[],
                                     int npieces)
{
   // assume the same vdim
   const int vdim = qf_array[0]->GetVDim();
   // assume the same quadrature rule
   QuadratureSpace *qspace = new QuadratureSpace(*mesh,
                                                 qf_array[0]->GetIntRule(0));
   SetQuadFunction(new QuadratureFunction(qspace, vdim));
   quad_f->SetOwnsSpace(true);
   real_t *g_data = quad_f->GetData();
   for (int p = 0; p < npieces; p++)
   {
      const real_t *l_data = qf_array[p]->GetData();
      const int l_size = qf_array[p]->Size();
      MFEM_ASSERT(g_data + l_size <= quad_f->GetData() + quad_f->Size(),
                  "Local parts do not fit to the global quadrature function!");
      memcpy(g_data, l_data, l_size * sizeof(real_t));
      g_data += l_size;
   }
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
      SetGridFunction(new GridFunction(cfes));
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

void StreamState::SetQuadSolution(QuadSolution type)
{
   // assume identical order
   const int order = quad_f->GetIntRule(0).GetOrder()/2; // <-- Gauss-Legendre
   // use the original mesh when available
   if (mesh_quad.get())
   {
      internal.mesh.swap(internal.mesh_quad);
      internal.mesh_quad.reset();
   }

   // check for tensor-product basis
   if (order > 0 && type != QuadSolution::HO_L2_projected)
   {
      Array<Geometry::Type> geoms;
      mesh->GetGeometries(mesh->Dimension(), geoms);
      for (auto geom : geoms)
         if (!Geometry::IsTensorProduct(geom))
         {
            cout << "High-order quadrature data are supported only for "
                 << "tensor-product finite elements with this representation"
                 << endl;
            SetQuadSolution(QuadSolution::HO_L2_projected);
            return;
         }
   }

   switch (type)
   {
      case QuadSolution::LOR_ClosedGL:
      {
         const int ref_factor = order + 1;
         if (ref_factor <= 1)
         {
            SetQuadSolution(QuadSolution::HO_L2_collocated); // low-order
            return;
         }
         Mesh *mesh_lor = new Mesh(Mesh::MakeRefined(*mesh, ref_factor,
                                                     BasisType::ClosedGL));
         FiniteElementCollection *fec = new L2_FECollection(0, mesh->Dimension());
         FiniteElementSpace *fes = new FiniteElementSpace(mesh_lor, fec,
                                                          quad_f->GetVDim(), Ordering::byVDIM);
         MFEM_ASSERT(quad_f->Size() == fes->GetVSize(), "Size mismatch");
         cout << "Representing quadrature data by piecewise-constant function on mesh refined "
              << ref_factor << " times" << endl;
         GridFunction *gf = new GridFunction(fes, *quad_f, 0);
         gf->MakeOwner(fec);
         internal.grid_f.reset(gf);
         internal.mesh.swap(internal.mesh_quad);
         internal.mesh.reset(mesh_lor);
      }
      break;
      case QuadSolution::HO_L2_collocated:
      {
         FiniteElementCollection *fec = new L2_FECollection(order, mesh->Dimension());
         FiniteElementSpace *fes = new FiniteElementSpace(mesh.get(), fec,
                                                          quad_f->GetVDim(), Ordering::byVDIM);
         MFEM_ASSERT(quad_f->Size() == fes->GetVSize(), "Size mismatch");
         cout << "Representing quadrature data by grid function of order "
              << order << endl;
         GridFunction *gf = new GridFunction(fes, *quad_f, 0);
         gf->MakeOwner(fec);
         internal.grid_f.reset(gf);
      }
      break;
      case QuadSolution::HO_L2_projected:
      {
         FiniteElementCollection *fec = new L2_FECollection(order, mesh->Dimension());
         FiniteElementSpace *fes = new FiniteElementSpace(mesh.get(), fec,
                                                          quad_f->GetVDim(), Ordering::byVDIM);
         cout << "Projecting quadrature data to grid function of order "
              << order << endl;
         GridFunction *gf = new GridFunction(fes);
         gf->MakeOwner(fec);

         VectorQuadratureFunctionCoefficient coeff_qf(*quad_f);
         VectorDomainLFIntegrator bi(coeff_qf);
         VectorMassIntegrator Mi;
         Mi.SetVDim(quad_f->GetVDim());

         DenseMatrix M_i;
         Vector x_i, b_i;
         DenseMatrixInverse M_i_inv;
         Array<int> vdofs;

         for (int i = 0; i < mesh->GetNE(); i++)
         {
            const FiniteElement *fe = fes->GetFE(i);
            ElementTransformation *Tr = mesh->GetElementTransformation(i);
            Mi.AssembleElementMatrix(*fe, *Tr, M_i);
            M_i_inv.Factor(M_i);

            bi.SetIntRule(&quad_f->GetIntRule(i));
            bi.AssembleRHSElementVect(*fe, *Tr, b_i);

            x_i.SetSize(b_i.Size());
            M_i_inv.Mult(b_i, x_i);

            fes->GetElementVDofs(i, vdofs);
            gf->SetSubVector(vdofs, x_i);
         }

         internal.grid_f.reset(gf);
      }
      break;
      default:
         cout << "Unknown quadrature data representation" << endl;
         return;
   }

   quad_sol = type;
}

void StreamState::SwitchQuadSolution(QuadSolution type, VisualizationScene *vs)
{
   unique_ptr<Mesh> old_mesh;
   // we must backup the refined mesh to prevent its deleting
   if (mesh_quad.get())
   {
      internal.mesh.swap(internal.mesh_quad);
      internal.mesh_quad.swap(old_mesh);
   }
   SetQuadSolution(type);
   ExtrudeMeshAndSolution();
   ResetMeshAndSolution(*this, vs);
}

// Read the content of an input stream (e.g. from socket/file)
StreamState::FieldType StreamState::ReadStream(istream &is,
                                               const string &data_type)
{
   FieldType field_type = FieldType::SCALAR;
   keys.clear();
   if (data_type == "fem2d_data")
   {
      SetMesh(new Mesh(is, 0, 0, fix_elem_orient));
      sol.Load(is, mesh->GetNV());
   }
   else if (data_type == "vfem2d_data" || data_type == "vfem2d_data_keys")
   {
      field_type = FieldType::VECTOR;
      SetMesh(new Mesh(is, 0, 0, fix_elem_orient));
      solu.Load(is, mesh->GetNV());
      solv.Load(is, mesh->GetNV());
      if (data_type == "vfem2d_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem3d_data")
   {
      SetMesh(new Mesh(is, 0, 0, fix_elem_orient));
      sol.Load(is, mesh->GetNV());
   }
   else if (data_type == "vfem3d_data" || data_type == "vfem3d_data_keys")
   {
      field_type = FieldType::VECTOR;
      SetMesh(new Mesh(is, 0, 0, fix_elem_orient));
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
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetGridFunction(new GridFunction(mesh.get(), is));
      if (data_type == "fem2d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "vfem2d_gf_data" || data_type == "vfem2d_gf_data_keys")
   {
      field_type = FieldType::VECTOR;
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetGridFunction(new GridFunction(mesh.get(), is));
      if (data_type == "vfem2d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem3d_gf_data" || data_type == "fem3d_gf_data_keys")
   {
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetGridFunction(new GridFunction(mesh.get(), is));
      if (data_type == "fem3d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "vfem3d_gf_data" || data_type == "vfem3d_gf_data_keys")
   {
      field_type = FieldType::VECTOR;
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetGridFunction(new GridFunction(mesh.get(), is));
      if (data_type == "vfem3d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "solution")
   {
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetGridFunction(new GridFunction(mesh.get(), is));
      field_type = (grid_f->VectorDim() == 1) ? FieldType::SCALAR : FieldType::VECTOR;
   }
   else if (data_type == "quadrature")
   {
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetQuadFunction(new QuadratureFunction(mesh.get(), is));
      SetQuadSolution();
      field_type = (quad_f->GetVDim() == 1) ? FieldType::SCALAR : FieldType::VECTOR;
   }
   else if (data_type == "mesh")
   {
      SetMesh(new Mesh(is, 1, 0, fix_elem_orient));
      SetMeshSolution();
      field_type = FieldType::MESH;
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

      SetMesh(new Mesh(2, tot_num_vert, tot_num_elem, 0));
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

      field_type = FieldType::SCALAR;
   }
   else
   {
      field_type = FieldType::UNKNOWN;
      cerr << "Unknown data format" << endl;
      cerr << data_type << endl;
   }

   if (field_type > FieldType::MIN && field_type < FieldType::MAX)
   {
      ExtrudeMeshAndSolution();
   }

   return field_type;
}

StreamState::FieldType StreamState::ReadStreams(const StreamCollection&
                                                input_streams)
{
   const int nproc = input_streams.size();
   Array<Mesh *> mesh_array(nproc);
   Array<GridFunction *> gf_array(nproc);
   Array<QuadratureFunction *> qf_array(nproc);

   std::string data_type;

   int gf_count = 0;
   int qf_count = 0;
   FieldType field_type = FieldType::SCALAR;

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
      mesh_array[p] = new Mesh(isock, 1, 0, fix_elem_orient);
      if (!keep_attr)
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

   SetMesh(new Mesh(mesh_array, nproc));
   if (gf_count > 0)
   {
      SetGridFunction(new GridFunction(mesh.get(), gf_array, nproc));
      field_type = (grid_f->VectorDim() == 1) ? FieldType::SCALAR : FieldType::VECTOR;
   }
   else if (qf_count > 0)
   {
      CollectQuadratures(qf_array, nproc);
      SetQuadSolution();
      field_type = (quad_f->GetVDim() == 1) ? FieldType::SCALAR : FieldType::VECTOR;
   }
   else
   {
      SetMeshSolution();
      field_type = FieldType::MESH;
   }

   for (int p = 0; p < nproc; p++)
   {
      delete mesh_array[nproc-1-p];
      delete gf_array[nproc-1-p];
   }

   ExtrudeMeshAndSolution();

   return field_type;
}

void StreamState::WriteStream(std::ostream &os)
{
   os.precision(8);
   if (quad_f)
   {
      os << "quadrature\n";
      if (mesh_quad.get())
      {
         mesh_quad->Print(os);
      }
      else
      {
         mesh->Print(os);
      }
      quad_f->Save(os);
   }
   else if (grid_f)
   {
      os << "solution\n";
      mesh->Print(os);
      grid_f->Save(os);
   }
}

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
std::unique_ptr<GridFunction>
StreamState::ProjectVectorFEGridFunction(std::unique_ptr<GridFunction> gf)
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
      ResetMeshAndSolution(new_state, vs);

      internal.grid_f = std::move(new_state.internal.grid_f);
      internal.mesh = std::move(new_state.internal.mesh);
      internal.quad_f = std::move(new_state.internal.quad_f);
      internal.mesh_quad = std::move(new_state.internal.mesh_quad);

      return true;
   }
   else
   {
      return false;
   }
}

void StreamState::ResetMeshAndSolution(StreamState &ss, VisualizationScene* vs)
{
   if (ss.mesh->SpaceDimension() == 2)
   {
      if (ss.grid_f->VectorDim() == 1)
      {
         VisualizationSceneSolution *vss =
            dynamic_cast<VisualizationSceneSolution *>(vs);
         ss.grid_f->GetNodalValues(ss.sol);
         vss->NewMeshAndSolution(ss.mesh.get(), ss.mesh_quad.get(), &ss.sol,
                                 ss.grid_f.get());
      }
      else
      {
         VisualizationSceneVector *vsv =
            dynamic_cast<VisualizationSceneVector *>(vs);
         vsv->NewMeshAndSolution(*ss.grid_f, ss.mesh_quad.get());
      }
   }
   else
   {
      if (ss.grid_f->VectorDim() == 1)
      {
         VisualizationSceneSolution3d *vss =
            dynamic_cast<VisualizationSceneSolution3d *>(vs);
         ss.grid_f->GetNodalValues(ss.sol);
         vss->NewMeshAndSolution(ss.mesh.get(), ss.mesh_quad.get(), &ss.sol,
                                 ss.grid_f.get());
      }
      else
      {
         ss.ProjectVectorFEGridFunction();

         VisualizationSceneVector3d *vss =
            dynamic_cast<VisualizationSceneVector3d *>(vs);
         vss->NewMeshAndSolution(ss.mesh.get(), ss.mesh_quad.get(), ss.grid_f.get());
      }
   }
}

QuadratureFunction *Extrude1DQuadFunction(Mesh *mesh, Mesh *mesh2d,
                                          QuadratureFunction *qf, int ny)
{
   // assume identical orders
   const int order = qf->GetIntRule(0).GetOrder();
   const int vdim = qf->GetVDim();
   QuadratureSpace *qspace2d = new QuadratureSpace(mesh2d, order);
   QuadratureFunction *qf2d = new QuadratureFunction(qspace2d);
   qf2d->SetOwnsSpace(true);

   DenseMatrix vals, vals2d;
   for (int ix = 0; ix < mesh->GetNE(); ix++)
   {
      qf->GetValues(ix, vals);
      const int nq = vals.Width();
      for (int iy = 0; iy < ny; iy++)
      {
         qf2d->GetValues(ix*ny+iy, vals2d);
         for (int qx = 0; qx < nq; qx++)
            for (int qy = 0; qy < nq; qy++)
               for (int d = 0; d < vdim; d++)
               {
                  vals2d(d, qy*nq+qx) = vals(d, qx);
               }
      }
   }

   return qf2d;
}
