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

#include "data_state.hpp"
#include "vsvector.hpp"
#include "vsvector3d.hpp"

using namespace std;
using namespace mfem;

/// Class used for extruding vector GridFunctions
class VectorExtrudeCoefficient : public VectorCoefficient
{
private:
   int n;
   Mesh *mesh_in;
   VectorCoefficient &sol_in;
public:
   VectorExtrudeCoefficient(Mesh *m, VectorCoefficient &s, int n_)
      : VectorCoefficient(s.GetVDim()), n(n_), mesh_in(m), sol_in(s) { }
   void Eval(Vector &v, ElementTransformation &T,
             const IntegrationPoint &ip) override;
   virtual ~VectorExtrudeCoefficient() { }
};

/// Helper function for extrusion of 1D quadrature functions to 2D
QuadratureFunction* Extrude1DQuadFunction(Mesh *mesh, Mesh *mesh2d,
                                          QuadratureFunction *qf, int ny);

DataState &DataState::operator=(DataState &&ss)
{
   internal = std::move(ss.internal);

   type = ss.type;
   quad_sol = ss.quad_sol;

   sol = std::move(ss.sol);
   solu = std::move(ss.solu);
   solv = std::move(ss.solv);
   solw = std::move(ss.solw);
   normals = std::move(ss.normals);
   keys = std::move(ss.keys);

   fix_elem_orient = ss.fix_elem_orient;
   save_coloring = ss.save_coloring;
   keep_attr = ss.keep_attr;

   return *this;
}

void DataState::SetMesh(Mesh *mesh_)
{
   internal.mesh.reset(mesh_);
   SetMesh(std::move(internal.mesh));
}

void DataState::SetMesh(std::unique_ptr<Mesh> &&pmesh)
{
   internal.mesh = std::move(pmesh);
   internal.mesh_quad.reset();
   if (grid_f && grid_f->FESpace()->GetMesh() != mesh.get()) { SetGridFunction(NULL); }
   if (quad_f && quad_f->GetSpace()->GetMesh() != mesh.get()) { SetQuadFunction(NULL); }
}

void DataState::SetGridFunction(GridFunction *gf, int component)
{
   internal.grid_f.reset(gf);
   SetGridFunction(std::move(internal.grid_f), component);
}

void DataState::SetGridFunction(std::unique_ptr<GridFunction> &&pgf,
                                int component)
{
   internal.grid_f = std::move(pgf);
   internal.quad_f.reset();
   quad_sol = QuadSolution::NONE;
   SetGridFunctionSolution(component);
}

void DataState::SetGridFunction(std::vector<mfem::GridFunction*> &gf_array,
                                int num_pieces)
{
   SetGridFunction(new GridFunction(mesh.get(), gf_array.data(), num_pieces));
   if (!keep_attr) { _ComputeDofsOffsets(gf_array); }
}

void DataState::SetQuadFunction(QuadratureFunction *qf, int component)
{
   if (quad_f.get() != qf)
   {
      internal.grid_f.reset();
      quad_sol = QuadSolution::NONE;
   }
   internal.quad_f.reset(qf);
   SetQuadFunctionSolution(component);
}

void DataState::SetQuadFunction(std::unique_ptr<QuadratureFunction> &&pqf,
                                int component)
{
   if (quad_f.get() != pqf.get())
   {
      internal.grid_f.reset();
      quad_sol = QuadSolution::NONE;
   }
   internal.quad_f = std::move(pqf);
   SetQuadFunctionSolution(component);
}

void DataState::SetQuadFunction(const std::vector<QuadratureFunction*>
                                &qf_array, int component)
{
   // assume the same vdim
   const int vdim = qf_array[0]->GetVDim();
   // assume the same quadrature rule
   QuadratureSpace *qspace = new QuadratureSpace(*mesh,
                                                 qf_array[0]->GetIntRule(0));
   QuadratureFunction *qf = new QuadratureFunction(qspace, vdim);
   qf->SetOwnsSpace(true);
   real_t *g_data = qf->GetData();
   for (QuadratureFunction *qf_piece : qf_array)
   {
      const real_t *l_data = qf_piece->GetData();
      const int l_size = qf_piece->Size();
      MFEM_ASSERT(g_data + l_size <= qf->GetData() + qf->Size(),
                  "Local parts do not fit to the global quadrature function!");
      memcpy(g_data, l_data, l_size * sizeof(real_t));
      g_data += l_size;
   }
   SetQuadFunction(qf, component);
}

void DataState::SetDataCollectionField(DataCollection *dc, int ti,
                                       const char *field, bool quad, int component)
{
   internal.data_coll.reset(dc);
   data_coll->Load(ti);
   Mesh *dc_mesh = data_coll->GetMesh();
   SetMesh(new Mesh(*dc_mesh, false));

   if (field)
   {
      if (!quad)
      {
         GridFunction *gf = data_coll->GetField(field);
         if (!gf)
         {
            std::cerr << "Field " << field << " not found in the collection!" << std::endl;
            return;
         }
         SetGridFunction(new GridFunction(gf->FESpace(), *gf), component);
      }
      else
      {
         QuadratureFunction *qf = data_coll->GetQField(field);
         if (!qf)
         {
            std::cerr << "Quadrature field " << field << " not found in the collection!" <<
                      std::endl;
            return;
         }
         SetQuadFunction(new QuadratureFunction(qf->GetSpace(), qf->GetData(),
                                                qf->GetVDim()), component);
      }
   }
   else
   {
      SetMeshSolution();
   }
}

void DataState::ExtrudeMeshAndSolution()
{
   Extrude1DMeshAndSolution();
   Extrude2D3VMeshAndSolution();
}

void DataState::Extrude1DMeshAndSolution()
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
      if (grid_f->VectorDim() > 1)
      {
         ProjectVectorFEGridFunction();
      }

      GridFunction *grid_f_2d = Extrude1DGridFunction(mesh.get(), mesh2d,
                                                      grid_f.get(), 1);

      if (grid_f_2d->VectorDim() < grid_f->VectorDim())
      {
         // workaround for older MFEM where Extrude1DGridFunction()
         // does not work for vector grid functions
         delete grid_f_2d;
         FiniteElementCollection *fec2d = new L2_FECollection(
            grid_f->FESpace()->FEColl()->GetOrder(), 2);
         FiniteElementSpace *fes2d = new FiniteElementSpace(mesh2d, fec2d,
                                                            grid_f->FESpace()->GetVDim());
         grid_f_2d = new GridFunction(fes2d);
         grid_f_2d->MakeOwner(fec2d);
         VectorGridFunctionCoefficient vcsol(grid_f.get());
         ::VectorExtrudeCoefficient vc2d(mesh.get(), vcsol, 1);
         grid_f_2d->ProjectCoefficient(vc2d);
      }

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

void DataState::Extrude2D3VMeshAndSolution()
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

void DataState::SetMeshSolution()
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
   type = FieldType::MESH;
}

void DataState::SetGridFunctionSolution(int gf_component)
{
   if (!grid_f)
   {
      type = (mesh)?(FieldType::MESH):(FieldType::UNKNOWN);
      return;
   }

   if (gf_component != -1)
   {
      if (gf_component < 0 || gf_component >= grid_f->FESpace()->GetVDim())
      {
         cerr << "Invalid component " << gf_component << '.' << endl;
         exit(1);
      }
      FiniteElementSpace *ofes = grid_f->FESpace();
      FiniteElementCollection *fec =
         FiniteElementCollection::New(ofes->FEColl()->Name());
      FiniteElementSpace *fes = new FiniteElementSpace(mesh.get(), fec);
      GridFunction *new_gf = new GridFunction(fes);
      new_gf->MakeOwner(fec);
      for (int i = 0; i < new_gf->Size(); i++)
      {
         (*new_gf)(i) = (*grid_f)(ofes->DofToVDof(i, gf_component));
      }
      SetGridFunction(new_gf);
      return;
   }

   if (grid_f->VectorDim() == 1)
   {
      grid_f->GetNodalValues(sol);
      type = FieldType::SCALAR;
   }
   else
   {
      type = FieldType::VECTOR;
   }
}

void DataState::SetQuadFunctionSolution(int qf_component)
{
   if (!quad_f)
   {
      type = (mesh)?(FieldType::MESH):(FieldType::UNKNOWN);
      return;
   }

   const int vdim = quad_f->GetVDim();
   if (qf_component != -1)
   {
      if (qf_component < 0 || qf_component >= vdim)
      {
         cerr << "Invalid component " << qf_component << '.' << endl;
         exit(1);
      }
      QuadratureSpaceBase *qspace = quad_f->GetSpace();
      QuadratureFunction *new_qf = new QuadratureFunction(qspace);
      for (int i = 0; i < new_qf->Size(); i++)
      {
         (*new_qf)(i) = (*quad_f)(i * vdim + qf_component);
      }
      quad_f->SetOwnsSpace(false);
      new_qf->SetOwnsSpace(true);
      SetQuadFunction(new_qf);
      return;
   }

   if (vdim == 1)
   {
      type = FieldType::SCALAR;
   }
   else
   {
      type = FieldType::VECTOR;
   }

   SetQuadSolution();
}

void DataState::SetQuadSolution(QuadSolution quad_type)
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
   if (order > 0 && quad_type != QuadSolution::HO_L2_projected)
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

   switch (quad_type)
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

   quad_sol = quad_type;
}

void DataState::SwitchQuadSolution(QuadSolution quad_type,
                                   VisualizationScene *vs)
{
   unique_ptr<Mesh> old_mesh;
   // we must backup the refined mesh to prevent its deleting
   if (mesh_quad.get())
   {
      internal.mesh.swap(internal.mesh_quad);
      internal.mesh_quad.swap(old_mesh);
   }
   SetQuadSolution(quad_type);
   ExtrudeMeshAndSolution();
   ResetMeshAndSolution(*this, vs);
}

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
std::unique_ptr<GridFunction>
DataState::ProjectVectorFEGridFunction(std::unique_ptr<GridFunction> gf)
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

bool DataState::SetNewMeshAndSolution(DataState new_state,
                                      VisualizationScene* vs)
{
   if (new_state.mesh->SpaceDimension() == mesh->SpaceDimension() &&
       new_state.grid_f->VectorDim() == grid_f->VectorDim())
   {
      ResetMeshAndSolution(new_state, vs);

      // do not move 'sol' vector as it is updated directly
      internal = std::move(new_state.internal);

      return true;
   }
   else
   {
      return false;
   }
}

void DataState::ResetMeshAndSolution(DataState &ss, VisualizationScene* vs)
{
   if (ss.mesh->SpaceDimension() == 2)
   {
      if (ss.grid_f->VectorDim() == 1)
      {
         auto *vss = dynamic_cast<VisualizationSceneSolution *>(vs);
         // use the local vector as pointer is invalid after the move
         ss.grid_f->GetNodalValues(sol);
         // update the offsets before the mesh and solution
         vss->SetDataOffsets(ss.offsets.get());
         vss->NewMeshAndSolution(ss.mesh.get(), ss.mesh_quad.get(), &sol,
                                 ss.grid_f.get());
      }
      else
      {
         auto *vsv = dynamic_cast<VisualizationSceneVector *>(vs);
         vsv->SetDataOffsets(ss.offsets.get());
         vsv->NewMeshAndSolution(*ss.grid_f, ss.mesh_quad.get());
      }
   }
   else
   {
      if (ss.grid_f->VectorDim() == 1)
      {
         auto *vss = dynamic_cast<VisualizationSceneSolution3d *>(vs);
         // use the local vector as pointer is invalid after the move
         ss.grid_f->GetNodalValues(sol);
         vss->NewMeshAndSolution(ss.mesh.get(), ss.mesh_quad.get(), &sol,
                                 ss.grid_f.get());
      }
      else
      {
         ss.ProjectVectorFEGridFunction();

         auto *vss = dynamic_cast<VisualizationSceneVector3d *>(vs);
         vss->NewMeshAndSolution(ss.mesh.get(), ss.mesh_quad.get(), ss.grid_f.get());
      }
   }
}

void DataState::_ComputeDofsOffsets(std::vector<mfem::GridFunction*> &gf_array)
{
   const int nprocs = static_cast<int>(gf_array.size());
   MFEM_VERIFY(!gf_array.empty(), "No grid functions provided for offsets");

   // only 2D meshes are supported for dofs offsets computation
   if (gf_array[0]->FESpace()->GetMesh()->Dimension() != 2) { return; }

   internal.offsets = std::make_unique<DataState::Offsets>(nprocs);

   DenseMatrix pointmat;
   Array<int> dofs, vertices;
   for (int i = 0, g_e = 0; i < nprocs; i++)
   {
      const GridFunction *gf = gf_array[i];
      const FiniteElementSpace *l_fes = gf->FESpace();
      Mesh *l_mesh = l_fes->GetMesh();
      // store the dofs numbers as they are fespace dependent
      auto &offset = (*offsets)[i];
      for (int l_e = 0; l_e < l_mesh->GetNE(); l_e++, g_e++)
      {
#ifdef GLVIS_DEBUG
         // Store elements centers
         l_mesh->GetPointMatrix(l_e, pointmat);
         const int nv = pointmat.Width();
         double xs = 0.0, ys = 0.0;
         for (int j = 0; j < nv; j++)
         {
            xs += pointmat(0,j), ys += pointmat(1,j);
         }
         xs /= nv, ys /= nv;
         offset.exy_map[ {l_e, i} ] = {xs, ys};
#endif // end GLVIS_DEBUG
         l_fes->GetElementDofs(l_e, dofs);
         for (int k = 0; k < dofs.Size(); k++)
         {
            offset[ {g_e, k} ] = dofs[k];
         }
      }
      if (i + 1 == nprocs) { continue; }
      auto &next = (*offsets)[i+1];
      // for NE, NV and NEdges, we accumulate the values
      next.nelems = offset.nelems + l_mesh->GetNE();
      next.nedges = offset.nedges + l_mesh->GetNEdges();
      next.nverts = offset.nverts + l_mesh->GetNV();
   }
}

void ::VectorExtrudeCoefficient::Eval(Vector &v, ElementTransformation &T,
                                      const IntegrationPoint &ip)
{
   ElementTransformation *T_in =
      mesh_in->GetElementTransformation(T.ElementNo / n);
   T_in->SetIntPoint(&ip);
   sol_in.Eval(v, *T_in, ip);
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
