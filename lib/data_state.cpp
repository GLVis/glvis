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
#include "visual.hpp"

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
   cmplx_sol = ss.cmplx_sol;
   quad_sol = ss.quad_sol;
   keys = std::move(ss.keys);

   fix_elem_orient = ss.fix_elem_orient;
   save_coloring = ss.save_coloring;
   keep_attr = ss.keep_attr;
   cmplx_phase = ss.cmplx_phase;

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
   if (cgrid_f && cgrid_f->FESpace()->GetMesh() != mesh.get()) { SetCmplxGridFunction(NULL); }
   if (quad_f && quad_f->GetSpace()->GetMesh() != mesh.get()) { SetQuadFunction(NULL); }
}

void DataState::SetScalarData(Vector new_sol)
{
   internal.grid_f.reset();
   internal.quad_f.reset();

   if (!sol)
   {
      internal.sol.reset(new Vector(std::move(new_sol)));
   }
   else
   {
      *sol = std::move(new_sol);
   }
   type = DataState::FieldType::SCALAR;
}

void DataState::SetNormals(Vector new_normals)
{
   if (!normals)
   {
      internal.sol.reset(new Vector(std::move(new_normals)));
   }
   else
   {
      *sol = std::move(new_normals);
   }
}

void DataState::SetVectorData(Vector new_solx, Vector new_soly)
{
   MFEM_VERIFY(mesh->SpaceDimension() == 2, "Incompatible space dimension");

   internal.grid_f.reset();
   internal.quad_f.reset();

   if (!solx || !soly)
   {
      internal.solx.reset(new Vector(std::move(new_solx)));
      internal.soly.reset(new Vector(std::move(new_soly)));
   }
   else
   {
      *solx = std::move(new_solx);
      *soly = std::move(new_soly);
   }
   type = DataState::FieldType::VECTOR;
}

void DataState::SetVectorData(Vector new_solx, Vector new_soly,
                              Vector new_solz)
{
   MFEM_VERIFY(mesh->SpaceDimension() == 3, "Incompatible space dimension");

   internal.grid_f.reset();
   internal.quad_f.reset();

   if (!solx || !soly || !solz)
   {
      internal.solx.reset(new Vector(std::move(new_solx)));
      internal.soly.reset(new Vector(std::move(new_soly)));
      internal.solz.reset(new Vector(std::move(new_solz)));
   }
   else
   {
      *solx = std::move(new_solx);
      *soly = std::move(new_soly);
      *solz = std::move(new_solz);
   }
   type = DataState::FieldType::VECTOR;
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
   internal.cgrid_f.reset();
   cmplx_sol = ComplexSolution::NONE;
   internal.quad_f.reset();
   quad_sol = QuadSolution::NONE;
   SetGridFunctionSolution(component);
}

void DataState::SetGridFunction(std::vector<GridFunction*> &gf_array,
                                int num_pieces, int component)
{
   SetGridFunction(new GridFunction(mesh.get(), gf_array.data(), num_pieces),
                   component);
   if (!keep_attr) { ComputeDofsOffsets(gf_array); }
}

void DataState::SetCmplxGridFunction(ComplexGridFunction *gf,
                                     int component)
{
   if (cgrid_f.get() != gf)
   {
      internal.grid_f.reset();
      cmplx_sol = ComplexSolution::NONE;
   }
   internal.cgrid_f.reset(gf);
   internal.quad_f.reset();
   quad_sol = QuadSolution::NONE;
   SetComplexFunctionSolution(component);
}

void DataState::SetCmplxGridFunction(
   std::unique_ptr<ComplexGridFunction> &&pgf, int component)
{
   if (cgrid_f.get() != pgf.get())
   {
      internal.grid_f.reset();
      cmplx_sol = ComplexSolution::NONE;
   }
   internal.cgrid_f = std::move(pgf);
   internal.quad_f.reset();
   quad_sol = QuadSolution::NONE;
   SetComplexFunctionSolution(component);
}

void DataState::SetCmplxGridFunction(
   const std::vector<ComplexGridFunction *> &cgf_array, int component)
{
   const int np = cgf_array.size();
   std::vector<GridFunction *> r_array(np), i_array(np);
   for (int p = 0; p < np; p++)
   {
      r_array[p] = &(cgf_array[p]->real());
      i_array[p] = &(cgf_array[p]->imag());
   }
   GridFunction *rgf = new GridFunction(mesh.get(), r_array.data(), np);
   GridFunction *igf = new GridFunction(mesh.get(), i_array.data(), np);
   ComplexGridFunction *cgf = new ComplexGridFunction(rgf->FESpace());
   // transfer ownership of the FES
   cgf->MakeOwner(rgf->OwnFEC());
   rgf->MakeOwner(NULL);
   cgf->real() = *rgf;
   cgf->imag() = *igf;
   delete rgf;
   delete igf;
   SetCmplxGridFunction(cgf, component);
}

void DataState::SetQuadFunction(QuadratureFunction *qf, int component)
{
   if (quad_f.get() != qf)
   {
      internal.grid_f.reset();
      quad_sol = QuadSolution::NONE;
   }
   internal.quad_f.reset(qf);
   internal.cgrid_f.reset();
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
   internal.cgrid_f.reset();
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
   else if (sol)
   {
      Vector sol2d(mesh2d->GetNV());
      for (int i = 0; i < mesh->GetNV(); i++)
      {
         sol2d(2*i+0) = sol2d(2*i+1) = (*sol)(i);
      }
      *sol = std::move(sol2d);
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
      internal.sol.reset(new Vector(mesh->GetNV()));
      *sol = 0.0;
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

   type = (grid_f->VectorDim() == 1) ? FieldType::SCALAR : FieldType::VECTOR;
}

void DataState::SetComplexFunctionSolution(int gf_component)
{
   if (!cgrid_f)
   {
      type = (mesh)?(FieldType::MESH):(FieldType::UNKNOWN);
      return;
   }

   if (gf_component != -1)
   {
      if (gf_component < 0 || gf_component >= cgrid_f->FESpace()->GetVDim())
      {
         cerr << "Invalid component " << gf_component << '.' << endl;
         exit(1);
      }
      FiniteElementSpace *ofes = cgrid_f->FESpace();
      FiniteElementCollection *fec =
         FiniteElementCollection::New(ofes->FEColl()->Name());
      FiniteElementSpace *fes = new FiniteElementSpace(mesh.get(), fec);
      ComplexGridFunction *new_gf = new ComplexGridFunction(fes);
      new_gf->MakeOwner(fec);
      for (int i = 0; i < new_gf->real().Size(); i++)
      {
         (new_gf->real())(i) = (cgrid_f->real())(ofes->DofToVDof(i, gf_component));
         (new_gf->imag())(i) = (cgrid_f->imag())(ofes->DofToVDof(i, gf_component));
      }
      SetCmplxGridFunction(new_gf);
      return;
   }

   if (cgrid_f->VectorDim() == 1)
   {
      type = FieldType::SCALAR;
   }
   else
   {
      type = FieldType::VECTOR;
   }

   SetComplexSolution();
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

void DataState::SetComplexSolution(ComplexSolution cmplx_type, bool print)
{
   double (*cmplx_2_scalar)(double, double);
   const char *str_cmplx_2_scalar;

   switch (cmplx_type)
   {
      case ComplexSolution::Magnitude:
         cmplx_2_scalar = [](double r, double i) { return hypot(r,i); };
         str_cmplx_2_scalar = "magnitude";
         break;
      case ComplexSolution::Phase:
         cmplx_2_scalar = [](double r, double i) { return atan2(r,i); };
         str_cmplx_2_scalar = "phase";
         break;
      case ComplexSolution::Real:
         cmplx_2_scalar = [](double r, double i) { return r; };
         str_cmplx_2_scalar = "real part";
         break;
      case ComplexSolution::Imag:
         cmplx_2_scalar = [](double r, double i) { return i; };
         str_cmplx_2_scalar = "imaginary part";
         break;
      default:
         cout << "Unknown complex data representation" << endl;
         return;
   }

   if (print)
   {
      cout << "Representing complex function by: " << str_cmplx_2_scalar;
      if (cmplx_type != ComplexSolution::Magnitude) { cout << " (+ animated harmonic phase)"; }
      cout << endl;
   }
   GridFunction *gf = new GridFunction(cgrid_f->FESpace());
   const FiniteElementSpace *fes = cgrid_f->FESpace();
   const Mesh *msh = fes->GetMesh();
   const int dim = msh->Dimension();

   const double cos_ph = cos(2. * M_PI * cmplx_phase);
   const double sin_ph = sin(2. * M_PI * cmplx_phase);
   auto rot_ph = [cos_ph, sin_ph](double &r, double &i)
   {
      double rph = +r * cos_ph - i * sin_ph;
      double iph = +r * sin_ph + i * cos_ph;
      r = rph, i = iph;
   };

   cmplx_mag_max = 0.;

   if (cgrid_f->FESpace()->FEColl()->GetMapType(dim) == FiniteElement::VALUE)
   {
      for (int i = 0; i < gf->Size(); i++)
      {
         double rval = cgrid_f->real()(i);
         double ival = cgrid_f->imag()(i);
         double mag = hypot(rval, ival);
         cmplx_mag_max = max(cmplx_mag_max, mag);
         rot_ph(rval, ival);
         (*gf)(i) = cmplx_2_scalar(rval, ival);
      }
   }
   else
   {
      const int vdim = fes->GetVDim();
      const int sdim = fes->GetVectorDim();
      Array<int> vdofs;
      ElementTransformation *Tr;
      Vector r_data, i_data, z_data;
      DenseMatrix r_vals, i_vals, z_vals;
      Vector r_vec, i_vec, z_vec;
      DenseMatrix vshape;
      Vector shape;
      for (int z = 0; z < msh->GetNE(); z++)
      {
         fes->GetElementVDofs(z, vdofs);
         cgrid_f->real().GetSubVector(vdofs, r_data);
         cgrid_f->imag().GetSubVector(vdofs, i_data);
         z_data.SetSize(vdofs.Size());

         Tr = fes->GetElementTransformation(z);
         const FiniteElement *fe = fes->GetFE(z);
         const IntegrationRule &ir = fe->GetNodes();
         const int nnp = ir.GetNPoints();
         r_vals.Reset(r_data.GetData(), nnp, vdim);
         i_vals.Reset(i_data.GetData(), nnp, vdim);
         z_vals.Reset(z_data.GetData(), nnp, vdim);

         if (fe->GetRangeType() == FiniteElement::SCALAR)
         {
            shape.SetSize(nnp);
         }
         else
         {
            vshape.SetSize(nnp, sdim);
         }
         for (int n = 0; n < nnp; n++)
         {
            const IntegrationPoint &ip = ir.IntPoint(n);
            Tr->SetIntPoint(&ip);
            double w;
            if (fe->GetRangeType() == FiniteElement::SCALAR)
            {
               fe->CalcPhysShape(*Tr, shape);
               w = shape(n);
            }
            else
            {
               fe->CalcPhysVShape(*Tr, vshape);
               Vector vec(sdim);
               vshape.GetRow(n, vec);
               w = vec.Norml2();
            }
            for (int d = 0; d < vdim; d++)
            {
               double rval = r_vals(n,d) * w;
               double ival = i_vals(n,d) * w;
               double mag = hypot(rval, ival);
               cmplx_mag_max = max(cmplx_mag_max, mag);
               rot_ph(rval, ival);
               const double zval = cmplx_2_scalar(rval, ival);
               z_vals(n,d) = (w != 0.)?(zval / w):(0.);
            }
         }

         gf->SetSubVector(vdofs, z_data);
      }
   }
   internal.grid_f.reset(gf);

   cmplx_sol = cmplx_type;
}

void DataState::SwitchComplexSolution(ComplexSolution cmplx_type, bool print)
{
   SetComplexSolution(cmplx_type, print);
   ExtrudeMeshAndSolution();
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

void DataState::SwitchQuadSolution(QuadSolution quad_type)
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

std::unique_ptr<ComplexGridFunction>
DataState::ProjectVectorFEGridFunction(std::unique_ptr<ComplexGridFunction> gf)
{
   if ((gf->VectorDim() == 3) && (gf->FESpace()->GetVDim() == 1))
   {
      int p = gf->FESpace()->GetOrder(0);
      cout << "Switching to order " << p
           << " discontinuous complex vector grid function..." << endl;
      int dim = gf->FESpace()->GetMesh()->Dimension();
      FiniteElementCollection *d_fec =
         (p == 1 && dim == 3) ?
         (FiniteElementCollection*)new LinearDiscont3DFECollection :
         (FiniteElementCollection*)new L2_FECollection(p, dim, 1);
      FiniteElementSpace *d_fespace =
         new FiniteElementSpace(gf->FESpace()->GetMesh(), d_fec, 3);
      ComplexGridFunction *d_gf = new ComplexGridFunction(d_fespace);
      d_gf->MakeOwner(d_fec);
      gf->real().ProjectVectorFieldOn(d_gf->real());
      gf->imag().ProjectVectorFieldOn(d_gf->imag());
      gf.reset(d_gf);
   }
   return gf;
}

void DataState::ProjectVectorFEGridFunction()
{
   if (cgrid_f)
   {
      internal.cgrid_f = ProjectVectorFEGridFunction(std::move(internal.cgrid_f));
      SwitchComplexSolution(GetComplexSolution(), false);
   }
   else
   {
      internal.grid_f = ProjectVectorFEGridFunction(std::move(internal.grid_f));
   }
}

void DataState::ComputeDofsOffsets(std::vector<GridFunction*> &gf_array)
{
   const int nprocs = static_cast<int>(gf_array.size());
   MFEM_VERIFY(!gf_array.empty(), "No grid functions provided for offsets");

   // only 2D meshes are supported for dofs offsets computation
   if (gf_array[0]->FESpace()->GetMesh()->Dimension() != 2) { return; }

   internal.offsets = std::make_unique<DataState::Offsets>(nprocs);

   DenseMatrix pointmat;
   Array<int> dofs, vertices;
   for (int rank = 0, g_e = 0; rank < nprocs; rank++)
   {
      const GridFunction *gf = gf_array[rank];
      const FiniteElementSpace *l_fes = gf->FESpace();
      Mesh *l_mesh = l_fes->GetMesh();
      // store the dofs numbers as they are fespace dependent
      auto &offset = (*offsets)[rank];
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
         offset.exy_map[ {g_e, rank} ] = {xs, ys};
#endif // end GLVIS_DEBUG
         l_fes->GetElementDofs(l_e, dofs);
         l_fes->AdjustVDofs(dofs);
         for (int k = 0; k < dofs.Size(); k++)
         {
            offset[ {g_e, k} ] = dofs[k];
         }
      }
      if (rank + 1 == nprocs) { continue; }
      auto &next = (*offsets)[rank+1];
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
