// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.


// GLVis - an OpenGL visualization server based on the MFEM library


#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <ctime>
#include <csignal>

#include <X11/keysym.h>
#include <unistd.h>

#include "mfem.hpp"
#include "lib/visual.hpp"

using namespace std;
using namespace mfem;

const char *string_none    = "(none)";
const char *string_default = "(default)";

// Global variables for command line arguments
const char *mesh_file       = string_none;
const char *sol_file        = string_none;
const char *vec_sol_file    = string_none;
const char *gfunc_file      = string_none;
const char *arg_keys        = string_none;
int         np              = 0;
int         pad_digits      = 6;
int         gf_component    = -1;
bool        fix_elem_orient = false;
bool        save_coloring   = false;
bool        keep_attr       = false;
int         window_x        = 0; // not a command line option
int         window_y        = 0; // not a command line option
int         window_w        = 400;
int         window_h        = 350;
const char *window_title    = string_default;
const char *c_plot_caption  = string_none;
string      plot_caption;
string      extra_caption;

// Global variables
int input = 1;
Mesh *mesh = NULL;
Vector sol, solu, solv, solw, normals;
GridFunction *grid_f = NULL;
int is_gf = 0;
string keys;
VisualizationSceneScalarData *vs = NULL;

GeometryRefiner GLVisGeometryRefiner;

const char *window_titles[] = { "GLVis [scalar data]",
                                "GLVis [vector data]", "GLVis [mesh]"
                              };
istream *script = NULL;
int scr_running = 0;
int scr_level = 0;
Vector *init_nodes = NULL;
double scr_min_val, scr_max_val;

Array<istream *> input_streams;

extern char **environ;

void Set_Palette(int); // defined in palettes.hpp

void PrintSampleUsage(ostream &out);

// read the mesh and the solution from a file
void ReadSerial();

// choose grid function component and set the input flag
void SetGridFunction();

// set a (checkerboard) solution when only the mesh is given
void SetMeshSolution(Mesh *mesh, GridFunction *&grid_f, bool save_coloring);

// read the mesh and the solution from multiple files
void ReadParallel();

int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                               const char *sol_prefix, Mesh **mesh_p,
                               GridFunction **sol_p, int keep_attr);

int ReadInputStreams();

void CloseInputStreams(bool);

GridFunction *ProjectVectorFEGridFunction(GridFunction*);

void Extrude1DMeshAndSolution(Mesh **mesh_p, GridFunction **grid_f_p,
                              Vector *sol);

// Read the content of an input stream (e.g. from socket/file)
int ReadStream(istream &is, const string &data_type)
{
   // 0 - scalar data, 1 - vector data, 2 - mesh only, (-1) - unknown
   int field_type = 0;

   delete mesh; mesh = NULL;
   delete grid_f; grid_f = NULL;
   keys.clear();
   if (data_type == "fem2d_data")
   {
      mesh = new Mesh(is, 0, 0, fix_elem_orient);
      sol.Load(is, mesh->GetNV());
   }
   else if (data_type == "vfem2d_data" || data_type == "vfem2d_data_keys")
   {
      field_type = 1;
      mesh = new Mesh(is, 0, 0, fix_elem_orient);
      solu.Load(is, mesh->GetNV());
      solv.Load(is, mesh->GetNV());
      if (data_type == "vfem2d_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem3d_data")
   {
      mesh = new Mesh(is, 0, 0, fix_elem_orient);
      sol.Load(is, mesh->GetNV());
   }
   else if (data_type == "vfem3d_data" || data_type == "vfem3d_data_keys")
   {
      field_type = 1;
      mesh = new Mesh(is, 0, 0, fix_elem_orient);
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
      mesh = new Mesh(is, 1, 0, fix_elem_orient);
      grid_f = new GridFunction(mesh, is);
      if (data_type == "fem2d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "vfem2d_gf_data" || data_type == "vfem2d_gf_data_keys")
   {
      field_type = 1;
      mesh = new Mesh(is, 1, 0, fix_elem_orient);
      grid_f = new GridFunction(mesh, is);
      if (data_type == "vfem2d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "fem3d_gf_data" || data_type == "fem3d_gf_data_keys")
   {
      mesh = new Mesh(is, 1, 0, fix_elem_orient);
      grid_f = new GridFunction(mesh, is);
      if (data_type == "fem3d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "vfem3d_gf_data" || data_type == "vfem3d_gf_data_keys")
   {
      field_type = 1;
      mesh = new Mesh(is, 1, 0, fix_elem_orient);
      grid_f = new GridFunction(mesh, is);
      if (data_type == "vfem3d_gf_data_keys")
      {
         is >> keys;
      }
   }
   else if (data_type == "solution")
   {
      mesh = new Mesh(is, 1, 0, fix_elem_orient);
      grid_f = new GridFunction(mesh, is);
      field_type = (grid_f->VectorDim() == 1) ? 0 : 1;
   }
   else if (data_type == "mesh")
   {
      mesh = new Mesh(is, 1, 0, fix_elem_orient);
      SetMeshSolution(mesh, grid_f, save_coloring);
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

      mesh = new Mesh(2, tot_num_vert, tot_num_elem, 0);
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
      if (grid_f)
      {
         Extrude1DMeshAndSolution(&mesh, &grid_f, NULL);
      }
      else
      {
         Extrude1DMeshAndSolution(&mesh, NULL, &sol);
      }
   }

   return field_type;
}

int InitVis(int t)
{
   const char *win_title =
      (window_title == string_default) ? window_titles[t] : window_title;

   return InitVisualization(win_title, window_x, window_y,
                            window_w, window_h);
}

// Visualize the data in the global variables mesh, sol/grid_f, etc
void StartVisualization(int field_type)
{
   if (field_type < 0 || field_type > 2)
   {
      return;
   }

   if (InitVis(field_type))
   {
      cerr << "Initializing the visualization failed." << endl;
      return;
   }

   communication_thread *comm_thread = NULL;

   if (input_streams.Size() > 0)
   {
      auxModKeyFunc(XK_space, ThreadsPauseFunc);
      glvis_command = new GLVisCommand(&vs, &mesh, &grid_f, &sol, &keep_attr,
                                       &fix_elem_orient);
      comm_thread = new communication_thread(input_streams);
   }

   double mesh_range = -1.0;
   if (field_type == 0 || field_type == 2)
   {
      if (grid_f)
      {
         grid_f->GetNodalValues(sol);
      }
      if (mesh->SpaceDimension() == 2)
      {
         VisualizationSceneSolution *vss;
         if (field_type == 2)
         {
            Set_Palette(4);
         }
         if (normals.Size() > 0)
         {
            vs = vss = new VisualizationSceneSolution(*mesh, sol, &normals);
         }
         else
         {
            vs = vss = new VisualizationSceneSolution(*mesh, sol);
         }
         if (grid_f)
         {
            vss->SetGridFunction(*grid_f);
         }
         if (field_type == 2)
         {
            vs->OrthogonalProjection = 1;
            vs->light = 0;
            vs->Zoom(1.8);
         }
      }
      else if (mesh->SpaceDimension() == 3)
      {
         VisualizationSceneSolution3d *vss;
         vs = vss = new VisualizationSceneSolution3d(*mesh, sol);
         if (grid_f)
         {
            vss->SetGridFunction(grid_f);
         }
         if (field_type == 2)
         {
            if (mesh->Dimension() == 3)
            {
               Set_Palette(4);
               // Set_Palette(11);
               // Set_Material_And_Light(4,3);
            }
            else
            {
               Set_Palette(4);
            }
            vss->ToggleDrawAxes();
            vss->ToggleDrawMesh();
         }
      }
      if (field_type == 2)
      {
         if (grid_f)
         {
            mesh_range = grid_f->Max() + 1.0;
         }
         else
         {
            mesh_range = sol.Max() + 1.0;
         }
      }
   }
   else if (field_type == 1)
   {
      if (mesh->SpaceDimension() == 2)
      {
         if (grid_f)
         {
            vs = new VisualizationSceneVector(*grid_f);
         }
         else
         {
            vs = new VisualizationSceneVector(*mesh, solu, solv);
         }
      }
      else if (mesh->SpaceDimension() == 3)
      {
         if (grid_f)
         {
            grid_f = ProjectVectorFEGridFunction(grid_f);
            vs = new VisualizationSceneVector3d(*grid_f);
         }
         else
         {
            vs = new VisualizationSceneVector3d(*mesh, solu, solv, solw);
         }
      }
   }

   if (vs)
   {
      // increase the refinement factors if visualizing a GridFunction
      if (grid_f)
      {
         vs->AutoRefine();
         vs->SetShading(2, true);
      }
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
         vs->SetAutoscale(0);
      }
      if (mesh->SpaceDimension() == 2 && field_type == 2)
      {
         SetVisualizationScene(vs, 2, keys.c_str());
      }
      else
      {
         SetVisualizationScene(vs, 3, keys.c_str());
      }
   }

   KillVisualization(); // deletes vs
   vs = NULL;
   if (input_streams.Size() > 0)
   {
      glvis_command->Terminate();
      delete comm_thread;
      delete glvis_command;
      glvis_command = NULL;
   }
   delete grid_f; grid_f = NULL;
   delete mesh; mesh = NULL;
   cout << "GLVis window closed." << endl;
}

int ScriptReadSolution(istream &scr, Mesh **mp, GridFunction **sp)
{
   string word;

   cout << "Script: solution: " << flush;
   // read the mesh
   scr >> ws >> word; // mesh filename (can't contain spaces)
   {
      named_ifgzstream imesh(word.c_str());
      if (!imesh)
      {
         cout << "Can not open mesh file: " << word << endl;
         return 1;
      }
      *mp = new Mesh(imesh, 1, 0, fix_elem_orient);
      cout << "mesh: " << word << "; " << flush;
   }

   // read the solution (GridFunction)
   scr >> ws >> word;
   {
      ifgzstream isol(word.c_str());
      if (!isol)
      {
         cout << "Can not open solution file: " << word << endl;
         delete *mp; *mp = NULL;
         return 2;
      }
      *sp = new GridFunction(*mp, isol);
      cout << "solution: " << word << endl;
   }

   Extrude1DMeshAndSolution(mp, sp, NULL);

   return 0;
}

int ScriptReadParSolution(istream &scr, Mesh **mp, GridFunction **sp)
{
   int np, keep_attr, err;
   string mesh_prefix, sol_prefix;

   cout << "Script: psolution: " << flush;
   // read number of processors
   scr >> np;
   cout << "# processors: " << np << "; " << flush;
   // read the mesh prefix
   scr >> ws >> mesh_prefix; // mesh prefix (can't contain spaces)
   cout << "mesh prefix: " << mesh_prefix << "; " << flush;
   scr >> ws >> keep_attr;
   if (keep_attr)
   {
      cout << "(real attributes); " << flush;
   }
   else
   {
      cout << "(processor attributes); " << flush;
   }
   // read the solution prefix
   scr >> ws >> sol_prefix;
   cout << "solution prefix: " << sol_prefix << endl;

   err = ReadParMeshAndGridFunction(np, mesh_prefix.c_str(),
                                    sol_prefix.c_str(), mp, sp, keep_attr);
   if (!err)
   {
      Extrude1DMeshAndSolution(mp, sp, NULL);
   }
   return err;
}

int ScriptReadDisplMesh(istream &scr, Mesh **mp, GridFunction **sp)
{
   Mesh *m;
   string word;

   cout << "Script: mesh: " << flush;
   scr >> ws >> word;
   {
      named_ifgzstream imesh(word.c_str());
      if (!imesh)
      {
         cout << "Can not open mesh file: " << word << endl;
         return 1;
      }
      cout << word << endl;
      m = new Mesh(imesh, 1, 0, fix_elem_orient);
   }
   Extrude1DMeshAndSolution(&m, NULL, NULL);
   if (init_nodes == NULL)
   {
      init_nodes = new Vector;
      m->GetNodes(*init_nodes);
      delete m;
      *mp = NULL;
      *sp = NULL;
   }
   else
   {
      FiniteElementCollection  *vfec = NULL;
      FiniteElementSpace *vfes;
      GridFunction *g;
      vfes = (FiniteElementSpace *)m->GetNodalFESpace();
      if (vfes == NULL)
      {
         vfec = new LinearFECollection;
         vfes = new FiniteElementSpace(m, vfec, m->SpaceDimension());
      }

      g = new GridFunction(vfes);
      if (vfec)
      {
         g->MakeOwner(vfec);
      }
      m->GetNodes(*g);
      if (g->Size() == init_nodes->Size())
      {
         subtract(*init_nodes, *g, *g);
      }
      else
      {
         cout << "Script: incompatible meshes!" << endl;
         *g = 0.0;
      }
      *mp = m;
      *sp = g;
   }

   return 0;
}

void ExecuteScriptCommand()
{
   if (!script)
   {
      cout << "No script stream defined! (Bug?)" << endl;
      return;
   }

   istream &scr = *script;
   string word;
   int done_one_command = 0;
   while (!done_one_command)
   {
      scr >> ws;
      if (!scr.good())
      {
         cout << "End of script." << endl;
         scr_level = 0;
         return;
      }
      if (scr.peek() == '#')
      {
         getline(scr, word);
         continue;
      }
      scr >> word;
      if (word == "{")
      {
         scr_level++;
      }
      else if (word == "}")
      {
         scr_level--;
         if (scr_level < 0)
         {
            scr_level = 0;
         }
      }
      else if (word == "solution" || word == "mesh" || word == "psolution")
      {
         Mesh *new_m = NULL;
         GridFunction *new_g;

         if (word == "solution")
         {
            if (ScriptReadSolution(scr, &new_m, &new_g))
            {
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "mesh")
         {
            if (ScriptReadDisplMesh(scr, &new_m, &new_g))
            {
               done_one_command = 1;
               continue;
            }
            if (new_m == NULL)
            {
               cout << "Script: unexpected 'mesh' command!" << endl;
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "psolution")
         {
            if (ScriptReadParSolution(scr, &new_m, &new_g))
            {
               done_one_command = 1;
               continue;
            }
         }

         if (new_m->SpaceDimension() == mesh->SpaceDimension() &&
             new_g->VectorDim() == grid_f->VectorDim())
         {
            if (new_m->SpaceDimension() == 2)
            {
               if (new_g->VectorDim() == 1)
               {
                  VisualizationSceneSolution *vss =
                     dynamic_cast<VisualizationSceneSolution *>(vs);
                  new_g->GetNodalValues(sol);
                  vss->NewMeshAndSolution(new_m, &sol, new_g);
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
                  vss->NewMeshAndSolution(new_m, &sol, new_g);
               }
               else
               {
                  new_g = ProjectVectorFEGridFunction(new_g);
                  VisualizationSceneVector3d *vss =
                     dynamic_cast<VisualizationSceneVector3d *>(vs);
                  vss->NewMeshAndSolution(new_m, new_g);
               }
            }
            delete grid_f; grid_f = new_g;
            delete mesh; mesh = new_m;

            vs->Draw();
         }
         else
         {
            cout << "Different type of mesh / solution." << endl;
            delete new_g;
            delete new_m;
         }
      }
      else if (word == "screenshot")
      {
         scr >> ws >> word;

         cout << "Script: screenshot: " << flush;

         if (Screenshot(word.c_str(), true))
         {
            cout << "Screenshot(" << word << ") failed." << endl;
            done_one_command = 1;
            continue;
         }
         cout << "-> " << word << endl;

         if (scr_min_val > vs->GetMinV())
         {
            scr_min_val = vs->GetMinV();
         }
         if (scr_max_val < vs->GetMaxV())
         {
            scr_max_val = vs->GetMaxV();
         }
      }
      else if (word == "viewcenter")
      {
         scr >> vs->ViewCenterX >> vs->ViewCenterY;
         cout << "Script: viewcenter: "
              << vs->ViewCenterX << ' ' << vs->ViewCenterY << endl;
         MyExpose();
      }
      else if (word ==  "perspective")
      {
         scr >> ws >> word;
         cout << "Script: perspective: " << word;
         if (word == "off")
         {
            vs->OrthogonalProjection = 1;
         }
         else if (word == "on")
         {
            vs->OrthogonalProjection = 0;
         }
         else
         {
            cout << '?';
         }
         cout << endl;
         MyExpose();
      }
      else if (word ==  "light")
      {
         scr >> ws >> word;
         cout << "Script: light: " << word;
         if (word == "off")
         {
            vs->light = 0;
         }
         else if (word == "on")
         {
            vs->light = 1;
         }
         else
         {
            cout << '?';
         }
         cout << endl;
         MyExpose();
      }
      else if (word == "view")
      {
         double theta, phi;
         scr >> theta >> phi;
         cout << "Script: view: " << theta << ' ' << phi << endl;
         vs->SetView(theta, phi);
         MyExpose();
      }
      else if (word == "zoom")
      {
         double factor;
         scr >> factor;
         cout << "Script: zoom: " << factor << endl;
         vs->Zoom(factor);
         MyExpose();
      }
      else if (word == "shading")
      {
         scr >> ws >> word;
         cout << "Script: shading: " << flush;
         int s = -1;
         if (word == "flat")
         {
            s = 0;
         }
         else if (word == "smooth")
         {
            s = 1;
         }
         else if (word == "cool")
         {
            s = 2;
         }
         if (s != -1)
         {
            vs->SetShading(s, false);
            cout << word << endl;
            MyExpose();
         }
         else
         {
            cout << word << " ?" << endl;
         }
      }
      else if (word == "subdivisions")
      {
         int t, b;
         scr >> t >> b;
         cout << "Script: subdivisions: " << flush;
         vs->SetRefineFactors(t, b);
         cout << t << ' ' << b << endl;
         MyExpose();
      }
      else if (word == "valuerange")
      {
         double min, max;
         scr >> min >> max;
         cout << "Script: valuerange: " << flush;
         vs->SetValueRange(min, max);
         cout << min << ' ' << max << endl;
         MyExpose();
      }
      else if (word == "autoscale")
      {
         scr >> ws >> word;
         cout << "Script: autoscale: " << word;
         if (word == "off")
         {
            vs->SetAutoscale(0);
         }
         else if (word == "on")
         {
            vs->SetAutoscale(1);
         }
         else if (word == "value")
         {
            vs->SetAutoscale(2);
         }
         else if (word == "mesh")
         {
            vs->SetAutoscale(3);
         }
         else
         {
            cout << '?';
         }
         cout << endl;
      }
      else if (word == "window")
      {
         scr >> window_x >> window_y >> window_w >> window_h;
         cout << "Script: window: " << window_x << ' ' << window_y
              << ' ' << window_w << ' ' << window_h << endl;
         MoveResizeWindow(window_x, window_y, window_w, window_h);
         MyExpose();
      }
      else if (word == "keys")
      {
         scr >> keys;
         cout << "Script: keys: '" << keys << "'" << endl;
         // SendKeySequence(keys.c_str());
         CallKeySequence(keys.c_str());
         MyExpose();
      }
      else if (word == "palette")
      {
         int pal;
         scr >> pal;
         cout << "Script: palette: " << pal << endl;
         Set_Palette(pal-1);
         if (!GetUseTexture())
         {
            vs->EventUpdateColors();
         }
         MyExpose();
      }
      else if (word == "toggle_attributes")
      {
         Array<int> attr_list;
         cout << "Script: toggle_attributes:";
         for (scr >> ws; scr.peek() != ';'; scr >> ws)
         {
            attr_list.Append(0);
            scr >> attr_list.Last();
            if (attr_list.Size() <= 256)
            {
               cout << ' ' << attr_list.Last();
            }
            else if (attr_list.Size() == 257)
            {
               cout << " ... " << flush;
            }
         }
         scr.get(); // read the end symbol: ';'
         cout << endl;
         vs->ToggleAttributes(attr_list);
         MyExpose();
      }
      else if (word == "rotmat")
      {
         cout << "Script: rotmat:";
         for (int i = 0; i < 16; i++)
         {
            scr >> vs->rotmat[i];
            cout << ' ' << vs->rotmat[i];
         }
         cout << endl;
         MyExpose();
      }
      else if (word == "camera")
      {
         double cam[9];
         cout << "Script: camera:";
         for (int i = 0; i < 9; i++)
         {
            scr >> cam[i];
            cout << ' ' << cam[i];
         }
         cout << endl;
         vs->cam.Set(cam);
         MyExpose();
      }
      else if (word == "scale")
      {
         double scale;
         cout << "Script: scale:";
         scr >> scale;
         cout << ' ' << scale;
         cout << endl;
         vs->Scale(scale);
         MyExpose();
      }
      else if (word == "translate")
      {
         double x, y, z;
         cout << "Script: translate:";
         scr >> x >> y >> z;
         cout << ' ' << x << ' ' << y << ' ' << z;
         cout << endl;
         vs->Translate(x, y, z);
         MyExpose();
      }
      else if (word == "plot_caption")
      {
         char delim;
         scr >> ws >> delim;
         getline(scr, plot_caption, delim);
         vs->UpdateCaption(); // turn on or off the caption
         MyExpose();
      }
      else
      {
         cout << "Unknown command in script: " << word << endl;
      }

      done_one_command = 1;
   }
}

void ScriptControl();

void ScriptIdleFunc()
{
   ExecuteScriptCommand();
   if (scr_level == 0)
   {
      ScriptControl();
   }
}

void ScriptControl()
{
   if (scr_running)
   {
      scr_running = 0;
      RemoveIdleFunc(ScriptIdleFunc);
   }
   else
   {
      scr_running = 1;
      AddIdleFunc(ScriptIdleFunc);
   }
}

void PlayScript(istream &scr)
{
   string word;

   scr_min_val = numeric_limits<double>::infinity();
   scr_max_val = -scr_min_val;

   // read initializing commands
   while (1)
   {
      scr >> ws;
      if (!scr.good())
      {
         cout << "Error in script" << endl;
         return;
      }
      if (scr.peek() == '#')
      {
         getline(scr, word);
         continue;
      }
      scr >> word;
      if (word == "window")
      {
         scr >> window_x >> window_y >> window_w >> window_h;
      }
      else if (word == "solution")
      {
         if (ScriptReadSolution(scr, &mesh, &grid_f))
         {
            return;
         }

         // start the visualization
         break;
      }
      else if (word == "psolution")
      {
         if (ScriptReadParSolution(scr, &mesh, &grid_f))
         {
            return;
         }

         // start the visualization
         break;
      }
      else if (word == "mesh")
      {
         if (ScriptReadDisplMesh(scr, &mesh, &grid_f))
         {
            return;
         }
         if (mesh)
         {
            break;
         }
      }
      else
      {
         cout << "Unknown command in script: " << word << endl;
      }
   }

   scr_level = scr_running = 0;
   auxKeyFunc(XK_space, ScriptControl);
   script = &scr;
   keys.clear();

   StartVisualization((grid_f->VectorDim() == 1) ? 0 : 1);

   delete init_nodes; init_nodes = NULL;

   cout << "Script: min_val = " << scr_min_val
        << ", max_val = " << scr_max_val << endl;

   script = NULL;
}


int main (int argc, char *argv[])
{
   // variables for command line arguments
   bool        multi_session = true;   // not added as option
   bool        mac           = false;
   const char *stream_file   = string_none;
   const char *script_file   = string_none;
   const char *font_name     = string_default;
   int         portnum       = 19916;
   bool        secure        = socketstream::secure_default;
   int         multisample   = GetMultisample();
   double      line_width    = Get_LineWidth();
   double      ms_line_width = Get_MS_LineWidth();
   int         geom_ref_type = 0;

   OptionsParser args(argc, argv);

   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to visualize.");
   args.AddOption(&gfunc_file, "-g", "--grid-function",
                  "Solution (GridFunction) file to visualize.");
   args.AddOption(&gf_component, "-gc", "--grid-function-component",
                  "Select a grid function component, [0-<num-comp>) or"
                  " -1 for all.");
   args.AddOption(&sol_file, "-s", "--scalar-solution",
                  "Scalar solution (vertex values) file to visualize.");
   args.AddOption(&vec_sol_file, "-v", "--vector-solution",
                  "Vector solution (vertex values) file to visualize.");
   args.AddOption(&np, "-np", "--num-proc",
                  "Load mesh/solution from multiple processors.");
   args.AddOption(&pad_digits, "-d", "--pad-digits",
                  "Number of digits used for processor ranks in file names.");
   args.AddOption(&script_file, "-run", "--run-script",
                  "Run a GLVis script file.");
   args.AddOption(&arg_keys, "-k", "--keys",
                  "Execute key shortcut commands in the GLVis window.");
   args.AddOption(&fix_elem_orient, "-fo", "--fix-orientations",
                  "-no-fo", "--dont-fix-orientations",
                  "Attempt to fix the orientations of inverted elements.");
   args.AddOption(&keep_attr, "-a", "--real-attributes",
                  "-ap", "--processor-attributes",
                  "When opening a parallel mesh, use the real mesh attributes"
                  " or replace them with the processor rank.");
   args.AddOption(&geom_ref_type, "-grt", "--geometry-refiner-type",
                  "Set of points to use when refining geometry:"
                  " 0 = uniform, 1 = Gauss-Lobatto.");
   args.AddOption(&save_coloring, "-sc", "--save-coloring",
                  "-no-sc", "--dont-save-coloring",
                  "Save the mesh coloring generated when opening only a mesh.");
   args.AddOption(&portnum, "-p", "--listen-port",
                  "Specify the port number on which to accept connections.");
   args.AddOption(&secure, "-sec", "--secure-sockets",
                  "-no-sec", "--standard-sockets",
                  "Enable or disable GnuTLS secure sockets.");
   args.AddOption(&mac, "-mac", "--save-stream",
                  "-no-mac", "--dont-save-stream",
                  "In server mode, save incoming data to a file before"
                  " visualization.");
   args.AddOption(&stream_file, "-saved", "--saved-stream",
                  "Load a GLVis stream saved to a file.");
   args.AddOption(&window_w, "-ww", "--window-width",
                  "Set the window width.");
   args.AddOption(&window_h, "-wh", "--window-height",
                  "Set the window height.");
   args.AddOption(&window_title, "-wt", "--window-title",
                  "Set the window title.");
   args.AddOption(&c_plot_caption, "-c", "--plot-caption",
                  "Set the plot caption (visible when colorbar is visible).");
   args.AddOption(&font_name, "-fn", "--font",
                  "Set the font: <font-name>[-<font-size>].");
   args.AddOption(&multisample, "-ms", "--multisample",
                  "Set the multisampling mode (toggled with the 'A' key).");
   args.AddOption(&line_width, "-lw", "--line-width",
                  "Set the line width (multisampling off).");
   args.AddOption(&ms_line_width, "-mslw", "--multisample-line-width",
                  "Set the line width (multisampling on).");

   cout << endl
        << "       _/_/_/  _/      _/      _/  _/"          << endl
        << "    _/        _/      _/      _/        _/_/_/" << endl
        << "   _/  _/_/  _/      _/      _/  _/  _/_/"      << endl
        << "  _/    _/  _/        _/  _/    _/      _/_/"   << endl
        << "   _/_/_/  _/_/_/_/    _/      _/  _/_/_/"      << endl
        << endl ;

   args.Parse();
   if (!args.Good())
   {
      if (!args.Help())
      {
         args.PrintError(cout);
         cout << endl;
      }
      PrintSampleUsage(cout);
      args.PrintHelp(cout);
      return 1;
   }

   // set options
   if (mesh_file != string_none)
   {
      input |= 2;
   }
   if (sol_file != string_none)
   {
      input |= 4;
   }
   if (vec_sol_file != string_none)
   {
      sol_file = vec_sol_file;
      input |= 8;
   }
   if (gfunc_file != string_none)
   {
      sol_file = gfunc_file;
      is_gf = 255;
   }
   if (np > 0)
   {
      input |= 256;
   }
   if (arg_keys != string_none)
   {
      keys = arg_keys;
   }
   if (font_name != string_default)
   {
      SetFont(font_name);
   }
   if (multisample != GetMultisample())
   {
      SetMultisample(multisample);
   }
   if (line_width != Get_LineWidth())
   {
      Set_LineWidth(line_width);
   }
   if (ms_line_width != Get_MS_LineWidth())
   {
      Set_MS_LineWidth(ms_line_width);
   }
   if (c_plot_caption != string_none)
   {
      plot_caption = c_plot_caption;
   }

   GLVisGeometryRefiner.SetType(geom_ref_type);

   string data_type;

   // check for saved stream file
   if (stream_file != string_none)
   {
      ifstream ifs(stream_file);
      if (!ifs)
      {
         cout << "Can not open stream file: " << stream_file << endl;
         return 1;
      }
      ifs >> data_type >> ws;
      int ft = ReadStream(ifs, data_type);
      input_streams.Append(&ifs);
      StartVisualization(ft);
      return 0;
   }

   // check for script file
   if (script_file != string_none)
   {
      ifstream scr(script_file);
      if (!scr)
      {
         cout << "Can not open script: " << script_file << endl;
         return 1;
      }
      PlayScript(scr);
      return 0;
   }

   // print help for wrong input
   if (!(input == 1 || input == 3 || input == 7 || input == 11 || input == 259 ||
         (is_gf && (input == 3 || input == 259))))
   {
      cout << "Invalid combination of mesh/solution options!\n\n";
      PrintSampleUsage(cout);
      args.PrintHelp(cout);
      return 1;
   }

#ifndef MFEM_USE_GNUTLS
   if (secure)
   {
      cout << "The secure option can only be used when MFEM is compiled with"
           " GnuTLS support." << endl;
      return 1;
   }
#endif

   int childPID, viscount = 0, nproc = 1, proc = 0;

   // server mode, read the mesh and the solution from a socket
   if (input == 1)
   {
      // get rid of zombies
      if (multi_session)
      {
         signal(SIGCHLD, SIG_IGN);
      }

#ifdef MFEM_USE_GNUTLS
      GnuTLS_global_state *state = NULL;
      // state->set_log_level(1000);
      GnuTLS_session_params *params = NULL;
      if (secure)
      {
         state = new GnuTLS_global_state;
         string home_dir(getenv("HOME"));
         string server_dir = home_dir + "/.config/glvis/server/";
         string pubkey  = server_dir + "pubring.gpg";
         string privkey = server_dir + "secring.gpg";
         string trustedkeys = server_dir + "trusted-clients.gpg";
         params = new GnuTLS_session_params(
            *state, pubkey.c_str(), privkey.c_str(),
            trustedkeys.c_str(), GNUTLS_SERVER);
         if (!params->status.good())
         {
            cout << "  public key   = " << pubkey << '\n'
                 << "  private key  = " << privkey << '\n'
                 << "  trusted keys = " << trustedkeys << endl;
            cout << "Error setting GLVis server parameters.\n"
                 "Generate your GLVis keys with:"
                 " bash glvis-keygen.sh [\"Your Name\"] [\"Your Email\"]"
                 << endl;
            delete params; delete state;
            return 3;
         }
      }
#endif

      const int backlog = 128;
      socketserver server(portnum, backlog);
      if (server.good())
      {
         cout << "Waiting for data on port " << portnum << " ..." << endl;
      }
      else
      {
         cout << "Server already running on port " << portnum << ".\n" << endl;
#ifdef MFEM_USE_GNUTLS
         delete params; delete state;
#endif
         return 2;
      }

      socketstream *isock;
#ifndef MFEM_USE_GNUTLS
      isock = new socketstream;
#else
      isock = secure ? new socketstream(*params) : new socketstream(false);
#endif
      while (1)
      {
         while (server.accept(*isock) < 0)
         {
#ifdef GLVIS_DEBUG
            cout << "GLVis: server.accept(...) failed." << endl;
#endif
         }

         *isock >> data_type >> ws;

         if (mac)
         {
            viscount++;
         }

         int par_data = 0;
         if (data_type == "parallel")
         {
            par_data = 1;
            np = 0;
            do
            {
               *isock >> nproc >> proc;
#ifdef GLVIS_DEBUG
               cout << "new connection: parallel " << nproc << ' ' << proc
                    << endl;
#endif
               if (np == 0)
               {
                  if (nproc <= 0)
                  {
                     cout << "Invalid number of processors: " << nproc << endl;
                     mfem_error();
                  }
                  input_streams.SetSize(nproc);
                  input_streams = NULL;
               }
               else
               {
                  if (nproc != input_streams.Size())
                  {
                     cout << "Unexpected number of processors: " << nproc
                          << ", expected: " << input_streams.Size() << endl;
                     mfem_error();
                  }
               }
               if (0 > proc || proc >= nproc)
               {
                  cout << "Invalid processor rank: " << proc
                       << ", number of processors: " << nproc << endl;
                  mfem_error();
               }
               if (input_streams[proc])
               {
                  cout << "Second connection attempt from processor rank: "
                       << proc << endl;
                  mfem_error();
               }

               input_streams[proc] = isock;
#ifndef MFEM_USE_GNUTLS
               isock = new socketstream;
#else
               isock = secure ? new socketstream(*params) :
                       new socketstream(false);
#endif
               np++;
               if (np == nproc)
               {
                  break;
               }
               // read next available socket stream
               while (server.accept(*isock) < 0)
               {
#ifdef GLVIS_DEBUG
                  cout << "GLVis: server.accept(...) failed." << endl;
#endif
               }
               *isock >> data_type >> ws; // "parallel"
               if (data_type != "parallel")
               {
                  cout << "Expected keyword \"parallel\", got \"" << data_type
                       << '"' << endl;
                  mfem_error();
               }
            }
            while (1);
         }

         char tmp_file[50];
         if (multi_session)
         {
            if (mac)
            {
               sprintf(tmp_file,"glvis-saved.%04d",viscount);
               ofstream ofs(tmp_file);
               if (!par_data)
               {
                  ofs << data_type << '\n';
                  ofs << isock->rdbuf();
                  isock->close();
               }
               else
               {
                  ReadInputStreams();
                  CloseInputStreams(false);
                  ofs.precision(8);
                  ofs << "solution\n";
                  mesh->Print(ofs);
                  grid_f->Save(ofs);
                  delete grid_f; grid_f = NULL;
                  delete mesh; mesh = NULL;
               }
               ofs.close();
               cout << "Data saved in " << tmp_file << endl;
            }
            childPID = fork();
         }
         else
         {
            childPID = 0;
         }

         switch (childPID)
         {
            case -1:
               cout << "The process couldn't fork. Exit." << endl;
               exit(1);

            case 0:                       // This is the child process
               server.close();
               if (mac)
               {
                  // exec ourself
                  const char *args[4] = { argv[0], "-saved", tmp_file, NULL };
                  execve(args[0], (char* const*)args, environ);
                  exit(0);
               }
               else
               {
                  if (multi_session)
                  {
                     signal(SIGINT, SIG_IGN);
                  }
                  int ft;
                  if (!par_data)
                  {
                     ft = ReadStream(*isock, data_type);
                     input_streams.Append(isock);
                  }
                  else
                  {
                     delete isock;
                     ft = ReadInputStreams();
                  }
                  StartVisualization(ft);
                  CloseInputStreams(false);
                  exit(0);
               }

            default :                     // This is the parent process
               if (!par_data)
               {
                  isock->rdbuf()->socketbuf::close();
               }
               else
               {
                  CloseInputStreams(true);
               }
         }
      }
#ifdef MFEM_USE_GNUTLS
      delete params; delete state;
#endif
   }
   else  // input != 1, non-server mode
   {
      if (input & 256)
      {
         ReadParallel();
      }
      else
      {
         ReadSerial();
      }

      int window_err;
      double mesh_range = -1.0;
      if (mesh->SpaceDimension() == 2)
      {
         if ((input & 8) == 0)
         {
            VisualizationSceneSolution *vss;
            window_err = InitVis((input & 4) ? 0 : 2);
            if (!window_err)
            {
               if ((input & 4) == 0)
               {
                  Set_Palette(4);   // Set_Palette(11);
               }
               vs = vss = new VisualizationSceneSolution(*mesh, sol);
               if (is_gf)
               {
                  vss->SetGridFunction(*grid_f);
               }
               if ((input & 4) == 0)
               {
                  vs->OrthogonalProjection = 1;
                  vs->light = 0;
                  vs->Zoom(1.8);
                  if (grid_f)
                  {
                     mesh_range = grid_f->Max() + 1.0;
                  }
                  else
                  {
                     mesh_range = sol.Max() + 1.0;
                  }
               }
            }
         }
         else
         {
            window_err = InitVis(1);
            if (!window_err)
            {
               if (is_gf)
               {
                  vs = new VisualizationSceneVector(*grid_f);
               }
               else
               {
                  vs = new VisualizationSceneVector(*mesh, solu, solv);
               }
            }
         }
      }
      else // 3D
      {
         if ((input & 8) == 0 && (input & 512) == 0)
         {
            VisualizationSceneSolution3d *vss;
            window_err = InitVis((input & 4) ? 0 : 2);
            if (!window_err)
            {
               vs = vss = new VisualizationSceneSolution3d(*mesh, sol);
               if (is_gf)
               {
                  vss->SetGridFunction(grid_f);
               }
               if ((input & 4) == 0)
               {
                  if (mesh->Dimension() == 3)
                  {
                     // Set_Palette(4);
                     Set_Palette(11);
                     Set_Material_And_Light(4,3);
                  }
                  else
                  {
                     Set_Palette(4);
                  }
                  vss->ToggleDrawAxes();
                  vss->ToggleDrawMesh();
                  if (grid_f)
                  {
                     mesh_range = grid_f->Max() + 1.0;
                  }
                  else
                  {
                     mesh_range = sol.Max() + 1.0;
                  }
               }
            }
         }
         else
         {
            window_err = InitVis(1);
            if (!window_err)
            {
               if (is_gf)
               {
                  grid_f = ProjectVectorFEGridFunction(grid_f);
                  vs = new VisualizationSceneVector3d(*grid_f);
               }
               else
               {
                  vs = new VisualizationSceneVector3d(*mesh, solu, solv, solw);
               }
            }
         }
      }
      if (!window_err)
      {
         // increase the refinement factors if visualizing a GridFunction
         if (is_gf)
         {
            vs->AutoRefine();
            vs->SetShading(2, true);
         }
         if (mesh_range > 0.0)
         {
            vs->SetValueRange(-mesh_range, mesh_range);
            vs->SetAutoscale(0);
         }
         if (mesh->SpaceDimension() == 2 && (input & 12) == 0)
         {
            SetVisualizationScene(vs, 2, keys.c_str());
         }
         else
         {
            SetVisualizationScene(vs, 3, keys.c_str());
         }
         KillVisualization(); // deletes vs
         if (is_gf) { delete grid_f; }
         delete mesh;
      }
      else
      {
         cerr << "Initializing the visualization failed." << endl;
         return 1;
      }
   }

   cout << "Thank you for using GLVis." << endl;

   return 0;
}


void PrintSampleUsage(ostream &out)
{
   out <<
       "Start a GLVis server:\n"
       "   glvis\n"
       "Visualize a mesh:\n"
       "   glvis -m <mesh_file>\n"
       "Visualize mesh and solution (grid function):\n"
       "   glvis -m <mesh_file> -g <grid_function_file> [-gc <component>]\n"
       "Visualize parallel mesh and solution (grid function):\n"
       "   glvis -np <#proc> -m <mesh_prefix> [-g <grid_function_prefix>]\n\n"
       "All Options:\n";
}


void ReadSerial()
{
   // get the mesh from a file
   {
      ifstream meshin(mesh_file);
      if (!meshin)
      {
         cerr << "Can not open mesh file " << mesh_file << ". Exit.\n";
         exit(1);
      }

      // mesh = new Mesh(meshin, is_gf?1:0, 0);
      mesh = new Mesh(mesh_file, 1, 0, fix_elem_orient);
   }

   if (is_gf || (input & 4) || (input & 8))
   {
      // get the solution from file
      ifgzstream solin(sol_file);
      if (!solin)
      {
         cerr << "Can not open solution file " << sol_file << ". Exit.\n";
         exit(1);
      }

      if (is_gf)
      {
         grid_f = new GridFunction(mesh, solin);
         SetGridFunction();
      }
      else if (input & 4)
      {
         // get rid of NetGen's info line
         char buff[128];
         solin.getline(buff,128);
         sol.Load(solin, mesh->GetNV());
      }
      else if (input & 8)
      {
         solu.Load(solin, mesh->GetNV());
         solv.Load(solin, mesh->GetNV());
         if (mesh->SpaceDimension() == 3)
         {
            solw.Load(solin, mesh->GetNV());
         }
      }
   }
   else
   {
      SetMeshSolution(mesh, grid_f, save_coloring);
   }

   Extrude1DMeshAndSolution(&mesh, &grid_f, &sol);
}


void SetGridFunction()
{
   if (gf_component != -1)
   {
      if (gf_component < 0 || gf_component >= grid_f->VectorDim())
      {
         cerr << "Invalid component " << gf_component << '.' << endl;
         exit(1);
      }
      FiniteElementSpace *ofes = grid_f->FESpace();
      FiniteElementCollection *fec =
         FiniteElementCollection::New(ofes->FEColl()->Name());
      FiniteElementSpace *fes = new FiniteElementSpace(mesh, fec);
      GridFunction *new_gf = new GridFunction(fes);
      new_gf->MakeOwner(fec);
      for (int i = 0; i < new_gf->Size(); i++)
      {
         (*new_gf)(i) = (*grid_f)(ofes->DofToVDof(i, gf_component));
      }
      delete grid_f;
      grid_f = new_gf;
   }
   if (grid_f->VectorDim() == 1)
   {
      grid_f->GetNodalValues(sol);
      input |= 4;
   }
   else
   {
      input |= 8;
   }
}


void SetMeshSolution(Mesh *mesh, GridFunction *&grid_f, bool save_coloring)
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


void ReadParallel()
{
   int err;

   if (is_gf)
   {
      err = ReadParMeshAndGridFunction(np, mesh_file, sol_file,
                                       &mesh, &grid_f, keep_attr);
      if (!err)
      {
         SetGridFunction();
      }
   }
   else
   {
      err = ReadParMeshAndGridFunction(np, mesh_file, NULL,
                                       &mesh, NULL, keep_attr);
      if (!err)
      {
         SetMeshSolution(mesh, grid_f, save_coloring);
      }
   }

   if (err)
   {
      exit(1);
   }

   Extrude1DMeshAndSolution(&mesh, &grid_f, &sol);
}

int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                               const char *sol_prefix, Mesh **mesh_p,
                               GridFunction **sol_p, int keep_attr)
{
   Array<Mesh *> mesh_array;

   mesh_array.SetSize(np);
   for (int p = 0; p < np; p++)
   {
      ostringstream fname;
      fname << mesh_prefix << '.' << setfill('0') << setw(pad_digits) << p;
      named_ifgzstream meshfile(fname.str().c_str());
      if (!meshfile)
      {
         cerr << "Can not open mesh file: " << fname.str().c_str()
              << '!' << endl;
         for (p--; p >= 0; p--)
         {
            delete mesh_array[p];
         }
         return 1;
      }
      mesh_array[p] = new Mesh(meshfile, 1, 0, fix_elem_orient);
      if (!keep_attr)
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
   }
   *mesh_p = new Mesh(mesh_array, np);

   if (sol_prefix && sol_p)
   {
      Array<GridFunction *> gf_array(np);
      for (int p = 0; p < np; p++)
      {
         ostringstream fname;
         fname << sol_prefix << '.' << setfill('0') << setw(pad_digits) << p;
         ifgzstream solfile(fname.str().c_str());
         if (!solfile)
         {
            cerr << "Can not open solution file " << fname.str().c_str()
                 << '!' << endl;
            for (p--; p >= 0; p--)
            {
               delete gf_array[p];
            }
            delete *mesh_p;
            *mesh_p = NULL;
            for (p = 0; p < np; p++)
            {
               delete mesh_array[np-1-p];
            }
            return 2;
         }
         gf_array[p] = new GridFunction(mesh_array[p], solfile);
      }
      *sol_p = new GridFunction(*mesh_p, gf_array, np);

      for (int p = 0; p < np; p++)
      {
         delete gf_array[np-1-p];
      }
   }

   for (int p = 0; p < np; p++)
   {
      delete mesh_array[np-1-p];
   }

   return 0;
}

int ReadInputStreams()
{
   int nproc = input_streams.Size();
   Array<Mesh *> mesh_array(nproc);
   Array<GridFunction *> gf_array(nproc);
   string data_type;

   int gf_count = 0;
   int field_type = 0;

   for (int p = 0; p < nproc; p++)
   {
#ifdef GLVIS_DEBUG
      cout << "connection[" << p << "]: reading initial data ... " << flush;
#endif
      istream &isock = *input_streams[p];
      // assuming the "parallel nproc p" part of the stream has been read
      isock >> ws >> data_type >> ws; // "*_data" / "mesh" / "solution"
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
      if (data_type != "mesh")
      {
         gf_array[p] = new GridFunction(mesh_array[p], isock);
         gf_count++;
      }
#ifdef GLVIS_DEBUG
      cout << "done." << endl;
#endif
   }

   if (gf_count > 0 && gf_count != nproc)
   {
      mfem_error("Input streams contain a mixture of data types!");
   }

   mesh = new Mesh(mesh_array, nproc);
   if (gf_count == 0)
   {
      SetMeshSolution(mesh, grid_f, save_coloring);
      field_type = 2;
   }
   else
   {
      grid_f = new GridFunction(mesh, gf_array, nproc);
      field_type = (grid_f->VectorDim() == 1) ? 0 : 1;
   }

   for (int p = 0; p < nproc; p++)
   {
      delete mesh_array[nproc-1-p];
      delete gf_array[nproc-1-p];
   }

   Extrude1DMeshAndSolution(&mesh, &grid_f, NULL);

   return field_type;
}

void CloseInputStreams(bool parent)
{
   for (int i = 0; i < input_streams.Size(); i++)
   {
      if (parent)
      {
         socketstream *sock = dynamic_cast<socketstream*>(input_streams[i]);
         if (sock) { sock->rdbuf()->socketbuf::close(); }
      }
      delete input_streams[i];
   }
   input_streams.DeleteAll();
}

// Replace a given VectorFiniteElement-based grid function (e.g. from a Nedelec
// or Raviart-Thomas space) with a discontinuous piece-wise polynomial Cartesian
// product vector grid function of the same order.
GridFunction *ProjectVectorFEGridFunction(GridFunction *gf)
{
   if ((gf->VectorDim() == 3) && (gf->FESpace()->GetVDim() == 1))
   {
      int p = gf->FESpace()->GetOrder(0);
      cout << "Switching to order " << p
           << " discontinuous vector grid function..." << endl;
      int dim = gf->FESpace()->GetMesh()->Dimension();
      FiniteElementCollection *d_fec = new L2_FECollection(p, dim, 1);
      FiniteElementSpace *d_fespace =
         new FiniteElementSpace(gf->FESpace()->GetMesh(), d_fec, 3);
      GridFunction *d_gf = new GridFunction(d_fespace);
      d_gf->MakeOwner(d_fec);
      gf->ProjectVectorFieldOn(*d_gf);
      delete gf;
      return d_gf;
   }
   return gf;
}

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
