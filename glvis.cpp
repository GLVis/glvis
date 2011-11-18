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


// GLVis - an OpenGL visualization server based on the MFEM library


#include <limits>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <signal.h>

#include "mfem.hpp"
#include "lib/visual.hpp"

#include <X11/keysym.h>

#include <unistd.h>
extern char ** environ;

void Set_Palette(int);

VisualizationSceneScalarData *vs = NULL;

int portnum = 19916, input = 1, np = 4;
char mesh_file[128], sol_file[128], keys[1000];
Mesh *mesh = NULL;
Vector sol, solu, solv, solw;
int is_gf = 0, gf_component = -1;
GridFunction *grid_f = NULL;
int mac = 0;
int viscount = 0;
int save_coloring = 0;
int keep_attr = 0;

int window_x = 0;
int window_y = 0;
int window_w = 400;
int window_h = 350;
const char *window_titles[2] ={ "GLVis [scalar data]",
                                "GLVis [vector data]" };
istream *script = NULL;
int scr_sol_autoscale = 0;
int scr_running = 0;
int scr_level = 0;
Vector *init_nodes = NULL;
double scr_min_val, scr_max_val;

Array<istream *> input_streams;

// read the mesh and the solution from a file
void ReadSerial();

// chose grid function component and set the input flag
void SetGridFunction();

// set a (checkerboard) solution when only the mesh is given
void SetMeshSolution();

// read the mesh and the solution from multiple files
void ReadParallel();

int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                               const char *sol_prefix, Mesh **mesh_p,
                               GridFunction **sol_p, int keep_attr);

void ReadInputStreams();

void CloseInputStreams();

GridFunction *ProjectVectorFEGridFunction(GridFunction*);

// Read the content of an input stream (e.g. from socket/file)
int ReadStream(istream &is, const char data_type[])
{
   int field_type = 0; // 0 - scalar data, 1 - vector data, (-1) - unknown

   delete mesh; mesh = NULL;
   delete grid_f; grid_f = NULL;
   keys[0] = '\0';
   if (strcmp(data_type,"fem2d_data") == 0)
   {
      mesh = new Mesh(is, 0, 0);
      sol.Load(is, mesh->GetNV());
   }
   else if (strcmp(data_type,"vfem2d_data") == 0 ||
            strcmp(data_type,"vfem2d_data_keys") == 0 )
   {
      field_type = 1;
      mesh = new Mesh(is, 0, 0);
      solu.Load(is, mesh->GetNV());
      solv.Load(is, mesh->GetNV());
      if (strcmp(data_type,"vfem2d_data_keys") == 0)
         is >> keys;
   }
   else if (strcmp(data_type,"fem3d_data") == 0)
   {
      mesh = new Mesh(is, 0, 0);
      sol.Load(is, mesh->GetNV());
   }
   else if (strcmp(data_type,"vfem3d_data") == 0 ||
            strcmp(data_type,"vfem3d_data_keys") == 0 )
   {
      field_type = 1;
      mesh = new Mesh(is, 0, 0);
      solu.Load(is, mesh->GetNV());
      solv.Load(is, mesh->GetNV());
      solw.Load(is, mesh->GetNV());
      if (strcmp(data_type,"vfem3d_data_keys") == 0)
         is >> keys;
   }
   else if (strcmp(data_type,"fem2d_gf_data") == 0 ||
            strcmp(data_type,"fem2d_gf_data_keys") == 0)
   {
      mesh = new Mesh(is, 1, 0);
      grid_f = new GridFunction(mesh, is);
      if (strcmp(data_type,"fem2d_gf_data_keys") == 0)
         is >> keys;
   }
   else if (strcmp(data_type,"vfem2d_gf_data") == 0 ||
            strcmp(data_type,"vfem2d_gf_data_keys") == 0 )
   {
      field_type = 1;
      mesh = new Mesh(is, 1, 0);
      grid_f = new GridFunction(mesh, is);
      if (strcmp(data_type,"vfem2d_gf_data_keys") == 0)
         is >> keys;
   }
   else if (strcmp(data_type,"fem3d_gf_data") == 0 ||
            strcmp(data_type,"fem3d_gf_data_keys") == 0)
   {
      mesh = new Mesh(is, 1, 0);
      grid_f = new GridFunction(mesh, is);
      if (strcmp(data_type,"fem3d_gf_data_keys") == 0)
         is >> keys;
   }
   else if (strcmp(data_type,"vfem3d_gf_data") == 0 ||
            strcmp(data_type,"vfem3d_gf_data_keys") == 0)
   {
      field_type = 1;
      mesh = new Mesh(is, 1, 0);
      grid_f = new GridFunction(mesh, is);
      if (strcmp(data_type,"vfem3d_gf_data_keys") == 0)
         is >> keys;
   }
   else if (strcmp(data_type,"solution") == 0)
   {
      mesh = new Mesh(is, 1, 0);
      grid_f = new GridFunction(mesh, is);
      field_type = (grid_f->VectorDim() == 1) ? 0 : 1;
   }
   else
   {
      field_type = -1;
      cerr << "Unknown data format" << endl;
      cerr << data_type << endl;
   }

   return field_type;
}

int InitVis(int t)
{
   return InitVisualization(window_titles[t], window_x, window_y,
                            window_w, window_h);
}

// Visualize the data in the global variables mesh, sol/grid_f, etc
void StartVisualization(int field_type)
{
   if (field_type != 0 && field_type != 1)
      return;

   if (InitVis(field_type))
   {
      cerr << "Initializing the visualization failed." << endl;
      return;
   }

   communication_thread *comm_thread;

   if (input_streams.Size() > 0)
   {
      auxKeyFunc(XK_space, ToggleThreads);
      glvis_command = new GLVisCommand(&vs, &mesh, &grid_f, &sol,
                                       &scr_sol_autoscale, &keep_attr);
      comm_thread = new communication_thread(input_streams);
   }

   if (field_type == 0)
   {
      if (grid_f)
         grid_f->GetNodalValues(sol);
      if (mesh->Dimension() == 2)
      {
         VisualizationSceneSolution *vss;
         vs = vss = new VisualizationSceneSolution(*mesh, sol);
         if (grid_f)
            vss->SetGridFunction(*grid_f);
      }
      else if (mesh->Dimension() == 3)
      {
         VisualizationSceneSolution3d *vss;
         vs = vss = new VisualizationSceneSolution3d(*mesh, sol);
         if (grid_f)
            vss->SetGridFunction(grid_f);
      }
   }
   else if (field_type == 1)
   {
      if (mesh->Dimension() == 2)
      {
         if (grid_f)
            vs = new VisualizationSceneVector(*grid_f);
         else
            vs = new VisualizationSceneVector(*mesh, solu, solv);
      }
      else if (mesh->Dimension() == 3)
      {
         if (grid_f)
         {
            grid_f = ProjectVectorFEGridFunction(grid_f);
            vs = new VisualizationSceneVector3d(*grid_f);
         }
         else
            vs = new VisualizationSceneVector3d(*mesh, solu, solv, solw);
      }
   }

   if (vs)
   {
      SetVisualizationScene(vs, 3, keys);
   }

   KillVisualization(); // deletes vs
   vs = NULL;
   if (input_streams.Size() > 0)
   {
      glvis_command->Terminate();
      delete comm_thread;
      delete glvis_command;
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
      ifstream imesh(word.c_str());
      if (!imesh)
      {
         cout << "Can not open mesh file: " << word << endl;
         return 1;
      }
      *mp = new Mesh(imesh, 1, 0);
      cout << "mesh: " << word << "; " << flush;
   }

   // read the solution (GridFunction)
   scr >> ws >> word;
   {
      ifstream isol(word.c_str());
      if (!isol)
      {
         cout << "Can not open solution file: " << word << endl;
         delete *mp; *mp = NULL;
         return 2;
      }
      *sp = new GridFunction(*mp, isol);
      cout << "solution: " << word << endl;
   }

   return 0;
}

int ScriptReadParSolution(istream &scr, Mesh **mp, GridFunction **sp)
{
   int np, keep_attr;
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
      cout << "(real attributes); " << flush;
   else
      cout << "(processor attributes); " << flush;
   // read the solution prefix
   scr >> ws >> sol_prefix;
   cout << "solution prefix: " << sol_prefix << endl;

   return ReadParMeshAndGridFunction(np, mesh_prefix.c_str(),
                                     sol_prefix.c_str(), mp, sp, keep_attr);
}

int ScriptReadDisplMesh(istream &scr, Mesh **mp, GridFunction **sp)
{
   Mesh *m;
   string word;

   cout << "Script: mesh: " << flush;
   scr >> ws >> word;
   {
      ifstream imesh(word.c_str());
      if (!imesh)
      {
         cout << "Can not open mesh file: " << word << endl;
         return 1;
      }
      cout << word << endl;
      m = new Mesh(imesh, 1, 0);
   }
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
         vfes = new FiniteElementSpace(m, vfec, mesh->Dimension());
      }

      g = new GridFunction(vfes);
      if (vfec)
         g->MakeOwner(vfec);
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

#ifdef GLVIS_USE_LIBTIFF
   const string ext = ".tif";
#else
   const string ext = ".xwd";
#endif

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
            scr_level = 0;
      }
      else if (word == "solution" || word == "mesh" || word == "psolution")
      {
         Mesh *new_m;
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

         if (new_m->Dimension() == mesh->Dimension() &&
             new_g->VectorDim() == grid_f->VectorDim())
         {
            if (new_m->Dimension() == 2)
            {
               if (new_g->VectorDim() == 1)
               {
                  VisualizationSceneSolution *vss =
                     dynamic_cast<VisualizationSceneSolution *>(vs);
                  new_g->GetNodalValues(sol);
                  vss->NewMeshAndSolution(new_m, &sol, new_g,
                                          scr_sol_autoscale);
               }
               else
               {
                  VisualizationSceneVector *vsv =
                     dynamic_cast<VisualizationSceneVector *>(vs);
                  vsv->NewMeshAndSolution(*new_g, scr_sol_autoscale);
               }
            }
            else
            {
               if (new_g->VectorDim() == 1)
               {
                  VisualizationSceneSolution3d *vss =
                     dynamic_cast<VisualizationSceneSolution3d *>(vs);
                  new_g->GetNodalValues(sol);
                  vss->NewMeshAndSolution(new_m, &sol, new_g,
                                          scr_sol_autoscale);
               }
               else
               {
                  new_g = ProjectVectorFEGridFunction(new_g);
                  VisualizationSceneVector3d *vss =
                     dynamic_cast<VisualizationSceneVector3d *>(vs);
                  vss->NewMeshAndSolution(new_m, new_g, scr_sol_autoscale);
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
         int err;
         scr >> ws >> word;

         cout << "Script: screenshot: " << flush;

         if (Screenshot(word.c_str()))
         {
            cout << "Screenshot(" << word << ") failed." << endl;
            done_one_command = 1;
            continue;
         }
         cout << "-> " << word << ext << flush;

         ostringstream cmd;
         cmd << "convert " << word << ext << " " << word;
         err = system(cmd.str().c_str());
         if (err)
            cout << "; convert failed." << endl;
         else
            cout << " -> " << word << endl;

         cmd.str("");
         cmd << word << ext;
         remove(cmd.str().c_str());

         if (scr_min_val > vs->GetMinV())
            scr_min_val = vs->GetMinV();
         if (scr_max_val < vs->GetMaxV())
            scr_max_val = vs->GetMaxV();
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
            s = 0;
         else if (word == "smooth")
            s = 1;
         else if (word == "cool")
            s = 2;
         if (s != -1)
         {
            vs->SetShading(s);
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
            scr_sol_autoscale = 0;
         }
         else if (word == "on")
         {
            scr_sol_autoscale = 1;
         }
         else if (word == "value")
         {
            scr_sol_autoscale = 2;
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
         SendKeySequence(keys);
         // MyExpose();
      }
      else if (word == "palette")
      {
         int pal;
         scr >> pal;
         cout << "Script: palette: " << pal << endl;
         Set_Palette(pal-1);
         vs->EventUpdateColors();
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
               cout << ' ' << attr_list.Last();
            else if (attr_list.Size() == 257)
               cout << " ... " << flush;
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
      ScriptControl();
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
            return;

         // start the visualization
         break;
      }
      else if (word == "psolution")
      {
         if (ScriptReadParSolution(scr, &mesh, &grid_f))
            return;

         // start the visualization
         break;
      }
      else if (word == "mesh")
      {
         if (ScriptReadDisplMesh(scr, &mesh, &grid_f))
            return;
         if (mesh)
            break;
      }
      else
      {
         cout << "Unknown command in script: " << word << endl;
      }
   }

   scr_level = scr_running = 0;
   auxKeyFunc(XK_space, ScriptControl);
   script = &scr;
   keys[0] = '\0';

   StartVisualization((grid_f->VectorDim() == 1) ? 0 : 1);

   delete init_nodes; init_nodes = NULL;

   cout << "Script: min_val = " << scr_min_val
        << ", max_val = " << scr_max_val << endl;

   script = NULL;
}


int main (int argc, char *argv[])
{
   int i, childPID;
   int multi_session = 1;
   int saved = 0;

   char data_type[32];
   keys[0] = '\0';

   // parse the program arguments
   if (argc != 1)
   {
      for (i = 1; i < argc; i++)
         if (strcmp("-m", argv[i])==0)
            strcpy(mesh_file, argv[++i]),        input |= 2;
         else if (strcmp("-s", argv[i])==0)
            strcpy(sol_file, argv[++i]),         input |= 4;
         else if (strcmp("-v", argv[i])==0)
            strcpy(sol_file, argv[++i]),         input |= 8;
         else if (strcmp("-g", argv[i])==0)
            strcpy(sol_file, argv[++i]),         is_gf = 255;
         else if (strcmp("-np", argv[i])==0)
            np  = atoi(argv[++i]),               input |= 256;
         else if (strcmp("-p", argv[i])==0)
            portnum = atoi(argv[++i]);
         else if (strcmp("-t", argv[i])==0)
            multi_session = 0;
         else if (strcmp("-k", argv[i])==0)
            strcpy(keys, argv[++i]);
         else if (strcmp("-saved", argv[i])==0)
            saved = 1;
         else if (strcmp("-mac", argv[i])==0)
            mac = 1;
         else if (strcmp("-sc", argv[i])==0)
            save_coloring = 1;
         else if (strcmp("-gc", argv[i])==0)
            gf_component = atoi(argv[++i]);
         else if (strcmp("-a", argv[i])==0)
            keep_attr = 1;
         else if (strcmp("-run", argv[i])==0 && i == 1 && argc == 3)
         {
            ifstream scr(argv[2]);
            if (!scr)
            {
               cout << "Can not open script: " << argv[2] << endl;
               return 1;
            }
            PlayScript(scr);
            return 0;
         }
         else // display help
            input |= 31415;
   }

   if (saved)
   {
      // This is the child process after exec()
      ifstream ifs(argv[2]);
      ifs >> data_type >> ws;
      int ft = ReadStream(ifs, data_type);
      ifs.close();
      StartVisualization(ft);
      exit(0);
   }

   cout << endl
        << "       _/_/_/  _/      _/      _/  _/"          << endl
        << "    _/        _/      _/      _/        _/_/_/" << endl
        << "   _/  _/_/  _/      _/      _/  _/  _/_/"      << endl
        << "  _/    _/  _/        _/  _/    _/      _/_/"   << endl
        << "   _/_/_/  _/_/_/_/    _/      _/  _/_/_/"      << endl
        << endl ;

   // print help for wrong input
   if (input != 1  && input != 3 && input != 7 && input != 11 && input != 259)
   {
      cerr <<
         "Usage:\n"
         "  glvis [-mac] [-t] [-p <port_number>] [-k keys] [-a]\n"
         "or\n"
         "  glvis -m <mesh_file> [-k keys] [-sc]\n"
         "or\n"
         "  glvis -m <mesh_file> -s <scalar_solution_file> [-k keys]\n"
         "or\n"
         "  glvis -m <mesh_file> -v <vector_solution_file> [-k keys]\n"
         "or\n"
         "  glvis -m <mesh_file> -g <grid_function_solution_file> [-gc <component>] [-k keys]\n"
         "or\n"
         "  glvis -np <#proc> -m <mesh_prefix> [-g <grid_function_prefix>] [-gc <component>] [-k keys] [-a]\n"
         "or\n"
         "  glvis -saved <glvis-saved-file>\n"
         "or\n"
         "  glvis -run <glvis-script>\n" << endl;
      exit(1);
   }

   int nproc = 1, proc = 0;

   // server mode, read the mesh and the solution from a socket
   if (input == 1)
   {
      // get rid of zombies
      if (multi_session)
         signal(SIGCHLD, SIG_IGN);

      socketserver server(portnum);
      if (server.good())
         cout << "Waiting for data on port "
              << portnum << " ...  " << endl;
      else
      {
         cout << "Server already running on port " << portnum << ".\n" << endl;
         return 2;
      }

      socketstream *isock = new socketstream;
      while (1)
      {
         while (server.accept(*isock) < 0);

         *isock >> data_type >> ws;

         if (mac)
            viscount++;

         int par_data = 0;
         if (strcmp(data_type,"parallel") == 0)
         {
            par_data = 1;
            np = 0;
            do
            {
               *isock >> nproc >> proc >> ws;
               input_streams.SetSize(nproc);
               input_streams[proc] = isock;
               isock = new socketstream;
               np++;
               if (np == nproc)
                  break;
               // read next available socket stream
               while (server.accept(*isock) < 0);
               *isock >> data_type >> ws; // "parallel"
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
                  CloseInputStreams();
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
                  signal(SIGINT, SIG_IGN);
               int ft;
               if (!par_data)
               {
                  ft = ReadStream(*isock, data_type);
                  input_streams.Append(isock);
               }
               else
               {
                  delete isock;
                  ReadInputStreams();
                  ft = (grid_f->VectorDim() == 1) ? 0 : 1;
               }
               StartVisualization(ft);
               CloseInputStreams();
               exit(0);
            }

         default :                     // This is the parent process
            if (!par_data)
               isock->close();
            else
               CloseInputStreams();
         }
      }
   }
   else  // input != 1, non-server mode
   {
      if (input & 256)
         ReadParallel();
      else
         ReadSerial();

      int window_err;
      if (mesh->Dimension() == 2)
      {
         if ((input & 8) == 0)
         {
            VisualizationSceneSolution *vss;
            window_err = InitVis(0);
            if (!window_err)
            {
               if ((input & 4) == 0)
                  Set_Palette(4); // Set_Palette(11);
               vs = vss = new VisualizationSceneSolution(*mesh, sol);
               if (is_gf)
               {
                  vss->SetGridFunction(*grid_f);
                  vss->ToggleShading();
               }
               if ((input & 4) == 0)
               {
                  vs->OrthogonalProjection = 1;
                  vs->light = 0;
                  vs->Zoom(1.8);
                  // increase the refinement factors if the mesh is "curved"
                  if (mesh->GetNodalFESpace())
                  {
                     cout << "Curvilinear mesh, subdivision factors = 4, 1"
                          << endl;
                     vs->SetRefineFactors(4, 1);
                  }
                  if (grid_f)
                     vs->SetValueRange(-1.0*(grid_f->Max()+1.0),
                                       1.0*(grid_f->Max()+1.0));
                  else
                     vs->SetValueRange(-1.0*(sol.Max()+1.0),
                                       1.0*(sol.Max()+1.0));
               }
            }
         }
         else
         {
            window_err = InitVis(1);
            if (!window_err)
            {
               if (is_gf)
                  vs = new VisualizationSceneVector(*grid_f);
               else
                  vs = new VisualizationSceneVector(*mesh, solu, solv);
            }
         }
      }
      else // 3D
      {
         if ((input & 8) == 0 && (input & 512) == 0)
         {
            VisualizationSceneSolution3d *vss;
            window_err = InitVis(0);
            if (!window_err)
            {
               vs = vss = new VisualizationSceneSolution3d(*mesh, sol);
               if (is_gf)
               {
                  vss->SetGridFunction (grid_f);
                  vss->ToggleShading();
                  vss->Prepare();
                  // vss -> Draw(); // ???
               }
               if ((input & 4) == 0)
               {
                  // Set_Palette(4);
                  Set_Palette(11);
                  Set_Material_And_Light(4,3);
                  vss->ToggleDrawAxes();
                  vss->ToggleDrawMesh();
                  vs->SetValueRange(-1.0*(grid_f->Max()+1.0),
                                    1.0*(grid_f->Max()+1.0));
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
                  vs = new VisualizationSceneVector3d(*mesh, solu, solv, solw);
            }
         }
      }
      if (!window_err)
      {
         if (mesh->Dimension() == 2 && (input & 12) == 0)
            SetVisualizationScene(vs, 2, keys);
         else
            SetVisualizationScene(vs, 3, keys);
         KillVisualization(); // deletes vs
         if (is_gf)  delete grid_f;
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
      mesh = new Mesh(meshin, 1, 0);
   }

   if (is_gf || (input & 4) || (input & 8))
   {
      ifstream solin;
      // get the solution from file
      solin.open(sol_file);
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
         if (mesh->Dimension() == 3)
            solw.Load(solin, mesh->GetNV());
      }
   }
   else
      SetMeshSolution();
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
         (*new_gf)(i) = (*grid_f)(ofes->DofToVDof(i, gf_component));
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


void SetMeshSolution()
{
   if (1) // checkerboard solution
   {
      FiniteElementCollection *cfec;
      if (mesh->Dimension() == 2)
         cfec = new Const2DFECollection;
      else
         cfec = new Const3DFECollection;
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
            (*grid_f)(i) = coloring[i];
         cout << "Number of colors: " << grid_f->Max() + 1 << endl;
      }
      grid_f->GetNodalValues(sol);
      is_gf = 1;
      if (save_coloring)
      {
         ofstream fgrid("GLVis_coloring.gf");
         cout << "Saving the coloring function -> " << flush;
         grid_f->Save(fgrid);
         cout << "done." << endl;
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
         SetGridFunction();
   }
   else
   {
      err = ReadParMeshAndGridFunction(np, mesh_file, NULL,
                                       &mesh, NULL, keep_attr);
      if (!err)
         SetMeshSolution();
   }

   if (err)
      exit(1);
}

int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                               const char *sol_prefix, Mesh **mesh_p,
                               GridFunction **sol_p, int keep_attr)
{
   Array<Mesh *> mesh_array;

   mesh_array.SetSize(np);
   ifstream meshfile;
   for (int p = 0; p < np; p++)
   {
      ostringstream fname;
      fname << mesh_prefix << '.' << setfill('0') << setw(6) << p;
      meshfile.open(fname.str().c_str());
      if (!meshfile)
      {
         cerr << "Can not open mesh file: " << fname.str().c_str()
              << '!' << endl;
         for (p--; p >= 0; p--)
            delete mesh_array[p];
         return 1;
      }
      mesh_array[p] = new Mesh(meshfile, 1, 0);
      if (!keep_attr)
      {
         // set element and boundary attributes to be the processor number + 1
         for (int i = 0; i < mesh_array[p]->GetNE(); i++)
            mesh_array[p]->GetElement(i)->SetAttribute(p+1);
         for (int i = 0; i < mesh_array[p]->GetNBE(); i++)
            mesh_array[p]->GetBdrElement(i)->SetAttribute(p+1);
      }
      meshfile.close();
   }
   *mesh_p = new Mesh(mesh_array, np);

   if (sol_prefix && sol_p)
   {
      Array<GridFunction *> gf_array(np);
      ifstream solfile;
      for (int p = 0; p < np; p++)
      {
         ostringstream fname;
         fname << sol_prefix << '.' << setfill('0') << setw(6) << p;
         solfile.open(fname.str().c_str());
         if (!solfile)
         {
            cerr << "Can not open solution file " << fname.str().c_str()
                 << '!' << endl;
            for (p--; p >= 0; p--)
               delete gf_array[p];
            delete *mesh_p;
            *mesh_p = NULL;
            for (p = 0; p < np; p++)
               delete mesh_array[np-1-p];
            return 2;
         }
         gf_array[p] = new GridFunction(mesh_array[p], solfile);
         solfile.close();
      }
      *sol_p = new GridFunction(*mesh_p, gf_array, np);

      for (int p = 0; p < np; p++)
         delete gf_array[np-1-p];
   }

   for (int p = 0; p < np; p++)
      delete mesh_array[np-1-p];

   return 0;
}

void ReadInputStreams()
{
   int nproc = input_streams.Size();
   Array<Mesh *> mesh_array(nproc);
   Array<GridFunction *> gf_array(nproc);
   string data_type;

   for (int p = 0; p < nproc; p++)
   {
      istream &isock = *input_streams[p];
      // assuming the "parallel nproc p" part of the stream has been read
      isock >> data_type >> ws; // "*_data" / "solution"
      mesh_array[p] = new Mesh(isock, 1, 0);
      if (!keep_attr)
      {
         // set element and boundary attributes to proc+1
         for (int i = 0; i < mesh_array[p]->GetNE(); i++)
            mesh_array[p]->GetElement(i)->SetAttribute(p+1);
         for (int i = 0; i < mesh_array[p]->GetNBE(); i++)
            mesh_array[p]->GetBdrElement(i)->SetAttribute(p+1);
      }
      gf_array[p] = new GridFunction(mesh_array[p], isock);
   }

   mesh = new Mesh(mesh_array, nproc);
   grid_f = new GridFunction(mesh, gf_array, nproc);

   for (int p = 0; p < nproc; p++)
   {
      delete gf_array[nproc-1-p];
      delete mesh_array[nproc-1-p];
   }
}

void CloseInputStreams()
{
   for (int i = 0; i < input_streams.Size(); i++)
      delete input_streams[i];
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
      FiniteElementCollection *d_fec = new L2_FECollection(p, 3);
      FiniteElementSpace *d_fespace = new FiniteElementSpace(mesh, d_fec, 3);
      GridFunction *d_gf = new GridFunction(d_fespace);
      d_gf->MakeOwner(d_fec);
      gf->ProjectVectorFieldOn(*d_gf);
      delete gf;
      return d_gf;
   }
   return gf;
}
