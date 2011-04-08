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

#include "mfem.hpp"
#include "lib/visual.hpp"

#include <X11/keysym.h>

#include <unistd.h>
extern char ** environ;

void Set_Palette(int);

VisualizationSceneScalarData *vs = NULL;

int portnum = 19916, input = 1, np = 4;
char mesh_file[128], sol_file[128], keys[1000];
Mesh * mesh = NULL;
Vector sol, solu, solv, solw, *data[3] = {NULL, NULL, NULL};
int is_gf = 0, gf_component = -1;
GridFunction *grid_f = NULL;
int mac = 0;
int viscount = 0;
int save_coloring = 0;

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

// read the mesh and the solution from a file
void ReadSerial();

// read the mesh and the solution from multiple files
void ReadParallel();

int InitVis(int t)
{
   return InitVisualization(window_titles[t], window_x, window_y,
                            window_w, window_h);
}

// Visualize the content of a string stream (e.g. from socket/file)
void VisStream(istream * iss, char * data_type)
{
   int window_err = 0;

   if (strcmp(data_type,"fem2d_data") == 0)
   {
      window_err = InitVis(0);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 0, 0);
         sol.Load(*iss, mesh->GetNV());

         vs = new VisualizationSceneSolution(*mesh, sol);
      }
   }
   else if (strcmp(data_type,"vfem2d_data") == 0 ||
            strcmp(data_type,"vfem2d_data_keys") == 0 )
   {
      window_err = InitVis(1);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 0, 0);
         solu.Load(*iss, mesh->GetNV());
         solv.Load(*iss, mesh->GetNV());

         vs = new VisualizationSceneVector(*mesh, solu, solv);

         if (strcmp(data_type,"vfem2d_data_keys") == 0)
            *iss >> keys;
      }
   }
   else if (strcmp(data_type,"fem3d_data") == 0)
   {
      window_err = InitVis(0);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 0, 0);
         sol.Load(*iss, mesh->GetNV());

         vs = new VisualizationSceneSolution3d(*mesh, sol);
      }
   }
   else if (strcmp(data_type,"vfem3d_data") == 0 ||
            strcmp(data_type,"vfem3d_data_keys") == 0 )
   {
      window_err = InitVis(1);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 0, 0);
         solu.Load(*iss, mesh->GetNV());
         solv.Load(*iss, mesh->GetNV());
         solw.Load(*iss, mesh->GetNV());

         vs = new VisualizationSceneVector3d(*mesh, solu, solv, solw);

         if (strcmp(data_type,"vfem3d_data_keys") == 0)
            *iss >> keys;
      }
   }
   else if (strcmp(data_type,"fem3d_paragrid") == 0)
   {
      window_err = InitVis(0);
      if (!window_err)
         vs = new VisualizationSceneSolution3d(*mesh, *data[0]);
   }
   else if (strcmp(data_type,"vfem3d_paragrid") == 0)
   {
      window_err = InitVis(1);
      if (!window_err)
         vs = new VisualizationSceneVector3d(*mesh, *data[0],
                                             *data[1], *data[2]);
   }
   else if (strcmp(data_type,"fem2d_gf_data") == 0 ||
            strcmp(data_type,"fem2d_gf_data_keys") == 0)
   {
      window_err = InitVis(0);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 1, 0);
         *iss >> ws;
         grid_f = new GridFunction (mesh, *iss);
         grid_f -> GetNodalValues (sol);
         VisualizationSceneSolution *vss;
         vs = vss = new VisualizationSceneSolution(*mesh, sol);
         vss -> SetGridFunction (*grid_f);
         if (strcmp(data_type,"fem2d_gf_data_keys") == 0)
            *iss >> keys;
      }
   }
   else if (strcmp(data_type,"vfem2d_gf_data") == 0 ||
            strcmp(data_type,"vfem2d_gf_data_keys") == 0 )
   {
      window_err = InitVis(1);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 1, 0);
         *iss >> ws;
         grid_f = new GridFunction (mesh, *iss);
         vs = new VisualizationSceneVector (*grid_f);
         if (strcmp(data_type,"vfem2d_gf_data_keys") == 0)
            *iss >> keys;
      }
   }
   else if (strcmp(data_type,"fem3d_gf_data") == 0 ||
            strcmp(data_type,"fem3d_gf_data_keys") == 0)
   {
      window_err = InitVis(0);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 1, 0);
         *iss >> ws;
         grid_f = new GridFunction (mesh, *iss);
         grid_f -> GetNodalValues (sol);
         VisualizationSceneSolution3d *vss;
         vs = vss = new VisualizationSceneSolution3d(*mesh, sol);
         vss -> SetGridFunction (grid_f);
         if (strcmp(data_type,"fem3d_gf_data_keys") == 0)
            *iss >> keys;
      }
   }
   else if (strcmp(data_type,"vfem3d_gf_data") == 0 ||
            strcmp(data_type,"vfem3d_gf_data_keys") == 0)
   {
      window_err = InitVis(1);
      if (!window_err)
      {
         mesh = new Mesh(*iss, 1, 0);
         *iss >> ws;
         grid_f = new GridFunction (mesh, *iss);
         vs = new VisualizationSceneVector3d (*grid_f);
         if (strcmp(data_type,"vfem3d_gf_data_keys") == 0)
            *iss >> keys;
      }
   }
   else
   {
      cerr << "Unknown data format" << endl;
      cerr << data_type << endl;
      return;
   }
   if (!window_err)
   {
      SetVisualizationScene (vs, 3, keys);
      KillVisualization(); // deletes vs
      vs = NULL;
      delete grid_f; grid_f = NULL;
      delete mesh; mesh = NULL;
   }
   else
   {
      cerr << "Initializing the visualization failed." << endl;
   }
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
      else if (word == "solution" || word == "mesh")
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
         else // word == "mesh"
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
                  // 3D vector field ...
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

   int window_err;
   if (grid_f->VectorDim() == 1)
      window_err = InitVis(0);
   else
      window_err = InitVis(1);
   if (window_err)
   {
      cerr << "Initializing the visualization failed." << endl;
      return;
   }

   if (mesh->Dimension() == 2)
   {
      if (grid_f->VectorDim() == 1)
      {
         grid_f->GetNodalValues(sol);
         VisualizationSceneSolution *vss;
         vs = vss = new VisualizationSceneSolution(*mesh, sol);
         vss->SetGridFunction(*grid_f);
      }
      else
      {
         vs = new VisualizationSceneVector(*grid_f);
      }
   }
   else if (mesh->Dimension() == 3)
   {
      if (grid_f->VectorDim() == 1)
      {
         grid_f->GetNodalValues(sol);
         VisualizationSceneSolution3d *vss;
         vs = vss = new VisualizationSceneSolution3d(*mesh, sol);
         vss->SetGridFunction(grid_f);
      }
      else
      {
         // ...
      }
   }

   scr_level = scr_running = 0;
   auxKeyFunc(XK_space, ScriptControl);
   script = &scr;
   scr_level = scr_running = 0;

   SetVisualizationScene (vs, 3);
   KillVisualization(); // deletes vs
   vs = NULL;
   delete grid_f; grid_f = NULL;
   delete mesh; mesh = NULL;
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

   isockstream *in;
   istringstream *iss = NULL;

   char data_type[32];

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
         else if (strcmp("-par3d", argv[i])==0)  input |= 128;
         else if (strcmp("-np", argv[i])==0)
            np  = atoi(argv[++i]),               input |= 256;
         else if (strcmp("-vpar3d", argv[i])==0) input |= 512;
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
      ifs >> data_type;
      VisStream(&ifs, data_type);
      ifs.close();
      exit(0);
   }

   cout << endl
        << "       _/_/_/  _/      _/      _/  _/"          << endl
        << "    _/        _/      _/      _/        _/_/_/" << endl
        << "   _/  _/_/  _/      _/      _/  _/  _/_/"      << endl
        << "  _/    _/  _/        _/  _/    _/      _/_/"   << endl
        << "   _/_/_/  _/_/_/_/    _/      _/  _/_/_/"      << endl
        << endl ;

   // get rid of zombies
   signal(SIGCHLD, SIG_IGN);

   // print help for wrong input
   if (input != 1  && input != 3 && input != 7 && input != 11 &&
       input != 391 && input != 775)
   {
      cerr << "Usage:" << endl
           << "  glvis [-mac] [-t] [-p <port_number>] [-k keys]" << endl
           << "or" << endl
           << "  glvis -m <mesh_file> [-k keys] [-sc]\n"
           << "or" << endl
           << "  glvis -m <mesh_file> -s <scalar_solution_file> [-k keys]\n"
           << "or" << endl
           << "  glvis -m <mesh_file> -v <vector_solution_file> [-k keys]\n"
           << "or" << endl
           << "  glvis -m <mesh_file> -g <grid_function_solution_file> [-gc <component>] [-k keys]\n"
           << "or" << endl
           << "  glvis {-par3d|-vpar3d} -np <#domains> -m <mesh_prefix> -s <solution_prefix> [-k keys]" << endl
           << "or" << endl
           << "  glvis -saved <glvis-saved-file>" << endl
           << "or" << endl
           << "  glvis -run <glvis-script>" << endl
           << endl;
      exit(1);
   }

   int nprocessors = 1, currentprocessor = 0;

   // server mode, read the mesh and the solution from a socket
   if (input == 1)
   {
      in = new isockstream(portnum);
      if (in->good())
         cout << "Waiting for data on port "
              << portnum << " ...  " << endl;
      else
         return 2;

      while (1)
      {
         do
         {
            in->receive(&iss);
         }
         while(!in->good());

         (*iss) >> setw(32) >> data_type;

         if (mac)
            viscount++;

         int np;
         if (strcmp(data_type,"fem3d_paragrid") == 0)
         {
            np = 0;
            do
            {
               *iss >> nprocessors >> currentprocessor;
               if (np == 0)
               {
                  if (mesh != NULL)
                     delete mesh;
                  data[0] = &sol;
                  mesh = new Mesh(*iss, data, nprocessors, currentprocessor);
               }
               else
                  mesh -> Load (*iss, data, nprocessors, currentprocessor);

               if (nprocessors - 1 == np)
                  break;

               in->receive(&iss);
               *iss >> data_type;
               np++;
            }
            while ( 1 );
         }

         if (strcmp(data_type,"vfem3d_paragrid") == 0)
         {
            np = 0;
            do
            {
               *iss >> nprocessors >> currentprocessor;
               if (np == 0)
               {
                  if (mesh != NULL)
                     delete mesh;
                  data[0] = &solu;
                  data[1] = &solv;
                  data[2] = &solw;
                  mesh = new Mesh ( *iss, data, nprocessors, currentprocessor);
               }
               else
                  mesh -> Load ( *iss, data, nprocessors, currentprocessor);

               if (nprocessors - 1 == np)
                  break;

               in->receive(&iss);
               *iss >> data_type;
               np++;
            }
            while ( 1 );
         }

         char tmp_file[50];
         if (multi_session)
         {
            if (mac)
            {
               sprintf(tmp_file,"glvis-saved.%04d",viscount);
               ofstream ofs(tmp_file);
               ofs << data_type << '\n' << iss -> rdbuf();
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
            delete in;
            if (mac)
            {
               // exec ourself
               const char *args[3] = { argv[0], "-saved", tmp_file };
               execve(args[0], (char* const*)args, environ);
               exit(0);
            }
            else
            {
               if (multi_session)
                  signal(SIGINT, SIG_IGN);
               VisStream(iss, data_type);
               delete iss;
               exit(0);
            }

         default :                     // This is the parent process
            data[0] = data[1] = data[2] = NULL;
         }
      }

      // delete in;
   }
   else
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
                  vs = new VisualizationSceneVector3d(*grid_f);
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
   ifstream solin;

   // get the mesh from a file
   {
      ifstream meshin(mesh_file);
      if (!meshin)
      {
         cerr << "Can not open mesh file " << mesh_file << ". Exit. \n";
         exit(1);
      }

      // mesh = new Mesh(meshin, is_gf?1:0, 0);
      mesh = new Mesh(meshin, 1, 0);
   }

   if (is_gf || (input & 4) || (input & 8))
   {
      // get the solution from file
      ifstream solin(sol_file);
      if (!solin)
      {
         cerr << "Can not open solution file " << sol_file << ". Exit.\n";
         exit(1);
      }
      if (is_gf)
      {
         grid_f = new GridFunction(mesh, solin);
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
   {
      if (1) // if not given one, use a checkerboard solution
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
      else // use a zero solution by default
      {
         sol.SetSize (mesh -> GetNV());
         sol = 0.0;
      }
   }
}


void ReadParallel()
{
   int i;
   char fname[100];
   istream ** parin = new istream* [np];

   int * dim = new int[np];
   for (i = 0; i < np; i++)
   {
      sprintf(fname,"%s_%d",mesh_file,i);
      parin[i] = new ifstream(fname);
      if (!(*parin[i]))
      {
         cerr << "Can not open mesh file " << fname << ". Exit. \n";
         exit(1);
      }
   }
   mesh = new Mesh(parin,np,dim);
   for (i = 0; i < np; i++) delete parin[i];

   ofstream meshout("mesh.out");
   meshout.setf(ios::scientific);
   meshout.precision(12);
   mesh -> Print(meshout);

   for (i = 0; i < np; i++)
   {
      sprintf(fname,"%s_%d",sol_file,i);
      parin[i] = new ifstream(fname);
      if (!(*parin[i]))
      {
         cerr << "Can not open solution file " << fname << ". Exit. \n";
         exit(1);
      }
   }

   if (input & 128)
   {
      // get rid of NetGen's info line
      for (i = 0; i < np; i++)
         parin[i] -> getline(fname,128);
      sol.Load(parin,np,dim);
   }
   else if (input & 512)
   {
      solu.Load(parin, np, dim);
      solv.Load(parin, np, dim);
      solw.Load(parin, np, dim);
   }

   for (i = 0; i < np; i++) delete parin[i];
   delete [] parin;
   delete [] dim;
}
