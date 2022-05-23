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


// GLVis - an OpenGL visualization server based on the MFEM library

#include <limits>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include <ctime>

// SDL may redefine main() as SDL_main() ostensibly to ease portability.
// (WinMain() instead of main() is used as the entry point in a non-console
// Windows program.)
//
// We must instead define SDL_MAIN_HANDLED so that SDL doesn't do this
// substitution, since we need a console to accept certain user input from
// stdin.
#ifdef _WIN32
#define SDL_MAIN_HANDLED
#endif

#include "mfem.hpp"
#include "lib/palettes.hpp"
#include "lib/visual.hpp"
#include "lib/stream_reader.hpp"

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
int         pad_digits      = 6;
int         gf_component    = -1;
bool        keep_attr       = false;
int         window_x        = 0; // not a command line option
int         window_y        = 0; // not a command line option
int         window_w        = 400;
int         window_h        = 350;
const char *window_title    = string_default;
const char *c_plot_caption  = string_none;
thread_local string      plot_caption;
thread_local string      extra_caption;
bool        secure          = socketstream::secure_default;

// Global variables
int input = 1;
thread_local StreamState stream_state;
thread_local VisualizationSceneScalarData *vs = NULL;
extern thread_local GLVisCommand* glvis_command;
thread_local communication_thread *comm_thread = NULL;

thread_local GeometryRefiner GLVisGeometryRefiner;

const char *window_titles[] = { "GLVis [scalar data]",
                                "GLVis [vector data]", "GLVis [mesh]"
                              };
istream *script = NULL;
int scr_running = 0;
int scr_level = 0;
Vector *init_nodes = NULL;
double scr_min_val, scr_max_val;

extern char **environ;

using StreamCollection = vector<unique_ptr<istream>>;

void PrintSampleUsage(ostream &out);

// read the mesh and the solution from a file
void ReadSerial(StreamState& state);

// choose grid function component and set the input flag
void SetGridFunction(StreamState& state);

// read the mesh and the solution from multiple files
void ReadParallel(int np, StreamState& state);

int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                               const char *sol_prefix, StreamState& state,
                               int keep_attr);

int ReadInputStreams(StreamState& state, const StreamCollection& input_streams);

// Visualize the data in the global variables mesh, sol/grid_f, etc
// 0 - scalar data, 1 - vector data, 2 - mesh only, (-1) - unknown
bool GLVisInitVis(int field_type, StreamCollection input_streams)
{
   if (field_type < 0 || field_type > 2)
   {
      return false;
   }

   const char *win_title = (window_title == string_default) ?
                           window_titles[field_type] : window_title;

   if (InitVisualization(win_title, window_x, window_y, window_w, window_h))
   {
      cerr << "Initializing the visualization failed." << endl;
      return false;
   }

   if (input_streams.size() > 0)
   {
      GetAppWindow()->setOnKeyDown(SDLK_SPACE, ThreadsPauseFunc);
      glvis_command = new GLVisCommand(&vs, stream_state, &keep_attr);
      comm_thread = new communication_thread(std::move(input_streams), glvis_command);
   }

   double mesh_range = -1.0;
   if (field_type == 0 || field_type == 2)
   {
      if (stream_state.grid_f)
      {
         stream_state.grid_f->GetNodalValues(stream_state.sol);
      }
      if (stream_state.mesh->SpaceDimension() == 2)
      {
         VisualizationSceneSolution * vss;
         if (stream_state.normals.Size() > 0)
         {
            vs = vss = new VisualizationSceneSolution(*stream_state.mesh, stream_state.sol,
                                                      &stream_state.normals);
         }
         else
         {
            vs = vss = new VisualizationSceneSolution(*stream_state.mesh, stream_state.sol);
         }
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(*stream_state.grid_f);
         }
         if (field_type == 2)
         {
            vs->OrthogonalProjection = 1;
            vs->SetLight(false);
            vs->Zoom(1.8);
            // Use the 'bone' palette when visualizing a 2D mesh only (otherwise
            // the 'jet-like' palette is used in 2D, see vssolution.cpp).
            vs->palette.SetIndex(4);
         }
      }
      else if (stream_state.mesh->SpaceDimension() == 3)
      {
         VisualizationSceneSolution3d * vss;
         vs = vss = new VisualizationSceneSolution3d(*stream_state.mesh,
                                                     stream_state.sol);
         if (stream_state.grid_f)
         {
            vss->SetGridFunction(stream_state.grid_f.get());
         }
         if (field_type == 2)
         {
            if (stream_state.mesh->Dimension() == 3)
            {
               // Use the 'white' palette when visualizing a 3D volume mesh only
               vss->palette.SetIndex(11);
               vss->SetLightMatIdx(4);
            }
            else
            {
               // Use the 'bone' palette when visualizing a surface mesh only
               vss->palette.SetIndex(4);
            }
            // Otherwise, the 'vivid' palette is used in 3D see vssolution3d.cpp
            vss->ToggleDrawAxes();
            vss->ToggleDrawMesh();
         }
      }
      if (field_type == 2)
      {
         if (stream_state.grid_f)
         {
            mesh_range = stream_state.grid_f->Max() + 1.0;
         }
         else
         {
            mesh_range = stream_state.sol.Max() + 1.0;
         }
      }
   }
   else if (field_type == 1)
   {
      if (stream_state.mesh->SpaceDimension() == 2)
      {
         if (stream_state.grid_f)
         {
            vs = new VisualizationSceneVector(*stream_state.grid_f);
         }
         else
         {
            vs = new VisualizationSceneVector(*stream_state.mesh, stream_state.solu,
                                              stream_state.solv);
         }
      }
      else if (stream_state.mesh->SpaceDimension() == 3)
      {
         if (stream_state.grid_f)
         {
            stream_state.grid_f
               = ProjectVectorFEGridFunction(std::move(stream_state.grid_f));
            vs = new VisualizationSceneVector3d(*stream_state.grid_f);
         }
         else
         {
            vs = new VisualizationSceneVector3d(*stream_state.mesh, stream_state.solu,
                                                stream_state.solv, stream_state.solw);
         }
      }
   }

   if (vs)
   {
      // increase the refinement factors if visualizing a GridFunction
      if (stream_state.grid_f)
      {
         vs->AutoRefine();
         vs->SetShading(2, true);
      }
      if (mesh_range > 0.0)
      {
         vs->SetValueRange(-mesh_range, mesh_range);
         vs->SetAutoscale(0);
      }
      if (stream_state.mesh->SpaceDimension() == 2 && field_type == 2)
      {
         SetVisualizationScene(vs, 2, stream_state.keys.c_str());
      }
      else
      {
         SetVisualizationScene(vs, 3, stream_state.keys.c_str());
      }
   }
   return true;
}

void GLVisStartVis()
{
   RunVisualization(); // deletes vs
   vs = NULL;
   if (glvis_command)
   {
      glvis_command->Terminate();
      delete comm_thread;
      delete glvis_command;
      glvis_command = NULL;
   }
   cout << "GLVis window closed." << endl;
}

int ScriptReadSolution(istream &scr, StreamState& state)
{
   string mword,sword;

   cout << "Script: solution: " << flush;
   // read the mesh
   scr >> ws >> mword; // mesh filename (can't contain spaces)
   cout << "mesh: " << mword << "; " << flush;
   named_ifgzstream imesh(mword.c_str());
   if (!imesh)
   {
      cout << "Can not open mesh file: " << mword << endl;
      return 1;
   }
   state.mesh.reset(new Mesh(imesh, 1, 0, state.fix_elem_orient));

   // read the solution (GridFunction)
   scr >> ws >> sword;
   if (sword == mword) // mesh and solution in the same file
   {
      cout << "solution: " << mword << endl;
      state.grid_f.reset(new GridFunction(state.mesh.get(), imesh));
   }
   else
   {
      cout << "solution: " << sword << endl;
      ifgzstream isol(sword.c_str());
      if (!isol)
      {
         cout << "Can not open solution file: " << sword << endl;
         return 2;
      }
      state.grid_f.reset(new GridFunction(state.mesh.get(), isol));
   }

   state.Extrude1DMeshAndSolution();

   return 0;
}

int ScriptReadParSolution(istream &scr, StreamState& state)
{
   int np, scr_keep_attr, err_read;
   string mesh_prefix, sol_prefix;

   cout << "Script: psolution: " << flush;
   // read number of processors
   scr >> np;
   cout << "# processors: " << np << "; " << flush;
   // read the mesh prefix
   scr >> ws >> mesh_prefix; // mesh prefix (can't contain spaces)
   cout << "mesh prefix: " << mesh_prefix << "; " << flush;
   scr >> ws >> scr_keep_attr;
   if (scr_keep_attr)
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

   err_read = ReadParMeshAndGridFunction(np, mesh_prefix.c_str(),
                                         sol_prefix.c_str(), state, scr_keep_attr);
   if (!err_read)
   {
      state.Extrude1DMeshAndSolution();
   }
   return err_read;
}

int ScriptReadDisplMesh(istream &scr, StreamState& state)
{
   StreamState meshstate;
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
      meshstate.mesh.reset(new Mesh(imesh, 1, 0, state.fix_elem_orient));
   }
   meshstate.Extrude1DMeshAndSolution();
   Mesh* const m = meshstate.mesh.get();
   if (init_nodes == NULL)
   {
      init_nodes = new Vector;
      meshstate.mesh->GetNodes(*init_nodes);
      state.mesh = NULL;
      state.grid_f = NULL;
   }
   else
   {
      FiniteElementCollection  *vfec = NULL;
      FiniteElementSpace *vfes;
      vfes = (FiniteElementSpace *)m->GetNodalFESpace();
      if (vfes == NULL)
      {
         vfec = new LinearFECollection;
         vfes = new FiniteElementSpace(m, vfec, m->SpaceDimension());
      }

      meshstate.grid_f.reset(new GridFunction(vfes));
      GridFunction * const g = meshstate.grid_f.get();
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

      state.mesh = std::move(meshstate.mesh);
      state.grid_f = std::move(meshstate.grid_f);
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
         StreamState new_state;

         if (word == "solution")
         {
            if (ScriptReadSolution(scr, new_state))
            {
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "mesh")
         {
            if (ScriptReadDisplMesh(scr, new_state))
            {
               done_one_command = 1;
               continue;
            }
            if (new_state.mesh == NULL)
            {
               cout << "Script: unexpected 'mesh' command!" << endl;
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "psolution")
         {
            if (ScriptReadParSolution(scr, new_state))
            {
               done_one_command = 1;
               continue;
            }
         }

         if (stream_state.SetNewMeshAndSolution(std::move(new_state), vs))
         {
            MyExpose();
         }
         else
         {
            cout << "Different type of mesh / solution." << endl;
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
            vs->SetLight(false);
         }
         else if (word == "on")
         {
            vs->SetLight(true);
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
         scr >> stream_state.keys;
         cout << "Script: keys: '" << stream_state.keys << "'" << endl;
         // SendKeySequence(keys.c_str());
         CallKeySequence(stream_state.keys.c_str());
         MyExpose();
      }
      else if (word == "palette")
      {
         int pal;
         scr >> pal;
         cout << "Script: palette: " << pal << endl;
         vs->palette.SetIndex(pal-1);
         MyExpose();
      }
      else if (word == "palette_repeat")
      {
         int rpt_times;
         scr >> rpt_times;
         cout << "Script: palette_repeat: " << rpt_times << endl;
         vs->palette.SetRepeatTimes(rpt_times);
         vs->palette.Init();
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
            scr >> vs->rotmat[i/4][i%4];
            cout << ' ' << vs->rotmat[i/4][i%4];
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
         vs->PrepareCaption(); // turn on or off the caption
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
         if (ScriptReadSolution(scr, stream_state))
         {
            return;
         }

         // start the visualization
         break;
      }
      else if (word == "psolution")
      {
         if (ScriptReadParSolution(scr, stream_state))
         {
            return;
         }

         // start the visualization
         break;
      }
      else if (word == "mesh")
      {
         if (ScriptReadDisplMesh(scr, stream_state))
         {
            return;
         }
         if (stream_state.mesh)
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
   script = &scr;
   stream_state.keys.clear();

   std::thread worker_thread
   {
      [&](StreamState local_state)
      {
         // set the thread-local StreamState
         stream_state = std::move(local_state);
         if (c_plot_caption != string_none)
         {
            plot_caption = c_plot_caption;
         }
         if (GLVisInitVis((stream_state.grid_f->VectorDim() == 1) ? 0 : 1, {}))
         {
            GetAppWindow()->setOnKeyDown(SDLK_SPACE, ScriptControl);
            GLVisStartVis();
         }
      },
      std::move(stream_state)
   };

   SDLMainLoop();
   worker_thread.join();

   delete init_nodes; init_nodes = NULL;

   cout << "Script: min_val = " << scr_min_val
        << ", max_val = " << scr_max_val << endl;

   script = NULL;
}

struct Session
{
   StreamCollection input_streams;
   StreamState state;
   int ft = -1;
   std::thread handler;

   Session(bool fix_elem_orient,
           bool save_coloring)
   {
      state.fix_elem_orient = fix_elem_orient;
      state.save_coloring = save_coloring;
   }

   Session(int other_ft, StreamState other_state)
      : state(std::move(other_state))
      , ft(other_ft)
   { }

   ~Session() = default;

   Session(Session&& from) = default;
   Session& operator= (Session&& from) = default;

   void StartSession()
   {
      auto funcThread =
         [](StreamState thread_state, int ftype, StreamCollection is)
      {
         // Set thread-local stream state
         stream_state = std::move(thread_state);
         if (c_plot_caption != string_none)
         {
            plot_caption = c_plot_caption;
         }

         if (GLVisInitVis(ftype, std::move(is)))
         {
            GLVisStartVis();
         }
      };
      handler = std::thread {funcThread,
                             std::move(state), ft, std::move(input_streams)};
      handler.detach();
   }

   bool StartSavedSession(std::string stream_file)
   {
      unique_ptr<ifstream> ifs(new ifstream(stream_file));
      if (!(*ifs))
      {
         cout << "Can not open stream file: " << stream_file << endl;
         return false;
      }
      string data_type;
      *ifs >> data_type >> ws;
      ft = state.ReadStream(*ifs, data_type);
      input_streams.emplace_back(std::move(ifs));

      StartSession();
      return true;
   }
};

void GLVisServer(int portnum, bool save_stream, bool fix_elem_orient,
                 bool save_coloring)
{
   std::vector<Session> current_sessions;
   string data_type;
   int viscount = 0;
   unsigned int nproc = 1, proc = 0;

#ifdef MFEM_USE_GNUTLS
   unique_ptr<GnuTLS_global_state> state;
   unique_ptr<GnuTLS_session_params> params;
   if (secure)
   {
      state.reset(new GnuTLS_global_state);
      // state->set_log_level(1000);
      string home_dir(getenv("HOME"));
      string server_dir = home_dir + "/.config/glvis/server/";
#ifndef MFEM_USE_GNUTLS_X509
      string pubkey  = server_dir + "pubring.gpg";
      string privkey = server_dir + "secring.gpg";
      string trustedkeys = server_dir + "trusted-clients.gpg";
#else
      string pubkey  = server_dir + "cert.pem";
      string privkey = server_dir + "key.pem";
      string trustedkeys = server_dir + "trusted-clients.pem";
#endif
      params.reset(new GnuTLS_session_params(
                      *state, pubkey.c_str(), privkey.c_str(),
                      trustedkeys.c_str(), GNUTLS_SERVER));
      if (!params->status.good())
      {
         cout << "  public key   = " << pubkey << '\n'
              << "  private key  = " << privkey << '\n'
              << "  trusted keys = " << trustedkeys << endl;
         cout << "Error setting GLVis server parameters.\n"
              "Generate your GLVis keys with:"
              " bash glvis-keygen.sh [\"Your Name\"] [\"Your Email\"]"
              << endl;
         return;
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
      exit(2);
   }
   while (1)
   {

      unique_ptr<socketstream> isock;
#ifndef MFEM_USE_GNUTLS
      isock.reset(new socketstream);
#else
      isock.reset(secure ? new socketstream(*params) : new socketstream(false));
#endif
      vector<unique_ptr<istream>> input_streams;
      while (server.accept(*isock) < 0)
      {
#ifdef GLVIS_DEBUG
         cout << "GLVis: server.accept(...) failed." << endl;
#endif
      }

      *isock >> data_type >> ws;

      if (save_stream)
      {
         viscount++;
      }

      int par_data = 0;
      if (data_type == "parallel")
      {
         par_data = 1;
         unsigned int np = 0;
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
               input_streams.resize(nproc);
            }
            else
            {
               if (nproc != input_streams.size())
               {
                  cout << "Unexpected number of processors: " << nproc
                       << ", expected: " << input_streams.size() << endl;
                  mfem_error();
               }
            }
            if (proc >= nproc)
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

            input_streams[proc] = std::move(isock);
#ifndef MFEM_USE_GNUTLS
            isock.reset(new socketstream);
#else
            isock.reset(secure ? new socketstream(*params) :
                        new socketstream(false));
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

      Session new_session(fix_elem_orient, save_coloring);

      char tmp_file[50];
      if (save_stream)
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
            ReadInputStreams(new_session.state, input_streams);
            ofs.precision(8);
            ofs << "solution\n";
            new_session.state.mesh->Print(ofs);
            new_session.state.grid_f->Save(ofs);
         }
         ofs.close();
         cout << "Data saved in " << tmp_file << endl;

         new_session.StartSavedSession(tmp_file);
      }
      else
      {
         if (!par_data)
         {
            new_session.ft = new_session.state.ReadStream(*isock, data_type);
            input_streams.emplace_back(std::move(isock));
         }
         else
         {
            new_session.ft = ReadInputStreams(new_session.state, input_streams);
         }
         // Pass ownership of input streams into session object
         new_session.input_streams = std::move(input_streams);
         new_session.StartSession();
      }
      current_sessions.emplace_back(std::move(new_session));
   }
}

int main (int argc, char *argv[])
{
#ifdef _WIN32
   // Call needed to avoid SDL_Init failure when not substituting main() for
   // SDL_main().
   SDL_SetMainReady();
#endif
   // variables for command line arguments
   int         np            = 0;
   bool        save_stream   = false;
   const char *stream_file   = string_none;
   const char *script_file   = string_none;
   const char *font_name     = string_default;
   int         portnum       = 19916;
   int         multisample   = GetMultisample();
   double      line_width    = 1.0;
   double      ms_line_width = gl3::LINE_WIDTH_AA;
   int         geom_ref_type = Quadrature1D::ClosedUniform;
   bool        legacy_gl_ctx = false;
   bool        enable_hidpi  = true;

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
   args.AddOption(&stream_state.fix_elem_orient, "-fo", "--fix-orientations",
                  "-no-fo", "--dont-fix-orientations",
                  "Attempt to fix the orientations of inverted elements.");
   args.AddOption(&keep_attr, "-a", "--real-attributes",
                  "-ap", "--processor-attributes",
                  "When opening a parallel mesh, use the real mesh attributes"
                  " or replace them with the processor rank.");
   args.AddOption(&geom_ref_type, "-grt", "--geometry-refiner-type",
                  "Set of points to use when refining geometry:"
                  " 3 = uniform, 1 = Gauss-Lobatto, (see mfem::Quadrature1D).");
   args.AddOption(&stream_state.save_coloring, "-sc", "--save-coloring",
                  "-no-sc", "--dont-save-coloring",
                  "Save the mesh coloring generated when opening only a mesh.");
   args.AddOption(&portnum, "-p", "--listen-port",
                  "Specify the port number on which to accept connections.");
   args.AddOption(&secure, "-sec", "--secure-sockets",
                  "-no-sec", "--standard-sockets",
                  "Enable or disable GnuTLS secure sockets.");
   args.AddOption(&save_stream, "-save", "--save-stream",
                  "-no-save", "--dont-save-stream",
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
                  "Set the font: [<font-name>[:style=<style>]][-<font-size>],"
                  " e.g. -fn \"Helvetica:style=Bold-16\".");
   args.AddOption(&multisample, "-ms", "--multisample",
                  "Set the multisampling mode (toggled with the 'A' key).");
   args.AddOption(&line_width, "-lw", "--line-width",
                  "Set the line width (multisampling off).");
   args.AddOption(&ms_line_width, "-mslw", "--multisample-line-width",
                  "Set the line width (multisampling on).");
   args.AddOption(&legacy_gl_ctx, "-oldgl", "--legacy-gl",
                  "-anygl", "--any-gl",
                  "Only try to create a legacy OpenGL (< 2.1) context.");
   args.AddOption(&enable_hidpi, "-hidpi", "--high-dpi",
                  "-nohidpi", "--no-high-dpi",
                  "Enable/disable support for HiDPI at runtime, if supported.");

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
      stream_state.is_gf = 255;
   }
   if (np > 0)
   {
      input |= 256;
   }
   if (arg_keys != string_none)
   {
      stream_state.keys = arg_keys;
   }
   if (font_name != string_default)
   {
      SetFont(font_name);
   }
   if (multisample != GetMultisample())
   {
      SetMultisample(multisample);
   }
   if (line_width != GetLineWidth())
   {
      SetLineWidth(line_width);
   }
   if (ms_line_width != GetLineWidthMS())
   {
      SetLineWidthMS(ms_line_width);
   }
   if (c_plot_caption != string_none)
   {
      plot_caption = c_plot_caption;
   }
   if (legacy_gl_ctx == true)
   {
      SetLegacyGLOnly(legacy_gl_ctx);
   }
   SetUseHiDPI(enable_hidpi);

   GLVisGeometryRefiner.SetType(geom_ref_type);

   string data_type;

   // check for saved stream file
   if (stream_file != string_none)
   {
      Session stream_session(stream_state.fix_elem_orient,
                             stream_state.save_coloring);

      if (!stream_session.StartSavedSession(stream_file))
      {
         return 1;
      }

      SDLMainLoop();
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
         (stream_state.is_gf && (input == 3 || input == 259))))
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

   // server mode, read the mesh and the solution from a socket
   if (input == 1)
   {
      // Run server in new thread
      std::thread serverThread{GLVisServer, portnum, save_stream,
                               stream_state.fix_elem_orient,
                               stream_state.save_coloring};

      // Start SDL in main thread
      SDLMainLoop(true);

      serverThread.detach();
   }
   else  // input != 1, non-server mode
   {
      if (input & 256)
      {
         ReadParallel(np, stream_state);
      }
      else
      {
         ReadSerial(stream_state);
      }

      bool use_vector_soln = (input & 8);
      bool use_soln = (input & 4);
      int field_type;
      if (use_vector_soln)
      {
         field_type = 1;
      }
      else
      {
         field_type = (use_soln) ? 0 : 2;
      }
      Session single_session(field_type, std::move(stream_state));
      single_session.StartSession();

      SDLMainLoop();
   }

   cout << "Thank you for using GLVis." << endl;

   return 0;
}


void PrintSampleUsage(ostream &os)
{
   os <<
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


void ReadSerial(StreamState& state)
{
   // get the mesh from a file
   named_ifgzstream meshin(mesh_file);
   if (!meshin)
   {
      cerr << "Can not open mesh file " << mesh_file << ". Exit.\n";
      exit(1);
   }

   state.mesh.reset(new Mesh(meshin, 1, 0, state.fix_elem_orient));

   if (state.is_gf || (input & 4) || (input & 8))
   {
      // get the solution from file
      bool freesolin = false;
      ifgzstream *solin = NULL;
      if (!strcmp(mesh_file,sol_file))
      {
         solin = &meshin;
      }
      else
      {
         solin = new ifgzstream(sol_file);
         freesolin = true;
         if (!(*solin))
         {
            cerr << "Can not open solution file " << sol_file << ". Exit.\n";
            exit(1);
         }
      }

      if (state.is_gf)
      {
         state.grid_f.reset(new GridFunction(state.mesh.get(), *solin));
         SetGridFunction(state);
      }
      else if (input & 4)
      {
         // get rid of NetGen's info line
         char buff[128];
         solin->getline(buff,128);
         state.sol.Load(*solin, state.mesh->GetNV());
      }
      else if (input & 8)
      {
         state.solu.Load(*solin, state.mesh->GetNV());
         state.solv.Load(*solin, state.mesh->GetNV());
         if (state.mesh->SpaceDimension() == 3)
         {
            state.solw.Load(*solin, state.mesh->GetNV());
         }
      }
      if (freesolin)
      {
         delete solin;
      }
   }
   else
   {
      state.SetMeshSolution();
   }

   state.Extrude1DMeshAndSolution();
}


void SetGridFunction(StreamState& state)
{
   if (gf_component != -1)
   {
      if (gf_component < 0 || gf_component >= state.grid_f->VectorDim())
      {
         cerr << "Invalid component " << gf_component << '.' << endl;
         exit(1);
      }
      FiniteElementSpace *ofes = state.grid_f->FESpace();
      FiniteElementCollection *fec =
         FiniteElementCollection::New(ofes->FEColl()->Name());
      FiniteElementSpace *fes = new FiniteElementSpace(state.mesh.get(), fec);
      GridFunction *new_gf = new GridFunction(fes);
      new_gf->MakeOwner(fec);
      for (int i = 0; i < new_gf->Size(); i++)
      {
         (*new_gf)(i) = (*state.grid_f)(ofes->DofToVDof(i, gf_component));
      }
      state.grid_f.reset(new_gf);
   }
   if (state.grid_f->VectorDim() == 1)
   {
      state.grid_f->GetNodalValues(state.sol);
      input |= 4;
   }
   else
   {
      input |= 8;
   }
}


void ReadParallel(int np, StreamState& state)
{
   int read_err;

   if (state.is_gf)
   {
      read_err = ReadParMeshAndGridFunction(np, mesh_file, sol_file,
                                            state, keep_attr);
      if (!read_err)
      {
         SetGridFunction(state);
      }
   }
   else
   {
      read_err = ReadParMeshAndGridFunction(np, mesh_file, NULL,
                                            state, keep_attr);
      if (!read_err)
      {
         state.SetMeshSolution();
      }
   }

   if (read_err)
   {
      exit(1);
   }

   state.Extrude1DMeshAndSolution();
}

int ReadParMeshAndGridFunction(int np, const char *mesh_prefix,
                               const char *sol_prefix,
                               StreamState& state, int keep_attribs)
{
   state.mesh = NULL;

   // are the solutions bundled together with the mesh files?
   bool same_file = false;
   if (sol_prefix)
   {
      same_file = !strcmp(sol_prefix, mesh_prefix);
      state.grid_f = NULL;
   }

   Array<Mesh *> mesh_array(np);
   Array<GridFunction *> gf_array(np);
   mesh_array = NULL;
   gf_array = NULL;

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

      mesh_array[p] = new Mesh(meshfile, 1, 0, state.fix_elem_orient);

      if (!keep_attribs)
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

            gf_array[p] = new GridFunction(mesh_array[p], solfile);
         }
         else  // mesh and solution in the same file
         {
            gf_array[p] = new GridFunction(mesh_array[p], meshfile);
         }
      }
   }

   if (!read_err)
   {
      // create the combined mesh and gf
      state.mesh.reset(new Mesh(mesh_array, np));
      if (sol_prefix)
      {
         state.grid_f.reset(new GridFunction(state.mesh.get(), gf_array, np));
      }
   }

   for (int p = 0; p < np; p++)
   {
      delete gf_array[np-1-p];
      delete mesh_array[np-1-p];
   }

   return read_err;
}

int ReadInputStreams(StreamState& state, const StreamCollection& input_streams)
{
   int nproc = input_streams.size();
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
#ifdef GLVIS_DEBUG
      cout << " type " << data_type << " ... " << flush;
#endif
      mesh_array[p] = new Mesh(isock, 1, 0, state.fix_elem_orient);
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

   state.mesh.reset(new Mesh(mesh_array, nproc));
   if (gf_count == 0)
   {
      state.SetMeshSolution();
      field_type = 2;
   }
   else
   {
      state.grid_f.reset(new GridFunction(state.mesh.get(), gf_array, nproc));
      field_type = (state.grid_f->VectorDim() == 1) ? 0 : 1;
   }

   for (int p = 0; p < nproc; p++)
   {
      delete mesh_array[nproc-1-p];
      delete gf_array[nproc-1-p];
   }

   state.Extrude1DMeshAndSolution();

   return field_type;
}
