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
#include "coll_reader.hpp"
#include "stream_reader.hpp"
#include "visual.hpp"

using namespace std;
using namespace mfem;

extern const char *string_none;
extern const char *string_default;

extern string      dc_protocol;
extern int         dc_cycle;
extern int         window_x;
extern int         window_y;
extern int         window_w;
extern int         window_h;
extern const char *c_plot_caption;

extern thread_local DataState stream_state;
extern thread_local VisualizationSceneScalarData *vs;

istream *script = NULL;
int scr_running = 0;
int scr_level = 0;
Vector *init_nodes = NULL;
double scr_min_val, scr_max_val;

bool GLVisInitVis(StreamCollection input_streams);
void GLVisStartVis();

int ScriptReadSolution(istream &scr, DataState& state)
{
   int err_read;
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
   state.SetMesh(new Mesh(imesh, 1, 0, state.fix_elem_orient));

   // read the solution (GridFunction)
   scr >> ws >> sword;
   cout << "solution: " << sword << endl;

   FileReader reader(state);
   err_read = reader.ReadSerial(FileReader::FileType::GRID_FUNC, mword.c_str(),
                                sword.c_str());

   return err_read;
}

int ScriptReadQuadrature(istream &scr, DataState& state)
{
   int err_read;
   string mword,sword;

   cout << "Script: quadrature: " << flush;
   // read the mesh
   scr >> ws >> mword; // mesh filename (can't contain spaces)
   cout << "mesh: " << mword << "; " << flush;
   named_ifgzstream imesh(mword.c_str());
   if (!imesh)
   {
      cout << "Can not open mesh file: " << mword << endl;
      return 1;
   }
   state.SetMesh(new Mesh(imesh, 1, 0, state.fix_elem_orient));

   // read the quadrature (QuadratureFunction)
   scr >> ws >> sword;
   cout << "quadrature: " << sword << endl;

   FileReader reader(state);
   err_read = reader.ReadSerial(FileReader::FileType::QUAD_FUNC, mword.c_str(),
                                sword.c_str());

   return err_read;
}

int ScriptReadParSolution(istream &scr, DataState& state)
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

   FileReader reader(state);
   err_read = reader.ReadParallel(np, FileReader::FileType::GRID_FUNC,
                                  mesh_prefix.c_str(), sol_prefix.c_str());
   return err_read;
}

int ScriptReadParQuadrature(istream &scr, DataState& state)
{
   int np, scr_keep_attr, err_read;
   string mesh_prefix, quad_prefix;

   cout << "Script: pquadrature: " << flush;
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
   // read the quadrature prefix
   scr >> ws >> quad_prefix;
   cout << "quadrature prefix: " << quad_prefix << endl;

   FileReader reader(state);
   err_read = reader.ReadParallel(np, FileReader::FileType::QUAD_FUNC,
                                  mesh_prefix.c_str(), quad_prefix.c_str());
   return err_read;
}

int ScriptReadDisplMesh(istream &scr, DataState& state)
{
   DataState meshstate;
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
      meshstate.SetMesh(new Mesh(imesh, 1, 0, state.fix_elem_orient));
   }
   meshstate.ExtrudeMeshAndSolution();
   Mesh* const m = meshstate.mesh.get();
   if (init_nodes == NULL)
   {
      init_nodes = new Vector;
      meshstate.mesh->GetNodes(*init_nodes);
      state.SetMesh(NULL);
      state.SetGridFunction(NULL);
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

      meshstate.SetGridFunction(new GridFunction(vfes));
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

      state = std::move(meshstate);
   }

   return 0;
}

int ScriptReadDataColl(istream &scr, DataState &state, bool mesh_only = true,
                       bool quad = false)
{
   int err_read;
   int type;
   string cword, fword;

   cout << "Script: data_collection: " << flush;
   // read the collection
   scr >> ws >> type; // collection type
   cout << "type: " << type << "; " << flush;
   scr >> ws >> cword; // collection filename (can't contain spaces)
   cout << "collection: " << cword << "; " << flush;

   if (!mesh_only)
   {
      // read the field
      scr >> ws >> fword;
      cout << "field: " << fword << endl;
   }

   DataCollectionReader reader(state);
   if (dc_protocol != string_default)
   {
      reader.SetProtocol(dc_protocol.c_str());
   }

   if (mesh_only)
      err_read = reader.ReadSerial((DataCollectionReader::CollType)type,
                                   cword.c_str(), dc_cycle);
   else
      err_read = reader.ReadSerial((DataCollectionReader::CollType)type,
                                   cword.c_str(), dc_cycle, fword.c_str(), quad);

   return err_read;
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
      else if (word == "data_coll_cycle")
      {
         scr >> dc_cycle;
      }
      else if (word == "data_coll_protocol")
      {
         scr >> dc_protocol;
      }
      else if (word == "solution" || word == "mesh" || word == "psolution"
               || word == "quadrature" || word == "pquadrature" || word == "data_coll_mesh"
               || word == "data_coll_field" || word == "data_coll_quad")
      {
         DataState new_state;

         if (word == "solution")
         {
            if (ScriptReadSolution(scr, new_state))
            {
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "quadrature")
         {
            if (ScriptReadQuadrature(scr, new_state))
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
         else if (word == "pquadrature")
         {
            if (ScriptReadParQuadrature(scr, new_state))
            {
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "data_coll_mesh")
         {
            if (ScriptReadDataColl(scr, new_state))
            {
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "data_coll_field")
         {
            if (ScriptReadDataColl(scr, new_state, false))
            {
               done_one_command = 1;
               continue;
            }
         }
         else if (word == "data_coll_quad")
         {
            if (ScriptReadDataColl(scr, new_state, false, true))
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
         VisualizationSceneScalarData::Shading s =
            VisualizationSceneScalarData::Shading::Invalid;
         if (word == "flat")
         {
            s = VisualizationSceneScalarData::Shading::Flat;
         }
         else if (word == "smooth")
         {
            s = VisualizationSceneScalarData::Shading::Smooth;
         }
         else if (word == "cool")
         {
            s = VisualizationSceneScalarData::Shading::Noncomforming;
         }
         if (s != VisualizationSceneScalarData::Shading::Invalid)
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
      else if (word == "levellines")
      {
         double min, max;
         int num;
         scr >> min >> max >> num;
         cout << "Script: levellines: " << flush;
         vs->SetLevelLines(min, max, num);
         vs->UpdateLevelLines();
         cout << min << ' ' << max << ' ' << num << endl;
         MyExpose();
      }
      else if (word == "axis_numberformat")
      {
         char delim;
         string axis_formatting;
         scr >> ws >> delim;
         getline(scr, axis_formatting, delim);
         cout << "Script: axis_numberformat: " << flush;
         vs->SetAxisNumberFormat(axis_formatting);
         cout << axis_formatting << endl;
         MyExpose();
      }
      else if (word == "colorbar_numberformat")
      {
         char delim;
         string colorbar_formatting;
         scr >> ws >> delim;
         getline(scr, colorbar_formatting, delim);
         cout << "Script: colorbar_numberformat: " << flush;
         vs->SetColorbarNumberFormat(colorbar_formatting);
         cout << colorbar_formatting << endl;
         MyExpose();
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
         vs->palette.GenerateTextures();
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
      else if (word == "data_coll_cycle")
      {
         scr >> dc_cycle;
      }
      else if (word == "data_coll_protocol")
      {
         scr >> dc_protocol;
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
      else if (word == "quadrature")
      {
         if (ScriptReadQuadrature(scr, stream_state))
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
      else if (word == "pquadrature")
      {
         if (ScriptReadParQuadrature(scr, stream_state))
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
      else if (word == "data_coll_mesh")
      {
         if (ScriptReadDataColl(scr, stream_state))
         {
            return;
         }

         // start the visualization
         break;
      }
      else if (word == "data_coll_field")
      {
         if (ScriptReadDataColl(scr, stream_state, false))
         {
            return;
         }

         // start the visualization
         break;
      }
      else if (word == "data_coll_quad")
      {
         if (ScriptReadDataColl(scr, stream_state, false, true))
         {
            return;
         }

         // start the visualization
         break;
      }
      else
      {
         cout << "Unknown command in script: " << word << endl;
      }
   }

   scr_level = scr_running = 0;
   script = &scr;
   stream_state.keys.clear();

   // Make sure the singleton object returned by GetMainThread() is
   // initialized from the main thread.
   GetMainThread();

   std::thread worker_thread
   {
      [&](DataState local_state)
      {
         // set the thread-local DataState
         stream_state = std::move(local_state);
         if (c_plot_caption != string_none)
         {
            plot_caption = c_plot_caption;
         }
         if (GLVisInitVis({}))
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
