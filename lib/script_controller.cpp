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

#include "script_controller.hpp"
#include "file_reader.hpp"
#include "coll_reader.hpp"
#include "stream_reader.hpp"
#include "visual.hpp"

#include <vector>
#include <algorithm>

using namespace std;
using namespace mfem;

enum class Command
{
   Mesh,
   Solution,
   ParSolution,
   Quadrature,
   ParQuadrature,
   DataCollMesh,
   DataCollField,
   DataCollQuad,
   DataCollCycle,
   DataCollProto,
   Screenshot,
   Viewcenter,
   Perspective,
   Light,
   View,
   Zoom,
   Shading,
   Subdivisions,
   Valuerange,
   Autoscale,
   Levellines,
   AxisNumberFormat,
   ColorbarNumberFormat,
   Window,
   Keys,
   Palette,
   PaletteRepeat,
   ToggleAttributes,
   Rotmat,
   Camera,
   Scale,
   Translate,
   PlotCaption,
   //----------
   Max
};

struct CmdItem
{
   const char *keyword;
   const char *params;
   const char *desc;

   bool operator==(const string &key) const { return key == keyword; }
};

static vector<CmdItem> commands;

int ScriptController::ScriptReadSolution(istream &scr, DataState &state)
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

int ScriptController::ScriptReadQuadrature(istream &scr, DataState &state)
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

int ScriptController::ScriptReadParSolution(istream &scr, DataState &state)
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
   state.keep_attr = scr_keep_attr;

   // read the solution prefix
   scr >> ws >> sol_prefix;
   cout << "solution prefix: " << sol_prefix << endl;

   FileReader reader(state);
   err_read = reader.ReadParallel(np, FileReader::FileType::GRID_FUNC,
                                  mesh_prefix.c_str(), sol_prefix.c_str());
   return err_read;
}

int ScriptController::ScriptReadParQuadrature(istream &scr, DataState &state)
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
   state.keep_attr = scr_keep_attr;

   // read the quadrature prefix
   scr >> ws >> quad_prefix;
   cout << "quadrature prefix: " << quad_prefix << endl;

   FileReader reader(state);
   err_read = reader.ReadParallel(np, FileReader::FileType::QUAD_FUNC,
                                  mesh_prefix.c_str(), quad_prefix.c_str());
   return err_read;
}

int ScriptController::ScriptReadDisplMesh(istream &scr, DataState &state)
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

int ScriptController::ScriptReadDataColl(istream &scr, DataState &state,
                                         bool mesh_only, bool quad)
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

void ScriptController::InitCommands()
{
   commands.resize((size_t)Command::Max);

   commands[(size_t)Command::Mesh]                 = {"mesh", "<file>", "Visualize the mesh."};
   commands[(size_t)Command::Solution]             = {"solution", "<mesh> <solution>", "Visualize the solution."};
   commands[(size_t)Command::ParSolution]          = {"psolution", "<np> <mesh prefix> <keep attributes> <solution prefix>", "Visualize the distributed solution."};
   commands[(size_t)Command::Quadrature]           = {"quadrature", "<mesh> <quadrature>", "Visualize the quadrature."};
   commands[(size_t)Command::ParQuadrature]        = {"pquadrature", "<np> <mesh prefix> <keep attributes> <quadrature prefix>", "Visualize the distributed quadrature."};
   commands[(size_t)Command::DataCollMesh]         = {"data_coll_mesh", "<type> <data coll>", "Visualize the mesh from data collection."};
   commands[(size_t)Command::DataCollField]        = {"data_coll_field", "<type> <data coll> <field>", "Visualize the field from data collection."};
   commands[(size_t)Command::DataCollQuad]         = {"data_coll_quad", "<type> <data coll> <quad>", "Visualize the Q-field from data collection."};
   commands[(size_t)Command::DataCollCycle]        = {"data_coll_cycle", "<cycle>", "Preset the cycle of the data collection."};
   commands[(size_t)Command::DataCollProto]        = {"data_coll_protocol", "<protocol>", "Preset the protocol of the data collection."};
   commands[(size_t)Command::Screenshot]           = {"screenshot", "<file>", "Take a screenshot, saving it to the file."};
   commands[(size_t)Command::Viewcenter]           = {"viewcenter", "<x> <y>", "Change the viewcenter."};
   commands[(size_t)Command::Perspective]          = {"perspective", "<on/off>", "Turn on or off perspective projection."};
   commands[(size_t)Command::Light]                = {"light", "<on/off>", "Turn on or off light."};
   commands[(size_t)Command::View]                 = {"view", "<theta> <phi>", "Change the solid angle of view."};
   commands[(size_t)Command::Zoom]                 = {"zoom", "<zoom>", "Change the zoom factor."};
   commands[(size_t)Command::Shading]              = {"shading", "<flat/smooth/cool>", "Change the shading algorithm."};
   commands[(size_t)Command::Subdivisions]         = {"subdivisions", "<times> <dummy>", "Change the refinement level."};
   commands[(size_t)Command::Valuerange]           = {"valuerange", "<min> <max>", "Change the value range."};
   commands[(size_t)Command::Autoscale]            = {"autoscale", "<off/on/value/mesh>", "Change the autoscale algorithm."};
   commands[(size_t)Command::Levellines]           = {"levellines", "<min> <max> <num>", "Set the level lines."};
   commands[(size_t)Command::AxisNumberFormat]     = {"axis_numberformat", "'<format>'", "Set the axis number format."};
   commands[(size_t)Command::ColorbarNumberFormat] = {"colorbar_numberformat", "'<format>'", "Set the colorbar number format."};
   commands[(size_t)Command::Window]               = {"window", "<x> <y> <w> <h>", "Set the position and size of the window."};
   commands[(size_t)Command::Keys]                 = {"keys", "<keys>", "Send the control key sequence."};
   commands[(size_t)Command::Palette]              = {"palette", "<index>", "Set the palette index."};
   commands[(size_t)Command::PaletteRepeat]        = {"palette_repeat", "<times>", "Set the repetition of the palette."};
   commands[(size_t)Command::ToggleAttributes]     = {"toggle_attributes", "<1/0> [[<1/0>] ...];", "Toggle visibility of the attributes."};
   commands[(size_t)Command::Rotmat]               = {"rotmat", "<[0,0]> <[1,0]> ... <[3,3]>", "Set the rotation matrix."};
   commands[(size_t)Command::Camera]               = {"camera", "<cam[0]> ... <cam[2]> <dir[0]> ... <dir[2]> <up[0]> ... <up[2]>", "Set the camera position, direction and upward vector."};
   commands[(size_t)Command::Scale]                = {"scale", "<scale>", "Set the scaling factor."};
   commands[(size_t)Command::Translate]            = {"translate", "<x> <y> <z>", "Set the translation coordinates."};
   commands[(size_t)Command::PlotCaption]          = {"plot_caption", "'<caption>'", "Set the plot caption."};
}

void ScriptController::PrintCommands()
{
   cout << "Available commands are:" << endl;

   for (const CmdItem &ci : commands)
   {
      cout << "\t" << ci.keyword << " " << ci.params << " - " << ci.desc << endl;
   }
}

void ScriptController::ExecuteScriptCommand()
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
         continue;
      }
      else if (word == "}")
      {
         scr_level--;
         if (scr_level < 0)
         {
            scr_level = 0;
         }
         continue;
      }

      auto it = find(commands.begin(), commands.end(), word);
      if (it == commands.end())
      {
         cout << "Unknown command in script: " << word << endl;
         PrintCommands();
         break;
      }

      const Command cmd = (Command)(it - commands.begin());
      switch (cmd)
      {
         case Command::Mesh:
         case Command::Solution:
         case Command::ParSolution:
         case Command::Quadrature:
         case Command::ParQuadrature:
         case Command::DataCollMesh:
         case Command::DataCollField:
         case Command::DataCollQuad:
         {
            DataState new_state;

            switch (cmd)
            {
               case Command::Solution:
                  if (ScriptReadSolution(scr, new_state))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               case Command::Quadrature:
                  if (ScriptReadQuadrature(scr, new_state))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               case Command::Mesh:
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
                  break;
               case Command::ParSolution:
                  if (ScriptReadParSolution(scr, new_state))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               case Command::ParQuadrature:
                  if (ScriptReadParQuadrature(scr, new_state))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               case Command::DataCollMesh:
                  if (ScriptReadDataColl(scr, new_state))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               case Command::DataCollField:
                  if (ScriptReadDataColl(scr, new_state, false))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               case Command::DataCollQuad:
                  if (ScriptReadDataColl(scr, new_state, false, true))
                  {
                     done_one_command = 1;
                     continue;
                  }
                  break;
               default:
                  break;
            }

            if (win.SetNewMeshAndSolution(std::move(new_state)))
            {
               MyExpose();
            }
            else
            {
               cout << "Different type of mesh / solution." << endl;
            }
         }
         break;
         case Command::DataCollCycle:
            scr >> dc_cycle;
            break;
         case Command::DataCollProto:
            scr >> dc_protocol;
            break;
         case Command::Screenshot:
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

            if (scr_min_val > win.vs->GetMinV())
            {
               scr_min_val = win.vs->GetMinV();
            }
            if (scr_max_val < win.vs->GetMaxV())
            {
               scr_max_val = win.vs->GetMaxV();
            }
         }
         break;
         case Command::Viewcenter:
         {
            scr >> win.vs->ViewCenterX >> win.vs->ViewCenterY;
            cout << "Script: viewcenter: "
                 << win.vs->ViewCenterX << ' ' << win.vs->ViewCenterY << endl;
            MyExpose();
         }
         break;
         case Command::Perspective:
         {
            scr >> ws >> word;
            cout << "Script: perspective: " << word;
            if (word == "off")
            {
               win.vs->OrthogonalProjection = 1;
            }
            else if (word == "on")
            {
               win.vs->OrthogonalProjection = 0;
            }
            else
            {
               cout << '?';
            }
            cout << endl;
            MyExpose();
         }
         break;
         case Command::Light:
         {
            scr >> ws >> word;
            cout << "Script: light: " << word;
            if (word == "off")
            {
               win.vs->SetLight(false);
            }
            else if (word == "on")
            {
               win.vs->SetLight(true);
            }
            else
            {
               cout << '?';
            }
            cout << endl;
            MyExpose();
         }
         break;
         case Command::View:
         {
            double theta, phi;
            scr >> theta >> phi;
            cout << "Script: view: " << theta << ' ' << phi << endl;
            win.vs->SetView(theta, phi);
            MyExpose();
         }
         break;
         case Command::Zoom:
         {
            double factor;
            scr >> factor;
            cout << "Script: zoom: " << factor << endl;
            win.vs->Zoom(factor);
            MyExpose();
         }
         break;
         case Command::Shading:
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
               win.vs->SetShading(s, false);
               cout << word << endl;
               MyExpose();
            }
            else
            {
               cout << word << " ?" << endl;
            }
         }
         break;
         case Command::Subdivisions:
         {
            int t, b;
            scr >> t >> b;
            cout << "Script: subdivisions: " << flush;
            win.vs->SetRefineFactors(t, b);
            cout << t << ' ' << b << endl;
            MyExpose();
         }
         break;
         case Command::Valuerange:
         {
            double min, max;
            scr >> min >> max;
            cout << "Script: valuerange: " << flush;
            win.vs->SetValueRange(min, max);
            cout << min << ' ' << max << endl;
            MyExpose();
         }
         break;
         case Command::Autoscale:
         {
            scr >> ws >> word;
            cout << "Script: autoscale: " << word;
            if (word == "off")
            {
               win.vs->SetAutoscale(0);
            }
            else if (word == "on")
            {
               win.vs->SetAutoscale(1);
            }
            else if (word == "value")
            {
               win.vs->SetAutoscale(2);
            }
            else if (word == "mesh")
            {
               win.vs->SetAutoscale(3);
            }
            else
            {
               cout << '?';
            }
            cout << endl;
         }
         break;
         case Command::Levellines:
         {
            double min, max;
            int num;
            scr >> min >> max >> num;
            cout << "Script: levellines: " << flush;
            win.vs->SetLevelLines(min, max, num);
            win.vs->UpdateLevelLines();
            cout << min << ' ' << max << ' ' << num << endl;
            MyExpose();
         }
         break;
         case Command::AxisNumberFormat:
         {
            char delim;
            string axis_formatting;
            scr >> ws >> delim;
            getline(scr, axis_formatting, delim);
            cout << "Script: axis_numberformat: " << flush;
            win.vs->SetAxisNumberFormat(axis_formatting);
            cout << axis_formatting << endl;
            MyExpose();
         }
         break;
         case Command::ColorbarNumberFormat:
         {
            char delim;
            string colorbar_formatting;
            scr >> ws >> delim;
            getline(scr, colorbar_formatting, delim);
            cout << "Script: colorbar_numberformat: " << flush;
            win.vs->SetColorbarNumberFormat(colorbar_formatting);
            cout << colorbar_formatting << endl;
            MyExpose();
         }
         break;
         case Command::Window:
         {
            scr >> win.window_x >> win.window_y >> win.window_w >> win.window_h;
            cout << "Script: window: " << win.window_x << ' ' << win.window_y
                 << ' ' << win.window_w << ' ' << win.window_h << endl;
            MoveResizeWindow(win.window_x, win.window_y, win.window_w, win.window_h);
            MyExpose();
         }
         break;
         case Command::Keys:
         {
            scr >> win.data_state.keys;
            cout << "Script: keys: '" << win.data_state.keys << "'" << endl;
            // SendKeySequence(keys.c_str());
            CallKeySequence(win.data_state.keys.c_str());
            MyExpose();
         }
         break;
         case Command::Palette:
         {
            int pal;
            scr >> pal;
            cout << "Script: palette: " << pal << endl;
            win.vs->palette.SetIndex(pal-1);
            MyExpose();
         }
         case Command::PaletteRepeat:
         {
            int rpt_times;
            scr >> rpt_times;
            cout << "Script: palette_repeat: " << rpt_times << endl;
            win.vs->palette.SetRepeatTimes(rpt_times);
            win.vs->palette.GenerateTextures();
            MyExpose();
         }
         break;
         case Command::ToggleAttributes:
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
            win.vs->ToggleAttributes(attr_list);
            MyExpose();
         }
         break;
         case Command::Rotmat:
         {
            cout << "Script: rotmat:";
            for (int i = 0; i < 16; i++)
            {
               scr >> win.vs->rotmat[i/4][i%4];
               cout << ' ' << win.vs->rotmat[i/4][i%4];
            }
            cout << endl;
            MyExpose();
         }
         break;
         case Command::Camera:
         {
            double cam[9];
            cout << "Script: camera:";
            for (int i = 0; i < 9; i++)
            {
               scr >> cam[i];
               cout << ' ' << cam[i];
            }
            cout << endl;
            win.vs->cam.Set(cam);
            MyExpose();
         }
         break;
         case Command::Scale:
         {
            double scale;
            cout << "Script: scale:";
            scr >> scale;
            cout << ' ' << scale;
            cout << endl;
            win.vs->Scale(scale);
            MyExpose();
         }
         break;
         case Command::Translate:
         {
            double x, y, z;
            cout << "Script: translate:";
            scr >> x >> y >> z;
            cout << ' ' << x << ' ' << y << ' ' << z;
            cout << endl;
            win.vs->Translate(x, y, z);
            MyExpose();
         }
         break;
         case Command::PlotCaption:
         {
            char delim;
            scr >> ws >> delim;
            getline(scr, win.plot_caption, delim);
            win.vs->PrepareCaption(); // turn on or off the caption
            MyExpose();
         }
         break;
         case Command::Max: //dummy
            break;
      }

      done_one_command = 1;
   }
}

thread_local ScriptController *ScriptController::script_ctrl = NULL;

void ScriptController::ScriptIdleFunc()
{
   script_ctrl->ExecuteScriptCommand();
   if (script_ctrl->scr_level == 0)
   {
      ScriptControl();
   }
}

void ScriptController::ScriptControl()
{
   if (script_ctrl->scr_running)
   {
      script_ctrl->scr_running = 0;
      RemoveIdleFunc(ScriptIdleFunc);
   }
   else
   {
      script_ctrl->scr_running = 1;
      AddIdleFunc(ScriptIdleFunc);
   }
}

ScriptController::ScriptController(Window &win_)
   : win(win_)
{
   if (commands.empty())
   {
      InitCommands();
   }
}

void ScriptController::PlayScript(istream &scr)
{
   string word;
   bool done = false;

   scr_min_val = numeric_limits<double>::infinity();
   scr_max_val = -scr_min_val;

   // read initializing commands
   while (!done)
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

      auto it = find(commands.begin(), commands.end(), word);
      if (it == commands.end())
      {
         cout << "Unknown command in script: " << word << endl;
         PrintCommands();
         return;
      }

      const Command cmd = (Command)(it - commands.begin());
      switch (cmd)
      {
         case Command::Window:
            scr >> win.window_x >> win.window_y >> win.window_w >> win.window_h;
            break;
         case Command::DataCollCycle:
            scr >> dc_cycle;
            break;
         case Command::DataCollProto:
            scr >> dc_protocol;
            break;
         case Command::Solution:
            if (ScriptReadSolution(scr, win.data_state))
            {
               return;
            }
            done = true; // start the visualization
            break;
         case Command::Quadrature:
            if (ScriptReadQuadrature(scr, win.data_state))
            {
               return;
            }
            done = true; // start the visualization
            break;
         case Command::ParSolution:
            if (ScriptReadParSolution(scr, win.data_state))
            {
               return;
            }
            done = true; // start the visualization
            break;
         case Command::ParQuadrature:
            if (ScriptReadParQuadrature(scr, win.data_state))
            {
               return;
            }
            done = true; // start the visualization
            break;
         case Command::Mesh:
            if (ScriptReadDisplMesh(scr, win.data_state))
            {
               return;
            }
            done = win.data_state.mesh != nullptr;
            break;
         case Command::DataCollMesh:
            if (ScriptReadDataColl(scr, win.data_state))
            {
               return;
            }
            done = true; // start the visualization
            break;
         case Command::DataCollField:
            if (ScriptReadDataColl(scr, win.data_state, false))
            {
               return;
            }
            done = true; // start the visualization
            break;
         case Command::DataCollQuad:
            if (ScriptReadDataColl(scr, win.data_state, false, true))
            {
               return;
            }
            done = true; // start the visualization
            break;
         default:
            break;
      }
   }

   scr_level = scr_running = 0;
   script = &scr;
   script_ctrl = this;
   win.data_state.keys.clear();

   // Make sure the singleton object returned by GetMainThread() is
   // initialized from the main thread.
   GetMainThread();

   std::thread worker_thread
   {
      [&](Window local_win)
      {
         if (local_win.GLVisInitVis({}))
         {
            GetAppWindow()->setOnKeyDown(SDLK_SPACE, ScriptControl);
            local_win.GLVisStartVis();
         }
      },
      std::move(win)
   };

   SDLMainLoop();
   worker_thread.join();

   delete init_nodes; init_nodes = NULL;

   cout << "Script: min_val = " << scr_min_val
        << ", max_val = " << scr_max_val << endl;

   script = NULL;
}
