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

#include <array>
#include <algorithm>
#include <thread>
#include <vector>

#include <mfem.hpp>

#include "threads.hpp"
#include "vsdata.hpp"
#include "palettes.hpp"

using namespace std;
using namespace mfem;

extern const char *strings_off_on[]; // defined in vsdata.cpp

GLVisCommand::GLVisCommand(Window &win_)
   : win(win_)
{
   // should be set in this thread by a call to InitVisualization()
   thread_wnd = win.wnd.get();

   num_waiting = 0;
   terminating = false;
   command = Command::NO_COMMAND;

   autopause = 0;
}

int GLVisCommand::lock()
{
   int my_id;
   unique_lock<mutex> scope_lock(glvis_mutex);
   if (terminating)
   {
      return -1;
   }
   my_id = num_waiting++;
   while (my_id > 0)
   {
      glvis_cond.wait(scope_lock);
      if (terminating)
      {
         num_waiting--;
         return -1;
      }
      my_id--;
   }
   return 0;
}

int GLVisCommand::signal()
{
   command_ready = true;

   if (thread_wnd)
   {
      thread_wnd->signalLoop();
   }

   return 0;
}

void GLVisCommand::unlock()
{
   command_ready = false;

   lock_guard<mutex> scope_lock(glvis_mutex);
   num_waiting--;
   if (num_waiting > 0)
   {
      glvis_cond.notify_all();
   }
}

int GLVisCommand::NewMeshAndSolution(DataState &&ss)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::NEW_MESH_AND_SOLUTION;
   new_state = std::move(ss);
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Screenshot(const char *filename)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::SCREENSHOT;
   screenshot_filename = filename;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::KeyCommands(const char *keys)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::KEY_COMMANDS;
   key_commands = keys;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::WindowSize(int w, int h)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::WINDOW_SIZE;
   window_w = w;
   window_h = h;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::WindowGeometry(int x, int y, int w, int h)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::WINDOW_GEOMETRY;
   window_x = x;
   window_y = y;
   window_w = w;
   window_h = h;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::WindowTitle(const char *title)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::WINDOW_TITLE;
   window_title = title;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::PlotCaption(const char *caption)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::PLOT_CAPTION;
   plot_caption = caption;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::AxisLabels(const char *a_x, const char *a_y, const char *a_z)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::AXIS_LABELS;
   axis_label_x = a_x;
   axis_label_y = a_y;
   axis_label_z = a_z;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::AxisNumberFormat(string formatting)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::AXIS_NUMBERFORMAT;
   axis_formatting = formatting;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::ColorbarNumberFormat(string formatting)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::COLORBAR_NUMBERFORMAT;
   colorbar_formatting = formatting;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Pause()
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::PAUSE;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::ViewAngles(double theta, double phi)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::VIEW_ANGLES;
   view_ang_theta = theta;
   view_ang_phi   = phi;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Zoom(double factor)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::ZOOM;
   zoom_factor = factor;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Subdivisions(int tot, int bdr)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::SUBDIVISIONS;
   subdiv_tot = tot;
   subdiv_bdr = bdr;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::ValueRange(double minv, double maxv)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::VALUE_RANGE;
   val_min = minv;
   val_max = maxv;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Levellines(double minv, double maxv, int number)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::LEVELLINES;
   lvl_min = minv;
   lvl_max = maxv;
   lvl_num = number;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::SetShading(const char *shd)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::SHADING;
   shading = shd;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::ViewCenter(double x, double y)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::VIEW_CENTER;
   view_center_x = x;
   view_center_y = y;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Autoscale(const char *mode)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::AUTOSCALE;
   autoscale_mode = mode;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Palette(int pal)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::PALETTE;
   palette = pal;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::PaletteName(std::string palname)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::PALETTE_NAME;
   palette_name = palname;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::PaletteFile(std::string filename)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::PALETTE_FILE;
   palette_file = filename;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::PaletteRepeat(int n)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::PALETTE_REPEAT;
   palette_repeat = n;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Camera(const double cam[])
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::CAMERA;
   for (int i = 0; i < 9; i++)
   {
      camera[i] = cam[i];
   }
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Autopause(const char *mode)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::AUTOPAUSE;
   autopause_mode = mode;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Quit()
{
   if (lock() < 0)
   {
      return -1;
   }
   command = Command::QUIT;
   if (signal() < 0)
   {
      return -2;
   }
   return 0;
}

int GLVisCommand::Execute()
{
   if (!command_ready)
   {
      return 1;
   }

   switch (command)
   {
      case Command::NO_COMMAND:
         break;

      case Command::NEW_MESH_AND_SOLUTION:
      {
         double mesh_range = -1.0;
         switch (new_state.GetType())
         {
            case DataState::FieldType::MESH:
               mesh_range = new_state.grid_f->Max() + 1.0;
               break;
            case DataState::FieldType::SCALAR:
            case DataState::FieldType::VECTOR:
               if (new_state.quad_f)
               {
                  auto qs = win.data_state.GetQuadSolution();
                  if (qs != DataState::QuadSolution::NONE)
                  {
                     new_state.SetQuadSolution(qs);
                  }
                  else
                  {
                     new_state.SetQuadSolution();
                  }
                  new_state.ExtrudeMeshAndSolution();
               }
               else if (new_state.cgrid_f)
               {
                  auto cs = win.data_state.GetComplexSolution();
                  if (cs != DataState::ComplexSolution::NONE)
                  {
                     new_state.SetComplexSolution(cs);
                  }
                  else
                  {
                     new_state.SetComplexSolution();
                  }
                  new_state.ExtrudeMeshAndSolution();
               }
               break;
            default:
               cerr << "Unknown field type" << endl;
               break;
         }
         if (win.SetNewMeshAndSolution(std::move(new_state)))
         {
            if (mesh_range > 0.0)
            {
               win.vs->SetValueRange(-mesh_range, mesh_range);
            }
            MyExpose();
         }
         else
         {
            cout << "Stream: field type does not match!" << endl;
         }
         if (autopause)
         {
            cout << "Autopause ..." << endl;
            ThreadsStop();
         }
         break;
      }

      case Command::SCREENSHOT:
      {
         cout << "Command: screenshot -> " << screenshot_filename << endl;
         // Allow GlWindow to handle the expose and screenshot action, in case
         // any actions need to be taken before MyExpose().
         thread_wnd->screenshot(screenshot_filename, true);
         break;
      }

      case Command::KEY_COMMANDS:
      {
         cout << "Command: keys: '" << key_commands << "'" << endl;
         // SendKeySequence(key_commands.c_str());
         CallKeySequence(key_commands.c_str());
         MyExpose();
         break;
      }

      case Command::WINDOW_SIZE:
      {
         cout << "Command: window_size: " << window_w << " x " << window_h << endl;
         ResizeWindow(window_w, window_h);
         break;
      }

      case Command::WINDOW_GEOMETRY:
      {
         cout << "Command: window_geometry: "
              << "@(" << window_x << "," << window_y << ") "
              << window_w << " x " << window_h << endl;
         MoveResizeWindow(window_x, window_y, window_w, window_h);
         break;
      }

      case Command::WINDOW_TITLE:
      {
         cout << "Command: window_title: " << window_title << endl;
         SetWindowTitle(window_title.c_str());
         break;
      }

      case Command::PLOT_CAPTION:
      {
         cout << "Command: plot_caption: " << plot_caption << endl;
         win.plot_caption = plot_caption;
         win.vs->PrepareCaption(); // turn on or off the caption
         MyExpose();
         break;
      }

      case Command::AXIS_LABELS:
      {
         cout << "Command: axis_labels: '" << axis_label_x << "' '"
              << axis_label_y << "' '" << axis_label_z << "'" << endl;
         win.vs->SetAxisLabels(axis_label_x.c_str(), axis_label_y.c_str(),
                               axis_label_z.c_str());
         MyExpose();
         break;
      }

      case Command::AXIS_NUMBERFORMAT:
      {
         cout << "Command: axis_numberformat: '"
              << axis_formatting << "'" << endl;
         win.vs->SetAxisNumberFormat(axis_formatting);
         MyExpose();
         break;
      }

      case Command::COLORBAR_NUMBERFORMAT:
      {
         cout << "Command: colorbar_numberformat: '"
              << colorbar_formatting << "'" << endl;
         win.vs->SetColorbarNumberFormat(colorbar_formatting);
         MyExpose();
         break;
      }

      case Command::PAUSE:
      {
         cout << "Command: pause: ";
         ToggleThreads();
         break;
      }

      case Command::VIEW_ANGLES:
      {
         cout << "Command: view: " << view_ang_theta << ' ' << view_ang_phi
              << endl;
         win.vs->SetView(view_ang_theta, view_ang_phi);
         MyExpose();
         break;
      }

      case Command::ZOOM:
      {
         cout << "Command: zoom: " << zoom_factor << endl;
         win.vs->Zoom(zoom_factor);
         MyExpose();
         break;
      }

      case Command::SUBDIVISIONS:
      {
         cout << "Command: subdivisions: " << flush;
         win.vs->SetRefineFactors(subdiv_tot, subdiv_bdr);
         cout << subdiv_tot << ' ' << subdiv_bdr << endl;
         MyExpose();
         break;
      }

      case Command::VALUE_RANGE:
      {
         cout << "Command: valuerange: " << flush;
         win.vs->SetValueRange(val_min, val_max);
         cout << val_min << ' ' << val_max << endl;
         MyExpose();
         break;
      }

      case Command::LEVELLINES:
      {
         cout << "Command: levellines: " << flush;
         win.vs->SetLevelLines(lvl_min, lvl_max, lvl_num);
         win.vs->UpdateLevelLines();
         cout << lvl_min << ' ' << lvl_max << ' ' << lvl_num << endl;
         MyExpose();
         break;
      }

      case Command::SHADING:
      {
         cout << "Command: shading: " << flush;
         VisualizationSceneScalarData::Shading s =
            VisualizationSceneScalarData::Shading::Invalid;
         if (shading == "flat")
         {
            s = VisualizationSceneScalarData::Shading::Flat;
         }
         else if (shading == "smooth")
         {
            s = VisualizationSceneScalarData::Shading::Smooth;
         }
         else if (shading == "cool")
         {
            s = VisualizationSceneScalarData::Shading::Noncomforming;
         }
         if (s != VisualizationSceneScalarData::Shading::Invalid)
         {
            win.vs->SetShading(s, false);
            cout << shading << endl;
            MyExpose();
         }
         else
         {
            cout << shading << " ?" << endl;
         }
         break;
      }

      case Command::VIEW_CENTER:
      {
         cout << "Command: viewcenter: "
              << view_center_x << ' ' << view_center_y << endl;
         win.vs->ViewCenterX = view_center_x;
         win.vs->ViewCenterY = view_center_y;
         MyExpose();
         break;
      }

      case Command::AUTOSCALE:
      {
         cout << "Command: autoscale: " << autoscale_mode;
         if (autoscale_mode == "off")
         {
            win.vs->SetAutoscale(VisualizationSceneScalarData::Autoscale::None);
         }
         else if (autoscale_mode == "on")
         {
            win.vs->SetAutoscale(VisualizationSceneScalarData::Autoscale::MeshAndValue);
         }
         else if (autoscale_mode == "value")
         {
            win.vs->SetAutoscale(VisualizationSceneScalarData::Autoscale::Value);
         }
         else if (autoscale_mode == "mesh")
         {
            win.vs->SetAutoscale(VisualizationSceneScalarData::Autoscale::Mesh);
         }
         else
         {
            cout << '?';
         }
         cout << endl;
         break;
      }

      case Command::PALETTE:
      {
         cout << "Command: palette: " << palette << endl;
         win.vs->palette.SetIndex(palette-1);
         if (!GetUseTexture())
         {
            win.vs->EventUpdateColors();
         }
         MyExpose();
         break;
      }

      case Command::PALETTE_NAME:
      {
         cout << "Command: palette_name: " << palette_name << endl;
         win.vs->palette.SetByName(palette_name);
         if (!GetUseTexture())
         {
            win.vs->EventUpdateColors();
         }
         MyExpose();
         break;
      }

      case Command::PALETTE_FILE:
      {
         cout << "Command: palette_file: " << palette_file << endl;
         BasePalettes.Load(palette_file);
         win.vs->palette.GenerateTextures(true); // need to reinitialize
         MyExpose();
         break;
      }

      case Command::PALETTE_REPEAT:
      {
         cout << "Command: palette_repeat: " << palette_repeat << endl;
         win.vs->palette.SetRepeatTimes(palette_repeat);
         win.vs->palette.GenerateTextures();

         if (!GetUseTexture())
         {
            win.vs->EventUpdateColors();
         }
         MyExpose();
         break;
      }

      case Command::CAMERA:
      {
         cout << "Command: camera: ";
         for (int i = 0; i < 9; i++)
         {
            cout << ' ' << camera[i];
         }
         cout << endl;
         win.vs->cam.Set(camera);
         MyExpose();
         break;
      }

      case Command::AUTOPAUSE:
      {
         if (autopause_mode == "off" || autopause_mode == "0")
         {
            autopause = 0;
         }
         else
         {
            autopause = 1;
         }
         cout << "Command: autopause: " << strings_off_on[autopause] << endl;
         if (autopause)
         {
            ThreadsStop();
         }
         else
         {
            ThreadsRun();   // probably not needed
         }
         break;
      }

      case Command::QUIT:
      {
         thread_wnd->signalQuit();
         break;
      }
   }

   command = Command::NO_COMMAND;
   unlock();
   return 0;
}

void GLVisCommand::Terminate()
{
   {
      lock_guard<mutex> scope_lock(glvis_mutex);
      terminating = true;
   }
   {
      lock_guard<mutex> scope_lock(glvis_mutex);
      if (num_waiting > 0)
      {
         glvis_cond.notify_all();
      }
   }
}

void GLVisCommand::ToggleAutopause()
{
   autopause = autopause ? 0 : 1;
   cout << "Autopause: " << strings_off_on[autopause] << endl;
   if (autopause)
   {
      ThreadsStop();
   }
   else
   {
      ThreadsRun();
   }
}

GLVisCommand::~GLVisCommand()
{
   if (num_waiting > 0)
   {
      cout << "\nGLVisCommand::~GLVisCommand() : num_waiting = "
           << num_waiting << '\n' << endl;
   }
}

enum class ThreadCommand
{
   Parallel,
   Screenshot,
   Viewcenter,
   View,
   Zoom,
   Shading,
   Subdivisions,
   Valuerange,
   Autoscale,
   Levellines,
   AxisNumberFormat,
   ColorbarNumberFormat,
   WindowSize,
   WindowGeometry,
   WindowTitle,
   Keys,
   Palette,
   PaletteName,
   PaletteFile,
   PaletteRepeat,
   Camera,
   PlotCaption,
   AxisLabels,
   Pause,
   Autopause,
   //----------
   Max
};
class ThreadCommands
{
   struct CmdItem
   {
      const char *keyword;
      const char *params;
      const char *desc;

      bool operator==(const string &key) const { return key == keyword; }
   };
   array<CmdItem,(size_t)ThreadCommand::Max> commands;

public:
   ThreadCommands();

   decltype(commands)::const_iterator begin() const { return commands.begin(); }
   decltype(commands)::const_iterator end() const { return commands.end(); }
   CmdItem& operator[](ThreadCommand cmd) { return commands[(size_t)cmd]; }
   const CmdItem& operator[](ThreadCommand cmd) const { return commands[(size_t)cmd]; }
};
static const ThreadCommands commands;

ThreadCommands::ThreadCommands()
{
   (*this)[ThreadCommand::Parallel]             = {"parallel", "<num proc> <proc>", "Prefix for distributed mesh/solution/quadrature."};
   (*this)[ThreadCommand::Screenshot]           = {"screenshot", "<file>", "Take a screenshot, saving it to the file."};
   (*this)[ThreadCommand::Viewcenter]           = {"viewcenter", "<x> <y>", "Change the viewcenter."};
   (*this)[ThreadCommand::View]                 = {"view", "<theta> <phi>", "Change the solid angle of view."};
   (*this)[ThreadCommand::Zoom]                 = {"zoom", "<zoom>", "Change the zoom factor."};
   (*this)[ThreadCommand::Shading]              = {"shading", "<flat/smooth/cool>", "Change the shading algorithm."};
   (*this)[ThreadCommand::Subdivisions]         = {"subdivisions", "<times> <dummy>", "Change the refinement level."};
   (*this)[ThreadCommand::Valuerange]           = {"valuerange", "<min> <max>", "Change the value range."};
   (*this)[ThreadCommand::Autoscale]            = {"autoscale", "<off/on/value/mesh>", "Change the autoscale algorithm."};
   (*this)[ThreadCommand::Levellines]           = {"levellines", "<min> <max> <num>", "Set the level lines."};
   (*this)[ThreadCommand::AxisNumberFormat]     = {"axis_numberformat", "'<format>'", "Set the axis number format."};
   (*this)[ThreadCommand::ColorbarNumberFormat] = {"colorbar_numberformat", "'<format>'", "Set the colorbar number format."};
   (*this)[ThreadCommand::WindowSize]           = {"window_size", "<w> <h>", "Set the size of the window."};
   (*this)[ThreadCommand::WindowGeometry]       = {"window_geometry", "<x> <y> <w> <h>", "Set the position and size of the window."};
   (*this)[ThreadCommand::WindowTitle]          = {"window_title", "'<title>'", "Set title of the window."};
   (*this)[ThreadCommand::Keys]                 = {"keys", "<keys>", "Send the control key sequence."};
   (*this)[ThreadCommand::Palette]              = {"palette", "<index>", "Set the palette index."};
   (*this)[ThreadCommand::PaletteFile]          = {"palette_file", "<filename>", "Load in a palette file."};
   (*this)[ThreadCommand::PaletteName]          = {"palette_name", "<palette_name>", "Use palette with given name."};
   (*this)[ThreadCommand::PaletteRepeat]        = {"palette_repeat", "<times>", "Set the repetition of the palette."};
   (*this)[ThreadCommand::Camera]               = {"camera", "<cam[0]> ... <cam[2]> <dir[0]> ... <dir[2]> <up[0]> ... <up[2]>", "Set the camera position, direction and upward vector."};
   (*this)[ThreadCommand::PlotCaption]          = {"plot_caption", "'<caption>'", "Set the plot caption."};
   (*this)[ThreadCommand::AxisLabels]           = {"axis_labels", "'<x label>' '<y label>' '<z label>'", "Set labels of the axes."};
   (*this)[ThreadCommand::Pause]                = {"pause", "", "Stop the stream until space is pressed."};
   (*this)[ThreadCommand::Autopause]            = {"autopause", "<0/off/1/on>", "Turns off or on autopause."};
}

communication_thread::communication_thread(StreamCollection _is,
                                           GLVisCommand* cmd,
                                           bool end_quit_)
   : is(std::move(_is)), glvis_command(cmd), end_quit(end_quit_)
{
   new_m = NULL;
   new_g = NULL;

   if (is.size() > 0)
   {
      tid = std::thread(&communication_thread::execute, this);
   }
}

communication_thread::~communication_thread()
{
   if (is.size() > 0)
   {
      terminate_thread = true;
      tid.join();
   }
}

void communication_thread::print_commands()
{
   StreamReader::PrintCommands();

   for (const auto &ci : commands)
   {
      cout << "\t" << ci.keyword << " " << ci.params << " - " << ci.desc << endl;
   }
}

void communication_thread::execute()
{
   while (1)
   {
      *is[0] >> ws;
      // thread cancellation point
      if (terminate_thread) { break; }

      *is[0] >> ident;
      if (!(*is[0]))
      {
         break;
      }

      // new solution handled by StreamReader
      if (StreamReader::SupportsDataType(ident))
      {
         DataState tmp;
         tmp.fix_elem_orient = glvis_command->FixElementOrientations();
         StreamReader reader(tmp);
         reader.ReadStream(*is[0], ident);
         if (!(*is[0]))
         {
            break;
         }

         // cout << "Stream: new solution" << endl;

         if (glvis_command->NewMeshAndSolution(std::move(tmp)))
         {
            goto comm_terminate;
         }
         if (!tmp.keys.empty())
         {
            if (glvis_command->KeyCommands(tmp.keys.c_str()))
            {
               goto comm_terminate;
            }
         }
         continue;
      }

      auto it = find(commands.begin(), commands.end(), ident);
      if (it == commands.end())
      {
         cout << "Stream: unknown command: " << ident << endl;
         print_commands();
         goto comm_terminate;
      }

      const ThreadCommand cmd = (ThreadCommand)(it - commands.begin());
      switch (cmd)
      {
         case ThreadCommand::Parallel:
         {
            unsigned int proc, nproc, np = 0;
            do
            {
               istream &isock = *is[np];
               isock >> nproc >> proc >> ws;
#ifdef GLVIS_DEBUG
               cout << "connection[" << np << "]: parallel " << nproc << ' '
                    << proc << endl;
#endif
               if (nproc != is.size())
               {
                  cout << "Unexpected number of processors: " << nproc
                       << ", expected: " << is.size() << endl;
                  mfem_error();
               }
               if (proc >= nproc)
               {
                  cout << "Invalid processor rank: " << proc
                       << ", number of processors: " << nproc << endl;
                  mfem_error();
               }
               np++;
               if (np == nproc)
               {
                  break;
               }
               *is[np] >> ident >> ws; // 'parallel'
               if (ident != "parallel")
               {
                  cout << "Expected keyword \"parallel\", got \"" << ident
                       << '"' << endl;
                  mfem_error();
               }
            }
            while (1);

            DataState tmp;
            tmp.fix_elem_orient = glvis_command->FixElementOrientations();
            tmp.keep_attr = glvis_command->KeepAttrib();
            StreamReader reader(tmp);
            reader.ReadStreams(is);

            // cout << "Stream: new solution" << endl;

            if (glvis_command->NewMeshAndSolution(std::move(tmp)))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Screenshot:
         {
            string filename;

            *is[0] >> ws >> filename;

            // all processors sent the screenshot command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'screenshot'
               *is[i] >> ws >> ident; // filename
            }

            if (glvis_command->Screenshot(filename.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Keys:
         {
            string keys;

            *is[0] >> ws >> keys;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'keys'
               *is[i] >> ws >> ident; // keys
            }

            if (glvis_command->KeyCommands(keys.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::WindowSize:
         {
            int w, h, t;

            *is[0] >> w >> h;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'window_size'
               *is[i] >> t >> t;
            }

            if (glvis_command->WindowSize(w, h))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::WindowGeometry:
         {
            int x, y, w, h, t;

            *is[0] >> x >> y >> w >> h;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'window_geometry'
               *is[i] >> t >> t >> t >> t;
            }

            if (glvis_command->WindowGeometry(x, y, w, h))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::WindowTitle:
         {
            char c;
            string title;

            // read the opening char
            *is[0] >> ws >> c;
            // use the opening char as termination as well
            getline(*is[0], title, c);

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'window_title'
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
            }

            if (glvis_command->WindowTitle(title.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::PlotCaption:
         {
            char c;
            string caption;

            // read the opening char
            *is[0] >> ws >> c;
            // use the opening char as termination as well
            getline(*is[0], caption, c);

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'plot_caption'
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
            }

            if (glvis_command->PlotCaption(caption.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::AxisLabels:
         {
            char c;
            string label_x, label_y, label_z;

            // read the opening char
            *is[0] >> ws >> c;
            // use the opening char as termination as well
            getline(*is[0], label_x, c);
            *is[0] >> ws >> c;
            getline(*is[0], label_y, c);
            *is[0] >> ws >> c;
            getline(*is[0], label_z, c);

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'axis_label'
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
            }

            if (glvis_command->AxisLabels(label_x.c_str(),
                                          label_y.c_str(),
                                          label_z.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Pause:
         {
            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'pause'
            }

            if (glvis_command->Pause())
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::View:
         {
            double theta, phi, a;

            *is[0] >> theta >> phi;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'view'
               *is[i] >> a >> a;
            }

            if (glvis_command->ViewAngles(theta, phi))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Zoom:
         {
            double factor, a;

            *is[0] >> factor;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'zoom'
               *is[i] >> a;
            }

            if (glvis_command->Zoom(factor))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Subdivisions:
         {
            int tot, bdr, a;

            *is[0] >> tot >> bdr;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'subdivisions'
               *is[i] >> a >> a;
            }

            if (glvis_command->Subdivisions(tot, bdr))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Valuerange:
         {
            double minv, maxv, a;

            *is[0] >> minv >> maxv;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'valuerange'
               *is[i] >> a >> a;
            }

            if (glvis_command->ValueRange(minv, maxv))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Levellines:
         {
            double minv, maxv, a;
            int num, b;

            *is[0] >> minv >> maxv >> num;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'levellines'
               *is[i] >> a >> a >> b;
            }

            if (glvis_command->Levellines(minv, maxv, num))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::AxisNumberFormat:
         {
            char c;
            string formatting;

            // read the opening char
            *is[0] >> ws >> c;
            // read formatting string & use c for termination
            getline(*is[0], formatting, c);

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'axis_numberformat'
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
            }

            if (glvis_command->AxisNumberFormat(formatting))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::ColorbarNumberFormat:
         {
            char c;
            string formatting;

            // read the opening char
            *is[0] >> ws >> c;
            // read formatting string & use c for termination
            getline(*is[0], formatting, c);

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'colorbar_numberformat'
               *is[i] >> ws >> c;
               getline(*is[i], ident, c);
            }

            if (glvis_command->ColorbarNumberFormat(formatting))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Shading:
         {
            string shd;

            *is[0] >> ws >> shd;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'shading'
               *is[i] >> ws >> ident;
            }

            if (glvis_command->SetShading(shd.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Viewcenter:
         {
            double x, y, a;

            *is[0] >> x >> y;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'viewcenter'
               *is[i] >> a >> a;
            }

            if (glvis_command->ViewCenter(x, y))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Autoscale:
         {
            string mode;

            *is[0] >> ws >> mode;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'autoscale'
               *is[i] >> ws >> ident;
            }

            if (glvis_command->Autoscale(mode.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Palette:
         {
            int pal, a;

            *is[0] >> pal;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'palette'
               *is[i] >> a;
            }

            if (glvis_command->Palette(pal))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::PaletteFile:
         {
            std::string filename, a;

            *is[0] >> ws;
            std::getline(*is[0], filename);

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'palette_file'
               *is[i] >> a;
            }

            if (glvis_command->PaletteFile(filename))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::PaletteName:
         {
            std::string palname, a;

            *is[0] >> palname;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'palette_name'
               *is[i] >> a;
            }

            if (glvis_command->PaletteName(palname))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::PaletteRepeat:
         {
            int n, a;

            *is[0] >> n;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'palette_repeat'
               *is[i] >> a;
            }

            if (glvis_command->PaletteRepeat(n))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Camera:
         {
            double cam[9], a;

            for (int i = 0; i < 9; i++)
            {
               *is[0] >> cam[i];
            }

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'camera'
               for (int j = 0; j < 9; j++)
               {
                  *is[i] >> a;
               }
            }

            if (glvis_command->Camera(cam))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Autopause:
         {
            string mode;

            *is[0] >> ws >> mode;

            // all processors sent the command
            for (size_t i = 1; i < is.size(); i++)
            {
               *is[i] >> ws >> ident; // 'autopause'
               *is[i] >> ws >> ident;
            }

            if (glvis_command->Autopause(mode.c_str()))
            {
               goto comm_terminate;
            }
         }
         break;
         case ThreadCommand::Max: //dummy
            break;
      }
   }

   cout << "Stream: end of input." << endl;

   if (end_quit)
   {
      glvis_command->Quit();
   }

comm_terminate:
   for (size_t i = 0; i < is.size(); i++)
   {
      socketstream *isock = dynamic_cast<socketstream *>(is[i].get());
      if (isock)
      {
         isock->close();
      }
   }
}
