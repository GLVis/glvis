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

#include "visual.hpp"
#include "palettes.hpp"

#include <array>
#include <algorithm>

using namespace std;

extern const char *strings_off_on[]; // defined in vsdata.cpp


GLVisCommand::GLVisCommand(Window &win_)
   : win(win_)
{
   // should be set in this thread by a call to InitVisualization()
   thread_wnd = win.wnd.get();

   num_waiting = 0;
   terminating = false;
   command = NO_COMMAND;

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
   command = NEW_MESH_AND_SOLUTION;
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
   command = SCREENSHOT;
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
   command = KEY_COMMANDS;
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
   command = WINDOW_SIZE;
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
   command = WINDOW_GEOMETRY;
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
   command = WINDOW_TITLE;
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
   command = PLOT_CAPTION;
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
   command = AXIS_LABELS;
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
   command = AXIS_NUMBERFORMAT;
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
   command = COLORBAR_NUMBERFORMAT;
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
   command = PAUSE;
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
   command = VIEW_ANGLES;
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
   command = ZOOM;
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
   command = SUBDIVISIONS;
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
   command = VALUE_RANGE;
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
   command = LEVELLINES;
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
   command = SHADING;
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
   command = VIEW_CENTER;
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
   command = AUTOSCALE;
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
   command = PALETTE;
   palette = pal;
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
   command = PALETTE_REPEAT;
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
   command = CAMERA;
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
   command = AUTOPAUSE;
   autopause_mode = mode;
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
      case NO_COMMAND:
         break;

      case NEW_MESH_AND_SOLUTION:
      {
         double mesh_range = -1.0;
         if (!new_state.grid_f)
         {
            if (!new_state.quad_f)
            {
               new_state.save_coloring = false;
               new_state.SetMeshSolution();
               mesh_range = new_state.grid_f->Max() + 1.0;
            }
            else
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

      case SCREENSHOT:
      {
         cout << "Command: screenshot -> " << screenshot_filename << endl;
         // Allow SdlWindow to handle the expose and screenshot action, in case
         // any actions need to be taken before MyExpose().
         thread_wnd->screenshot(screenshot_filename, true);
         break;
      }

      case KEY_COMMANDS:
      {
         cout << "Command: keys: '" << key_commands << "'" << endl;
         // SendKeySequence(key_commands.c_str());
         CallKeySequence(key_commands.c_str());
         MyExpose();
         break;
      }

      case WINDOW_SIZE:
      {
         cout << "Command: window_size: " << window_w << " x " << window_h << endl;
         ResizeWindow(window_w, window_h);
         break;
      }

      case WINDOW_GEOMETRY:
      {
         cout << "Command: window_geometry: "
              << "@(" << window_x << "," << window_y << ") "
              << window_w << " x " << window_h << endl;
         MoveResizeWindow(window_x, window_y, window_w, window_h);
         break;
      }

      case WINDOW_TITLE:
      {
         cout << "Command: window_title: " << window_title << endl;
         SetWindowTitle(window_title.c_str());
         break;
      }

      case PLOT_CAPTION:
      {
         cout << "Command: plot_caption: " << plot_caption << endl;
         win.plot_caption = plot_caption;
         win.vs->PrepareCaption(); // turn on or off the caption
         MyExpose();
         break;
      }

      case AXIS_LABELS:
      {
         cout << "Command: axis_labels: '" << axis_label_x << "' '"
              << axis_label_y << "' '" << axis_label_z << "'" << endl;
         win.vs->SetAxisLabels(axis_label_x.c_str(), axis_label_y.c_str(),
                               axis_label_z.c_str());
         MyExpose();
         break;
      }

      case AXIS_NUMBERFORMAT:
      {
         cout << "Command: axis_numberformat: '"
              << axis_formatting << "'" << endl;
         win.vs->SetAxisNumberFormat(axis_formatting);
         MyExpose();
         break;
      }

      case COLORBAR_NUMBERFORMAT:
      {
         cout << "Command: colorbar_numberformat: '"
              << colorbar_formatting << "'" << endl;
         win.vs->SetColorbarNumberFormat(colorbar_formatting);
         MyExpose();
         break;
      }

      case PAUSE:
      {
         cout << "Command: pause: ";
         ToggleThreads();
         break;
      }

      case VIEW_ANGLES:
      {
         cout << "Command: view: " << view_ang_theta << ' ' << view_ang_phi
              << endl;
         win.vs->SetView(view_ang_theta, view_ang_phi);
         MyExpose();
         break;
      }

      case ZOOM:
      {
         cout << "Command: zoom: " << zoom_factor << endl;
         win.vs->Zoom(zoom_factor);
         MyExpose();
         break;
      }

      case SUBDIVISIONS:
      {
         cout << "Command: subdivisions: " << flush;
         win.vs->SetRefineFactors(subdiv_tot, subdiv_bdr);
         cout << subdiv_tot << ' ' << subdiv_bdr << endl;
         MyExpose();
         break;
      }

      case VALUE_RANGE:
      {
         cout << "Command: valuerange: " << flush;
         win.vs->SetValueRange(val_min, val_max);
         cout << val_min << ' ' << val_max << endl;
         MyExpose();
         break;
      }

      case LEVELLINES:
      {
         cout << "Command: levellines: " << flush;
         win.vs->SetLevelLines(lvl_min, lvl_max, lvl_num);
         win.vs->UpdateLevelLines();
         cout << lvl_min << ' ' << lvl_max << ' ' << lvl_num << endl;
         MyExpose();
         break;
      }

      case SHADING:
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

      case VIEW_CENTER:
      {
         cout << "Command: viewcenter: "
              << view_center_x << ' ' << view_center_y << endl;
         win.vs->ViewCenterX = view_center_x;
         win.vs->ViewCenterY = view_center_y;
         MyExpose();
         break;
      }

      case AUTOSCALE:
      {
         cout << "Command: autoscale: " << autoscale_mode;
         if (autoscale_mode == "off")
         {
            win.vs->SetAutoscale(0);
         }
         else if (autoscale_mode == "on")
         {
            win.vs->SetAutoscale(1);
         }
         else if (autoscale_mode == "value")
         {
            win.vs->SetAutoscale(2);
         }
         else if (autoscale_mode == "mesh")
         {
            win.vs->SetAutoscale(3);
         }
         else
         {
            cout << '?';
         }
         cout << endl;
         break;
      }

      case PALETTE:
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

      case PALETTE_REPEAT:
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

      case CAMERA:
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

      case AUTOPAUSE:
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

   }

   command = NO_COMMAND;
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

enum class Command
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
   array<CmdItem,(size_t)Command::Max> commands;

public:
   ThreadCommands();

   decltype(commands)::const_iterator begin() const { return commands.begin(); }
   decltype(commands)::const_iterator end() const { return commands.end(); }
   CmdItem& operator[](Command cmd) { return commands[(size_t)cmd]; }
   const CmdItem& operator[](Command cmd) const { return commands[(size_t)cmd]; }
};
static const ThreadCommands commands;

ThreadCommands::ThreadCommands()
{
   (*this)[Command::Parallel]             = {"parallel", "<num proc> <proc>", "Prefix for distributed mesh/solution/quadrature."};
   (*this)[Command::Screenshot]           = {"screenshot", "<file>", "Take a screenshot, saving it to the file."};
   (*this)[Command::Viewcenter]           = {"viewcenter", "<x> <y>", "Change the viewcenter."};
   (*this)[Command::View]                 = {"view", "<theta> <phi>", "Change the solid angle of view."};
   (*this)[Command::Zoom]                 = {"zoom", "<zoom>", "Change the zoom factor."};
   (*this)[Command::Shading]              = {"shading", "<flat/smooth/cool>", "Change the shading algorithm."};
   (*this)[Command::Subdivisions]         = {"subdivisions", "<times> <dummy>", "Change the refinement level."};
   (*this)[Command::Valuerange]           = {"valuerange", "<min> <max>", "Change the value range."};
   (*this)[Command::Autoscale]            = {"autoscale", "<off/on/value/mesh>", "Change the autoscale algorithm."};
   (*this)[Command::Levellines]           = {"levellines", "<min> <max> <num>", "Set the level lines."};
   (*this)[Command::AxisNumberFormat]     = {"axis_numberformat", "'<format>'", "Set the axis number format."};
   (*this)[Command::ColorbarNumberFormat] = {"colorbar_numberformat", "'<format>'", "Set the colorbar number format."};
   (*this)[Command::WindowSize]           = {"window_size", "<w> <h>", "Set the size of the window."};
   (*this)[Command::WindowGeometry]       = {"window_geometry", "<x> <y> <w> <h>", "Set the position and size of the window."};
   (*this)[Command::WindowTitle]          = {"window_title", "'<title>'", "Set title of the window."};
   (*this)[Command::Keys]                 = {"keys", "<keys>", "Send the control key sequence."};
   (*this)[Command::Palette]              = {"palette", "<index>", "Set the palette index."};
   (*this)[Command::PaletteRepeat]        = {"palette_repeat", "<times>", "Set the repetition of the palette."};
   (*this)[Command::Camera]               = {"camera", "<cam[0]> ... <cam[2]> <dir[0]> ... <dir[2]> <up[0]> ... <up[2]>", "Set the camera position, direction and upward vector."};
   (*this)[Command::PlotCaption]          = {"plot_caption", "'<caption>'", "Set the plot caption."};
   (*this)[Command::AxisLabels]           = {"axis_labels", "'<x label>' '<y label>' '<z label>'", "Set labels of the axes."};
   (*this)[Command::Pause]                = {"pause", "", "Stop the stream until space is pressed."};
   (*this)[Command::Autopause]            = {"autopause", "<0/off/1/on>", "Turns off or on autopause."};
}

communication_thread::communication_thread(StreamCollection _is,
                                           GLVisCommand* cmd, bool multithread)
   : is(std::move(_is)), glvis_command(cmd), is_multithread(multithread)
{
   if (is_multithread && is.size() > 0)
   {
      tid = std::thread(&communication_thread::execute, this);
   }
}

bool communication_thread::process_one()
{
   std::string word;
   *is[0] >> ws;
   *is[0] >> word;
   return execute_one(word);
}

communication_thread::~communication_thread()
{
   if (is_multithread && is.size() > 0)
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

bool communication_thread::execute_one(std::string ident)
{
   // new solution handled by StreamReader
   if (StreamReader::SupportsDataType(ident))
   {
      DataState tmp;
      tmp.fix_elem_orient = glvis_command->FixElementOrientations();
      StreamReader reader(tmp);
      reader.ReadStream(*is[0], ident);

      // cout << "Stream: new solution" << endl;

      if (glvis_command->NewMeshAndSolution(std::move(tmp)))
      {
         return false;
      }
      if (!tmp.keys.empty())
      {
         if (glvis_command->KeyCommands(tmp.keys.c_str()))
         {
            return false;
         }
      }
      return true;
   }

   auto it = find(commands.begin(), commands.end(), ident);
   if (it == commands.end())
   {
      cout << "Stream: unknown command: " << ident << endl;
      print_commands();
      return false;
   }

   const Command cmd = (Command)(it - commands.begin());
   switch (cmd)
   {
      case Command::Parallel:
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
            *is[np] >> ident >> ws; // "parallel"
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
            return false;
         }
      }
      break;
      case Command::Screenshot:
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
            return false;
         }
      }
      break;
      case Command::Keys:
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
            return false;
         }
      }
      break;
      case Command::WindowSize:
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
            return false;
         }
      }
      break;
      case Command::WindowGeometry:
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
            return false;
         }
      }
      break;
      case Command::WindowTitle:
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
            return false;
         }
      }
      break;
      case Command::PlotCaption:
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
            return false;
         }
      }
      break;
      case Command::AxisLabels:
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
            return false;
         }
      }
      break;
      case Command::Pause:
      {
         // all processors sent the command
         for (size_t i = 1; i < is.size(); i++)
         {
            *is[i] >> ws >> ident; // 'pause'
         }

         if (glvis_command->Pause())
         {
            return false;
         }
      }
      break;
      case Command::View:
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
            return false;
         }
      }
      break;
      case Command::Zoom:
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
            return false;
         }
      }
      break;
      case Command::Subdivisions:
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
            return false;
         }
      }
      break;
      case Command::Valuerange:
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
            return false;
         }
      }
      break;
      case Command::Levellines:
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
            return false;
         }
      }
      break;
      case Command::AxisNumberFormat:
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
            return false;
         }
      }
      break;
      case Command::ColorbarNumberFormat:
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
            return false;
         }
      }
      break;
      case Command::Shading:
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
            return false;
         }
      }
      break;
      case Command::Viewcenter:
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
            return false;
         }
      }
      break;
      case Command::Autoscale:
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
            return false;
         }
      }
      break;
      case Command::Palette:
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
            return false;
         }
      }
      break;
      case Command::PaletteRepeat:
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
            return false;
         }
      }
      break;
      case Command::Camera:
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
            return false;
         }
      }
      break;
      case Command::Autopause:
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
            return false;
         }
      }
      break;
      case Command::Max: //dummy
         break;
   }

   return true;
}

void communication_thread::execute()
{
   std::string ident;
   bool status = true;

   while (status)
   {
      *is[0] >> ws;
      // thread cancellation point
      if (terminate_thread) { break; }

      *is[0] >> ident;
      if (!(*is[0]))
      {
         break;
      }

      status = execute_one(ident);
   }

   if (status)
   {
      cout << "Stream: end of input." << endl;
   }

   for (size_t i = 0; i < is.size(); i++)
   {
      socketstream *isock = dynamic_cast<socketstream *>(is[i].get());
      if (isock)
      {
         isock->close();
      }
   }
}
