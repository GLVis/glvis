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

#include "visual.hpp"
#include "palettes.hpp"

using namespace std;

extern const char *strings_off_on[]; // defined in vsdata.cpp


GLVisCommand::GLVisCommand(
   VisualizationSceneScalarData **_vs, StreamState& state, bool *_keep_attr)
   : curr_state(state)
{
   vs        = _vs;
   keep_attr = _keep_attr;
   // should be set in this thread by a call to InitVisualization()
   thread_wnd = GetAppWindow();

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

int GLVisCommand::NewMeshAndSolution(std::unique_ptr<Mesh> _new_m,
                                     std::unique_ptr<GridFunction> _new_g)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = NEW_MESH_AND_SOLUTION;
   new_state.mesh = std::move(_new_m);
   new_state.grid_f = std::move(_new_g);
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
            new_state.save_coloring = false;
            new_state.SetMeshSolution();
            mesh_range = new_state.grid_f->Max() + 1.0;
         }
         if (curr_state.SetNewMeshAndSolution(std::move(new_state), *vs))
         {
            if (mesh_range > 0.0)
            {
               (*vs)->SetValueRange(-mesh_range, mesh_range);
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
         GetAppWindow()->screenshot(screenshot_filename, true);
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
         ::plot_caption = plot_caption;
         (*vs)->PrepareCaption(); // turn on or off the caption
         MyExpose();
         break;
      }

      case AXIS_LABELS:
      {
         cout << "Command: axis_labels: '" << axis_label_x << "' '"
              << axis_label_y << "' '" << axis_label_z << "'" << endl;
         (*vs)->SetAxisLabels(axis_label_x.c_str(), axis_label_y.c_str(),
                              axis_label_z.c_str());
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
         (*vs)->SetView(view_ang_theta, view_ang_phi);
         MyExpose();
         break;
      }

      case ZOOM:
      {
         cout << "Command: zoom: " << zoom_factor << endl;
         (*vs)->Zoom(zoom_factor);
         MyExpose();
         break;
      }

      case SUBDIVISIONS:
      {
         cout << "Command: subdivisions: " << flush;
         (*vs)->SetRefineFactors(subdiv_tot, subdiv_bdr);
         cout << subdiv_tot << ' ' << subdiv_bdr << endl;
         MyExpose();
         break;
      }

      case VALUE_RANGE:
      {
         cout << "Command: valuerange: " << flush;
         (*vs)->SetValueRange(val_min, val_max);
         cout << val_min << ' ' << val_max << endl;
         MyExpose();
         break;
      }

      case SHADING:
      {
         cout << "Command: shading: " << flush;
         int s = -1;
         if (shading == "flat")
         {
            s = 0;
         }
         else if (shading == "smooth")
         {
            s = 1;
         }
         else if (shading == "cool")
         {
            s = 2;
         }
         if (s != -1)
         {
            (*vs)->SetShading(s, false);
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
         (*vs)->ViewCenterX = view_center_x;
         (*vs)->ViewCenterY = view_center_y;
         MyExpose();
         break;
      }

      case AUTOSCALE:
      {
         cout << "Command: autoscale: " << autoscale_mode;
         if (autoscale_mode == "off")
         {
            (*vs)->SetAutoscale(0);
         }
         else if (autoscale_mode == "on")
         {
            (*vs)->SetAutoscale(1);
         }
         else if (autoscale_mode == "value")
         {
            (*vs)->SetAutoscale(2);
         }
         else if (autoscale_mode == "mesh")
         {
            (*vs)->SetAutoscale(3);
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
         (*vs)->palette.SetIndex(palette-1);
         if (!GetUseTexture())
         {
            (*vs)->EventUpdateColors();
         }
         MyExpose();
         break;
      }

      case PALETTE_REPEAT:
      {
         cout << "Command: palette_repeat: " << palette_repeat << endl;
         (*vs)->palette.SetRepeatTimes(palette_repeat);
         (*vs)->palette.Init();

         if (!GetUseTexture())
         {
            (*vs)->EventUpdateColors();
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
         (*vs)->cam.Set(camera);
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

communication_thread::communication_thread(StreamCollection _is,
                                           GLVisCommand* cmd)
   : is(std::move(_is)), glvis_command(cmd)
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

      if (ident == "mesh" || ident == "solution" ||
          ident == "parallel")
      {
         bool fix_elem_orient = glvis_command->FixElementOrientations();
         StreamState tmp;
         if (ident == "mesh")
         {
            tmp.mesh.reset(new Mesh(*is[0], 1, 0, fix_elem_orient));
            if (!(*is[0]))
            {
               break;
            }
            tmp.grid_f = NULL;
         }
         else if (ident == "solution")
         {
            tmp.mesh.reset(new Mesh(*is[0], 1, 0, fix_elem_orient));
            if (!(*is[0]))
            {
               break;
            }
            tmp.grid_f.reset(new GridFunction(tmp.mesh.get(), *is[0]));
            if (!(*is[0]))
            {
               break;
            }
         }
         else if (ident == "parallel")
         {
            Array<Mesh *> mesh_array;
            Array<GridFunction *> gf_array;
            int proc, nproc, np = 0;
            bool keep_attr = glvis_command->KeepAttrib();
            do
            {
               istream &isock = *is[np];
               isock >> nproc >> proc >> ws;
#ifdef GLVIS_DEBUG
               cout << "connection[" << np << "]: parallel " << nproc << ' '
                    << proc << endl;
#endif
               isock >> ident >> ws; // "solution"
               mesh_array.SetSize(nproc);
               gf_array.SetSize(nproc);
               mesh_array[proc] = new Mesh(isock, 1, 0, fix_elem_orient);
               if (!keep_attr)
               {
                  // set element and boundary attributes to proc+1
                  for (int i = 0; i < mesh_array[proc]->GetNE(); i++)
                  {
                     mesh_array[proc]->GetElement(i)->SetAttribute(proc+1);
                  }
                  for (int i = 0; i < mesh_array[proc]->GetNBE(); i++)
                  {
                     mesh_array[proc]->GetBdrElement(i)->SetAttribute(proc+1);
                  }
               }
               gf_array[proc] = new GridFunction(mesh_array[proc], isock);
               np++;
               if (np == nproc)
               {
                  break;
               }
               *is[np] >> ident >> ws; // "parallel"
            }
            while (1);
            tmp.mesh.reset(new Mesh(mesh_array, nproc));
            tmp.grid_f.reset(new GridFunction(tmp.mesh.get(), gf_array, nproc));

            for (int p = 0; p < nproc; p++)
            {
               delete gf_array[nproc-1-p];
               delete mesh_array[nproc-1-p];
            }
            gf_array.DeleteAll();
            mesh_array.DeleteAll();
         }

         // cout << "Stream: new solution" << endl;

         tmp.Extrude1DMeshAndSolution();

         if (glvis_command->NewMeshAndSolution(std::move(tmp.mesh),
                                               std::move(tmp.grid_f)))
         {
            goto comm_terminate;
         }
      }
      else if (ident == "screenshot")
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
      else if (ident == "keys")
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
      else if (ident == "window_size")
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
      else if (ident == "window_geometry")
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
      else if (ident == "window_title")
      {
         char c;
         string title;

         *is[0] >> ws >> c; // read the opening char
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
      else if (ident == "plot_caption")
      {
         char c;
         string caption;

         *is[0] >> ws >> c; // read the opening char
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
      else if (ident == "axis_labels")
      {
         char c;
         string label_x, label_y, label_z;

         *is[0] >> ws >> c; // read the opening char
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
      else if (ident == "pause")
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
      else if (ident == "view")
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
      else if (ident == "zoom")
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
      else if (ident == "subdivisions")
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
      else if (ident == "valuerange")
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
      else if (ident == "shading")
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
      else if (ident == "viewcenter")
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
      else if (ident == "autoscale")
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
      else if (ident == "palette")
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
      else if (ident == "palette_repeat")
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
      else if (ident == "camera")
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
      else if (ident == "autopause")
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
      else
      {
         cout << "Stream: unknown command: " << ident << endl;
      }
   }

   cout << "Stream: end of input." << endl;

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
