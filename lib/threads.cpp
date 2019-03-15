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

#include <unistd.h>    // pipe, fcntl, write
#include <fcntl.h>     // fcntl
#include <cerrno>      // errno, EAGAIN
#include <cstdio>      // perror

#include "palettes.hpp"
#include "visual.hpp"

using namespace std;

void SetMeshSolution(Mesh *mesh, GridFunction *&grid_f, bool save_coloring);
extern const char *strings_off_on[]; // defined in vsdata.cpp

GLVisCommand *glvis_command = NULL;

GLVisCommand::GLVisCommand(
   VisualizationSceneScalarData **_vs, Mesh **_mesh, GridFunction **_grid_f,
   Vector *_sol, bool *_keep_attr, bool *_fix_elem_orient)
{
   vs        = _vs;
   mesh      = _mesh;
   grid_f    = _grid_f;
   sol       = _sol;
   keep_attr = _keep_attr;
   fix_elem_orient = _fix_elem_orient;

   pthread_mutex_init(&glvis_mutex, NULL);
   pthread_cond_init(&glvis_cond, NULL);
   num_waiting = 0;
   terminating = false;
   if (pipe(pfd) == -1)
   {
      perror("pipe()");
      exit(EXIT_FAILURE);
   }
   int flag = fcntl(pfd[0], F_GETFL);
   fcntl(pfd[0], F_SETFL, flag | O_NONBLOCK);

   command = NO_COMMAND;

   autopause = 0;
}

int GLVisCommand::lock()
{
   int my_id;
   pthread_mutex_lock(&glvis_mutex);
   if (terminating)
   {
      pthread_mutex_unlock(&glvis_mutex);
      return -1;
   }
   my_id = num_waiting++;
   while (my_id > 0)
   {
      pthread_cond_wait(&glvis_cond, &glvis_mutex);
      if (terminating)
      {
         num_waiting--;
         pthread_mutex_unlock(&glvis_mutex);
         return -1;
      }
      my_id--;
   }
   pthread_mutex_unlock(&glvis_mutex);
   return 0;
}

int GLVisCommand::signal()
{
   char c = 's';
   if (write(pfd[1], &c, 1) != 1)
   {
      return -1;
   }
   return 0;
}

void GLVisCommand::unlock()
{
   pthread_mutex_lock(&glvis_mutex);
   num_waiting--;
   if (num_waiting > 0)
   {
      pthread_cond_broadcast(&glvis_cond);
   }
   pthread_mutex_unlock(&glvis_mutex);
}

int GLVisCommand::NewMeshAndSolution(Mesh *_new_m, GridFunction *_new_g)
{
   if (lock() < 0)
   {
      return -1;
   }
   command = NEW_MESH_AND_SOLUTION;
   new_m = _new_m;
   new_g = _new_g;
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

extern GridFunction *ProjectVectorFEGridFunction(GridFunction*);

int GLVisCommand::Execute()
{
   char c;
   int n = read(pfd[0], &c, 1);

   if (n == -1 && errno == EAGAIN)
   {
      return 1;
   }
   if (n != 1 || c != 's')
   {
      return -1;
   }

   switch (command)
   {
      case NO_COMMAND:
         break;

      case NEW_MESH_AND_SOLUTION:
      {
         double mesh_range = -1.0;
         if (new_g == NULL)
         {
            SetMeshSolution(new_m, new_g, false);
            mesh_range = new_g->Max() + 1.0;
         }
         if (new_m->SpaceDimension() == (*mesh)->SpaceDimension() &&
             new_g->VectorDim() == (*grid_f)->VectorDim())
         {
            if (new_m->SpaceDimension() == 2)
            {
               if (new_g->VectorDim() == 1)
               {
                  VisualizationSceneSolution *vss =
                     dynamic_cast<VisualizationSceneSolution *>(*vs);
                  new_g->GetNodalValues(*sol);
                  vss->NewMeshAndSolution(new_m, sol, new_g);
               }
               else
               {
                  VisualizationSceneVector *vsv =
                     dynamic_cast<VisualizationSceneVector *>(*vs);
                  vsv->NewMeshAndSolution(*new_g);
               }
            }
            else
            {
               if (new_g->VectorDim() == 1)
               {
                  VisualizationSceneSolution3d *vss =
                     dynamic_cast<VisualizationSceneSolution3d *>(*vs);
                  new_g->GetNodalValues(*sol);
                  vss->NewMeshAndSolution(new_m, sol, new_g);
               }
               else
               {
                  new_g = ProjectVectorFEGridFunction(new_g);
                  VisualizationSceneVector3d *vss =
                     dynamic_cast<VisualizationSceneVector3d *>(*vs);
                  vss->NewMeshAndSolution(new_m, new_g);
               }
            }
            if (mesh_range > 0.0)
            {
               (*vs)->SetValueRange(-mesh_range, mesh_range);
            }
            delete (*grid_f);
            *grid_f = new_g;
            delete (*mesh);
            *mesh = new_m;

            (*vs)->Draw();
         }
         else
         {
            cout << "Stream: field type does not match!" << endl;
            delete new_g;
            delete new_m;
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
         cout << "Command: screenshot: " << flush;
         if (::Screenshot(screenshot_filename.c_str(), true))
         {
            cout << "Screenshot(" << screenshot_filename << ") failed." << endl;
         }
         else
         {
            cout << "-> " << screenshot_filename << endl;
         }
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
         (*vs)->UpdateCaption(); // turn on or off the caption
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
         Set_Palette(palette-1);
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
         RepeatPaletteTimes = palette_repeat;
         Set_Texture_Image();
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
   char c;
   int n = read(pfd[0], &c, 1);

   pthread_mutex_lock(&glvis_mutex);
   terminating = true;
   pthread_mutex_unlock(&glvis_mutex);
   if (n == 1 && c == 's')
   {
      switch (command)
      {
         case NEW_MESH_AND_SOLUTION:
            delete new_g;
            delete new_m;
            break;
      }
      unlock();
   }
   pthread_mutex_lock(&glvis_mutex);
   if (num_waiting > 0)
   {
      pthread_cond_broadcast(&glvis_cond);
   }
   pthread_mutex_unlock(&glvis_mutex);
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
      cout << "\nGLVisCommand::~GLVisCommand() : num_waiting = "
           << num_waiting << '\n' << endl;
   close(pfd[0]);
   close(pfd[1]);
   pthread_cond_destroy(&glvis_cond);
   pthread_mutex_destroy(&glvis_mutex);
}

communication_thread::communication_thread(Array<istream *> &_is)
   : is(_is)
{
   new_m = NULL;
   new_g = NULL;

   if (is.Size() > 0)
   {
      pthread_create(&tid, NULL, communication_thread::execute, this);
   }
}

communication_thread::~communication_thread()
{
   if (is.Size() > 0)
   {
      pthread_cancel(tid);
      pthread_join(tid, NULL);
   }

   delete new_g;
   delete new_m;
}

// defined in glvis.cpp
extern void Extrude1DMeshAndSolution(Mesh **, GridFunction **, Vector *);

void *communication_thread::execute(void *p)
{
   communication_thread *_this = (communication_thread *)p;

   while (1)
   {
      *_this->is[0] >> ws; // thread cancellation point

      _this->cancel_off();

      *_this->is[0] >> _this->ident;
      if (!(*_this->is[0]))
      {
         break;
      }

      if (_this->ident == "mesh" || _this->ident == "solution" ||
          _this->ident == "parallel")
      {
         bool fix_elem_orient = glvis_command->FixElementOrientations();
         if (_this->ident == "mesh")
         {
            _this->new_m = new Mesh(*_this->is[0], 1, 0, fix_elem_orient);
            if (!(*_this->is[0]))
            {
               break;
            }
            _this->new_g = NULL;
         }
         else if (_this->ident == "solution")
         {
            _this->new_m = new Mesh(*_this->is[0], 1, 0, fix_elem_orient);
            if (!(*_this->is[0]))
            {
               break;
            }
            _this->new_g = new GridFunction(_this->new_m, *_this->is[0]);
            if (!(*_this->is[0]))
            {
               break;
            }
         }
         else if (_this->ident == "parallel")
         {
            Array<Mesh *> mesh_array;
            Array<GridFunction *> gf_array;
            int proc, nproc, np = 0;
            bool keep_attr = glvis_command->KeepAttrib();
            do
            {
               istream &isock = *_this->is[np];
               isock >> nproc >> proc >> ws;
#ifdef GLVIS_DEBUG
               cout << "connection[" << np << "]: parallel " << nproc << ' '
                    << proc << endl;
#endif
               isock >> _this->ident >> ws; // "solution"
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
               *_this->is[np] >> _this->ident >> ws; // "parallel"
            }
            while (1);
            _this->new_m = new Mesh(mesh_array, nproc);
            _this->new_g = new GridFunction(_this->new_m, gf_array, nproc);

            for (int p = 0; p < nproc; p++)
            {
               delete gf_array[nproc-1-p];
               delete mesh_array[nproc-1-p];
            }
            gf_array.DeleteAll();
            mesh_array.DeleteAll();
         }

         // cout << "Stream: new solution" << endl;

         Extrude1DMeshAndSolution(&_this->new_m, &_this->new_g, NULL);

         if (glvis_command->NewMeshAndSolution(_this->new_m, _this->new_g))
         {
            goto comm_terminate;
         }

         _this->new_m = NULL;
         _this->new_g = NULL;
      }
      else if (_this->ident == "screenshot")
      {
         string filename;

         *_this->is[0] >> ws >> filename;

         // all processors sent the screenshot command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'screenshot'
            *_this->is[i] >> ws >> _this->ident; // filename
         }

         if (glvis_command->Screenshot(filename.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "keys")
      {
         string keys;

         *_this->is[0] >> ws >> keys;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'keys'
            *_this->is[i] >> ws >> _this->ident; // keys
         }

         if (glvis_command->KeyCommands(keys.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "window_size")
      {
         int w, h, t;

         *_this->is[0] >> w >> h;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'window_size'
            *_this->is[i] >> t >> t;
         }

         if (glvis_command->WindowSize(w, h))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "window_geometry")
      {
         int x, y, w, h, t;

         *_this->is[0] >> x >> y >> w >> h;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'window_geometry'
            *_this->is[i] >> t >> t >> t >> t;
         }

         if (glvis_command->WindowGeometry(x, y, w, h))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "window_title")
      {
         char c;
         string title;

         *_this->is[0] >> ws >> c; // read the opening char
         // use the opening char as termination as well
         getline(*_this->is[0], title, c);

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'window_title'
            *_this->is[i] >> ws >> c;
            getline(*_this->is[i], _this->ident, c);
         }

         if (glvis_command->WindowTitle(title.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "plot_caption")
      {
         char c;
         string caption;

         *_this->is[0] >> ws >> c; // read the opening char
         // use the opening char as termination as well
         getline(*_this->is[0], caption, c);

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'plot_caption'
            *_this->is[i] >> ws >> c;
            getline(*_this->is[i], _this->ident, c);
         }

         if (glvis_command->PlotCaption(caption.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "axis_labels")
      {
         char c;
         string label_x, label_y, label_z;

         *_this->is[0] >> ws >> c; // read the opening char
         // use the opening char as termination as well
         getline(*_this->is[0], label_x, c);
         *_this->is[0] >> ws >> c;
         getline(*_this->is[0], label_y, c);
         *_this->is[0] >> ws >> c;
         getline(*_this->is[0], label_z, c);

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'axis_label'
            *_this->is[i] >> ws >> c;
            getline(*_this->is[i], _this->ident, c);
            *_this->is[i] >> ws >> c;
            getline(*_this->is[i], _this->ident, c);
            *_this->is[i] >> ws >> c;
            getline(*_this->is[i], _this->ident, c);
         }

         if (glvis_command->AxisLabels(label_x.c_str(),
                                       label_y.c_str(),
                                       label_z.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "pause")
      {
         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'pause'
         }

         if (glvis_command->Pause())
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "view")
      {
         double theta, phi, a;

         *_this->is[0] >> theta >> phi;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'view'
            *_this->is[i] >> a >> a;
         }

         if (glvis_command->ViewAngles(theta, phi))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "zoom")
      {
         double factor, a;

         *_this->is[0] >> factor;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'zoom'
            *_this->is[i] >> a;
         }

         if (glvis_command->Zoom(factor))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "subdivisions")
      {
         int tot, bdr, a;

         *_this->is[0] >> tot >> bdr;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'subdivisions'
            *_this->is[i] >> a >> a;
         }

         if (glvis_command->Subdivisions(tot, bdr))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "valuerange")
      {
         double minv, maxv, a;

         *_this->is[0] >> minv >> maxv;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'valuerange'
            *_this->is[i] >> a >> a;
         }

         if (glvis_command->ValueRange(minv, maxv))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "shading")
      {
         string shd;

         *_this->is[0] >> ws >> shd;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'shading'
            *_this->is[i] >> ws >> _this->ident;
         }

         if (glvis_command->SetShading(shd.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "viewcenter")
      {
         double x, y, a;

         *_this->is[0] >> x >> y;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'viewcenter'
            *_this->is[i] >> a >> a;
         }

         if (glvis_command->ViewCenter(x, y))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "autoscale")
      {
         string mode;

         *_this->is[0] >> ws >> mode;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'autoscale'
            *_this->is[i] >> ws >> _this->ident;
         }

         if (glvis_command->Autoscale(mode.c_str()))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "palette")
      {
         int pal, a;

         *_this->is[0] >> pal;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'palette'
            *_this->is[i] >> a;
         }

         if (glvis_command->Palette(pal))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "palette_repeat")
      {
         int n, a;

         *_this->is[0] >> n;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'palette_repeat'
            *_this->is[i] >> a;
         }

         if (glvis_command->PaletteRepeat(n))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "camera")
      {
         double cam[9], a;

         for (int i = 0; i < 9; i++)
         {
            *_this->is[0] >> cam[i];
         }

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'camera'
            for (int j = 0; j < 9; j++)
            {
               *_this->is[i] >> a;
            }
         }

         if (glvis_command->Camera(cam))
         {
            goto comm_terminate;
         }
      }
      else if (_this->ident == "autopause")
      {
         string mode;

         *_this->is[0] >> ws >> mode;

         // all processors sent the command
         for (int i = 1; i < _this->is.Size(); i++)
         {
            *_this->is[i] >> ws >> _this->ident; // 'autopause'
            *_this->is[i] >> ws >> _this->ident;
         }

         if (glvis_command->Autopause(mode.c_str()))
         {
            goto comm_terminate;
         }
      }
      else
      {
         cout << "Stream: unknown command: " << _this->ident << endl;
      }

      _this->cancel_on();
   }

   cout << "Stream: end of input." << endl;

comm_terminate:
   for (int i = 0; i < _this->is.Size(); i++)
   {
      socketstream *isock = dynamic_cast<socketstream *>(_this->is[i]);
      if (isock)
      {
         isock->close();
      }
   }
   _this->cancel_on();

   return p;
}
