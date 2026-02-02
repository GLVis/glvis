// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_THREADS_HPP
#define GLVIS_THREADS_HPP

#include <mfem.hpp>
#include <thread>
#include <atomic>
#include <condition_variable>

#include "window.hpp"
#include "vsdata.hpp"
#include "data_state.hpp"

class GLVisCommand
{
private:
   // Pointers to global GLVis data
   Window      &win;
   GLWindow    *thread_wnd;

   std::mutex glvis_mutex;
   std::condition_variable glvis_cond;

   int num_waiting;
   bool terminating;

   enum class Command
   {
      NO_COMMAND,
      NEW_MESH_AND_SOLUTION,
      SCREENSHOT,
      KEY_COMMANDS,
      WINDOW_SIZE,
      WINDOW_TITLE,
      PAUSE,
      VIEW_ANGLES,
      ZOOM,
      SUBDIVISIONS,
      VALUE_RANGE,
      SHADING,
      VIEW_CENTER,
      AUTOSCALE,
      PALETTE,
      PALETTE_NAME,
      PALETTE_FILE,
      CAMERA,
      AUTOPAUSE,
      WINDOW_GEOMETRY,
      PLOT_CAPTION,
      AXIS_LABELS,
      PALETTE_REPEAT,
      LEVELLINES,
      AXIS_NUMBERFORMAT,
      COLORBAR_NUMBERFORMAT,
      QUIT
   };

   std::atomic<bool> command_ready{false};

   // command to be executed
   Command command;

   // command arguments
   DataState     new_state;
   std::string   screenshot_filename;
   std::string   key_commands;
   int           window_x, window_y;
   int           window_w, window_h;
   std::string   window_title;
   std::string   plot_caption;
   std::string   axis_label_x;
   std::string   axis_label_y;
   std::string   axis_label_z;
   double        view_ang_theta, view_ang_phi;
   double        zoom_factor;
   int           subdiv_tot, subdiv_bdr;
   double        val_min, val_max;
   std::string   shading;
   double        view_center_x, view_center_y;
   std::string   autoscale_mode;
   int           palette, palette_repeat;
   std::string   palette_name;
   std::string   palette_file;
   double        lvl_min, lvl_max;
   int           lvl_num;
   double        camera[9];
   std::string   autopause_mode;
   std::string   axis_formatting;
   std::string   colorbar_formatting;

   // internal variables
   int autopause;

   int lock();
   int signal();
   void unlock();

public:
   // called by the main execution thread
   GLVisCommand(Window &win);

   // to be used by worker threads
   bool KeepAttrib() { return win.data_state.keep_attr; } // may need to sync this
   bool FixElementOrientations() { return win.data_state.fix_elem_orient; }

   // called by worker threads
   int NewMeshAndSolution(DataState &&ss);
   int Screenshot(const char *filename);
   int KeyCommands(const char *keys);
   int WindowSize(int w, int h);
   int WindowGeometry(int x, int y, int w, int h);
   int WindowTitle(const char *title);
   int PlotCaption(const char *caption);
   int AxisLabels(const char *a_x, const char *a_y, const char *a_z);
   int Pause();
   int ViewAngles(double theta, double phi);
   int Zoom(double factor);
   int Subdivisions(int tot, int bdr);
   int ValueRange(double minv, double maxv);
   int SetShading(const char *shd);
   int ViewCenter(double x, double y);
   int Autoscale(const char *mode);
   int Palette(int pal);
   int PaletteName(std::string palname);
   int PaletteFile(std::string filename);
   int PaletteRepeat(int n);
   int Levellines(double minv, double maxv, int number);
   int AxisNumberFormat(std::string formatting);
   int ColorbarNumberFormat(std::string formatting);
   int Camera(const double cam[]);
   int Autopause(const char *mode);
   int Quit();

   // called by the main execution thread
   int Execute();

   // called by the main execution thread
   void Terminate();

   void ToggleAutopause();

   // called by the main execution thread
   ~GLVisCommand();
};

class communication_thread
{
public:
   using StreamCollection = std::vector<std::unique_ptr<std::istream>>;

private:
   // streams to read data from
   StreamCollection is;

   GLVisCommand* glvis_command;

   // data that may be dynamically allocated by the thread
   std::unique_ptr<mfem::Mesh> new_m;
   std::unique_ptr<mfem::GridFunction> new_g;
   std::string ident;

   // thread object
   std::thread tid;
   // signal for thread cancellation
   std::atomic<bool> terminate_thread {false};

   // flag for closing the window at the end of stream
   bool end_quit;

   static void print_commands();
   void execute();

public:
   communication_thread(StreamCollection _is, GLVisCommand* cmd,
                        bool end_quit = false);

   ~communication_thread();
};

#endif // GLVIS_THREADS_HPP
