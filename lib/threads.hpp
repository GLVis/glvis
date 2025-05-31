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

#ifndef GLVIS_THREADS_HPP
#define GLVIS_THREADS_HPP

#include "window.hpp"
#include <mfem.hpp>
#include <thread>
#include <atomic>
#include <condition_variable>

class GLVisCommand
{
private:
   // Pointers to global GLVis data
   Window      &win;
   SdlWindow   *thread_wnd;

   std::mutex glvis_mutex;
   std::condition_variable glvis_cond;

   int num_waiting;
   bool terminating;

   enum
   {
      NO_COMMAND = 0,
      NEW_MESH_AND_SOLUTION = 1,
      SCREENSHOT = 2,
      KEY_COMMANDS = 3,
      WINDOW_SIZE = 4,
      WINDOW_TITLE = 5,
      PAUSE = 6,
      VIEW_ANGLES = 7,
      ZOOM = 8,
      SUBDIVISIONS = 9,
      VALUE_RANGE = 10,
      SHADING = 11,
      VIEW_CENTER = 12,
      AUTOSCALE = 13,
      PALETTE = 14,
      CAMERA = 15,
      AUTOPAUSE = 16,
      WINDOW_GEOMETRY = 17,
      PLOT_CAPTION = 18,
      AXIS_LABELS = 19,
      PALETTE_REPEAT = 20,
      LEVELLINES = 21,
      AXIS_NUMBERFORMAT = 22,
      COLORBAR_NUMBERFORMAT = 23
   };

   std::atomic<bool> command_ready{false};

   // command to be executed
   int command;

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
   int PaletteRepeat(int n);
   int Levellines(double minv, double maxv, int number);
   int AxisNumberFormat(string formatting);
   int ColorbarNumberFormat(string formatting);
   int Camera(const double cam[]);
   int Autopause(const char *mode);

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
   bool is_multithread;

   // thread object
   std::thread tid;
   // signal for thread cancellation
   std::atomic<bool> terminate_thread {false};

   static void print_commands();
   bool execute_one(std::string word);
   void execute();

public:
   communication_thread(StreamCollection _is, GLVisCommand* cmd,
                        bool mulithread = true);

   bool process_one();

   ~communication_thread();
};

#endif
