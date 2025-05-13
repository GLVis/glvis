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
#include "lib/window.hpp"
#include "lib/script_controller.hpp"
#include "lib/stream_reader.hpp"
#include "lib/file_reader.hpp"
#include "lib/coll_reader.hpp"

using namespace std;
using namespace mfem;

const char *string_none    = "(none)";
const char *string_default = "(default)";

// Global variables
enum InputOptions
{
   INPUT_SERVER_MODE = 1 << 0,
   INPUT_MESH        = 1 << 1,
   INPUT_SCALAR_SOL  = 1 << 2,
   INPUT_VECTOR_SOL  = 1 << 3,
   INPUT_GRID_FUNC   = 1 << 4,
   INPUT_QUAD_FUNC   = 1 << 5,
   INPUT_DATA_COLL   = 1 << 6,
   //...
   INPUT_PARALLEL    = 1 << 8,
};
int input = INPUT_SERVER_MODE;

thread_local GeometryRefiner GLVisGeometryRefiner;

void PrintSampleUsage(ostream &out);

class Session
{
   StreamCollection input_streams;
   Window win;
   std::thread handler;

public:
   Session(bool fix_elem_orient,
           bool save_coloring,
           string plot_caption)
   {
      win.data_state.fix_elem_orient = fix_elem_orient;
      win.data_state.save_coloring = save_coloring;
      win.plot_caption = plot_caption;
   }

   Session(Window other_win)
      : win(std::move(other_win))
   { }

   ~Session() = default;

   Session(Session&& from) = default;
   Session& operator= (Session&& from) = default;

   inline DataState& GetState() { return win.data_state; }
   inline const DataState& GetState() const { return win.data_state; }

   void StartSession(bool detached = true)
   {
      auto funcThread = [](Window w, StreamCollection is)
      {
         if (w.GLVisInitVis(std::move(is)))
         {
            w.GLVisStartVis();
         }
      };
      handler = std::thread {funcThread,
                             std::move(win), std::move(input_streams)};
      if (detached)
      {
         handler.detach();
      }
   }

   bool StartSavedSession(std::string stream_file, bool detached = true)
   {
      unique_ptr<ifstream> ifs(new ifstream(stream_file));
      if (!(*ifs))
      {
         cout << "Can not open stream file: " << stream_file << endl;
         return false;
      }
      string data_type;
      *ifs >> data_type >> ws;
      StreamReader reader(win.data_state);
      reader.ReadStream(*ifs, data_type);
      input_streams.emplace_back(std::move(ifs));

      StartSession(detached);
      return true;
   }

   int StartStreamSession(std::unique_ptr<mfem::socketstream> &&stream,
                          const std::string &data_type)
   {
      StreamReader reader(win.data_state);
      int ierr = reader.ReadStream(*stream, data_type);
      if (ierr) { return ierr; }
      input_streams.emplace_back(std::move(stream));

      StartSession();
      return 0;
   }

   int StartStreamSession(StreamCollection &&streams)
   {
      StreamReader reader(win.data_state);
      int ierr = reader.ReadStreams(streams);
      if (ierr) { return ierr; }
      input_streams = std::move(streams);

      StartSession();
      return 0;
   }

   void WaitForSession()
   {
      handler.join();
   }
};

void GLVisServer(int portnum, bool save_stream, bool fix_elem_orient,
                 bool save_coloring, string plot_caption)
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

      Session new_session(fix_elem_orient, save_coloring, plot_caption);

      constexpr int tmp_filename_size = 50;
      char tmp_file[tmp_filename_size];
      if (save_stream)
      {
         snprintf(tmp_file, tmp_filename_size, "glvis-saved.%04d", viscount);
         ofstream ofs(tmp_file);
         if (!par_data)
         {
            ofs << data_type << '\n';
            ofs << isock->rdbuf();
            isock->close();
         }
         else
         {
            StreamReader reader(new_session.GetState());
            reader.ReadStreams(input_streams);
            reader.WriteStream(ofs);
         }
         ofs.close();
         cout << "Data saved in " << tmp_file << endl;

         new_session.StartSavedSession(tmp_file);
      }
      else
      {
         if (!par_data)
         {
            new_session.StartStreamSession(std::move(isock), data_type);
         }
         else
         {
            new_session.StartStreamSession(std::move(input_streams));
         }
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
   // main Window structure
   Window win;

   // variables for command line arguments
   const char *mesh_file       = string_none;
   const char *sol_file        = string_none;
   const char *vec_sol_file    = string_none;
   const char *gfunc_file      = string_none;
   const char *qfunc_file      = string_none;
   string      dc_protocol     = string_default;
   int         dc_cycle        = 0;
   const char *arg_keys        = string_none;
   int         pad_digits      = 6;
   int         gf_component    = -1;
   int         qf_component    = -1;
   const char *c_plot_caption  = string_none;
   bool        secure          = socketstream::secure_default;
   const char *visit_coll    = string_none;
   const char *sidre_coll    = string_none;
   const char *fms_coll      = string_none;
   const char *conduit_coll  = string_none;
   int         np            = 0;
   bool        save_stream   = false;
   const char *stream_file   = string_none;
   const char *script_file   = string_none;
   const char *palette_file  = string_none;
   const char *font_name     = string_default;
   int         portnum       = 19916;
   int         multisample   = GetMultisample();
   double      line_width    = GetLineWidth();
   double      ms_line_width = GetLineWidthMS();
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
   args.AddOption(&qfunc_file, "-q", "--quadrature-function",
                  "Quadrature function file to visualize.");
   args.AddOption(&qf_component, "-qc", "--quadrature-function-component",
                  "Select a quadrature function component, [0-<num-comp>) or"
                  " -1 for all.");
   args.AddOption(&sol_file, "-s", "--scalar-solution",
                  "Scalar solution (vertex values) file to visualize.");
   args.AddOption(&vec_sol_file, "-v", "--vector-solution",
                  "Vector solution (vertex values) file to visualize.");
   args.AddOption(&visit_coll, "-visit", "--visit-datafiles",
                  "VisIt collection to load");
#ifdef MFEM_USE_SIDRE
   args.AddOption(&sidre_coll, "-sidre", "--sidre-datafiles",
                  "Sidre collection to load");
#endif // MFEM_USE_SIDRE
#ifdef MFEM_USE_FMS
   args.AddOption(&fms_coll, "-fms", "--fms-datafiles",
                  "FMS collection to load");
#endif // MFEM_USE_FMS
#ifdef MFEM_USE_CONDUIT
   args.AddOption(&conduit_coll, "-conduit", "--conduit-datafiles",
                  "Conduit collection to load");
#endif // MFEM_USE_CONDUIT
   args.AddOption(&dc_protocol, "-dc-prot", "--data-collection-protocol",
                  "Protocol of the data collection to load");
   args.AddOption(&dc_cycle, "-dc-cycle", "--data-collection-cycle",
                  "Cycle of the data collection to load");
   args.AddOption(&np, "-np", "--num-proc",
                  "Load mesh/solution from multiple processors.");
   args.AddOption(&pad_digits, "-d", "--pad-digits",
                  "Number of digits used for processor ranks in file names.");
   args.AddOption(&script_file, "-run", "--run-script",
                  "Run a GLVis script file.");
   args.AddOption(&palette_file, "-pal", "--palettes",
                  "Palette file.");
   args.AddOption(&arg_keys, "-k", "--keys",
                  "Execute key shortcut commands in the GLVis window.");
   args.AddOption(&win.data_state.fix_elem_orient, "-fo", "--fix-orientations",
                  "-no-fo", "--dont-fix-orientations",
                  "Attempt to fix the orientations of inverted elements.");
   args.AddOption(&win.data_state.keep_attr, "-a", "--real-attributes",
                  "-ap", "--processor-attributes",
                  "When opening a parallel mesh, use the real mesh attributes"
                  " or replace them with the processor rank.");
   args.AddOption(&geom_ref_type, "-grt", "--geometry-refiner-type",
                  "Set of points to use when refining geometry:"
                  " 3 = uniform, 1 = Gauss-Lobatto, (see mfem::Quadrature1D).");
   args.AddOption(&win.data_state.save_coloring, "-sc", "--save-coloring",
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
   args.AddOption(&win.window_w, "-ww", "--window-width",
                  "Set the window width.");
   args.AddOption(&win.window_h, "-wh", "--window-height",
                  "Set the window height.");
   args.AddOption(&win.window_title, "-wt", "--window-title",
                  "Set the window title.");
   args.AddOption(&win.headless, "-hl", "--headless",
                  "-no-hl", "--no-headless",
                  "Start headless (no GUI) visualization.");
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
      input |= INPUT_MESH;
   }
   if (sol_file != string_none)
   {
      input |= INPUT_SCALAR_SOL;
   }
   if (vec_sol_file != string_none)
   {
      sol_file = vec_sol_file;
      input |= INPUT_VECTOR_SOL;
   }
   if (gfunc_file != string_none)
   {
      sol_file = gfunc_file;
      input |= INPUT_GRID_FUNC;
   }
   if (qfunc_file != string_none)
   {
      sol_file = qfunc_file;
      input |= INPUT_QUAD_FUNC;
   }

   if (visit_coll != string_none)
   {
      mesh_file = visit_coll;
      input |= INPUT_DATA_COLL;
   }
   else if (sidre_coll != string_none)
   {
      mesh_file = sidre_coll;
      input |= INPUT_DATA_COLL;
   }
   else if (fms_coll != string_none)
   {
      mesh_file = fms_coll;
      input |= INPUT_DATA_COLL;
   }
   else if (conduit_coll != string_none)
   {
      mesh_file = conduit_coll;
      input |= INPUT_DATA_COLL;
   }

   if (np > 0)
   {
      input |= INPUT_PARALLEL;
   }
   if (arg_keys != string_none)
   {
      win.data_state.keys = arg_keys;
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
      win.plot_caption = c_plot_caption;
   }
   if (legacy_gl_ctx == true)
   {
      SetLegacyGLOnly(legacy_gl_ctx);
   }
   SetUseHiDPI(enable_hidpi);

   // Load in palette file, if specified
   if (palette_file != string_none)
   {
      BasePalettes.Load(palette_file);
   }

   GLVisGeometryRefiner.SetType(geom_ref_type);

   string data_type;

   // check for saved stream file
   if (stream_file != string_none)
   {
      // backup the headless flag as the window is moved
      const bool headless = win.headless;

      if (!headless)
      {
         // Make sure the singleton object returned by GetMainThread() is
         // initialized from the main thread.
         GetMainThread();
      }

      Session stream_session(std::move(win));

      if (!stream_session.StartSavedSession(stream_file, !headless))
      {
         return 1;
      }

      if (!headless)
      {
         SDLMainLoop();
      }
      else
      {
         stream_session.WaitForSession();
      }
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
      cout << "Running script from file: " << script_file << endl;
      cout << "You may need to press <space> to execute the script steps." << endl;
      ScriptController::PlayScript(std::move(win), scr);
      return 0;
   }

   // turn off the server mode if other options are present
   if (input & ~INPUT_SERVER_MODE) { input &= ~INPUT_SERVER_MODE; }

   // print help for wrong input
   if (!(input == INPUT_SERVER_MODE
         || input == (INPUT_MESH)
         || input == (INPUT_DATA_COLL)
         || input == (INPUT_MESH | INPUT_SCALAR_SOL)
         || input == (INPUT_MESH | INPUT_VECTOR_SOL)
         || input == (INPUT_MESH | INPUT_PARALLEL)
         || input == (INPUT_MESH | INPUT_GRID_FUNC)
         || input == (INPUT_MESH | INPUT_GRID_FUNC | INPUT_PARALLEL)
         || input == (INPUT_DATA_COLL | INPUT_GRID_FUNC)
         || input == (INPUT_MESH | INPUT_QUAD_FUNC)
         || input == (INPUT_MESH | INPUT_QUAD_FUNC | INPUT_PARALLEL)
         || input == (INPUT_DATA_COLL | INPUT_QUAD_FUNC)))
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
   if (input == INPUT_SERVER_MODE)
   {
      // Make sure the singleton object returned by GetMainThread() is
      // initialized from the main thread.
      GetMainThread();

      // Run server in new thread
      std::thread serverThread{GLVisServer, portnum, save_stream,
                               win.data_state.fix_elem_orient,
                               win.data_state.save_coloring,
                               win.plot_caption};

      // Start SDL in main thread
      SDLMainLoop(true);

      serverThread.detach();
   }
   else  // input != 1, non-server mode
   {
      if (!(input & INPUT_DATA_COLL))
      {
         FileReader::FileType file_type;
         int component = -1;
         if (input & INPUT_SCALAR_SOL)
         {
            file_type = FileReader::FileType::SCALAR_SOL;
         }
         else if (input & INPUT_VECTOR_SOL)
         {
            file_type = FileReader::FileType::VECTOR_SOL;
         }
         else if (input & INPUT_GRID_FUNC)
         {
            file_type = FileReader::FileType::GRID_FUNC;
            component = gf_component;
         }
         else if (input & INPUT_QUAD_FUNC)
         {
            file_type = FileReader::FileType::QUAD_FUNC;
            component = qf_component;
         }
         else if (input & INPUT_MESH)
         {
            file_type = FileReader::FileType::MESH;
         }
         else
         {
            cerr << "Unknown input type" << endl;
            return 1;
         }

         FileReader reader(win.data_state, pad_digits);
         int ierr;
         if (input & INPUT_PARALLEL)
         {
            ierr = reader.ReadParallel(np, file_type, mesh_file, sol_file, component);
         }
         else
         {
            ierr = reader.ReadSerial(file_type, mesh_file, sol_file, component);
         }
         if (ierr) { exit(ierr); }
      }
      else
      {
         DataCollectionReader::CollType coll_type;
         int component;
         bool quad_func, mesh_only = false;
         if (visit_coll != string_none)
         {
            coll_type = DataCollectionReader::CollType::VISIT;
         }
         else if (sidre_coll != string_none)
         {
            coll_type = DataCollectionReader::CollType::SIDRE;
         }
         else if (fms_coll != string_none)
         {
            coll_type = DataCollectionReader::CollType::FMS;
         }
         else if (conduit_coll != string_none)
         {
            coll_type = DataCollectionReader::CollType::CONDUIT;
         }

         if (input & INPUT_GRID_FUNC)
         {
            quad_func = false;
            component = gf_component;
         }
         else if (input & INPUT_QUAD_FUNC)
         {
            quad_func = true;
            component = qf_component;
         }
         else
         {
            mesh_only = true;
         }

         DataCollectionReader reader(win.data_state);
         reader.SetPadDigits(pad_digits);
         if (dc_protocol != string_default)
         {
            reader.SetProtocol(dc_protocol.c_str());
         }

         int ierr;
         if (mesh_only)
         {
            ierr = reader.ReadSerial(coll_type, mesh_file, dc_cycle);
         }
         else
         {
            ierr = reader.ReadSerial(coll_type, mesh_file, dc_cycle, sol_file,
                                     quad_func, component);
         }
         if (ierr) { exit(ierr); }
      }

      // Make sure the singleton object returned by GetMainThread() is
      // initialized from the main thread.
      GetMainThread();

      Session single_session(std::move(win));
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
      "   glvis -np <#proc> -m <mesh_prefix> [-g <grid_function_prefix>]\n"
      "Visualize mesh and quadrature function:\n"
      "   glvis -m <mesh_file> -q <quadrature_function_file> [-qc <component>]\n"
      "Visualize parallel mesh and quadrature function:\n"
      "   glvis -np <#proc> -m <mesh_prefix> [-q <quadrature_function_prefix>]\n\n"
      "All Options:\n";
}
