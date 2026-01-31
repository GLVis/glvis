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

#ifndef GLVIS_SCRIPT_CONTROLLER_HPP
#define GLVIS_SCRIPT_CONTROLLER_HPP

#include <iostream>
#include <memory>

#include "window.hpp"

extern const char *string_none;
extern const char *string_default;

class ScriptController
{
   Window win;

   std::string dc_protocol = string_default;
   int dc_cycle = 0;

   std::istream *script = NULL;
   int scr_running = 0;
   int scr_level = 0;
   std::unique_ptr<mfem::Vector> init_nodes;
   double scr_min_val, scr_max_val;

   static int ScriptReadSolution(std::istream &scr, DataState &state);
   static int ScriptReadQuadrature(std::istream &scr, DataState &state);
   static int ScriptReadParSolution(std::istream &scr, DataState &state);
   static int ScriptReadParQuadrature(std::istream &scr, DataState &state);
   int ScriptReadDisplMesh(std::istream &scr, DataState &state);
   int ScriptReadDataColl(std::istream &scr, DataState &state,
                          bool mesh_only = true, bool quad = false);

   //key handlers using thread-local singleton
   static thread_local ScriptController *script_ctrl;
   static void ScriptIdleFunc();
   static void ScriptControl();

   static void PrintCommands();
   bool ExecuteScriptCommand();

public:
   ScriptController(Window win_) : win(std::move(win_)) { }

   static void PlayScript(Window win, std::istream &scr);
};

#endif // GLVIS_SCRIPT_CONTROLLER_HPP
