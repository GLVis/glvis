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

#include "threads.hpp" // IWYU pragma: keep for GLVisCommand
#include "aux_vis.hpp" // SetUseHiDPI

#include "session.hpp"

#include <fem/geom.hpp>

extern thread_local mfem::GeometryRefiner GLVisGeometryRefiner;

int GLVisStreamSession(const bool fix_elem_orient,
                       const bool save_coloring,
                       const bool headless,
                       const std::string &plot_caption,
                       StreamCollection &&streams)
{
   const int geom_ref_type = mfem::Quadrature1D::ClosedUniform;
   const bool enable_hidpi = true;

   Window win;

   std::cout << std::endl
             << "       _/_/_/  _/      _/      _/  _/"          << std::endl
             << "    _/        _/      _/      _/        _/_/_/" << std::endl
             << "   _/  _/_/  _/      _/      _/  _/  _/_/"      << std::endl
             << "  _/    _/  _/        _/  _/    _/      _/_/"   << std::endl
             << "   _/_/_/  _/_/_/_/    _/      _/  _/_/_/"      << std::endl
             << std::endl;

   SetUseHiDPI(enable_hidpi);
   GLVisGeometryRefiner.SetType(geom_ref_type);

   GetMainThread(false);

   Session new_session(fix_elem_orient, save_coloring, plot_caption, headless);
   new_session.StartStreamSession(std::move(streams));

   std::vector<Session> current_sessions;
   current_sessions.emplace_back(std::move(new_session));

   MainThreadLoop(false, false);

   return EXIT_SUCCESS;
}
