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

#ifndef GLVIS_EGL_MAIN_HPP
#define GLVIS_EGL_MAIN_HPP

#if defined(GLVIS_USE_EGL) || defined(GLVIS_USE_CGL)

#include <list>
#include <future>

#include "egl.hpp"

class EglMainThread : public MainThread
{
   using Handle = EglWindow::Handle;
#ifdef GLVIS_USE_EGL
   EGLDisplay disp {EGL_NO_DISPLAY};
#endif
   bool server_mode {false};

   std::list<EglWindow*> windows;
   int num_windows {-1};

   struct CreateWndCmd;
   struct ResizeWndCmd;
   struct DeleteWndCmd;

   enum class CtrlCmdType
   {
      Create,
      Resize,
      Delete,
      Terminate,
   };

   struct CtrlCmd
   {
      CtrlCmdType type;
      union
      {
         CreateWndCmd *create_cmd;
         ResizeWndCmd *resize_cmd;
         DeleteWndCmd *delete_cmd;
      };

      std::promise<void> finished;
   };

   std::condition_variable events_available;
   std::mutex window_cmd_mtx;
   std::deque<CtrlCmd> window_cmds;

   bool CreateWndImpl(CreateWndCmd &cmd);
   bool ResizeWndImpl(ResizeWndCmd &cmd);
   bool DeleteWndImpl(DeleteWndCmd &cmd);
   void QueueWndCmd(CtrlCmd cmd, bool sync);

   static void InterruptHandler(int param);

public:

   EglMainThread();
   ~EglMainThread();

   static EglMainThread& Get();
#ifdef GLVIS_USE_EGL
   EGLDisplay GetDisplay() const { return disp; }
#endif

   Handle CreateWindow(EglWindow *caller, int w, int h, bool legacy_gl);
   void ResizeWindow(Handle &hnd, int w, int h);
   void DeleteWindow(EglWindow *caller, Handle &hnd);
   void Terminate();

   void MainLoop(bool server = false) override;
};

#endif // GLVIS_USE_EGL
#endif // GLVIS_EGL_MAIN_HPP
