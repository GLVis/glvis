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

#ifndef GLVIS_SDL_MAIN_HPP
#define GLVIS_SDL_MAIN_HPP

#include <set>
#include <condition_variable>

#include "sdl.hpp"

class SdlMainThread
{
private:
   using Handle = SdlWindow::Handle;
public:
   SdlMainThread();
   ~SdlMainThread();

   // If called, all SDL window operations are run immediately; this is for
   // running in a single-threaded mode.
   void SetSingleThread() { sdl_multithread = false; }

   bool SdlInitialized() const { return sdl_init; }

   // Handles all SDL operations that are expected to be handled on the main
   // SDL thread (i.e. events and window creation)
   void MainLoop(bool server_mode);

   // Dequeues all incoming events from SDL, and queues them up to their
   // matching windows. Intended to be called only in single-threaded mode.
   void DispatchSDLEvents();

   // Executed from a window worker thread. Returns a handle to a new window
   // and associated OpenGL context.
   Handle GetHandle(SdlWindow* wnd, const std::string& title,
                    int x, int y, int w, int h, bool legacyGlOnly);

   // Executed from a window worker thread. Deletes a handle to a window and
   // the corresponding OpenGL context.
   void DeleteHandle(Handle to_delete);

   // Issues a command on the main thread to set the window title.
   void SetWindowTitle(const Handle& handle, std::string title);

   // Issues a command on the main thread to set the window size.
   void SetWindowSize(const Handle& handle, int w, int h);

   // Issues a command on the main thread to set the window position.
   void SetWindowPosition(const Handle& handle, int x, int y);

   SdlNativePlatform* GetPlatform() const { return platform.get(); }

private:
   struct CreateWindowCmd;
   struct SdlCtrlCommand;

   enum class SdlCmdType
   {
      None,
      Create,
      Delete,
      SetTitle,
      SetSize,
      SetPosition
   };

   // Wakes up the main thread, if sleeping.
   void SendEvent()
   {
      std::unique_lock<std::mutex> platform_lk{event_mtx};
      event_cv.wait(platform_lk, [this]() { return try_create_platform; });
      if (platform)
      {
         platform->SendEvent();
      }
   }

   void queueWindowEvent(SdlCtrlCommand cmd, bool sync = false);

   template<typename T>
   Uint32 getWindowID(const T& eventStruct)
   {
      return eventStruct.windowID;
   }

   void setWindowIcon(SDL_Window* hwnd);

   void handleWindowCmdImpl(SdlCtrlCommand& cmd);

   // Setup the correct OpenGL context flags in SDL for when we actually open the
   // window.
   void probeGLContextSupport(bool legacyGlOnly);

   void getDpi(const Handle& handle, int& wdpi, int& hdpi);

   void createWindowImpl(CreateWindowCmd& cmd);

   void handleBackgroundWindowEvent(SDL_WindowEvent e);

   bool sdl_init {false};
   bool sdl_multithread {true};

   bool server_mode {false};

   // A flag indicating whether the main loop will *begin* terminating
   bool terminating {false};
   unique_ptr<SDL_Window, decltype(&SDL_DestroyWindow)> bg_wnd{nullptr, SDL_DestroyWindow};

   // -------------------------------------------------------------------------
   // Objects for handling passing of window control commands to the main event
   // loop.

   mutex window_cmd_mtx;
   vector<SdlCtrlCommand> window_cmds;

   int num_windows {-1}; // -1: waiting for window to be created

   // -------------------------------------------------------------------------
   // Objects for handling dispatching events from the main event loop to
   // worker threads.

   unordered_map<int, SdlWindow*> hwnd_to_window;
   unordered_map<int, vector<SDL_Event>> wnd_events;
   std::set<SDL_FingerID> fingers;
   bool disable_mouse {false};

   mutex gl_ctx_mtx;

   mutex event_mtx;
   condition_variable event_cv;
   bool try_create_platform{false};
   unique_ptr<SdlNativePlatform> platform;

   int title_height_offset {0};
};

#endif
