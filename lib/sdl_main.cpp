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

#include <future>
#include <thread>
#include <climits>

#include "sdl_main.hpp"
#include "sdl_helper.hpp"
#include "logo.hpp"

#ifdef SDL_VIDEO_DRIVER_COCOA
#include "sdl_mac.hpp"
#endif

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

extern int GetMultisample();
extern bool wndUseHiDPI;

struct SdlMainThread::CreateWindowCmd
{
   SdlWindow* wnd;
   std::string title;
   int x, y, w, h;
   bool legacy_gl_only;
   promise<Handle> out_handle;
};

struct SdlMainThread::SdlCtrlCommand
{
   SdlCmdType type {SdlCmdType::None};

   const Handle*            handle {nullptr};
   CreateWindowCmd*         cmd_create;
   Handle                   cmd_delete;
   string                   cmd_title;
   pair<int, int>           cmd_set_size;
   pair<int, int>           cmd_set_position;
   // Promise object for signaling completion of the command on the main thread
   promise<void>            finished;
};

SdlMainThread::SdlMainThread()
{
   SDL_version sdl_ver;
   SDL_GetVersion(&sdl_ver);
   PRINT_DEBUG("Using SDL " << (int)sdl_ver.major << "." << (int)sdl_ver.minor
               << "." << (int)sdl_ver.patch << std::endl);

#ifdef SDL_VIDEO_DRIVER_COCOA
   if (SDL_VERSIONNUM(sdl_ver.major, sdl_ver.minor, sdl_ver.patch)
       < SDL_VERSIONNUM(2, 0, 14))
   {
      std::cerr << "Warning: your current version of SDL ("
                << (int)sdl_ver.major << "." << (int)sdl_ver.minor << "."
                << (int)sdl_ver.patch
                << ") may be unsupported in a future version of GLVis on macOS."
                << std::endl;
      std::cerr << "If possible, upgrade to SDL version 2.0.14 or newer."
                << std::endl;
   }
#endif

   if (!SDL_WasInit(SDL_INIT_VIDEO | SDL_INIT_EVENTS))
   {
      if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0)
      {
         std::cerr << "FATAL: Failed to initialize SDL: " << SDL_GetError() << endl;
         return;
      }
      SDL_EnableScreenSaver();
   }

   sdl_init = true;
#ifdef __EMSCRIPTEN__
   SetSingleThread();
#endif
}

SdlMainThread::~SdlMainThread()
{
   SDL_Quit();
}

void SdlMainThread::MainLoop(bool server)
{
   if (!sdl_init) { return; }

   this->server_mode = server;
   if (server_mode)
   {
      // Create a hidden window for dock icon and suppressing SDL_QUIT on last
      // vis window closing
      int flags = SDL_WINDOW_HIDDEN
                  | SDL_WINDOW_BORDERLESS
                  | SDL_WINDOW_ALLOW_HIGHDPI;
      SDL_Window* bg_wnd_handle = SDL_CreateWindow("GLVis",
                                                   SDL_WINDOWPOS_CENTERED,
                                                   SDL_WINDOWPOS_CENTERED,
                                                   1,
                                                   1,
                                                   flags);
      setWindowIcon(bg_wnd_handle);
      bg_wnd.reset(bg_wnd_handle);
   }
   while (1)
   {
      // Process all pending window commands
      vector<SdlCtrlCommand> pending_cmds;
      {
         lock_guard<mutex> cmd_lock{window_cmd_mtx};
         pending_cmds = std::move(window_cmds);
      }

      for (SdlCtrlCommand& cmd : pending_cmds)
      {
         handleWindowCmdImpl(cmd);
      }
      // Check if all windows have been closed, and we're not in server mode.
      // Exiting via shortcut 'q' won't generate an SDL_QUIT event that can be
      // handled in DispatchSDLEvents().
      if (num_windows == 0 && !server_mode)
      {
         terminating = true;
      }

      DispatchSDLEvents();

      // At this point, we're ready to start cleaning up
      if (terminating)
      {
         // Signal worker threads of all remaining open windows to terminate
         for (auto wnds : hwnd_to_window)
         {
            SDL_Event sdl_quit;
            sdl_quit.type = SDL_QUIT;
            wnds.second->queueEvents({sdl_quit});
         }
         while (num_windows > 0)
         {
            // Process pending destroy window commands
            vector<SdlCtrlCommand> pending_destroys;
            {
               lock_guard<mutex> cmd_lock{window_cmd_mtx};
               pending_destroys = std::move(window_cmds);
            }
            for (SdlCtrlCommand& cmd : pending_destroys)
            {
               if (cmd.type == SdlCmdType::Delete)
               {
                  handleWindowCmdImpl(cmd);
               }
            }
         }
         // Exit main loop
         break;
      }

      // Wait for the next event (without consuming CPU cycles, if possible)
      // See also: SdlWindow::signalLoop()
      if (platform)
      {
         platform->WaitEvent();
      }
      else
      {
         if (!SDL_PollEvent(nullptr))
         {
            std::this_thread::sleep_for(std::chrono::milliseconds(8));
         }
      }
   }
}

void SdlMainThread::DispatchSDLEvents()
{
   SDL_Event e;
   while (SDL_PollEvent(&e))
   {
      unsigned int destWindow = UINT_MAX;
      bool sendToAll = false;
      switch (e.type)
      {
         case SDL_QUIT:
            terminating = true;
            break;
         case SDL_WINDOWEVENT:
            destWindow = getWindowID(e.window);
            if (server_mode && (destWindow == SDL_GetWindowID(bg_wnd.get())))
            {
               handleBackgroundWindowEvent(e.window);
            }
            break;
         case SDL_FINGERDOWN:
            fingers.insert(e.tfinger.fingerId);
            if (fingers.size() >= 2) { disable_mouse = true; }
            break;
         case SDL_FINGERUP:
            fingers.erase(e.tfinger.fingerId);
            if (fingers.size() < 2) { disable_mouse = false; }
            break;
         case SDL_MULTIGESTURE:
            // No window ID is provided with multigesture events. We'll
            // have to use focus status within the windows themselves in
            // order to resolve the correct target.
            sendToAll = true;
            break;
         case SDL_KEYDOWN:
         case SDL_KEYUP:
            destWindow = getWindowID(e.key);
            break;
         case SDL_TEXTINPUT:
            destWindow = getWindowID(e.text);
            break;
         case SDL_MOUSEMOTION:
            if (!disable_mouse) { destWindow = getWindowID(e.motion); }
            break;
         case SDL_MOUSEBUTTONDOWN:
         case SDL_MOUSEBUTTONUP:
            destWindow = getWindowID(e.button);
            break;
      }
      if (sendToAll == true)
      {
         for (auto wnds : hwnd_to_window)
         {
            int windowId = wnds.first;
            wnd_events[windowId].emplace_back(e);
         }
      }
      else if (destWindow != UINT_MAX)
      {
         wnd_events[destWindow].emplace_back(e);
      }
   }
   // Send events to window worker threads
   for (auto wnds : hwnd_to_window)
   {
      int windowId = wnds.first;
      SdlWindow* wnd = wnds.second;
      if (!wnd_events[windowId].empty())
      {
         wnd->queueEvents(std::move(wnd_events[windowId]));
      }
   }
}

void SdlMainThread::handleBackgroundWindowEvent(SDL_WindowEvent e)
{
   // Background window should always stay hidden.
   if (e.event == SDL_WINDOWEVENT_SHOWN)
   {
      SDL_HideWindow(bg_wnd.get());
   }
}

SdlMainThread::Handle SdlMainThread::GetHandle(SdlWindow* wnd,
                                               const std::string& title,
                                               int x, int y, int w, int h, bool legacyGlOnly)
{
   if (!sdl_init)
   {
      return Handle{};
   }

   CreateWindowCmd cmd_create = { wnd, title, x, y, w, h, legacyGlOnly, {} };
   future<Handle> res_handle = cmd_create.out_handle.get_future();

   SdlCtrlCommand main_thread_cmd;
   main_thread_cmd.type = SdlCmdType::Create;
   main_thread_cmd.cmd_create = &cmd_create;

   queueWindowEvent(std::move(main_thread_cmd));

   Handle out_hnd = res_handle.get();
   if (out_hnd.isInitialized())
   {
      // Make the OpenGL context current on the worker thread.
      // Since SDL calls aren't guaranteed to be thread-safe, we guard
      // the call to SDL_GL_MakeCurrent.
      lock_guard<mutex> ctx_lock{gl_ctx_mtx};
#ifdef SDL_VIDEO_DRIVER_COCOA
      // TODO: Temporary workaround - after merge, everyone should update to
      // latest SDL
      SdlCocoaPlatform* mac_platform
         = dynamic_cast<SdlCocoaPlatform*>(platform.get());
      if (mac_platform && mac_platform->UseThreadWorkaround())
      {
         int wnd_id = SDL_GetWindowID(out_hnd.hwnd);
         mac_platform->SetCurrentContext(wnd_id);
      }
      else
      {
         SDL_GL_MakeCurrent(out_hnd.hwnd, out_hnd.gl_ctx);
      }
#else
      SDL_GL_MakeCurrent(out_hnd.hwnd, out_hnd.gl_ctx);
#endif
   }
   return out_hnd;
}

void SdlMainThread::DeleteHandle(Handle to_delete)
{
   if (to_delete.isInitialized())
   {
      {
         lock_guard<mutex> ctx_lock{gl_ctx_mtx};
         // Unbinding the context before deletion seems to resolve some issues
         // with thread-local storage on Wayland/EGL.
         SDL_GL_MakeCurrent(to_delete.hwnd, nullptr);
      }
      SdlCtrlCommand main_thread_cmd;
      main_thread_cmd.type = SdlCmdType::Delete;
      main_thread_cmd.cmd_delete = std::move(to_delete);

      queueWindowEvent(std::move(main_thread_cmd));
   }
}

void SdlMainThread::SetWindowTitle(const Handle& handle, std::string title)
{
   if (handle.isInitialized())
   {
      SdlCtrlCommand main_thread_cmd;
      main_thread_cmd.type = SdlCmdType::SetTitle;
      main_thread_cmd.handle = &handle;
      main_thread_cmd.cmd_title = std::move(title);

      queueWindowEvent(std::move(main_thread_cmd));
   }
}

void SdlMainThread::SetWindowSize(const Handle& handle, int w, int h)
{
   if (handle.isInitialized())
   {
      SdlCtrlCommand main_thread_cmd;
      main_thread_cmd.type = SdlCmdType::SetSize;
      main_thread_cmd.handle = &handle;
      main_thread_cmd.cmd_set_size = {w, h};

      queueWindowEvent(std::move(main_thread_cmd), true);
   }
}

void SdlMainThread::SetWindowPosition(const Handle& handle, int x, int y)
{
   if (handle.isInitialized())
   {
      SdlCtrlCommand main_thread_cmd;
      main_thread_cmd.type = SdlCmdType::SetPosition;
      main_thread_cmd.handle = &handle;
      main_thread_cmd.cmd_set_position = {x, y};

      queueWindowEvent(std::move(main_thread_cmd), true);
   }
}

void SdlMainThread::queueWindowEvent(SdlCtrlCommand cmd, bool sync)
{
   future<void> wait_complete;
   if (sdl_multithread)
   {
      if (sync)
      {
         wait_complete = cmd.finished.get_future();
      }
      // queue up our event
      {
         lock_guard<mutex> req_lock{window_cmd_mtx};
         window_cmds.emplace_back(std::move(cmd));
      }
      // wake up the main thread to handle our event
      SendEvent();

      if (sync) { wait_complete.get(); }
   }
   else
   {
      // call the underlying SDL command immediately
      handleWindowCmdImpl(cmd);
   }
}

void SdlMainThread::handleWindowCmdImpl(SdlCtrlCommand& cmd)
{
   switch (cmd.type)
   {
      case SdlCmdType::Create:
         createWindowImpl(*cmd.cmd_create);
         break;
      case SdlCmdType::Delete:
         if (cmd.cmd_delete.isInitialized())
         {
            Handle to_delete = std::move(cmd.cmd_delete);
            if (platform)
            {
               platform->UnregisterWindow(to_delete.hwnd);
            }
            int wnd_id = SDL_GetWindowID(to_delete.hwnd);
            hwnd_to_window.erase(wnd_id);
            wnd_events.erase(wnd_id);
            num_windows--;
         }
         break;
      case SdlCmdType::SetTitle:
         SDL_SetWindowTitle(cmd.handle->hwnd, cmd.cmd_title.c_str());
         break;
      case SdlCmdType::SetSize:
         SDL_SetWindowSize(cmd.handle->hwnd,
                           cmd.cmd_set_size.first,
                           cmd.cmd_set_size.second);
         break;
      case SdlCmdType::SetPosition:
         SDL_SetWindowPosition(cmd.handle->hwnd,
                               cmd.cmd_set_position.first,
                               cmd.cmd_set_position.second + title_height_offset);
         break;
      default:
         cerr << "Error in main thread: unknown window control command.\n";
         break;
   }
   // Signal completion of the command, in case worker thread is waiting.
   cmd.finished.set_value();
}

void SdlMainThread::setWindowIcon(SDL_Window* hwnd)
{
#ifdef GLVIS_USE_LOGO
   const int PixelStride = 4;
   int stride = (int) sqrt(logo_rgba_len / PixelStride);
   if (unsigned(stride * stride * PixelStride) != logo_rgba_len)
   {
      cerr << "Unable to set window logo: icon size not square" << endl;
   }
   else
   {
      SDL_Surface* iconSurf =
         SDL_CreateRGBSurfaceFrom(logo_rgba,
                                  stride, stride,       // height, width
                                  8 * PixelStride,      // depth
                                  stride * PixelStride, // pitch
                                  0x000000FF,
                                  0x0000FF00,
                                  0x00FF0000,
                                  0xFF000000);
      if (iconSurf)
      {
         SDL_SetWindowIcon(hwnd, iconSurf);
         SDL_FreeSurface(iconSurf);
      }
      else
      {
         PRINT_DEBUG("Unable to set window logo: " << SDL_GetError() << endl);
      }
   }
#endif // GLVIS_USE_LOGO
}

void SdlMainThread::probeGLContextSupport(bool legacyGlOnly)
{
   Uint32 win_flags_hidden = SDL_WINDOW_OPENGL | SDL_WINDOW_HIDDEN;
#ifndef __EMSCRIPTEN__
   if (!legacyGlOnly)
   {
      // Try and probe for a core/compatibility context. Needed for Mac OS X,
      // which will only support OpenGL 2.1 if you don't create a core context.
      PRINT_DEBUG("Testing if OpenGL core profile window can be created..." << flush);
      // Try to create core profile context first
      SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
      {
         Handle testCore("",
                         SDL_WINDOWPOS_UNDEFINED,
                         SDL_WINDOWPOS_UNDEFINED,
                         100, 100, win_flags_hidden);
         if (testCore.isInitialized())
         {
            PRINT_DEBUG("success!" << endl);
            return;
         }
      }

      PRINT_DEBUG("failed." << endl);
      PRINT_DEBUG("Testing if OpenGL compatibility profile window can be created..."
                  << flush);
      SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK,
                          SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);
      {
         Handle testCompat("",
                           SDL_WINDOWPOS_UNDEFINED,
                           SDL_WINDOWPOS_UNDEFINED,
                           100, 100, win_flags_hidden);
         if (testCompat.isInitialized())
         {
            PRINT_DEBUG("success!" << endl);
            return;
         }
      }
      PRINT_DEBUG("failed." << endl);
   }
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, 0);
   if (GetMultisample() > 0)
   {
      SDL_GL_SetAttribute( SDL_GL_MULTISAMPLEBUFFERS, 1);
      SDL_GL_SetAttribute( SDL_GL_MULTISAMPLESAMPLES, GetMultisample());
      PRINT_DEBUG("Testing OpenGL with default profile and antialiasing..." << endl);
      {
         Handle testMsaa("",
                         SDL_WINDOWPOS_UNDEFINED,
                         SDL_WINDOWPOS_UNDEFINED,
                         100, 100, win_flags_hidden);
         if (testMsaa.isInitialized())
         {
            PRINT_DEBUG("success!" << endl);
            return;
         }
      }
      SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 0);
      SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 0);
      PRINT_DEBUG("failed." << endl);
   }
   PRINT_DEBUG("Opening OpenGL window with no flags.");
#else
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK,
                       SDL_GL_CONTEXT_PROFILE_ES);
   PRINT_DEBUG("Testing for WebGL 2.0 context availability..." << flush);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
   {
      Handle testWebGl2("",
                        SDL_WINDOWPOS_UNDEFINED,
                        SDL_WINDOWPOS_UNDEFINED,
                        100, 100, win_flags_hidden);
      if (testWebGl2.isInitialized())
      {
         PRINT_DEBUG("success!" << endl);
         return;
      }
   }
   PRINT_DEBUG("failed." << endl);
   PRINT_DEBUG("Testing for WebGL 1.0 context availability..." << flush);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
   // Note that WebGL 2 support controllable antialiasing, while WebGL 1 only
   // supports requesting multisampling at context creation time. The browser
   // is free to ignore this flag.
   SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
   SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, GetMultisample());
   {
      Handle testWebGl("",
                       SDL_WINDOWPOS_UNDEFINED,
                       SDL_WINDOWPOS_UNDEFINED,
                       100, 100, win_flags_hidden);
      if (testWebGl.isInitialized())
      {
         PRINT_DEBUG("success!" << endl);
         return;
      }
   }
#endif
}

const int default_dpi = 72;
void SdlMainThread::getDpi(const Handle& handle, int& w, int& h)
{
   w = default_dpi;
   h = default_dpi;
   if (!handle.isInitialized())
   {
      PRINT_DEBUG("warning: unable to get dpi: handle is null" << endl);
      return;
   }
   int disp = SDL_GetWindowDisplayIndex(handle.hwnd);
   if (disp < 0)
   {
      PRINT_DEBUG("warning: problem getting display index: " << SDL_GetError()
                  << endl);
      PRINT_DEBUG("returning default dpi of " << default_dpi << endl);
      return;
   }

   float f_w, f_h;
   if (SDL_GetDisplayDPI(disp, NULL, &f_w, &f_h))
   {
      PRINT_DEBUG("warning: problem getting dpi, setting to " << default_dpi
                  << ": " << SDL_GetError() << endl);
   }
   else
   {
      PRINT_DEBUG("Screen DPI: w = " << f_w << " ppi, h = " << f_h << " ppi");
      w = f_w + 0.5f;
      h = f_h + 0.5f;
      PRINT_DEBUG(" (" << w << " x " << h << ')' << endl);
   }
}

void SdlMainThread::createWindowImpl(CreateWindowCmd& cmd)
{
   Uint32 win_flags = SDL_WINDOW_OPENGL;
   // Hide window until we adjust its size for high-dpi displays
   win_flags |= SDL_WINDOW_HIDDEN;
#ifndef __EMSCRIPTEN__
   win_flags |= SDL_WINDOW_RESIZABLE;
#endif
   if (wndUseHiDPI)
   {
      win_flags |= SDL_WINDOW_ALLOW_HIGHDPI;
   }
   SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1);
   SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24);
   probeGLContextSupport(cmd.legacy_gl_only);
   Handle new_handle(cmd.title,
                     cmd.x, cmd.y + title_height_offset,
                     cmd.w, cmd.h, win_flags);

   // at this point, window should be up
   if (!new_handle.isInitialized())
   {
      PRINT_DEBUG("failed." << endl);
      cerr << "FATAL: window and/or OpenGL context creation failed." << endl;
      // If not in server mode, this was the only window we were going to open.
      if (!server_mode) { terminating = true; }
      return;
   }
   else
   {
      PRINT_DEBUG("Handle for window created." << endl);
   }

#ifndef __EMSCRIPTEN__
   SDL_GL_SetSwapInterval(0);
   glEnable(GL_DEBUG_OUTPUT);
#endif

   // Register window internally in the main thread so it can receive events
   int wnd_id = SDL_GetWindowID(new_handle.hwnd);
   hwnd_to_window[wnd_id] = cmd.wnd;

   {
      lock_guard<mutex> platform_lk{event_mtx};
      if (!platform)
      {
         platform = SdlNativePlatform::Create(new_handle.hwnd);
         try_create_platform = true;
      }
      event_cv.notify_all();
   }
   if (platform)
   {
      platform->RegisterWindow(new_handle.hwnd);
   }

   setWindowIcon(new_handle.hwnd);

   if (num_windows == -1)
   {
      // Get the window titlebar height after creating the first window.
      // Window coordinates are based on the client area, so placing a window
      // at (0, 0) will hide the title bar on Windows.
      SDL_GetWindowBordersSize(new_handle.hwnd, &title_height_offset, nullptr,
                               nullptr, nullptr);
      SDL_SetWindowPosition(new_handle.hwnd, cmd.x, cmd.y + title_height_offset);
   }

   // Detect if we are using a high-dpi display and resize the window unless it
   // was already resized by SDL's underlying backend.
   {
      SdlWindow* wnd = cmd.wnd;
      int scr_w, scr_h, pix_w, pix_h, wdpi, hdpi;
      // SDL_GetWindowSize() -- size in "screen coordinates"
      SDL_GetWindowSize(new_handle.hwnd, &scr_w, &scr_h);
      // SDL_GL_GetDrawableSize() -- size in pixels
      SDL_GL_GetDrawableSize(new_handle.hwnd, &pix_w, &pix_h);
      wnd->high_dpi = false;
      wnd->pixel_scale_x = wnd->pixel_scale_y = 1.0f;
      float sdl_pixel_scale_x = 1.0f, sdl_pixel_scale_y = 1.0f;
      // If "screen" and "pixel" sizes are different, assume high-dpi and no
      // need to scale the window.
      if (scr_w == pix_w && scr_h == pix_h)
      {
         getDpi(new_handle, wdpi, hdpi);
         if (std::max(wdpi, hdpi) >= SdlWindow::high_dpi_threshold)
         {
            wnd->high_dpi = true;
            wnd->pixel_scale_x = wnd->pixel_scale_y = 2.0f;
            // the following two calls use 'pixel_scale_*'
            SDL_SetWindowSize(new_handle.hwnd,
                              wnd->pixel_scale_x*cmd.w,
                              wnd->pixel_scale_y*cmd.h);
            bool uc_x = SDL_WINDOWPOS_ISUNDEFINED(cmd.x) ||
                        SDL_WINDOWPOS_ISCENTERED(cmd.x);
            bool uc_y = SDL_WINDOWPOS_ISUNDEFINED(cmd.y) ||
                        SDL_WINDOWPOS_ISCENTERED(cmd.y);
            SDL_SetWindowPosition(new_handle.hwnd,
                                  uc_x ? cmd.x : wnd->pixel_scale_x*cmd.x,
                                  uc_y ? cmd.y : wnd->pixel_scale_y*cmd.y);
         }
      }
      else
      {
         wnd->high_dpi = true;
         // keep 'pixel_scale_*' = 1, scaling is done inside SDL
         sdl_pixel_scale_x = float(pix_w)/scr_w;
         sdl_pixel_scale_y = float(pix_h)/scr_h;
      }
      if (wnd->high_dpi)
      {
         cout << "High-dpi display detected: using window scaling: "
              << sdl_pixel_scale_x*wnd->pixel_scale_x << " x "
              << sdl_pixel_scale_y*wnd->pixel_scale_y << endl;
      }
   }

   // Unset GL context in this thread
   {
      lock_guard<mutex> ctx_lock{gl_ctx_mtx};
#ifdef SDL_VIDEO_DRIVER_COCOA
      // TODO: Temporary workaround - after merge, everyone should update to
      // latest SDL
      SdlCocoaPlatform* mac_platform
         = dynamic_cast<SdlCocoaPlatform*>(platform.get());
      if (mac_platform && mac_platform->UseThreadWorkaround())
      {
         mac_platform->ClearCurrentContext(wnd_id);
      }
      else
      {
         SDL_GL_MakeCurrent(new_handle.hwnd, nullptr);
      }
#else
      SDL_GL_MakeCurrent(new_handle.hwnd, nullptr);
#endif
   }

   SDL_ShowWindow(new_handle.hwnd);
   SDL_RaiseWindow(new_handle.hwnd);
   if (num_windows == -1)
   {
      num_windows = 0;
   }
   num_windows++;
   cmd.out_handle.set_value(std::move(new_handle));
}
