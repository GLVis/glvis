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

#ifndef GLVIS_EGL_HPP
#define GLVIS_EGL_HPP
#ifdef GLVIS_USE_EGL

#include "../glwindow.hpp"

#include <EGL/egl.h>

#include <condition_variable>
#include <mutex>
#include <deque>
#include <list>

class EglWindow : public GLWindow
{
public:
   struct Handle
   {
      EGLSurface surf{EGL_NO_SURFACE};
      EGLContext ctx{EGL_NO_CONTEXT};
      EGLConfig eglCfg{};

      bool isInitialized()
      {
         return surf != EGL_NO_SURFACE && ctx != EGL_NO_CONTEXT;
      }
   };

private:
   Handle handle;

   bool running{false};

   bool is_multithreaded{true};
   bool call_idle_func{false};

   enum class EventType
   {
      Keydown,
      Screenshot,
      Quit,
   };
   struct Event
   {
      EventType type;
      union Events
      {
         struct Keydown
         {
            SDL_Keycode k;
            SDL_Keymod m;
         } keydown;

         struct Screenshot
         {
            bool convert;
         } screenshot;

         struct Quit { } quit;
      } event;
   };

   std::string screenshot_filename;

   std::condition_variable events_available;
   std::mutex event_mutex;
   std::deque<Event> waiting_events;

   void queueEvents(std::vector<Event> events);

public:
   EglWindow();
   ~EglWindow();

   /** @brief Creates a new OpenGL window. Returns false if EGL or OpenGL
       initialization fails. */
   bool createWindow(const char *title, int x, int y, int w, int h,
                     bool legacyGlOnly) override;

   /// Runs the window loop.
   void mainLoop() override;
   void mainIter() override;

   void signalLoop() override;

   void getWindowSize(int& w, int& h) const override { getGLDrawSize(w, h); }
   void getGLDrawSize(int& w, int& h) const override;

   // use the default values of SdlWindow, LoadFont() ignores them anyway
   void getDpi(int& wdpi, int& hdpi) const override { wdpi = hdpi = 72; }
   bool isHighDpi() const override;

   void setWindowSize(int w, int h) override;

   bool isWindowInitialized() const override { return handle.surf != EGL_NO_SURFACE; }
   bool isGlInitialized() const override { return handle.ctx != EGL_NO_CONTEXT; }

   void signalKeyDown(SDL_Keycode k, SDL_Keymod m = KMOD_NONE) override;
   void signalQuit() override;
   // as there is no swap, switch to updated state right away
   void signalSwap() override { wnd_state = RenderState::Updated; }

   // used in Screenshot, as there is no swapping, the single buffer is always
   // up to date and can be read directly
   bool isSwapPending() const override { return true; }

   void screenshot(std::string filename, bool convert = false) override;
};

#endif //GLVIS_USE_EGL
#endif //GLVIS_EGL_HPP