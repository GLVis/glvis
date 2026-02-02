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

#ifndef GLVIS_EGL_HPP
#define GLVIS_EGL_HPP

#if defined(GLVIS_USE_EGL) || defined(GLVIS_USE_CGL)

#include "../glwindow.hpp"

#ifdef GLVIS_USE_EGL
#include <EGL/egl.h>
#endif

#ifdef GLVIS_USE_CGL
#define GL_SILENCE_DEPRECATION // CGL has been deprecated since MacOS 10.14
#include <OpenGL/OpenGL.h>
#endif

#include <condition_variable>
#include <memory>
#include <mutex>
#include <deque>

class EglWindow : public GLWindow
{
public:
   struct Handle
   {
#ifdef GLVIS_USE_EGL
      EGLSurface surf {EGL_NO_SURFACE};
      EGLContext ctx{EGL_NO_CONTEXT};
      EGLConfig eglCfg{};
#endif
#ifdef GLVIS_USE_CGL
      GLuint buf_frame {}, buf_color {}, buf_depth {};
      CGLPixelFormatObj pix {};

      using CGLContextDeleter = void(*)(CGLContextObj);
      std::unique_ptr<_CGLContextObject, CGLContextDeleter> ctx
      {nullptr, [](CGLContextObj ctx) { CGLDestroyContext(ctx); }};
#endif

      bool isInitialized()
      {
#ifdef GLVIS_USE_EGL
         return surf != EGL_NO_SURFACE && ctx != EGL_NO_CONTEXT;
#endif
#ifdef GLVIS_USE_CGL
         return ctx != nullptr;
#endif
      }
   };

private:
   Handle handle;

   bool running {false};
   bool is_multithreaded {true};
   bool call_idle_func {false};

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
   EglWindow() = default;
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

#if defined(GLVIS_USE_EGL)
   bool isWindowInitialized() const override { return handle.surf != EGL_NO_SURFACE; }
   bool isGlInitialized() const override { return handle.ctx != EGL_NO_CONTEXT; }
#elif defined(GLVIS_USE_CGL)
   bool isWindowInitialized() const override { return isGlInitialized(); }
   bool isGlInitialized() const override { return handle.ctx != nullptr; }
#endif

   void signalKeyDown(SDL_Keycode k, SDL_Keymod m = KMOD_NONE) override;
   void signalQuit() override;
   // as there is no swap, switch to updated state right away
   void signalSwap() override { wnd_state = RenderState::Updated; }

   // used in Screenshot, as there is no swapping, the single buffer is always
   // up to date and can be read directly
   bool isSwapPending() const override { return true; }

   void screenshot(std::string filename, bool convert = false) override;
};

#endif // GLVIS_USE_EGL || GLVIS_USE_CGL
#endif // GLVIS_EGL_HPP
