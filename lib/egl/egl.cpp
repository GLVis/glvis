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

#if defined(GLVIS_USE_EGL) || defined(GLVIS_USE_CGL)

#include "egl.hpp"
#include "egl_main.hpp"
#include "../aux_vis.hpp"
#include <iostream>
#include <vector>

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

using namespace std;

EglWindow::EglWindow()
{
}

EglWindow::~EglWindow()
{
   EglMainThread::Get().DeleteWindow(this, handle);
}

#ifdef GLVIS_USE_CGL
bool EglWindow::initGLEW(bool legacyGlOnly)
{
   return GLWindow::initGLEW(legacyGlOnly);
}
#endif


bool EglWindow::createWindow(const char *, int, int, int w, int h,
                             bool legacyGlOnly)
{
   handle = EglMainThread::Get().CreateWindow(this, w, h, legacyGlOnly);
   assert(glGetError() == GL_NO_ERROR);

   if (!handle.isInitialized())
   {
      return false;
   }

#ifdef GLVIS_USE_CGL
   return true;
#endif

#ifndef __EMSCRIPTEN__
   glEnable(GL_DEBUG_OUTPUT);
#endif

   PRINT_DEBUG("EGL context is ready" << endl);

   return initGLEW(legacyGlOnly);
}

void EglWindow::queueEvents(vector<Event> events)
{
   {
      lock_guard<mutex> evt_guard{event_mutex};
      waiting_events.insert(waiting_events.end(), events.begin(), events.end());
   }
   if (is_multithreaded)
   {
      events_available.notify_all();
   }
}

void EglWindow::mainLoop()
{
   running = true;
   while (running)
   {
      mainIter();
   }
}

void EglWindow::mainIter()
{
   bool sleep = false;
   bool events_pending = false;
   {
      lock_guard<mutex> evt_guard{event_mutex};
      events_pending = !waiting_events.empty();
   }
   if (events_pending)
   {
      do
      {
         Event e;
         // Fetch next event from the queue
         {
            lock_guard<mutex> evt_guard{event_mutex};
            e = waiting_events.front();
            waiting_events.pop_front();
            events_pending = !waiting_events.empty();
         }

         switch (e.type)
         {
            case EventType::Keydown:
               if (onKeyDown[e.event.keydown.k])
               {
                  onKeyDown[e.event.keydown.k](e.event.keydown.m);
                  recordKey(e.event.keydown.k, e.event.keydown.m);
               }
               break;
            case EventType::Screenshot:
               if (isExposePending())
               {
                  onExpose();
                  wnd_state = RenderState::Updated;
               }
               Screenshot(screenshot_filename.c_str(), e.event.screenshot.convert);
               break;
            case EventType::Quit:
               running = false;
               break;
         }
      }
      while (events_pending);
   }
   else if (onIdle)
   {
      sleep = onIdle();
   }
   else
   {
      // No actions performed this iteration.
      sleep = true;
   }

   if (isExposePending())
   {
      onExpose();
      wnd_state = RenderState::Updated;
   }
   else if (sleep)
   {
      unique_lock<mutex> event_lock{event_mutex};
      events_available.wait(event_lock, [this]()
      {
         // Sleep until events from WM or glvis_command can be handled
         return !waiting_events.empty() || call_idle_func;
      });
   }
}

void EglWindow::signalLoop()
{
   // Note: not executed from the main thread
   {
      lock_guard<mutex> evt_guard{event_mutex};
      call_idle_func = true;
   }
   events_available.notify_all();
}

void EglWindow::getGLDrawSize(int& w, int& h) const
{
#ifdef GLVIS_USE_EGL
   EGLint egl_w, egl_h;

   EGLDisplay disp = EglMainThread::Get().GetDisplay();
   eglQuerySurface(disp, handle.surf, EGL_WIDTH, &egl_w);
   eglQuerySurface(disp, handle.surf, EGL_HEIGHT, &egl_h);
   w = egl_w;
   h = egl_h;
#endif
#ifdef GLVIS_USE_CGL

   const auto err = glGetError();
   std::cout << "glGetError:" << err << std::endl;

   assert(glGetError() == GL_NO_ERROR);
   GLint cgl_w, cgl_h;
   static_assert(sizeof(GLint) == sizeof(int),
                 "CGL width size must be 4 bytes");
   // CGLGetRenderbufferParameter(handle, cgl_w, cgl_h);
   {
      // dbg("glGetError:{}", glGetError());
      assert(glGetError() == GL_NO_ERROR);
      // assert(handle.isInitialized());

      // Optional: Ensure context is current if not guaranteed elsewhere
      CGLError ctx_err = CGLSetCurrentContext(handle.ctx.get());
      if (ctx_err != kCGLNoError)
      {
         PRINT_DEBUG("❌ Cannot set CGL context as current");
      }
      // Bind the color renderbuffer (assuming that's what we want to query)
      glBindRenderbuffer(GL_RENDERBUFFER, handle.buf_color);
      if (glGetError() != GL_NO_ERROR)
      {
         return;
      }
      // Query width and height
      glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_WIDTH, &cgl_w);
      glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_HEIGHT, &cgl_h);
      // Check for query errors
      GLenum gl_err = glGetError();
      if (gl_err != GL_NO_ERROR)
      {
         PRINT_DEBUG("❌ GL error during parameter query");
         // Unbind to clean up
         glBindRenderbuffer(GL_RENDERBUFFER, 0);
         return ;
      }
      // Unbind to restore state
      glBindRenderbuffer(GL_RENDERBUFFER, 0);
      assert(glGetError() == GL_NO_ERROR);
   }

   w = cgl_w;
   h = cgl_h;
#endif
}

bool EglWindow::isHighDpi() const
{
   return GetUseHiDPI();
}

void EglWindow::setWindowSize(int w, int h)
{
   EglMainThread::Get().ResizeWindow(handle, w, h);
}

void EglWindow::signalKeyDown(SDL_Keycode k, SDL_Keymod m)
{
   Event::Events e;
   e.keydown = {k, m};
   queueEvents({{EventType::Keydown, e}});
}

void EglWindow::signalQuit()
{
   queueEvents({{EventType::Quit, {}}});
}

void EglWindow::screenshot(string filename, bool convert)
{
   screenshot_filename = filename;
   Event::Events e;
   e.screenshot = {convert};
   queueEvents({{EventType::Screenshot, e}});
   // Queue up an expose, so Screenshot() can pull image from updated buffer
   signalExpose();
}

#endif // GLVIS_USE_EGL
