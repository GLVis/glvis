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

#ifdef GLVIS_USE_CGL

#include <iostream>
#include <vector>

#include "cgl.hpp"
#include "cgl_main.hpp"

#include "../aux_vis.hpp"

#define NVTX_COLOR ::nvtx::kMagenta
#include "../../nvtx.hpp" // IWYU pragma: keep

using namespace std;

CGLWindow::CGLWindow()
{
   dbg("Could create window");
}

CGLWindow::~CGLWindow()
{
   dbg();
   CGLMainThread::Get().DeleteWindow(this, handle);
}
// return initGLEW(legacyGlOnly);

bool CGLWindow::initGLEW(bool legacyGlOnly)
{
   return GLWindow::initGLEW(legacyGlOnly);
}

bool CGLWindow::createWindow(const char *title, int x, int y, int w, int h,
                             bool legacyGlOnly)
{
   dbg("title: '{}' x:{} y:{} {}x{} legacyGlOnly:{}",
       title, x, y, w, h, legacyGlOnly);
   handle = CGLMainThread::Get().CreateWindow(this, w, h, legacyGlOnly);
   assert(glGetError() == GL_NO_ERROR);

   if (!handle.isInitialized())
   {
      dbg("❌ Cannot create CGL window");
      return false;
   }

   // dbg("CGL context is ready");
   // return initGLEW(legacyGlOnly);
   return true;
}

void CGLWindow::queueEvents(vector<Event> events)
{
   dbg();
   {
      lock_guard<mutex> evt_guard{event_mutex};
      waiting_events.insert(waiting_events.end(), events.begin(), events.end());
   }
   if (is_multithreaded)
   {
      events_available.notify_all();
   }
}

void CGLWindow::mainLoop()
{
   dbg();
   running = true;
   while (running)
   {
      mainIter();
   }
}

void CGLWindow::mainIter()
{
   dbg();
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

void CGLWindow::signalLoop()
{
   dbg();
   // Note: not executed from the main thread
   {
      lock_guard<mutex> evt_guard{event_mutex};
      call_idle_func = true;
   }
   events_available.notify_all();
}

void CGLWindow::getGLDrawSize(int& w, int& h) const
{
   const auto err = glGetError();
   dbg("glGetError: {:x}", err);
   std::cout << "glGetError:" << err << std::endl;

   assert(glGetError() == GL_NO_ERROR);
   GLint cgl_w, cgl_h;
   static_assert(sizeof(GLint) == sizeof(int),
                 "CGL width size must be 4 bytes");
   // CGLGetRenderbufferParameter(handle, cgl_w, cgl_h);
   {
      dbg("glGetError:{}", glGetError());
      assert(glGetError() == GL_NO_ERROR);
      assert(handle.isInitialized());

      // Optional: Ensure context is current if not guaranteed elsewhere
      CGLError err = CGLSetCurrentContext(handle.ctx.get());
      if (err != kCGLNoError) { dbg("❌ Cannot set CGL context as current");}
      // Bind the color renderbuffer (assuming that's what we want to query)
      glBindRenderbuffer(GL_RENDERBUFFER, handle.buf_color);
      if (glGetError() != GL_NO_ERROR)
      {
         dbg("❌ Failed to bind renderbuffer");
         // return EXIT_FAILURE;
      }
      // Query width and height
      glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_WIDTH, &cgl_w);
      glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_HEIGHT, &cgl_h);
      // Check for query errors
      GLenum gl_err = glGetError();
      if (gl_err != GL_NO_ERROR)
      {
         dbg("❌ GL error during parameter query: {}", gl_err);
         // Unbind to clean up
         glBindRenderbuffer(GL_RENDERBUFFER, 0);
         return ;
      }
      dbg("Renderbuffer parameters - width: {}, height: {}", cgl_w, cgl_h);
      // Unbind to restore state
      glBindRenderbuffer(GL_RENDERBUFFER, 0);
      assert(glGetError() == GL_NO_ERROR);
   }

   // cgl_w = 1024, cgl_h = 768; // TEMP
   dbg("CGL draw size: {}x{}", (int) cgl_w, (int)cgl_h);
   w = cgl_w;
   h = cgl_h;
}

bool CGLWindow::isHighDpi() const
{
   dbg("{}", GetUseHiDPI());
   return GetUseHiDPI();
}

void CGLWindow::setWindowSize(int w, int h)
{
   dbg("w:{} h:{}", w, h);
   CGLMainThread::Get().ResizeWindow(handle, w, h);
}

void CGLWindow::signalKeyDown(SDL_Keycode k, SDL_Keymod m)
{
   dbg("k:{} m:{}", (int)k, (int)m);
   Event::Events e;
   e.keydown = {k, m};
   queueEvents({{EventType::Keydown, e}});
}

void CGLWindow::signalQuit()
{
   dbg("Quit");
   queueEvents({{EventType::Quit, {}}});
}

void CGLWindow::screenshot(string filename, bool convert)
{
   dbg("filename:{} convert:{}", filename, convert);
   screenshot_filename = filename;
   Event::Events e;
   e.screenshot = {convert};
   queueEvents({{EventType::Screenshot, e}});
   // Queue up an expose, so Screenshot() can pull image from updated buffer
   signalExpose();
}

#endif // GLVIS_USE_CGL
