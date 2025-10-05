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
#include <csignal>

#include "cgl_main.hpp"

#include "../aux_vis.hpp"

#define NVTX_COLOR ::nvtx::kYellow
#include "../../nvtx.hpp" // IWYU pragma: keep

using namespace std;

///////////////////////////////////////////////////////////////////////////////
struct CGLMainThread::CreateWndCmd
{
   CGLWindow *wnd;
   int w, h;
   bool legacy_gl;
   promise<Handle> handle;
};

struct CGLMainThread::ResizeWndCmd
{
   Handle *handle;
   int w, h;
};

struct CGLMainThread::DeleteWndCmd
{
   CGLWindow *wnd;
   Handle *handle;
};

///////////////////////////////////////////////////////////////////////////////
bool CGLMainThread::CreateWndImpl(CreateWndCmd &cmd)
{
   dbg();
   Handle handle;

   const int multisamples = GetMultisample();
   dbg("multisamples:{}", multisamples);

   vector<CGLPixelFormatAttribute> pixAttribs =
   {
      kCGLPFASampleBuffers, CGLPixelFormatAttribute((multisamples > 0)?(1):(0)), // must be first
      kCGLPFASamples,       CGLPixelFormatAttribute(multisamples),
      kCGLPFAColorSize,     CGLPixelFormatAttribute(24),
      kCGLPFAAlphaSize,     CGLPixelFormatAttribute(8),
      kCGLPFADepthSize,     CGLPixelFormatAttribute(24),
      kCGLPFAAccelerated, // must be last
      CGLPixelFormatAttribute(0)
   };

   if (cmd.legacy_gl)
   {
      dbg("Inserting legacy OpenGL profile");
      // insert legacy OpenGL compatibility requirement
      auto it = (pixAttribs.end() -= 2);
      pixAttribs.insert(it, {kCGLPFAOpenGLProfile, CGLPixelFormatAttribute(kCGLOGLPVersion_Legacy)});
   }

   dbg("Choosing CGL pixel format");
   GLint numConfigs;
   CGLError err = CGLChoosePixelFormat(pixAttribs.data(), &handle.pix,
                                       &numConfigs);
   if (multisamples > 0 && (err != kCGLNoError || numConfigs < 1))
   {
      std::cerr << "CGL with multisampling is not supported, turning it off" <<
                std::endl;
      pixAttribs[1] = CGLPixelFormatAttribute(0);
      err = CGLChoosePixelFormat(pixAttribs.data(), &handle.pix, &numConfigs);
   }
   if (err != kCGLNoError || numConfigs < 1)
   {
      std::cerr << "CGL with hardware acceleration not supported, turning it off" <<
                std::endl;
      pixAttribs.pop_back();
      pixAttribs.back() = CGLPixelFormatAttribute(0);
      err = CGLChoosePixelFormat(pixAttribs.data(), &handle.pix, &numConfigs);
   }
   if (err != kCGLNoError || numConfigs < 1)
   {
      std::cerr << "Cannot set up CGL pixel format configuration, error: "
                << CGLErrorString(err) << std::endl;
      return false;
   }

   dbg("Creating CGL context");
   CGLContextObj ctx;
   err = CGLCreateContext(handle.pix, nullptr, &ctx);
   assert(err == kCGLNoError);
   if (err != kCGLNoError)
   {
      std::cerr << "Cannot create an OpenGL context, error: "
                << CGLErrorString(err) << std::endl;
      return false;
   }

   handle.ctx.reset(ctx);
   assert(handle.ctx != nullptr);

   windows.push_back(cmd.wnd);
   if (num_windows < 0)
   {
      num_windows = 1;
   }
   else
   {
      num_windows++;
   }
   cmd.handle.set_value(std::move(handle));

   return true;
}

///////////////////////////////////////////////////////////////////////////////
bool CGLMainThread::ResizeWndImpl(ResizeWndCmd &cmd)
{
   dbg();
   return true;
}

///////////////////////////////////////////////////////////////////////////////
bool CGLMainThread::DeleteWndImpl(DeleteWndCmd &cmd)
{
   dbg();
   if (cmd.handle->isInitialized())
   {
      dbg("Destroying buffers");
      // CGLDeleteFrameBuffers(cmd.handle);

      dbg();
      assert(glGetError() == GL_NO_ERROR);
      assert(cmd.handle);
      static_assert(sizeof(cmd.handle->buf_frame) == sizeof(GLuint),
                    "Buffer frame size must be 4 bytes");
      dbg("Deleting frame buffers");
      glDeleteFramebuffers(1, &cmd.handle->buf_frame);
      dbg("Deleting color buffers");
      glDeleteRenderbuffers(1, &cmd.handle->buf_color);
      dbg("Deleting depth buffers");
      glDeleteRenderbuffers(1, &cmd.handle->buf_depth);
      dbg("✅");
      assert(glGetError() == GL_NO_ERROR);
   }

   if (cmd.handle->pix)
   {
      dbg("Destroying pixel format");
      CGLError err = CGLDestroyPixelFormat(cmd.handle->pix);
      if (err != kCGLNoError)
      {
         std::cerr << "Cannot destroy pixel format, error: "
                   << CGLErrorString(err) << std::endl;
         return false;
      }
      cmd.handle->pix = nullptr;
   }

   dbg("Deleting window");
   windows.remove(cmd.wnd);
   num_windows--;

   dbg("✅");
   return true;
}

///////////////////////////////////////////////////////////////////////////////
void CGLMainThread::QueueWndCmd(CtrlCmd cmd, bool sync)
{
   dbg();
   future<void> wait_complete;
   if (sync)
   {
      dbg("Waiting for command to complete");
      wait_complete = cmd.finished.get_future();
   }
   // queue up our event
   {
      dbg("Queueing command");
      lock_guard<mutex> req_lock{window_cmd_mtx};
      window_cmds.emplace_back(std::move(cmd));
   }
   // wake up the main thread to handle our event
   events_available.notify_all();

   if (sync) { dbg("Waiting for command to complete"); wait_complete.get(); }
}

///////////////////////////////////////////////////////////////////////////////
void CGLMainThread::InterruptHandler(int param)
{
   dbg();
   CGLMainThread::Get().Terminate();
}

///////////////////////////////////////////////////////////////////////////////
CGLMainThread::CGLMainThread()
{
   dbg();
   GLint major, minor;
   CGLGetVersion(&major, &minor);
   std::cout << "Using CGL " << major << "." << minor << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
CGLMainThread::~CGLMainThread()
{
}

///////////////////////////////////////////////////////////////////////////////
CGLMainThread &CGLMainThread::Get()
{
   dbg();
   static CGLMainThread singleton;
   return singleton;
}

///////////////////////////////////////////////////////////////////////////////
CGLWindow::CGLHandle CGLMainThread::CreateWindow(CGLWindow *caller,
                                                 int w,
                                                 int h, bool legacy_gl)
{
   dbg();
   CtrlCmd cmd;
   cmd.type = CtrlCmdType::Create;

   CreateWndCmd create_window_cmd;
   cmd.create_cmd = &create_window_cmd;
   create_window_cmd.wnd = caller;
   create_window_cmd.w = w;
   create_window_cmd.h = h;
   create_window_cmd.legacy_gl = legacy_gl;

   auto result_handle = create_window_cmd.handle.get_future();

   dbg("Queueing create window command");
   QueueWndCmd(std::move(cmd), false);

   dbg("Waiting for window creation to complete");
   Handle handle = result_handle.get();
   assert(handle.isInitialized());

   if (!handle.isInitialized())
   {
      dbg("❌ Failed to create window");
      return handle;
   }

   dbg("Setting CGL context");
   CGLError err = CGLSetCurrentContext(handle.ctx.get());
   if (err != kCGLNoError)
   {
      dbg("❌ Failed to set context");
      return handle;
   }

   {
      // need to initialize GLEW before using any OpenGL functions
      const auto status = caller->initGLEW(legacy_gl);
      assert(status && glGetError() == GL_NO_ERROR);

      dbg("GL error: {}", glGetError());
      const GLubyte* renderer = glGetString(GL_RENDERER);
      dbg("OpenGL renderer: {}", (const char*)renderer);
      if (strstr((const char*)renderer, "Software") != nullptr)
      {
         dbg("⚠️ Using software renderer; expect limited functionality");
      }

      const GLubyte* version = glGetString(GL_VERSION);
      dbg("OpenGL version: {}", (const char*)version);
      dbg("GL error: {}", glGetError());
   }

   glGenFramebuffers(1, &handle.buf_frame);
   assert(glGetError() == GL_NO_ERROR);
   glCheckFramebufferStatus(GL_FRAMEBUFFER);
   assert(glGetError() == GL_NO_ERROR);

   dbg("GL error: {}", glGetError());
   dbg("new_handle.buf_frame:{}", (int)handle.buf_frame);
   assert(handle.buf_frame != 0);
   dbg("Creating buf_color and buf_depth renderbuffers on main thread");
   glGenRenderbuffers(1, &handle.buf_color);
   dbg("new_handle.buf_color:{}", (int)handle.buf_color);
   glGenRenderbuffers(1, &handle.buf_depth);
   dbg("new_handle.buf_depth:{}", (int)handle.buf_depth);
   ResizeWindow(handle, w, h);  // Perform initial resize here

   // Release current context
   CGLSetCurrentContext(nullptr);
   assert(glGetError() == GL_NO_ERROR);

   // dbg("Making context current on this thread");
   // CGLInitializeFrameBuffers(handle, cmd.create_cmd->w, cmd.create_cmd->h);

   err = CGLSetCurrentContext(handle.ctx.get());
   if (err != kCGLNoError)
   {
      dbg("❌ Cannot set CGL context as current");
   }
   dbg("✅");
   return handle;
}

///////////////////////////////////////////////////////////////////////////////
void CGLMainThread::ResizeWindow(Handle &handle, int w, int h)
{
   dbg("Resizing window to {}x{}", w, h);
   // CGLResizeWindow(handle, w, h);
   dbg("Resizing window to {}x{}", w, h);
   assert(glGetError() == GL_NO_ERROR);

   if (!handle.isInitialized())
   {
      dbg("❌ Handle not initialized, cannot resize");
      return;
   }

   dbg("glBindRenderbuffer buf_color");
   glBindRenderbuffer(GL_RENDERBUFFER, handle.buf_color);

   dbg("glRenderbufferStorage");
   glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, w, h);

   dbg("glBindRenderbuffer buf_depth");
   glBindRenderbuffer(GL_RENDERBUFFER, handle.buf_depth);
   glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);

   dbg("glBindFramebuffer buf_frame");
   glBindFramebuffer(GL_FRAMEBUFFER, handle.buf_frame);
   glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                             GL_RENDERBUFFER, handle.buf_color);
   glFramebufferRenderbuffer(GL_DRAW_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                             GL_RENDERBUFFER, handle.buf_depth);

   if (glGetError() != GL_NO_ERROR ||
       glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
   {
      dbg("❌ Cannot resize the framebuffer!");
   }

   assert(glGetError() == GL_NO_ERROR);
   dbg("✅");
}

///////////////////////////////////////////////////////////////////////////////
void CGLMainThread::DeleteWindow(CGLWindow *caller, Handle &handle)
{
   dbg();
   if (!handle.isInitialized()) { return; }

   CtrlCmd cmd;
   cmd.type = CtrlCmdType::Delete;

   DeleteWndCmd del_cmd;
   cmd.delete_cmd = &del_cmd;

   del_cmd.wnd = caller;
   del_cmd.handle = &handle;

   QueueWndCmd(std::move(cmd), true);
   dbg("✅");
}

///////////////////////////////////////////////////////////////////////////////
void CGLMainThread::Terminate()
{
   dbg();
   CtrlCmd cmd;
   cmd.type = CtrlCmdType::Terminate;

   // queue to the front to get priority
   {
      lock_guard<mutex> req_lock{window_cmd_mtx};
      window_cmds.emplace_front(std::move(cmd));
   }
   // wake up the main thread to handle our event
   events_available.notify_all();
}

///////////////////////////////////////////////////////////////////////////////
void CGLMainThread::MainLoop(bool server)
{
   dbg();
   server_mode = server;
   bool terminating = false;

   // set up interrupt handler for graceful closing
   signal(SIGINT, InterruptHandler);

   while (true)
   {
      bool events_pending = false;
      {
         lock_guard<mutex> evt_guard{window_cmd_mtx};
         events_pending = !window_cmds.empty();
      }
      if (events_pending)
      {
         do
         {
            CtrlCmd cmd;
            // Fetch next event from the queue
            {
               lock_guard<mutex> evt_guard{window_cmd_mtx};
               cmd = std::move(window_cmds.front());
               window_cmds.pop_front();
               events_pending = !window_cmds.empty();
               // Skip non-delete events if terminating
               if (terminating && cmd.type != CtrlCmdType::Delete)
               {
                  continue;
               }
            }

            switch (cmd.type)
            {
               case CtrlCmdType::Create:
                  if (!CreateWndImpl(*cmd.create_cmd))
                  {
                     terminating = true;
                  }
                  break;
               case CtrlCmdType::Resize:
                  if (!ResizeWndImpl(*cmd.resize_cmd))
                  {
                     terminating = true;
                  }
                  break;
               case CtrlCmdType::Delete:
                  if (!DeleteWndImpl(*cmd.delete_cmd))
                  {
                     terminating = true;
                  }
                  break;
               case CtrlCmdType::Terminate:
                  // do not wait for windows to open
                  if (num_windows < 0) { num_windows = 0; }
                  terminating = true;
                  break;
            }

            // Signal completion of the command, in case worker thread is waiting.
            cmd.finished.set_value();
         }
         while (events_pending);
      }

      if (num_windows == 0)
      {
         if (!server_mode || terminating)
         {
            break;
         }
      }

      if (terminating)
      {
         for (CGLWindow *wnd : windows)
         {
            wnd->signalQuit();
         }
      }

      {
         unique_lock<mutex> event_lock{window_cmd_mtx};
         events_available.wait(event_lock, [this]()
         {
            // Sleep until events from windows can be handled
            return !window_cmds.empty();
         });
      }
   }
}

#endif // GLVIS_USE_CGL
