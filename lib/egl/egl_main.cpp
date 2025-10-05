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

#include <iostream>
#include <csignal>

#include "egl_main.hpp"
#include "../aux_vis.hpp"

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

#define NVTX_COLOR ::nvtx::kYellow
#include "../../nvtx.hpp" // IWYU pragma: keep

using namespace std;

struct EglMainThread::CreateWndCmd
{
   EglWindow *wnd;
   int w, h;
   bool legacy_gl;
   promise<Handle> out_handle;
};

struct EglMainThread::ResizeWndCmd
{
   Handle *handle;
   int w, h;
};

struct EglMainThread::DeleteWndCmd
{
   EglWindow *wnd;
   Handle *handle;
};

bool EglMainThread::CreateWndImpl(CreateWndCmd &cmd)
{
   Handle new_handle;

#ifdef GLVIS_USE_EGL
   // 1. Select an appropriate configuration

   const int multisamples = GetMultisample();

   EGLint configAttribs[] =
   {
      EGL_SAMPLE_BUFFERS, (multisamples > 0)?(1):(0), // must be first
      EGL_SAMPLES, multisamples,
      EGL_SURFACE_TYPE, EGL_PBUFFER_BIT,
      EGL_COLOR_BUFFER_TYPE, EGL_RGB_BUFFER,
      EGL_BLUE_SIZE, 8,
      EGL_GREEN_SIZE, 8,
      EGL_RED_SIZE, 8,
      EGL_ALPHA_SIZE, 8,
      EGL_DEPTH_SIZE, 24,
      EGL_CONFORMANT, EGL_OPENGL_BIT,
      EGL_RENDERABLE_TYPE, EGL_OPENGL_BIT,
      EGL_NONE
   };

   EGLint numConfigs;

   if (multisamples > 0)
   {
      if (!eglChooseConfig(disp, configAttribs, NULL, 0, &numConfigs) ||
          numConfigs < 1)
      {
         std::cerr << "EGL with multisampling is not supported, turning it off" <<
                   std::endl;
         // turn off multisampling
         configAttribs[1] = 0;
      }
   }

   if (!eglChooseConfig(disp, configAttribs, &new_handle.eglCfg, 1, &numConfigs) ||
       numConfigs < 1)
   {
      std::cerr << "Cannot find working EGL configuration!" << std::endl;
      return false;
   }

   // 2. Create a surface
   const EGLint pbufferAttribs[] =
   {
      EGL_WIDTH, cmd.w,
      EGL_HEIGHT, cmd.h,
      EGL_NONE
   };

   new_handle.surf = eglCreatePbufferSurface(disp, new_handle.eglCfg,
                                             pbufferAttribs);
   if (new_handle.surf == EGL_NO_SURFACE)
   {
      std::cerr << "Cannot create a pixel buffer, error: " << eglGetError() <<
                std::endl;
      return false;
   }
#endif
#ifdef GLVIS_USE_CGL
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
   CGLError err = CGLChoosePixelFormat(pixAttribs.data(), &new_handle.pix,
                                       &numConfigs);
   if (multisamples > 0 && (err != kCGLNoError || numConfigs < 1))
   {
      std::cerr << "CGL with multisampling is not supported, turning it off" <<
                std::endl;
      pixAttribs[1] = CGLPixelFormatAttribute(0);
      err = CGLChoosePixelFormat(pixAttribs.data(), &new_handle.pix, &numConfigs);
   }
   if (err != kCGLNoError || numConfigs < 1)
   {
      std::cerr << "CGL with hardware acceleration not supported, turning it off" <<
                std::endl;
      pixAttribs.pop_back();
      pixAttribs.back() = CGLPixelFormatAttribute(0);
      err = CGLChoosePixelFormat(pixAttribs.data(), &new_handle.pix, &numConfigs);
   }
   if (err != kCGLNoError || numConfigs < 1)
   {
      std::cerr << "Cannot set up CGL pixel format configuration, error: "
                << CGLErrorString(err) << std::endl;
      return false;
   }

   dbg("Creating CGL context");
   CGLContextObj ctx;
   err = CGLCreateContext(new_handle.pix, nullptr, &ctx);
   assert(err == kCGLNoError);
   if (err != kCGLNoError)
   {
      std::cerr << "Cannot create an OpenGL context, error: "
                << CGLErrorString(err) << std::endl;
      return false;
   }

   new_handle.ctx.reset(ctx);
   assert(new_handle.ctx != nullptr);
#endif

   windows.push_back(cmd.wnd);
   if (num_windows < 0)
   {
      num_windows = 1;
   }
   else
   {
      num_windows++;
   }
   cmd.out_handle.set_value(std::move(new_handle));

   return true;
}

bool EglMainThread::ResizeWndImpl(ResizeWndCmd &cmd)
{
#ifdef GLVIS_USE_EGL
   const EGLint pbufferAttribs[] =
   {
      EGL_WIDTH, cmd.w,
      EGL_HEIGHT, cmd.h,
      EGL_NONE
   };

   EGLSurface surf_new = eglCreatePbufferSurface(disp, cmd.handle->eglCfg,
                                                 pbufferAttribs);
   if (surf_new == EGL_NO_SURFACE)
   {
      std::cerr << "Cannot create a pixel buffer, error: " << eglGetError() <<
                std::endl;
      return false;
   }

   if (!eglDestroySurface(disp, cmd.handle->surf))
   {
      std::cerr << "Cannot destroy surface, error: " << eglGetError() << std::endl;
      return false;
   }

   cmd.handle->surf = surf_new;
#endif
   return true;
}

bool EglMainThread::DeleteWndImpl(DeleteWndCmd &cmd)
{
#ifdef GLVIS_USE_EGL
   if (cmd.handle->ctx != EGL_NO_CONTEXT)
   {
      if (!eglDestroyContext(disp, cmd.handle->ctx))
      {
         std::cerr << "Cannot destroy context, error: " << eglGetError() << std::endl;
         return false;
      }
      cmd.handle->ctx = EGL_NO_CONTEXT;
   }

   if (cmd.handle->surf != EGL_NO_SURFACE)
   {
      if (!eglDestroySurface(disp, cmd.handle->surf))
      {
         std::cerr << "Cannot destroy surface, error: " << eglGetError() << std::endl;
         return false;
      }
      cmd.handle->surf = EGL_NO_SURFACE;
   }
#endif
#ifdef GLVIS_USE_CGL

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
#endif
   windows.remove(cmd.wnd);
   num_windows--;

   return true;
}

void EglMainThread::QueueWndCmd(CtrlCmd cmd, bool sync)
{
   future<void> wait_complete;
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
   events_available.notify_all();

   if (sync) { wait_complete.get(); }
}

void EglMainThread::InterruptHandler(int param)
{
   EglMainThread::Get().Terminate();
}

EglMainThread::EglMainThread()
{
#ifdef GLVIS_USE_EGL
   if (disp != EGL_NO_DISPLAY)
   {
      return;
   }

   // Initialize EGL
   disp = eglGetDisplay(EGL_DEFAULT_DISPLAY);
   if (disp == EGL_NO_DISPLAY)
   {
      std::cerr << "FATAL: Failed to get an EGL display: " << eglGetError() <<
                std::endl;
      return;
   }

   EGLint major, minor;

   if (!eglInitialize(disp, &major, &minor))
   {
      std::cerr << "FATAL: Failed to initialize EGL: " << eglGetError() << std::endl;
      return;
   }

   std::cout << "Using EGL " << major << "." << minor << std::endl;
#endif
#ifdef GLVIS_USE_CGL
   GLint major, minor;
   CGLGetVersion(&major, &minor);
   std::cout << "Using CGL " << major << "." << minor << std::endl;
#endif
}

EglMainThread::~EglMainThread()
{
#ifdef GLVIS_USE_EGL
   eglTerminate(disp);
#endif
}

EglMainThread &EglMainThread::Get()
{
   static EglMainThread singleton;
   return singleton;
}

EglMainThread::Handle EglMainThread::CreateWindow(EglWindow *caller, int w,
                                                  int h, bool legacy_gl)
{
   CtrlCmd cmd;
   cmd.type = CtrlCmdType::Create;

   CreateWndCmd crt_cmd;
   cmd.create_cmd = &crt_cmd;

   crt_cmd.wnd = caller;
   crt_cmd.w = w;
   crt_cmd.h = h;
   crt_cmd.legacy_gl = legacy_gl;

   auto res_handle = crt_cmd.out_handle.get_future();

   QueueWndCmd(std::move(cmd), false);

   Handle out_hnd = res_handle.get();

#ifdef GLVIS_USE_EGL
   if (out_hnd.surf == EGL_NO_SURFACE)
   {
      return out_hnd;
   }

   // 3. Bind the API
   if (!eglBindAPI(EGL_OPENGL_API))
   {
      std::cerr << "Cannot bind OpenGL API, error: " << eglGetError() << std::endl;
      return out_hnd;
   }

   // 4. Create a context and make it current
   if (legacy_gl)
   {
      // Try and probe for a core/compatibility context. Needed for Mac OS X,
      // which will only support OpenGL 2.1 if you don't create a core context.
      PRINT_DEBUG("Opening OpenGL core profile context..." << std::flush);
      const EGLint attrListCore[] =
      {
         EGL_CONTEXT_OPENGL_PROFILE_MASK, EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT,
         EGL_NONE
      };
      out_hnd.ctx = eglCreateContext(disp, out_hnd.eglCfg, EGL_NO_CONTEXT,
                                     attrListCore);
      if (out_hnd.ctx == EGL_NO_CONTEXT)
      {
         PRINT_DEBUG("failed." << std::endl);
         PRINT_DEBUG("Opening OpenGL compatibility profile context..." << std::flush);
         const EGLint attrListCompat[] =
         {
            EGL_CONTEXT_OPENGL_PROFILE_MASK, EGL_CONTEXT_OPENGL_COMPATIBILITY_PROFILE_BIT,
            EGL_NONE
         };
         out_hnd.ctx = eglCreateContext(disp, out_hnd.eglCfg, EGL_NO_CONTEXT,
                                        attrListCompat);
         if (out_hnd.ctx == EGL_NO_CONTEXT)
         {
            PRINT_DEBUG("failed." << std::endl);
         }
         else
         {
            PRINT_DEBUG("success!" << std::endl);
         }
      }
      else
      {
         PRINT_DEBUG("success!" << std::endl);
      }
   }

   if (out_hnd.ctx == EGL_NO_CONTEXT)
   {
      PRINT_DEBUG("Opening OpenGL context with no flags..." << std::flush);
      out_hnd.ctx = eglCreateContext(disp, out_hnd.eglCfg, EGL_NO_CONTEXT,
                                     NULL);
      if (out_hnd.ctx == EGL_NO_CONTEXT)
      {
         PRINT_DEBUG("failed." << std::endl);
         std::cerr << "Cannot create an EGL context, error: " << eglGetError() <<
                   std::endl;

         return out_hnd;
      }
      else
      {
         PRINT_DEBUG("success!" << std::endl);
      }
   }

   if (!eglMakeCurrent(disp, out_hnd.surf, out_hnd.surf, out_hnd.ctx))
   {
      std::cerr << "Cannot set the EGL context as current, error: " << eglGetError()
                << std::endl;

      return out_hnd;
   }
#endif //GLVIS_USE_EGL
#ifdef GLVIS_USE_CGL

   if (!out_hnd.isInitialized())
   {
      dbg("❌ Failed to create window");
      return out_hnd;
   }

   // dbg("Setting CGL context");
   CGLError err = CGLSetCurrentContext(out_hnd.ctx.get());
   if (err != kCGLNoError)
   {
      dbg("❌ Failed to set context");
      return out_hnd;
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

   glGenFramebuffers(1, &out_hnd.buf_frame);
   assert(glGetError() == GL_NO_ERROR);
   glCheckFramebufferStatus(GL_FRAMEBUFFER);
   assert(glGetError() == GL_NO_ERROR);

   dbg("GL error: {}", glGetError());
   dbg("new_handle.buf_frame:{}", (int)out_hnd.buf_frame);
   assert(out_hnd.buf_frame != 0);
   dbg("Creating buf_color and buf_depth renderbuffers on main thread");
   glGenRenderbuffers(1, &out_hnd.buf_color);
   dbg("new_handle.buf_color:{}", (int)out_hnd.buf_color);
   glGenRenderbuffers(1, &out_hnd.buf_depth);
   dbg("new_handle.buf_depth:{}", (int)out_hnd.buf_depth);
   ResizeWindow(out_hnd, w, h);  // Perform initial resize here

   // Release current context
   CGLSetCurrentContext(nullptr);
   assert(glGetError() == GL_NO_ERROR);

   // dbg("Making context current on this thread");
   // CGLInitializeFrameBuffers(handle, cmd.create_cmd->w, cmd.create_cmd->h);

   err = CGLSetCurrentContext(out_hnd.ctx.get());
   if (err != kCGLNoError)
   {
      dbg("❌ Cannot set CGL context as current");
   }
   dbg("✅");
#endif

   return out_hnd;
}

void EglMainThread::ResizeWindow(Handle &handle, int w, int h)
{
   if (!handle.isInitialized()) { return; }

#ifdef GLVIS_USE_EGL
   CtrlCmd cmd;
   cmd.type = CtrlCmdType::Resize;

   ResizeWndCmd res_cmd;
   cmd.resize_cmd = &res_cmd;

   res_cmd.handle = &handle;
   res_cmd.w = w;
   res_cmd.h = h;

   QueueWndCmd(std::move(cmd), true);

   if (!eglMakeCurrent(disp, handle.surf, handle.surf, handle.ctx))
   {
      std::cerr << "Cannot set the EGL context as current, error: " << eglGetError()
                << std::endl;
      return;
   }
#endif
#ifdef GLVIS_USE_CGL

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
#endif
}

void EglMainThread::DeleteWindow(EglWindow *caller, Handle &handle)
{
   if (!handle.isInitialized()) { return; }

   CtrlCmd cmd;
   cmd.type = CtrlCmdType::Delete;

   DeleteWndCmd del_cmd;
   cmd.delete_cmd = &del_cmd;

   del_cmd.wnd = caller;
   del_cmd.handle = &handle;

   QueueWndCmd(std::move(cmd), true);
}

void EglMainThread::Terminate()
{
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

void EglMainThread::MainLoop(bool server)
{
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
         for (EglWindow *wnd : windows)
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

#endif // GLVIS_USE_EGL || GLVIS_USE_CGL
