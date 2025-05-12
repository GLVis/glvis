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

#include "egl.hpp"
#include "aux_vis.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

EglWindow::EglWindow()
{
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
}

EglWindow::~EglWindow()
{
   if (ctx != EGL_NO_CONTEXT) { eglDestroyContext(disp, ctx); }
   if (surf != EGL_NO_SURFACE) { eglDestroySurface(disp, surf); }
   eglTerminate(disp);
}

bool EglWindow::createWindow(int w, int h, bool legacyGlOnly)
{
   // 1. Select an appropriate configuration

   const int multisamples = GetMultisample();

   std::vector<EGLint> configAttribs =
   {
      EGL_SURFACE_TYPE, EGL_PBUFFER_BIT,
      EGL_COLOR_BUFFER_TYPE, EGL_RGB_BUFFER,
      EGL_BLUE_SIZE, 8,
      EGL_GREEN_SIZE, 8,
      EGL_RED_SIZE, 8,
      EGL_ALPHA_SIZE, 8,
      EGL_DEPTH_SIZE, 24,
      EGL_CONFORMANT, EGL_OPENGL_BIT,
      EGL_RENDERABLE_TYPE, EGL_OPENGL_BIT,
      EGL_SAMPLE_BUFFERS, (multisamples > 0)?(1):(0),
      EGL_SAMPLES, multisamples,
      EGL_NONE
   };

   EGLint numConfigs;

   if (multisamples > 0)
   {
      if (!eglChooseConfig(disp, configAttribs.data(), NULL, 0, &numConfigs) ||
          numConfigs < 1)
      {
         std::cerr << "EGL with multisampling is not supported, turning it off" <<
                   std::endl;
         // turn off multisampling
         auto it = std::find(configAttribs.begin(), configAttribs.end(),
                             EGL_SAMPLE_BUFFERS);
         *(++it) = 0;
      }
   }

   EGLConfig eglCfg;

   if (!eglChooseConfig(disp, configAttribs.data(), &eglCfg, 1, &numConfigs) ||
       numConfigs < 1)
   {
      std::cerr << "Cannot find working EGL configuration!" << std::endl;
      return false;
   }

   // 2. Create a surface
   const EGLint pbufferAttribs[] =
   {
      EGL_WIDTH, w,
      EGL_HEIGHT, h,
      EGL_NONE
   };

   surf = eglCreatePbufferSurface(disp, eglCfg, pbufferAttribs);
   if (surf == EGL_NO_SURFACE)
   {
      std::cerr << "Cannot create a pixel buffer, error: " << eglGetError() <<
                std::endl;
      return false;
   }

   // 3. Bind the API
   if (!eglBindAPI(EGL_OPENGL_API))
   {
      std::cerr << "Cannot bind OpenGL API, error: " << eglGetError() << std::endl;
      return false;
   }

   // 4. Create a context and make it current
   if (legacyGlOnly)
   {
      // Try and probe for a core/compatibility context. Needed for Mac OS X,
      // which will only support OpenGL 2.1 if you don't create a core context.
      PRINT_DEBUG("Opening OpenGL core profile context..." << std::flush);
      const EGLint attrListCore[] =
      {
         EGL_CONTEXT_OPENGL_PROFILE_MASK, EGL_CONTEXT_OPENGL_CORE_PROFILE_BIT,
         EGL_NONE
      };
      ctx = eglCreateContext(disp, eglCfg, EGL_NO_CONTEXT, attrListCore);
      if (ctx == EGL_NO_CONTEXT)
      {
         PRINT_DEBUG("failed." << std::endl);
         PRINT_DEBUG("Opening OpenGL core profile context..." << std::flush);
         const EGLint attrListCompat[] =
         {
            EGL_CONTEXT_OPENGL_PROFILE_MASK, EGL_CONTEXT_OPENGL_COMPATIBILITY_PROFILE_BIT,
            EGL_NONE
         };
         ctx = eglCreateContext(disp, eglCfg, EGL_NO_CONTEXT, attrListCompat);
         if (ctx == EGL_NO_CONTEXT)
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

   if (ctx == EGL_NO_CONTEXT)
   {
      PRINT_DEBUG("Opening OpenGL context with no flags..." << std::flush);
      ctx = eglCreateContext(disp, eglCfg, EGL_NO_CONTEXT, NULL);
      if (ctx == EGL_NO_CONTEXT)
      {
         PRINT_DEBUG("failed." << std::endl);
         std::cerr << "Cannot create an EGL context, error: " << eglGetError() <<
                   std::endl;
         return false;
      }
      else
      {
         PRINT_DEBUG("success!" << std::endl);
      }
   }

   if (!eglMakeCurrent(disp, surf, surf, ctx))
   {
      std::cerr << "Cannot set the EGL context as current, error: " << eglGetError()
                << std::endl;
      return false;
   }

#ifndef __EMSCRIPTEN__
   glEnable(GL_DEBUG_OUTPUT);
#endif

   PRINT_DEBUG("EGL context is ready" << std::endl);

   return initGLEW(legacyGlOnly);
}

void EglWindow::getGLDrawSize(int& w, int& h)
{
   EGLint egl_w, egl_h;
   eglQuerySurface(disp, surf, EGL_WIDTH, &egl_w);
   eglQuerySurface(disp, surf, EGL_HEIGHT, &egl_h);
   w = egl_w;
   h = egl_h;
}

void EglWindow::signalExpose()
{
   MyExpose();
}

void EglWindow::screenshot(std::string filename, bool convert)
{
   Screenshot(filename.c_str(), convert);
}
