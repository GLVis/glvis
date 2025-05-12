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

#include "glwindow.hpp"

#include <EGL/egl.h>

class EglWindow : public GLWindow
{
   EGLDisplay disp{EGL_NO_DISPLAY};
   EGLSurface surf{EGL_NO_SURFACE};
   EGLContext ctx{EGL_NO_CONTEXT};
   EGLConfig eglCfg{};

public:
   EglWindow();
   ~EglWindow();

   /** @brief Creates a new OpenGL window. Returns false if EGL or OpenGL
       initialization fails. */
   bool createWindow(const char *title, int x, int y, int w, int h,
                     bool legacyGlOnly) override;

   void getWindowSize(int& w, int& h) const override { getGLDrawSize(w, h); }
   void getGLDrawSize(int& w, int& h) const override;

   // use the default values of SdlWindow, LoadFont() ignores them anyway
   void getDpi(int& wdpi, int& hdpi) const override { wdpi = hdpi = 72; }
   bool isHighDpi() const override { return true; }

   void setWindowSize(int w, int h) override;

   bool isWindowInitialized() const override { return surf != EGL_NO_SURFACE; }
   bool isGlInitialized() const override { return ctx != EGL_NO_CONTEXT; }

   void signalExpose() override;
   void signalSwap() override { }

   bool isExposePending() const override { return false; }
   // used in Screenshot, as there is no swapping, the single buffer is always
   // up to date and can be read directly
   bool isSwapPending() const override { return true; }

   void screenshot(std::string filename, bool convert = false) override;
};

#endif //GLVIS_EGL_HPP