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

#ifndef GLVIS_GLWINDOW_HPP
#define GLVIS_GLWINDOW_HPP

#include "gl/renderer.hpp"
#include <memory>

class GLWindow
{
protected:
   std::unique_ptr<gl3::MeshRenderer> renderer;

   bool initGLEW(bool legacyGlOnly);
public:
   virtual ~GLWindow() = default;

   /// Creates a new OpenGL window
   virtual bool createWindow(const char *title, int x, int y, int w, int h,
                             bool legacyGlOnly) = 0;

   /// Returns the renderer object
   inline gl3::MeshRenderer& getRenderer() { return *renderer.get(); }

   /// Returns size of the window
   virtual void getWindowSize(int& w, int& h) const { w = h = 0; }

   /// Returns the drawable size
   virtual void getGLDrawSize(int& w, int& h) const { w = h = 0; }

   /// Returns the resolution (DPI) of the display
   virtual void getDpi(int& wdpi, int& hdpi) const { wdpi = hdpi = 0; }

   /// Checks if the display has high resolution (DPI)
   virtual bool isHighDpi() const { return false; }

   /// Set title of the window (string version)
   virtual void setWindowTitle(const std::string& title) { }

   /// Set title of the window (C-string version)
   virtual void setWindowTitle(const char* title) { }

   /// Set window size
   virtual void setWindowSize(int w, int h) { }

   /// Set window position
   virtual void setWindowPos(int x, int y) { }

   /// Returns true if the window has been succesfully initialized
   virtual bool isWindowInitialized() const { return false; }

   /// Returns true if the OpenGL context was successfully initialized
   virtual bool isGlInitialized() const { return false; }

   /// Signals expose event when objects have been updated
   virtual void signalExpose() = 0;

   /// Signals swap event when the back buffer is ready for swapping
   virtual void signalSwap() = 0;

   /// Checks if the swap event is pending
   virtual bool isSwapPending() const { return false; }

   /// Checks if the expose event is pending
   virtual bool isExposePending() const { return false; }

   /// Saves a screenshot ot the file, performing conversion optionally
   virtual void screenshot(std::string filename, bool convert = false) { }
};

#endif //GLVIS_GLWINDOW_HPP
