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

#ifndef GLVIS_SDL_HELPER_HPP
#define GLVIS_SDL_HELPER_HPP

#include "gl/platform_gl.hpp"
#include <memory>

class SdlNativePlatform
{
public:
   static std::unique_ptr<SdlNativePlatform> Create(SDL_Window* window);

   virtual ~SdlNativePlatform() = default;

   // Registers the window handle with this platform handler, in order to
   // wait for events. This is needed for X11, which has one event pipe
   // per window.
   virtual void RegisterWindow(SDL_Window* window) { }
   // Unregisters the window handle.
   virtual void UnregisterWindow(SDL_Window* window) { }
   // SDL_WaitEvent only polls for events and sleeps, instead of actually
   // blocking. This method calls the system-native blocking event pump.
   virtual void WaitEvent() = 0;
   // This method sends a system-native event, which will wake the blocking
   // event pump.
   virtual void SendEvent() = 0;
};

#endif
