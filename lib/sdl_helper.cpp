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

#include <iostream>

#include "sdl_helper.hpp"

#ifndef __EMSCRIPTEN__
#include <SDL2/SDL_syswm.h>
#else
#include <SDL_syswm.h>
#endif

#include "sdl_x11.hpp"
#include "sdl_mac.hpp"
#include "sdl_windows.hpp"

using namespace std;

std::unique_ptr<SdlNativePlatform>
SdlNativePlatform::Create(SDL_Window* window)
{
   SDL_SysWMinfo sysinfo;
   SDL_VERSION(&sysinfo.version);
   if (!SDL_GetWindowWMInfo(window, &sysinfo))
   {
      cerr << "Error: unable to get window manager information for the "
           << "current window." << endl;
      return {};
   }
   switch (sysinfo.subsystem)
   {
#if defined(SDL_VIDEO_DRIVER_WINDOWS)
      case SDL_SYSWM_WINDOWS:
         return std::unique_ptr<SdlNativePlatform> {new SdlWindowsPlatform};
#endif
#if defined(SDL_VIDEO_DRIVER_COCOA)
      case SDL_SYSWM_COCOA:
         return std::unique_ptr<SdlNativePlatform> {new SdlCocoaPlatform};
#endif
#if defined(SDL_VIDEO_DRIVER_X11)
      case SDL_SYSWM_X11:
         return std::unique_ptr<SdlNativePlatform> {new SdlX11Platform};
#endif
      case SDL_SYSWM_UNKNOWN:
      default:
         cerr << "Error: unrecognized window manager system." << endl;
         return {};
   }
}
