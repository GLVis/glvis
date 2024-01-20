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

#ifndef GLVIS_SDL_WINDOWS_HPP
#define GLVIS_SDL_WINDOWS_HPP

#include "sdl_helper.hpp"

#ifdef SDL_VIDEO_DRIVER_WINDOWS

class SdlWindowsPlatform final : public SdlNativePlatform
{
public:
   SdlWindowsPlatform();

   ~SdlWindowsPlatform();

   void WaitEvent();

   void SendEvent();
private:
   struct Impl;
   std::unique_ptr<Impl> m_impl;
};

#endif

#endif // GLVIS_SDL_WINDOWS_HPP
