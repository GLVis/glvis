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

#ifndef GLVIS_SDL_MAC_HPP
#define GLVIS_SDL_MAC_HPP
#include "sdl_helper.hpp"

class SdlCocoaPlatform final : public SdlNativePlatform
{
public:
   void WaitEvent();
   void SendEvent();

   void ContextUpdate();

   // TODO: Workaround methods for SDL 2.0.12, delete in follow-up PR
   bool UseThreadWorkaround() const;
   void ClearCurrentContext(int wnd_id);
   void SetCurrentContext(int wnd_id);
   void SwapWindow();
};

#endif
