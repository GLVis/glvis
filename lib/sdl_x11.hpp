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

#ifndef GLVIS_SDL_X11_HPP
#define GLVIS_SDL_X11_HPP

#include "sdl_helper.hpp"

#ifdef SDL_VIDEO_DRIVER_X11

#include <unordered_map>

#include <poll.h>

#include <unistd.h>    // pipe, fcntl, write
#include <fcntl.h>     // fcntl
#include <cerrno>      // errno, EAGAIN
#include <cstdio>      // perror

class SdlX11Platform final : public SdlNativePlatform
{
public:
   SdlX11Platform()
   {
      // Create pipe for external events
      if (pipe(event_pfd) == -1)
      {
         perror("pipe()");
         exit(EXIT_FAILURE);
      }
      int flag = fcntl(event_pfd[0], F_GETFL);
      fcntl(event_pfd[0], F_SETFL, flag | O_NONBLOCK);
   }

   ~SdlX11Platform()
   {
      close(event_pfd[0]);
      close(event_pfd[1]);
   }

   void RegisterWindow(SDL_Window* window);

   void UnregisterWindow(SDL_Window* window);

   void WaitEvent();

   void SendEvent()
   {
      char c = 's';
      if (write(event_pfd[1], &c, sizeof(c)) != 1)
      {
         perror("write()");
         exit(EXIT_FAILURE);
      }
   }

private:
   int event_pfd[2]; // pfd[0] -- reading, pfd[1] -- writing
   std::unordered_map<SDL_Window*, int> display_fds;
};

#endif // SDL_VIDEO_DRIVER_X11

#endif // GLVIS_SDL_X11_HPP
