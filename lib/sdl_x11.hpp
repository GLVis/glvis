// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
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

#ifdef SDL_VIDEO_DRIVER_X11

#include "sdl_helper.hpp"
#include "gl/platform_gl.hpp"
#include <poll.h>

#include <unistd.h>    // pipe, fcntl, write
#include <fcntl.h>     // fcntl
#include <cerrno>      // errno, EAGAIN
#include <cstdio>      // perror
#ifdef SDL_VIDEO_DRIVER_X11_XINPUT2
#include <X11/extensions/XInput2.h>
#endif // SDL_VIDEO_DRIVER_X11_XINPUT2

class SdlX11Platform final : public SdlNativePlatform
{
public:
   SdlX11Platform(Display* xdisplay, Window xwindow)
      : disp(xdisplay), wnd(xwindow)
   {
#ifdef SDL_VIDEO_DRIVER_X11_XINPUT2
      // Disable XInput extension events since they are generated even outside
      // the GLVis window.
      Window root_win = DefaultRootWindow(disp);
      unsigned char mask[4] = {0,0,0,0};
      XIEventMask event_mask;
      event_mask.deviceid = XIAllMasterDevices;
      event_mask.mask_len = sizeof(mask);
      event_mask.mask = mask;
#ifdef SDL_VIDEO_DRIVER_X11_DYNAMIC_XINPUT2
      const char Xi_lib[] = SDL_VIDEO_DRIVER_X11_DYNAMIC_XINPUT2;
#else
      const char Xi_lib[] = "libXi.so";
#endif
      typedef int (*XISelectEvents_ptr)(Display *, Window, XIEventMask *, int);
      XISelectEvents_ptr XISelectEvents_ = NULL;
      void *lib = SDL_LoadObject(Xi_lib);
      if (lib != NULL)
      {
         XISelectEvents_ =
            (XISelectEvents_ptr)SDL_LoadFunction(lib, "XISelectEvents");
      }
      if (XISelectEvents_ == NULL)
      {
         cerr << "Error accessing XISelectEvents!" << endl;
         exit(EXIT_FAILURE);
      }
      if (XISelectEvents_(disp, root_win, &event_mask, 1) != Success)
      {
         cerr << "Failed to disable XInput on the default root window!" << endl;
      }
      if (XISelectEvents_(disp, wnd, &event_mask, 1) != Success)
      {
         cerr << "Failed to disable XInput on the current window!" << endl;
      }
#ifndef SDL_VIDEO_DRIVER_X11_DYNAMIC_XINPUT2
      SDL_UnloadObject(lib);
#endif
#endif // SDL_VIDEO_DRIVER_X11_XINPUT2

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

   void WaitEvent()
   {
      int nstr, nfd = 1;
      struct pollfd pfd[2];

      pfd[0].fd     = ConnectionNumber(disp);
      pfd[0].events = POLLIN;
      pfd[0].revents = 0;

      pfd[1].fd     = event_pfd[0];
      pfd[1].events = POLLIN;
      pfd[1].revents = 0;

      do
      {
         nstr = poll(pfd, nfd, -1);
      }
      while (nstr == -1 && errno == EINTR);

      if (nstr == -1) { perror("poll()"); }

      int n = 0;
      // Read out the pending GLVisCommand-sent events, if any
      do
      {
          char c[10];
          n = read(event_pfd[0], c, 10);
      }
      while (n > 0);
   }
   void SendEvent()
   {
       char c = 's';
       if (write(pfd[1], &c, 1) != 1)
       {
         perror("write()");
         exit(EXIT_FAILURE);
       }
   }

private:
   Display* disp;
   Window wnd;
   int event_pfd[2]; // pfd[0] -- reading, pfd[1] -- writing
};

#endif // SDL_VIDEO_DRIVER_X11

#endif // GLVIS_SDL_X11_HPP
