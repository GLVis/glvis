// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
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
#include "threads.hpp"
#include <poll.h>
#ifdef SDL_VIDEO_DRIVER_X11_XINPUT2
#include <X11/extensions/XInput2.h>
#endif // SDL_VIDEO_DRIVER_X11_XINPUT2

extern int visualize;

class SdlX11Platform final : public SdlNativePlatform
{
public:
   SdlX11Platform(Display* xdisplay, Window xwindow)
      : disp(xdisplay), wnd(xwindow)
   {
      // Disable XInput extension events since they are generated even outside
      // the GLVis window.
      Window root_win = DefaultRootWindow(disp);
      unsigned char mask[4] = {0,0,0,0};
      XIEventMask event_mask;
      event_mask.deviceid = XIAllMasterDevices;
      event_mask.mask_len = sizeof(mask);
      event_mask.mask = mask;
#ifdef SDL_VIDEO_DRIVER_X11_DYNAMIC_XINPUT2
      typedef int (*XISelectEvents_ptr)(Display *, Window, XIEventMask *, int);
      static XISelectEvents_ptr XISelectEvents_ = NULL;
      if (XISelectEvents_ == NULL)
      {
         void *lib = SDL_LoadObject(SDL_VIDEO_DRIVER_X11_DYNAMIC_XINPUT2);
         if (lib != NULL)
         {
            XISelectEvents_ =
               (XISelectEvents_ptr)SDL_LoadFunction(lib, "XISelectEvents");
         }
      }
      if (XISelectEvents_ == NULL)
      {
         cerr << "Error accessing XISelectEvents!" << endl;
         exit(EXIT_FAILURE);
      }
#else
#define XISelectEvents_ XISelectEvents
#endif
      if (XISelectEvents_(disp, root_win, &event_mask, 1) != Success)
      {
         cerr << "Failed to disable XInput on the default root window!" << endl;
      }
      if (XISelectEvents_(disp, wnd, &event_mask, 1) != Success)
      {
         cerr << "Failed to disable XInput on the current window!" << endl;
      }
   }
   void WaitEvent()
   {
      int nstr, nfd = 1;
      struct pollfd pfd[2];

      pfd[0].fd     = ConnectionNumber(disp);
      pfd[0].events = POLLIN;
      pfd[0].revents = 0;
      if (glvis_command && visualize == 1)
      {
         pfd[1].fd     = glvis_command->ReadFD();
         pfd[1].events = POLLIN;
         pfd[1].revents = 0;
         nfd = 2;
      }
      do
      {
         nstr = poll(pfd, nfd, -1);
      }
      while (nstr == -1 && errno == EINTR);

      if (nstr == -1) { perror("poll()"); }
   }
   void SendEvent() {}

private:
   Display* disp;
   Window wnd;
};

#endif // SDL_VIDEO_DRIVER_X11

#endif // GLVIS_SDL_X11_HPP
