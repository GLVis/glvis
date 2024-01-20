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

#include "sdl_x11.hpp"

#ifdef SDL_VIDEO_DRIVER_X11

#include <vector>
#include <array>
#include <iostream>

#ifndef __EMSCRIPTEN__
#include <SDL2/SDL_syswm.h>
#else
#include <SDL_syswm.h>
#endif

#ifdef SDL_VIDEO_DRIVER_X11_XINPUT2
#include <X11/extensions/XInput2.h>
#endif // SDL_VIDEO_DRIVER_X11_XINPUT2

using namespace std;

void SdlX11Platform::RegisterWindow(SDL_Window* window)
{
   SDL_SysWMinfo sysinfo;
   SDL_VERSION(&sysinfo.version);
   if (!SDL_GetWindowWMInfo(window, &sysinfo))
   {
      cerr << "Error: unable to get window manager information for the "
           << "current window." << endl;
      return;
   }
   if (sysinfo.subsystem != SDL_SYSWM_X11)
   {
      cerr << "Error: created SDL window is not an X11 window." << endl;
      return;
   }
   Display *disp = sysinfo.info.x11.display;
   Window wnd = sysinfo.info.x11.window;

   display_fds.emplace(window, ConnectionNumber(disp));

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
}

void SdlX11Platform::UnregisterWindow(SDL_Window* window)
{
   display_fds.erase(window);
}

void SdlX11Platform::WaitEvent()
{
   vector<pollfd> pfd;
   pfd.emplace_back( pollfd{ event_pfd[0], POLLIN, 0 } );

   for (auto wnd : display_fds)
   {
      // add file descriptors for the window event queue
      pfd.emplace_back( pollfd{ wnd.second, POLLIN, 0 } );
   }

   int nstr;
   do
   {
      // We timeout the poll call after 500ms, just in case a pending SDL_QUIT
      // needs to be pumped into the SDL event queue
      nstr = poll(pfd.data(), pfd.size(), 500);
   }
   while (nstr == -1 && errno == EINTR);

   if (nstr == -1) { perror("poll()"); }

   int n = 0;
   // Read out the pending GLVisCommand-sent events, if any
   do
   {
      array<char, 16> buf;
      n = read(event_pfd[0], buf.data(), buf.size());
   }
   while (n > 0);
}

#endif // SDL_VIDEO_DRIVER_X11
