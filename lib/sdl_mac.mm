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

#include "sdl_mac.hpp"
#import <Cocoa/Cocoa.h>

#include <unordered_map>
#include <mutex>

#ifdef __clang__
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

std::unordered_map<int, NSOpenGLContext*>& GetContextMap()
{
   static std::unordered_map<int, NSOpenGLContext*> inst;
   return inst;
}

inline bool AtLeastVersion(SDL_version sdl_ver, int major, int minor, int patch)
{
   return ((sdl_ver.major > major) ||
           (sdl_ver.major == major && sdl_ver.minor > minor) ||
           (sdl_ver.major == major && sdl_ver.minor == minor && sdl_ver.patch >= patch));
}

void SdlCocoaPlatform::WaitEvent()
{
   @autoreleasepool
   {
      // We timeout the poll call after 500ms, just in case a pending SDL_QUIT
      // needs to be pumped into the SDL event queue
      [NSApp nextEventMatchingMask:NSEventMaskAny
             untilDate:[NSDate dateWithTimeIntervalSinceNow:0.500]
             inMode:NSDefaultRunLoopMode
             dequeue:NO];
   }
}

void SdlCocoaPlatform::SendEvent()
{
   @autoreleasepool
   {
      NSPoint loc = {0., 0.};
      [NSApp postEvent:[NSEvent otherEventWithType:NSEventTypeApplicationDefined
                        location:loc
                        modifierFlags:0
                        timestamp:0.0
                        windowNumber:0
                        context:nil
                        subtype:0
                        data1:0
                        data2:0]
             atStart:NO];
   }
}

bool SdlCocoaPlatform::UseThreadWorkaround() const
{
   static bool first_call = true;
   static bool value = false;
   if (first_call)
   {
      SDL_version sdl_ver;
      SDL_GetVersion(&sdl_ver);
      value = !AtLeastVersion(sdl_ver, 2, 0, 14);
      first_call = false;
   }
   return value;
}

void SdlCocoaPlatform::ContextUpdate()
{
   @autoreleasepool
   {
      NSOpenGLContext* ctx = [NSOpenGLContext currentContext];
      // This calls [SDLOpenGLContext update] defined in SDL_cocoaopengl.m
      dispatch_sync(dispatch_get_main_queue(), ^{ [ctx update]; });
   }
}

void SdlCocoaPlatform::ClearCurrentContext(int wnd_id)
{
   @autoreleasepool
   {
      NSOpenGLContext* ctx = [NSOpenGLContext currentContext];
      GetContextMap().emplace(wnd_id, ctx);
      [NSOpenGLContext clearCurrentContext];
   }
}

void SdlCocoaPlatform::SetCurrentContext(int wnd_id)
{
   @autoreleasepool
   {
      NSOpenGLContext* ctx = GetContextMap()[wnd_id];
      [ctx makeCurrentContext];
   }
}

std::mutex swap_mtx;

void SdlCocoaPlatform::SwapWindow()
{
   @autoreleasepool
   {
      std::lock_guard<std::mutex> lk{swap_mtx};
      [[NSOpenGLContext currentContext] flushBuffer];
   }
}
