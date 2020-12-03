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

#import <Cocoa/Cocoa.h>
#include "sdl_mac.hpp"

void SdlCocoaPlatform::WaitEvent()
{ @autoreleasepool
{
   [NSApp nextEventMatchingMask:NSEventMaskAny
          untilDate:[NSDate distantFuture]
          inMode:NSDefaultRunLoopMode
          dequeue:NO];
}}

void SdlCocoaPlatform::SendEvent()
{ @autoreleasepool
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
}}
