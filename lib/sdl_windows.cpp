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

#include "sdl_windows.hpp"

#ifdef SDL_VIDEO_DRIVER_WINDOWS

#include <SDL2/SDL_syswm.h>

#include <iostream>

using namespace std;

struct SdlWindowsPlatform::Impl
{
   HANDLE signal_evt;
};

SdlWindowsPlatform::SdlWindowsPlatform()
    : m_impl(new Impl)
{
   m_impl->signal_evt = CreateEventA(nullptr, false, false, "");
   if (m_impl->signal_evt == NULL)
   {
      cerr << "Error: CreateEventA() failed with code " << GetLastError() << endl;
   }
}

SdlWindowsPlatform::~SdlWindowsPlatform()
{
   CloseHandle(m_impl->signal_evt);
}

void SdlWindowsPlatform::WaitEvent()
{
   // This call either waits for an event from Windows to be posted to the main
   // thread's message queue, or for the event object to be set by one of the
   // worker threads.
   // A timeout of 500ms is set to allow for interrupts to be pumped into the
   // SDL event queue.
   DWORD retval = MsgWaitForMultipleObjectsEx(1,
                                              &(m_impl->signal_evt),
                                              500,
                                              QS_ALLINPUT,
                                              MWMO_INPUTAVAILABLE);
   if (retval == WAIT_FAILED)
   {
      cerr << "Error: MsgWaitForMultipleObjectsEx() failed with code "
           << GetLastError() << endl;
   }
   if (ResetEvent(m_impl->signal_evt) == NULL)
   {
      cerr << "Error: ResetEvent() failed with code " << GetLastError() << endl;
   }
}

void SdlWindowsPlatform::SendEvent()
{
   if (SetEvent(m_impl->signal_evt) == NULL)
   {
      cerr << "Error: SetEvent() failed with code " << GetLastError() << endl;
   }
}

#endif // SDL_VIDEO_DRIVER_WINDOWS
