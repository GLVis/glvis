// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
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
#include "aux_vis.hpp"
#include "gl/renderer_core.hpp"
#include "gl/renderer_ff.hpp"
#include "sdl.hpp"
#include "sdl_main.hpp"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/html5.h>
#endif

#ifdef SDL_VIDEO_DRIVER_COCOA
#include "sdl_mac.hpp"
#endif

using std::cerr;
using std::endl;

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

extern int GetMultisample();

SdlWindow::Handle::Handle(const std::string& title, int x, int y, int w, int h,
                          Uint32 wndflags)
   : hwnd(nullptr)
   , gl_ctx(0)
{
   hwnd = SDL_CreateWindow(title.c_str(), x, y, w, h, wndflags);
   if (!hwnd)
   {
      PRINT_DEBUG("SDL window creation failed with error: "
                  << SDL_GetError() << endl);
      return;
   }
   gl_ctx = SDL_GL_CreateContext(hwnd);
   if (!gl_ctx)
   {
      PRINT_DEBUG("OpenGL context creation failed with error: "
                  << SDL_GetError() << endl);
   }
}

SdlWindow::Handle::~Handle()
{
   if (gl_ctx)
   {
      SDL_GL_DeleteContext(gl_ctx);
   }
   if (hwnd)
   {
      SDL_DestroyWindow(hwnd);
   }
}

SdlMainThread& GetMainThread()
{
   static SdlMainThread inst;
   return inst;
}

bool SdlWindow::isGlInitialized()
{
   return (handle.gl_ctx != 0);
}

SdlWindow::SdlWindow() {}

void SdlWindow::StartSDL(bool server_mode)
{
   GetMainThread().MainLoop(server_mode);
}

const int default_dpi = 72;

bool SdlWindow::createWindow(const char* title, int x, int y, int w, int h,
                             bool legacyGlOnly)
{
#ifdef __EMSCRIPTEN__
   is_multithreaded = false;
#endif
   // create a new SDL window
   handle = GetMainThread().GetHandle(this, title, x, y, w, h, legacyGlOnly);

   // at this point, window should be up
   if (!handle.isInitialized())
   {
      return false;
   }

   window_id = SDL_GetWindowID(handle.hwnd);

   GLenum err = glewInit();
#ifdef GLEW_ERROR_NO_GLX_DISPLAY
   // NOTE: Hacky workaround for Wayland initialization failure
   // See https://github.com/nigels-com/glew/issues/172
   if (err == GLEW_ERROR_NO_GLX_DISPLAY)
   {
      cerr << "GLEW: No GLX display found. If you are using Wayland this can "
           << "be ignored." << endl;
      err = GLEW_OK;
   }
#endif
   if (err != GLEW_OK)
   {
      cerr << "FATAL: Failed to initialize GLEW: "
           << glewGetErrorString(err) << endl;
      return false;
   }

   // print versions
   PRINT_DEBUG("Using GLEW " << glewGetString(GLEW_VERSION) << std::endl);
   PRINT_DEBUG("Using GL " << glGetString(GL_VERSION) << std::endl);

   renderer.reset(new gl3::MeshRenderer);
   renderer->setSamplesMSAA(GetMultisample());
#ifndef __EMSCRIPTEN__
   if (!GLEW_VERSION_1_1)
   {
      cerr << "FATAL: Minimum of OpenGL 1.1 is required." << endl;
      return false;
   }
   if (!GLEW_VERSION_1_3)
   {
      // Multitexturing was introduced into the core OpenGL specification in
      // version 1.3; for versions before, we need to load the functions from
      // the ARB_multitexture extension.
      if (GLEW_ARB_multitexture)
      {
         glActiveTexture = glActiveTextureARB;
         glClientActiveTexture = glClientActiveTextureARB;
         glMultiTexCoord2f = glMultiTexCoord2fARB;
      }
      else
      {
         cerr << "FATAL: Missing OpenGL multitexture support." << endl;
         return false;
      }
   }
   if (!GLEW_VERSION_3_0 && GLEW_EXT_transform_feedback)
   {
      glBindBufferBase            = glBindBufferBaseEXT;
      // Use an explicit typecast to suppress an error from inconsistent types
      // that are present in older versions of GLEW.
      glTransformFeedbackVaryings =
         (PFNGLTRANSFORMFEEDBACKVARYINGSPROC)glTransformFeedbackVaryingsEXT;
      glBeginTransformFeedback    = glBeginTransformFeedbackEXT;
      glEndTransformFeedback      = glEndTransformFeedbackEXT;
   }
   if (!legacyGlOnly && (GLEW_VERSION_3_0
                         || (GLEW_VERSION_2_0 && GLEW_EXT_transform_feedback)))
   {
      // We require both shaders and transform feedback EXT_transform_feedback
      // was made core in OpenGL 3.0
      PRINT_DEBUG("Loading CoreGLDevice..." << endl);
      renderer->setDevice<gl3::CoreGLDevice>();
   }
   else
   {
      PRINT_DEBUG("Shader support missing, loading FFGLDevice..." << endl);
      renderer->setDevice<gl3::FFGLDevice>();
   }

#else
   renderer->setDevice<gl3::CoreGLDevice>();
#endif

   return true;
}

SdlWindow::~SdlWindow()
{
   // Let the main SDL thread delete the handles
   GetMainThread().DeleteHandle(std::move(handle));
}

void SdlWindow::windowEvent(SDL_WindowEvent& ew)
{
   switch (ew.event)
   {
      case SDL_WINDOWEVENT_EXPOSED:
      case SDL_WINDOWEVENT_RESIZED:
         update_before_expose = true;
         if (onExpose)
         {
            wnd_state = RenderState::ExposePending;
         }
         break;
      case SDL_WINDOWEVENT_CLOSE:
         running = false;
         break;
      case SDL_WINDOWEVENT_MOVED:
         update_before_expose = true;
         break;
      default:
         break;
   }
}

void SdlWindow::motionEvent(SDL_MouseMotionEvent& em)
{
   EventInfo info =
   {
      em.x, em.y,
      SDL_GetModState()
   };
   if (em.state & SDL_BUTTON_LMASK)
   {
      if (onMouseMove[SDL_BUTTON_LEFT])
      {
         onMouseMove[SDL_BUTTON_LEFT](&info);
      }
   }
   else if (em.state & SDL_BUTTON_RMASK)
   {
      if (onMouseMove[SDL_BUTTON_RIGHT])
      {
         onMouseMove[SDL_BUTTON_RIGHT](&info);
      }
   }
   else if (em.state & SDL_BUTTON_MMASK)
   {
      if (onMouseMove[SDL_BUTTON_MIDDLE])
      {
         onMouseMove[SDL_BUTTON_MIDDLE](&info);
      }
   }
}

void SdlWindow::mouseEventDown(SDL_MouseButtonEvent& eb)
{
   if (onMouseDown[eb.button])
   {
      EventInfo info =
      {
         eb.x, eb.y,
         SDL_GetModState()
      };
      onMouseDown[eb.button](&info);
   }
}

void SdlWindow::mouseEventUp(SDL_MouseButtonEvent& eb)
{
   if (onMouseUp[eb.button])
   {
      EventInfo info =
      {
         eb.x, eb.y,
         SDL_GetModState()
      };
      onMouseUp[eb.button](&info);
   }
}

void SdlWindow::keyDownEvent(SDL_Keysym& ks)
{
   // Some keyDown events will be followed by a textInput event which will
   // handle key translation due to Shift, CapsLock or AltGr, so we leave such
   // events to be processed there.
   // Note: the same condition has to be used in signalKeyDown().
   const char *scan_name = SDL_GetScancodeName(ks.scancode);
   if ((scan_name[0] >= 32 && scan_name[0] < 127) && scan_name[1] == '\0'
       && (ks.mod & (KMOD_CTRL | KMOD_LALT | KMOD_GUI)) == 0)
   {
      lastKeyDownProcessed = false;
      lastKeyDownMods = ks.mod;
      lastKeyDownChar = ks.sym;
      return;
   }
   // If any 'mod' key other than KMOD_SHIFT, KMOD_CAPS or KMOD_RALT is
   // pressed, or the key is not in the range [32,127) then we processed the
   // event here.
   lastKeyDownProcessed = true;
   if (onKeyDown[ks.sym])
   {
      onKeyDown[ks.sym](ks.mod);

      // Record the key in 'saved_keys':
      bool isAlt = ks.mod & (KMOD_ALT);
      bool isCtrl = ks.mod & (KMOD_CTRL);
      saved_keys += "[";
      if (isCtrl) { saved_keys += "C-"; }
      if (isAlt) { saved_keys += "Alt-"; }
      if (ks.sym >= 32 && ks.sym < 127)
      {
         saved_keys += (char)(ks.sym);
      }
      else
      {
         saved_keys += SDL_GetKeyName(ks.sym);
      }
      saved_keys += "]";
   }
}

void SdlWindow::textInputEvent(const SDL_TextInputEvent &tie)
{
   // This event follows a keyDown event where we've recorded if the event was
   // processed in keyDownEvent(). If it was not processed, we do it here.
   if (lastKeyDownProcessed) { return; }
   char c = tie.text[0];
   if (!onKeyDown[c])
   {
      // If the key was translated to something that is not handled, return to
      // the physical key passed in the keyDown event.
      c = lastKeyDownChar;
   }
   if (onKeyDown[c])
   {
      onKeyDown[c](lastKeyDownMods & ~(KMOD_CAPS | KMOD_LSHIFT | KMOD_RSHIFT));

      // Record the key in 'saved_keys':
      bool isAlt = lastKeyDownMods & (KMOD_ALT);
      bool isCtrl = lastKeyDownMods & (KMOD_CTRL);
      if (isAlt || isCtrl)
      {
         saved_keys += "[";
      }
      if (isCtrl) { saved_keys += "C-"; }
      if (isAlt) { saved_keys += "Alt-"; }
      saved_keys += c;
      if (isAlt || isCtrl)
      {
         saved_keys += "]";
      }
   }
}

void SdlWindow::multiGestureEvent(SDL_MultiGestureEvent & e)
{
   if (e.numFingers == 2)
   {
      if (onTouchPinch && fabs(e.dDist) > 0.00002)
      {
         onTouchPinch(e);
      }

      if (onTouchRotate)
      {
         onTouchRotate(e);
      }
   }
}

void SdlWindow::mainIter()
{
   if (!is_multithreaded)
   {
      // Pull events from GetMainThread() object
      GetMainThread().DispatchSDLEvents();
   }
   bool events_pending = false;
   bool sleep = false;
   {
      lock_guard<mutex> evt_guard{event_mutex};
      events_pending = !waiting_events.empty();
   }
   if (events_pending)
   {
      bool keep_going;
      do
      {
         SDL_Event e;
         // Fetch next event from the queue
         {
            lock_guard<mutex> evt_guard{event_mutex};
            e = waiting_events.front();
            waiting_events.pop_front();
            events_pending = !waiting_events.empty();
         }
         keep_going = false;
         switch (e.type)
         {
            case SDL_QUIT:
               running = false;
               break;
            case SDL_WINDOWEVENT:
               windowEvent(e.window);
               if (wnd_state != RenderState::ExposePending)
               {
                  keep_going = true;
               }
               break;
            case SDL_MULTIGESTURE:
               multiGestureEvent(e.mgesture);
               keep_going = true;
               break;
            case SDL_KEYDOWN:
               // For debugging: uncomment the next line to track key events.
               // #define TRACK_KEY_EVENTS
#ifdef TRACK_KEY_EVENTS
               cout << "Event: SDL_KEYDOWN sym=" << e.key.keysym.sym
                    << " mod=" << e.key.keysym.mod << endl;
#endif
               keyDownEvent(e.key.keysym);
               break;
            case SDL_KEYUP:
#ifdef TRACK_KEY_EVENTS
               cout << "Event: SDL_KEYUP sym=" << e.key.keysym.sym
                    << " mod=" << e.key.keysym.mod << endl;
#endif
               break;
            case SDL_TEXTINPUT:
#ifdef TRACK_KEY_EVENTS
               cout << "Event: SDL_TEXTINPUT text[0..3]="
                    << (int)(unsigned char)(e.text.text[0])
                    << ' ' << (int)(unsigned char)(e.text.text[1])
                    << ' ' << (int)(unsigned char)(e.text.text[2])
                    << ' ' << (int)(unsigned char)(e.text.text[3])
                    << " (as codes 0-255)" << endl;
#endif
               textInputEvent(e.text);
               break;
            case SDL_MOUSEMOTION:
               motionEvent(e.motion);
               // continue processing events
               keep_going = true;
               break;
            case SDL_MOUSEBUTTONDOWN:
               mouseEventDown(e.button);
               break;
            case SDL_MOUSEBUTTONUP:
               mouseEventUp(e.button);
               break;
         }
      }
      while (keep_going && events_pending);
   }
   else if (onIdle)
   {
      {
         unique_lock<mutex> event_lock{event_mutex};
         call_idle_func = false;
      }
      sleep = onIdle();
   }
   else
   {
      // No actions performed this iteration.
      sleep = true;
   }
   if (wnd_state == RenderState::ExposePending)
   {
#ifdef SDL_VIDEO_DRIVER_COCOA
      if (update_before_expose)
      {
         // On SDL, when the OpenGL context is on a separate thread from the
         // main thread, the call to [NSOpenGLContext update] after a resize or
         // move event is only scheduled for after the next swap event. Any
         // rendering/OpenGL commands occurring before this update will be
         // corrupted.
         //
         // To avoid this issue, we just call [NSOpenGLContext update]
         // immediately before the expose event.
         SdlCocoaPlatform* platform =
            dynamic_cast<SdlCocoaPlatform*>(GetMainThread().GetPlatform());
         if (platform)
         {
            platform->ContextUpdate();
         }
         update_before_expose = false;
      }
#endif
      onExpose();
      wnd_state = RenderState::SwapPending;
   }
   else if (is_multithreaded && sleep)
   {
      // No updates to vis, wait for next wakeup event from glvis_command or WM
      unique_lock<mutex> event_lock{event_mutex};
      events_available.wait(event_lock, [this]()
      {
         // Sleep until events from WM or glvis_command can be handled
         return !waiting_events.empty() || call_idle_func;
      });
   }
}

void SdlWindow::mainLoop()
{
   running = true;
#ifdef __EMSCRIPTEN__
   emscripten_set_main_loop_arg([](void* arg)
   {
      ((SdlWindow*) arg)->mainIter();
   }, this, 0, true);
#else
   while (running)
   {
      mainIter();
      if (takeScreenshot)
      {
         Screenshot(screenshot_file.c_str(), screenshot_convert);
         takeScreenshot = false;
      }
      if (wnd_state == RenderState::SwapPending)
      {
#ifdef SDL_VIDEO_DRIVER_COCOA
         // TODO: Temporary workaround - after merge, everyone should update to
         // latest SDL
         SdlCocoaPlatform* mac_platform
            = dynamic_cast<SdlCocoaPlatform*>(GetMainThread().GetPlatform());
         if (mac_platform && mac_platform->UseThreadWorkaround())
         {
            mac_platform->SwapWindow();
         }
         else
         {
            SDL_GL_SwapWindow(handle.hwnd);
         }
#else
         SDL_GL_SwapWindow(handle.hwnd);
#endif
         wnd_state = RenderState::Updated;
      }
   }
#endif
}

void SdlWindow::signalLoop()
{
   // Note: not executed from the main thread
   {
      lock_guard<mutex> evt_guard{event_mutex};
      call_idle_func = true;
   }
   if (is_multithreaded)
   {
      events_available.notify_all();
   }
}

void SdlWindow::getWindowSize(int& w, int& h)
{
   w = 0;
   h = 0;
   if (handle.isInitialized())
   {
#ifdef __EMSCRIPTEN__
      if (canvas_id_.empty())
      {
         std::cerr << "error: id is undefined: " << canvas_id_ << std::endl;
         return;
      }
      double dw, dh;
      auto err = emscripten_get_element_css_size(canvas_id_.c_str(), &dw, &dh);
      w = int(dw);
      h = int(dh);
      if (err != EMSCRIPTEN_RESULT_SUCCESS)
      {
         std::cerr << "error (emscripten_get_element_css_size): " << err << std::endl;
         return;
      }
#else
      SDL_GetWindowSize(handle.hwnd, &w, &h);
#endif
      w /= pixel_scale_x;
      h /= pixel_scale_y;
   }
}

void SdlWindow::getGLDrawSize(int& w, int& h)
{
   SDL_GL_GetDrawableSize(handle.hwnd, &w, &h);
}

void SdlWindow::getDpi(int& w, int& h)
{
   w = default_dpi;
   h = default_dpi;
   if (!handle.isInitialized())
   {
      PRINT_DEBUG("warning: unable to get dpi: handle is null" << endl);
      return;
   }
   int disp = SDL_GetWindowDisplayIndex(handle.hwnd);
   if (disp < 0)
   {
      PRINT_DEBUG("warning: problem getting display index: " << SDL_GetError()
                  << endl);
      PRINT_DEBUG("returning default dpi of " << default_dpi << endl);
      return;
   }

   float f_w, f_h;
   if (SDL_GetDisplayDPI(disp, NULL, &f_w, &f_h))
   {
      PRINT_DEBUG("warning: problem getting dpi, setting to " << default_dpi
                  << ": " << SDL_GetError() << endl);
   }
   else
   {
      PRINT_DEBUG("Screen DPI: w = " << f_w << " ppi, h = " << f_h << " ppi");
      w = f_w + 0.5f;
      h = f_h + 0.5f;
      PRINT_DEBUG(" (" << w << " x " << h << ')' << endl);
   }
}

void SdlWindow::setWindowTitle(std::string& title)
{
   setWindowTitle(title.c_str());
}

void SdlWindow::setWindowTitle(const char * title)
{
   GetMainThread().SetWindowTitle(handle, title);
}

void SdlWindow::setWindowSize(int w, int h)
{
   GetMainThread().SetWindowSize(handle, pixel_scale_x*w, pixel_scale_y*h);
   update_before_expose = true;

}

void SdlWindow::setWindowPos(int x, int y)
{
   bool uc_x = SDL_WINDOWPOS_ISUNDEFINED(x) ||
               SDL_WINDOWPOS_ISCENTERED(x);
   bool uc_y = SDL_WINDOWPOS_ISUNDEFINED(y) ||
               SDL_WINDOWPOS_ISCENTERED(y);
   GetMainThread().SetWindowPosition(handle,
                                     uc_x ? x : pixel_scale_x*x,
                                     uc_y ? y : pixel_scale_y*y);
   update_before_expose = true;
}

void SdlWindow::signalKeyDown(SDL_Keycode k, SDL_Keymod m)
{
   SDL_Event event;

   event.type = SDL_KEYDOWN;
   event.key.windowID = window_id;
   event.key.keysym.sym = k;
   event.key.keysym.scancode = SDL_GetScancodeFromKey(tolower(k));
   event.key.keysym.mod = m;
   queueEvents({ event });

   // The same condition as in keyDownEvent().
   if ((k >= 32 && k < 127) && (m & (KMOD_CTRL | KMOD_LALT | KMOD_GUI)) == 0)
   {
      event.type = SDL_TEXTINPUT;
      event.text.windowID = window_id;
      event.text.text[0] = k;
      queueEvents({ event });
   }
}

void SdlWindow::swapBuffer()
{
   SDL_GL_SwapWindow(handle.hwnd);
}
