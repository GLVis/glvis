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

#include <iostream>
#include <chrono>
#include "sdl.hpp"
#include "visual.hpp"
#include "logo.hpp"
#include "gl/renderer_core.hpp"
#include "gl/renderer_ff.hpp"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

using std::cerr;
using std::endl;

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

extern int GetMultisample();
extern int visualize;

struct SdlWindow::Handle
{
   SDL_Window * hwnd;
   SDL_GLContext gl_ctx;
   Handle(const std::string& title, int x, int y, int w, int h,
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

   ~Handle()
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

   bool isInitialized()
   {
      return (hwnd != nullptr && gl_ctx != 0);
   }
};

bool SdlWindow::isGlInitialized()
{
   return (handle->gl_ctx != 0);
}

SdlWindow::SdlWindow()
   : onIdle(nullptr)
   , onExpose(nullptr)
   , onReshape(nullptr)
   , ctrlDown(false)
   , requiresExpose(false)
   , takeScreenshot(false)
   , saved_keys("")
{
}

int SdlWindow::probeGLContextSupport()
{
   Uint32 win_flags_hidden = SDL_WINDOW_OPENGL | SDL_WINDOW_HIDDEN;
   PRINT_DEBUG("Testing if OpenGL core profile window can be created..." << flush);
   // Try to create core profile context first
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
   {
      Handle testCore("",
                      SDL_WINDOWPOS_UNDEFINED,
                      SDL_WINDOWPOS_UNDEFINED,
                      100, 100, win_flags_hidden);
      if (testCore.isInitialized())
      {
         PRINT_DEBUG("success!" << endl);
         return SDL_GL_CONTEXT_PROFILE_CORE;
      }
   }

   PRINT_DEBUG("failed." << endl);
   PRINT_DEBUG("Testing if OpenGL compatibility profile window can be created..."
               << flush);
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK,
                       SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);
   {
      Handle testCompat("",
                        SDL_WINDOWPOS_UNDEFINED,
                        SDL_WINDOWPOS_UNDEFINED,
                        100, 100, win_flags_hidden);
      if (testCompat.isInitialized())
      {
         PRINT_DEBUG("success!" << endl);
         return SDL_GL_CONTEXT_PROFILE_COMPATIBILITY;
      }
   }
   PRINT_DEBUG("failed." << endl);
   PRINT_DEBUG("No profile flags were accepted." << endl);
   return 0;
}

bool SdlWindow::createWindow(const char * title, int x, int y, int w, int h,
                             bool legacyGlOnly)
{
   if (!SDL_WasInit(SDL_INIT_VIDEO | SDL_INIT_EVENTS))
   {
      if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0)
      {
         cerr << "FATAL: Failed to initialize SDL: " << SDL_GetError() << endl;
         return false;
      }
   }

   // destroy any existing SDL window
   handle.reset();

   // If we want to use WebGL 2:
   // #ifdef __EMSCRIPTEN__
   // SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_ES);
   // SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
   // SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
   // #endif

   Uint32 win_flags = SDL_WINDOW_OPENGL;
   // Hide window until we adjust its size for high-dpi displays
   win_flags |= SDL_WINDOW_HIDDEN;
#ifndef __EMSCRIPTEN__
   win_flags |= SDL_WINDOW_ALLOW_HIGHDPI | SDL_WINDOW_RESIZABLE;
#endif
   SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1);
   SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 8);
   SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24);
   Uint32 win_gl_ctx = 0;
#ifndef __EMSCRIPTEN__
   if (!legacyGlOnly)
   {
      // Try and probe for a core/compatibility context. Needed for Mac OS X,
      // which will only support OpenGL 2.1 if you don't create a core context.
      win_gl_ctx = probeGLContextSupport();
   }
#endif
   SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, win_gl_ctx);
   if (GetMultisample() > 0)
   {
      SDL_GL_SetAttribute( SDL_GL_MULTISAMPLEBUFFERS, 1);
      SDL_GL_SetAttribute( SDL_GL_MULTISAMPLESAMPLES, GetMultisample());
      cerr << "Creating window..." << flush;
   }
   handle.reset(new Handle(title, x, y, w, h, win_flags));
   if (GetMultisample() > 0 && !handle->isInitialized())
   {
      // Antialiasing support might not be available on all platforms.
      PRINT_DEBUG("failed." << endl);
      PRINT_DEBUG("Disabling antialiasing and trying again..." << flush);
      SetMultisample(0);
      SDL_GL_SetAttribute( SDL_GL_MULTISAMPLEBUFFERS, 0);
      SDL_GL_SetAttribute( SDL_GL_MULTISAMPLESAMPLES, 0);
      handle.reset(new Handle(title, x, y, w, h, win_flags));
   }

   // at this point, window should be up
   if (!handle->isInitialized())
   {
      PRINT_DEBUG("failed." << endl);
      cerr << "FATAL: window and/or OpenGL context creation failed." << endl;
      return false;
   }
   else
   {
      PRINT_DEBUG("done." << endl);
   }

   SDL_Surface* iconSurf = SDL_CreateRGBSurfaceFrom(icon_pixels,
                                                    64, 64, // height, width
                                                    32,     // depth
                                                    64 * 4, // pitch
                                                    0xFF000000,
                                                    0x00FF0000,
                                                    0x0000FF00,
                                                    0x000000FF);
   if (iconSurf)
   {
      SDL_SetWindowIcon(handle->hwnd, iconSurf);
      SDL_FreeSurface(iconSurf);
   }

#ifndef __EMSCRIPTEN__
   SDL_GL_SetSwapInterval(0);
   glEnable(GL_DEBUG_OUTPUT);
#endif

   GLenum err = glewInit();
   if (err != GLEW_OK)
   {
      cerr << "FATAL: Failed to initialize GLEW: "
           << glewGetErrorString(err) << endl;
      return false;
   }

   // print versions
   PRINT_DEBUG("Using GLEW " << glewGetString(GLEW_VERSION) << std::endl);
   PRINT_DEBUG("Using GL " << glGetString(GL_VERSION) << std::endl);

   SDL_version sdl_ver;
   SDL_GetVersion(&sdl_ver);
   PRINT_DEBUG("Using SDL " << (int)sdl_ver.major << "." << (int)sdl_ver.minor
               << "." << (int)sdl_ver.patch << std::endl);

   renderer.reset(new gl3::MeshRenderer);
#ifndef __EMSCRIPTEN__
   if (GLEW_EXT_transform_feedback)
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

   // Detect if we are using a high-dpi display and resize the window unless it
   // was already resized by SDL's underlying backend.
   {
      int scr_w, scr_h, pix_w, pix_h, wdpi, hdpi;
      // SDL_GetWindowSize() -- size in "screen coordinates"
      SDL_GetWindowSize(handle->hwnd, &scr_w, &scr_h);
      // SDL_GL_GetDrawableSize() -- size in pixels
      SDL_GL_GetDrawableSize(handle->hwnd, &pix_w, &pix_h);
      high_dpi = false;
      pixel_scale_x = pixel_scale_y = 1.0f;
      float sdl_pixel_scale_x = 1.0f, sdl_pixel_scale_y = 1.0f;
      // If "screen" and "pixel" sizes are different, assume high-dpi and no
      // need to scale the window.
      if (scr_w == pix_w && scr_h == pix_h)
      {
         getDpi(wdpi, hdpi);
         if (std::max(wdpi, hdpi) >= high_dpi_threshold)
         {
            high_dpi = true;
            pixel_scale_x = pixel_scale_y = 2.0f;
            // the following two calls use 'pixel_scale_*'
            setWindowSize(w, h);
            setWindowPos(x, y);
         }
      }
      else
      {
         high_dpi = true;
         // keep 'pixel_scale_*' = 1, scaling is done inside SDL
         sdl_pixel_scale_x = float(pix_w)/scr_w;
         sdl_pixel_scale_y = float(pix_h)/scr_h;
      }
      if (high_dpi)
      {
         cout << "High-dpi display detected: using window scaling: "
              << sdl_pixel_scale_x*pixel_scale_x << " x "
              << sdl_pixel_scale_y*pixel_scale_y << endl;
      }
   }
   SDL_ShowWindow(handle->hwnd);

   return true;
}

// defined here because the Handle destructor needs to be visible
SdlWindow::~SdlWindow() {}

void SdlWindow::windowEvent(SDL_WindowEvent& ew)
{
   switch (ew.event)
   {
      case SDL_WINDOWEVENT_SIZE_CHANGED:
         if (onReshape)
         {
            onReshape(ew.data1, ew.data2);
         }
         break;
      case SDL_WINDOWEVENT_EXPOSED:
         if (onExpose)
         {
            requiresExpose = true;
         }
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

bool SdlWindow::keyEvent(SDL_Keysym& ks)
{
   bool handled = false;
   if (ks.sym >= 128 || ks.sym < 32)
   {
      if (onKeyDown[ks.sym])
      {
         onKeyDown[ks.sym](ks.mod);
         handled = true;
      }
   }
   else if (ctrlDown)
   {
      if (onKeyDown[ks.sym])
      {
         onKeyDown[ks.sym](ks.mod);
         handled = true;
      }
   }
   if (ks.sym == SDLK_RCTRL || ks.sym == SDLK_LCTRL)
   {
      ctrlDown = true;
   }
   if (handled)
   {
      bool isAlt = ks.mod & (KMOD_ALT);
      bool isCtrl = ks.mod & (KMOD_CTRL);
      saved_keys += "[";
      if (isCtrl) { saved_keys += "C-"; }
      if (isAlt) { saved_keys += "Alt-"; }
      if (std::isalpha(ks.sym))
      {
         //key with corresponding text output
         char c = ks.sym;
         if (!(ks.mod & KMOD_SHIFT)) { c = std::tolower(c); }
         saved_keys += c;
      }
      else
      {
         saved_keys += SDL_GetKeyName(ks.sym);
      }
      saved_keys += "]";
   }
   return handled;
}

bool SdlWindow::keyEvent(char c)
{
   if (onKeyDown[c])
   {
      SDL_Keymod mods = SDL_GetModState();
      bool isAlt = mods & (KMOD_ALT);
      bool isCtrl = mods & (KMOD_CTRL);
      onKeyDown[c](mods);
      if (isAlt || isCtrl)
      {
         saved_keys += "[";
         if (isCtrl) { saved_keys += "C-"; }
         if (isAlt) { saved_keys += "Alt-"; }
      }
      saved_keys += c;
      if (isAlt || isCtrl)
      {
         saved_keys += "]";
      }
      return true;
   }
   return false;
}

bool SdlWindow::mainIter()
{
   SDL_Event e;
   static bool useIdle = false;
   bool needsSwap = false;
   while (SDL_PollEvent(&e))
   {
      bool renderKeyEvent = false;
      switch (e.type)
      {
         case SDL_QUIT:
            running = false;
            break;
         case SDL_WINDOWEVENT:
            windowEvent(e.window);
            break;
         case SDL_KEYDOWN:
            renderKeyEvent = keyEvent(e.key.keysym);
            break;
         case SDL_KEYUP:
            if (e.key.keysym.sym == SDLK_LCTRL
                || e.key.keysym.sym == SDLK_RCTRL)
            {
               ctrlDown = false;
            }
            break;
         case SDL_TEXTINPUT:
            renderKeyEvent = keyEvent(e.text.text[0]);
            break;
         case SDL_MOUSEMOTION:
            motionEvent(e.motion);
            break;
         case SDL_MOUSEBUTTONDOWN:
            mouseEventDown(e.button);
            break;
         case SDL_MOUSEBUTTONUP:
            mouseEventUp(e.button);
            break;
      }
      if (renderKeyEvent)
      {
         break;
      }
   }
#ifndef __EMSCRIPTEN__
   if (onIdle)
   {
      if (glvis_command == NULL || visualize == 2 || useIdle)
      {
         onIdle();
         needsSwap = true;
      }
      else
      {
         if (glvis_command->Execute() < 0)
         {
            running = false;
         }
      }
      useIdle = !useIdle;
   }
   else if (glvis_command && visualize == 1)
   {
      if (glvis_command->Execute() < 0) { running = false; }
   }
#else
   if (onIdle)
   {
      onIdle();
      needsSwap = true;
   }
#endif
   if (requiresExpose)
   {
      onExpose();
      requiresExpose = false;
      needsSwap = true;
   }
   return needsSwap;
}

void SdlWindow::mainLoop()
{
   running = true;
#ifdef __EMSCRIPTEN__
   emscripten_set_main_loop_arg([](void* arg)
   {
      ((SdlWindow*) arg)->mainIter();
   }, this, 0, 1);
#else
   visualize = 1;
   while (running)
   {
      bool glSwap = mainIter();
      if (glSwap)
      {
         SDL_GL_SwapWindow(handle->hwnd);
      }
      if (takeScreenshot)
      {
         Screenshot(screenshot_file.c_str());
         takeScreenshot = false;
      }
      // sleep for n milliseconds to avoid pegging CPU at 100%
      SDL_WaitEventTimeout(NULL, 10);
   }
#endif
}

void SdlWindow::getWindowSize(int& w, int& h)
{
   if (handle)
   {
#ifdef __EMSCRIPTEN__
      int is_fullscreen;
      emscripten_get_canvas_size(&w, &h, &is_fullscreen);
      // TODO: ^ is deprecated but we need to store the id somewhere
      /*
      EMSCRIPTEN_RESULT r = emscripten_get_canvas_element_size("#canvas", &w, &h);
      if (r != EMSCRIPTEN_RESULT_SUCCESS) {
        std::cerr << "emscripten error" << std::endl;
      }
      */
#else
      SDL_GetWindowSize(handle->hwnd, &w, &h);
#endif
      w /= pixel_scale_x;
      h /= pixel_scale_y;
   }
}

void SdlWindow::getGLDrawSize(int& w, int& h)
{
   SDL_GL_GetDrawableSize(handle->hwnd, &w, &h);
}

const int default_dpi = 72;
void SdlWindow::getDpi(int& w, int& h)
{
   w = default_dpi;
   h = default_dpi;
   if (!handle)
   {
      PRINT_DEBUG("warning: unable to get dpi: handle is null" << endl);
      return;
   }
   int disp = SDL_GetWindowDisplayIndex(handle->hwnd);
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
   if (handle)
   {
      SDL_SetWindowTitle(handle->hwnd, title);
   }
}

void SdlWindow::setWindowSize(int w, int h)
{
   if (handle)
   {
      SDL_SetWindowSize(handle->hwnd, pixel_scale_x*w, pixel_scale_y*h);
   }
}

void SdlWindow::setWindowPos(int x, int y)
{
   if (handle)
   {
      bool uc_x = SDL_WINDOWPOS_ISUNDEFINED(x) ||
                  SDL_WINDOWPOS_ISCENTERED(x);
      bool uc_y = SDL_WINDOWPOS_ISUNDEFINED(y) ||
                  SDL_WINDOWPOS_ISCENTERED(y);
      SDL_SetWindowPosition(handle->hwnd,
                            uc_x ? x : pixel_scale_x*x,
                            uc_y ? y : pixel_scale_y*y);
   }
}

void SdlWindow::signalKeyDown(SDL_Keycode k, SDL_Keymod m)
{
   SDL_Event event;
   if (k >= 32 && k < 128)
   {
      event.type = SDL_TEXTINPUT;
      event.text.text[0] = k;
   }
   else
   {
      event.type = SDL_KEYDOWN;
      event.key.keysym.sym = k;
      event.key.keysym.mod = m;
   }
   SDL_PushEvent(&event);
}

void SdlWindow::swapBuffer()
{
   SDL_GL_SwapWindow(handle->hwnd);
}
