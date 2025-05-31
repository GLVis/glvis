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
#include <sstream>
#include <cmath>
#include <chrono>
#include <regex>

#include "gl/types.hpp"
#include "gl2ps.h"
#include "palettes.hpp"
#include "sdl.hpp"
#include "threads.hpp"

#include <mfem.hpp>

#if defined(GLVIS_USE_LIBTIFF)
#include "tiffio.h"
#elif defined(GLVIS_USE_LIBPNG)
#include <png.h>
#endif

#include "font.hpp"
#ifndef __EMSCRIPTEN__
#include <fontconfig/fontconfig.h>
#endif

thread_local int visualize = 0;
thread_local VisualizationScene * locscene;
thread_local GLVisCommand *glvis_command = nullptr;

#ifdef GLVIS_MULTISAMPLE
static int glvis_multisample = GLVIS_MULTISAMPLE;
#else
static int glvis_multisample = -1;
#endif

float line_w = 1.f;
float line_w_aa = gl3::LINE_WIDTH_AA;

thread_local SdlWindow * wnd = nullptr;
bool wndLegacyGl = false;
bool wndUseHiDPI = true;
void SDLMainLoop(bool server_mode)
{
   SdlWindow::StartSDL(server_mode);
}

SdlWindow * GetAppWindow()
{
   return wnd;
}

VisualizationScene * GetVisualizationScene()
{
   return locscene;
}

void SetLegacyGLOnly(bool status)
{
   wndLegacyGl = true;
}

void SetUseHiDPI(bool status)
{
   wndUseHiDPI = status;
}

void MyExpose(GLsizei w, GLsizei h);
void MyExpose();

int InitVisualization (const char name[], int x, int y, int w, int h)
{

#ifdef GLVIS_DEBUG
   cout << "OpenGL Visualization" << endl;
#endif
   if (!wnd)
   {
      wnd = new SdlWindow();
      if (!wnd->createWindow(name, x, y, w, h, wndLegacyGl))
      {
         return 1;
      }
   }
   else
   {
      wnd->clearEvents();
   }

#ifdef GLVIS_DEBUG
   cout << "Window should be up" << endl;
#endif
   InitFont();
   wnd->getRenderer().setLineWidth(line_w);
   wnd->getRenderer().setLineWidthMS(line_w_aa);

   // auxReshapeFunc (MyReshape); // not needed, MyExpose calls it
   // auxReshapeFunc (NULL);
   void (*exposeFunc)(void) = MyExpose;
   wnd->setOnExpose(exposeFunc);

   wnd->setOnMouseDown(SDL_BUTTON_LEFT, LeftButtonDown);
   wnd->setOnMouseUp(SDL_BUTTON_LEFT, LeftButtonUp);
   wnd->setOnMouseMove(SDL_BUTTON_LEFT, LeftButtonLoc);
   wnd->setOnMouseDown(SDL_BUTTON_MIDDLE, MiddleButtonDown);
   wnd->setOnMouseUp(SDL_BUTTON_MIDDLE, MiddleButtonUp);
   wnd->setOnMouseMove(SDL_BUTTON_MIDDLE, MiddleButtonLoc);
   wnd->setOnMouseDown(SDL_BUTTON_RIGHT, RightButtonDown);
   wnd->setOnMouseUp(SDL_BUTTON_RIGHT, RightButtonUp);
   wnd->setOnMouseMove(SDL_BUTTON_RIGHT, RightButtonLoc);

   wnd->setTouchPinchCallback(TouchPinch);

   // auxKeyFunc (AUX_p, KeyCtrlP); // handled in vsdata.cpp
   wnd->setOnKeyDown (SDLK_s, KeyS);
   wnd->setOnKeyDown ('S', KeyS);

   wnd->setOnKeyDown (SDLK_q, KeyQPressed);
   // wnd->setOnKeyDown (SDLK_Q, KeyQPressed);

   wnd->setOnKeyDown (SDLK_LEFT, KeyLeftPressed);
   wnd->setOnKeyDown (SDLK_RIGHT, KeyRightPressed);
   wnd->setOnKeyDown (SDLK_UP, KeyUpPressed);
   wnd->setOnKeyDown (SDLK_DOWN, KeyDownPressed);

   wnd->setOnKeyDown (SDLK_KP_0, Key0Pressed);
   wnd->setOnKeyDown (SDLK_KP_1, Key1Pressed);
   wnd->setOnKeyDown (SDLK_KP_2, Key2Pressed);
   wnd->setOnKeyDown (SDLK_KP_3, Key3Pressed);
   wnd->setOnKeyDown (SDLK_KP_4, Key4Pressed);
   wnd->setOnKeyDown (SDLK_KP_5, Key5Pressed);
   wnd->setOnKeyDown (SDLK_KP_6, Key6Pressed);
   wnd->setOnKeyDown (SDLK_KP_7, Key7Pressed);
   wnd->setOnKeyDown (SDLK_KP_8, Key8Pressed);
   wnd->setOnKeyDown (SDLK_KP_9, Key9Pressed);

   wnd->setOnKeyDown (SDLK_KP_MEMSUBTRACT, KeyMinusPressed);
   wnd->setOnKeyDown (SDLK_KP_MEMADD, KeyPlusPressed);

   wnd->setOnKeyDown (SDLK_KP_DECIMAL, KeyDeletePressed);
   wnd->setOnKeyDown (SDLK_KP_ENTER, KeyEnterPressed);

   wnd->setOnKeyDown (SDLK_PERIOD, KeyDeletePressed);
   wnd->setOnKeyDown (SDLK_RETURN, KeyEnterPressed);

   wnd->setOnKeyDown (SDLK_0, Key0Pressed);
   wnd->setOnKeyDown (SDLK_1, Key1Pressed);
   wnd->setOnKeyDown (SDLK_2, Key2Pressed);
   wnd->setOnKeyDown (SDLK_3, Key3Pressed);
   wnd->setOnKeyDown (SDLK_4, Key4Pressed);
   wnd->setOnKeyDown (SDLK_5, Key5Pressed);
   wnd->setOnKeyDown (SDLK_6, Key6Pressed);
   wnd->setOnKeyDown (SDLK_7, Key7Pressed);
   wnd->setOnKeyDown (SDLK_8, Key8Pressed);
   wnd->setOnKeyDown (SDLK_9, Key9Pressed);

   wnd->setOnKeyDown (SDLK_MINUS, KeyMinusPressed);
   wnd->setOnKeyDown (SDLK_PLUS, KeyPlusPressed);
   wnd->setOnKeyDown (SDLK_EQUALS, KeyPlusPressed);

   wnd->setOnKeyDown (SDLK_j, KeyJPressed);
   // wnd->setOnKeyDown (AUX_J, KeyJPressed);

   wnd->setOnKeyDown (SDLK_KP_MULTIPLY, ZoomIn);
   wnd->setOnKeyDown (SDLK_KP_DIVIDE, ZoomOut);

   wnd->setOnKeyDown (SDLK_ASTERISK, ZoomIn);
   wnd->setOnKeyDown (SDLK_SLASH, ZoomOut);

   wnd->setOnKeyDown (SDLK_LEFTBRACKET, ScaleDown);
   wnd->setOnKeyDown (SDLK_RIGHTBRACKET, ScaleUp);
   wnd->setOnKeyDown (SDLK_AT, LookAt);

#ifndef __EMSCRIPTEN__
   wnd->setOnKeyDown(SDLK_LEFTPAREN, ShrinkWindow);
   wnd->setOnKeyDown(SDLK_RIGHTPAREN, EnlargeWindow);

   if (locscene)
   {
      delete locscene;
   }
#endif
   locscene = nullptr;

   return 0;
}

void SendKeySequence(const char *seq)
{
   for (const char* key = seq; *key != '\0'; key++)
   {
      if (*key == '~')
      {
         key++;
         switch (*key)
         {
            case 'e': // expose event
               SendExposeEvent();
               break;
            case 'l': // left arrow
               wnd->signalKeyDown(SDLK_LEFT);
               break;
            case 'r': // right arrow
               wnd->signalKeyDown(SDLK_RIGHT);
               break;
            case 'u': // up arrow
               wnd->signalKeyDown(SDLK_UP);
               break;
            case 'd': // down arrow
               wnd->signalKeyDown(SDLK_DOWN);
               break;
            case '3': // F3
               wnd->signalKeyDown(SDLK_F3);
               break;
            case '5': // F5
               wnd->signalKeyDown(SDLK_F5);
               break;
            case '6': // F6
               wnd->signalKeyDown(SDLK_F6);
               break;
            case '7': // F7
               wnd->signalKeyDown(SDLK_F7);
               break;
            case '1': // F11 or F12
               key++;
               switch (*key)
               {
                  case '1': // F11
                     wnd->signalKeyDown(SDLK_F11);
                     break;
                  case '2': // F12
                     wnd->callKeyDown(SDLK_F12);
                     break;
               }
               break;
            case '.': // Keypad ./Del
               wnd->signalKeyDown(SDLK_PERIOD);
               break;
            case 'E': // Keypad Enter
               wnd->signalKeyDown(SDLK_RETURN);
               break;
         }
         continue;
      }
      else
      {
         wnd->signalKeyDown(*key);
      }
   }
}


thread_local bool disableSendExposeEvent = false;

void CallKeySequence(const char *seq)
{
   const char *key = seq;

   disableSendExposeEvent = true;
   for ( ; *key != '\0'; key++ ) // see /usr/include/X11/keysymdef.h
   {
      if (*key != '~')
      {
         wnd->callKeyDown(*key);
      }
      else
      {
         key++;
         switch (*key)
         {
            case 'l': // left arrow
               wnd->callKeyDown(SDLK_LEFT);
               break;
            case 'r': // right arrow
               wnd->callKeyDown(SDLK_RIGHT);
               break;
            case 'u': // up arrow
               wnd->callKeyDown(SDLK_UP);
               break;
            case 'd': // down arrow
               wnd->callKeyDown(SDLK_DOWN);
               break;
            case '3': // F3
               wnd->callKeyDown(SDLK_F3);
               break;
            case '5': // F5
               wnd->callKeyDown(SDLK_F5);
               break;
            case '6': // F6
               wnd->callKeyDown(SDLK_F6);
               break;
            case '7': // F7
               wnd->callKeyDown(SDLK_F7);
               break;
            case '1': // F11 or F12
               key++;
               switch (*key)
               {
                  case '1': // F11
                     wnd->callKeyDown(SDLK_F11);
                     break;
                  case '2': // F12
                     wnd->callKeyDown(SDLK_F12);
                     break;
               }
               break;
            case '.': // Keypad ./Del
               wnd->callKeyDown(SDLK_PERIOD);
               break;
            case 'E': // Keypad Enter
               wnd->callKeyDown(SDLK_RETURN);
               break;
         }
      }
   }
   disableSendExposeEvent = false;
}

void InitIdleFuncs();

void SetVisualizationScene(VisualizationScene * scene, int view,
                           const char *keys)
{
   locscene = scene;
   locscene -> view = view;
   if (view == 2)
   {
      scene -> CenterObject2D();
   }
   else
   {
      scene -> CenterObject();
   }

   InitIdleFuncs();
   if (scene -> spinning)
   {
      AddIdleFunc(MainLoop);
   }

   if (keys)
   {
      SendKeySequence(keys);
   }
   wnd->getRenderer().setPalette(&locscene->palette);
}

void RunVisualization()
{
   visualize = 1;
#ifndef __EMSCRIPTEN__
   wnd->mainLoop();
#endif
   InitIdleFuncs();
   delete locscene;
   delete wnd;
   wnd = nullptr;
}

void SendExposeEvent()
{
   if (disableSendExposeEvent) { return; }
   wnd->signalExpose();
}

void MyReshape(GLsizei w, GLsizei h)
{
   wnd->getRenderer().setViewport(w, h);

   gl3::GlMatrix projmtx;
   projmtx.identity();

   double ViewCenterX = locscene->ViewCenterX;
   double ViewCenterY = locscene->ViewCenterY;

   if (locscene->OrthogonalProjection)
   {
      double scale = locscene->ViewScale;
      if (w <= h)
      {
         projmtx.ortho(-1.0, 1.0, -double(h)/w, double(h)/w, -10, 10);
      }
      else
      {
         projmtx.ortho(-double(w)/h, double(w)/h, -1, 1, -10, 10);
      }
      projmtx.scale(scale, scale, 1.0);
   }
   else
   {
      double ViewAngle = locscene->ViewAngle;

      if (w < h)
         ViewAngle = (360.0/M_PI)*atan( tan( ViewAngle*(M_PI/360.0) ) *
                                        double(h)/w );

      projmtx.perspective(ViewAngle, double(w)/h, 0.1, 5.0);
   }
   projmtx.translate(-ViewCenterX, -ViewCenterY, 0.0);
   locscene->SetProjectionMtx(projmtx.mtx);
}

void MyExpose(GLsizei w, GLsizei h)
{
   MyReshape (w, h);
   GLuint color_tex = locscene->palette.GetColorTexture();
   GLuint alpha_tex = locscene->palette.GetAlphaTexture();
   wnd->getRenderer().setColorTexture(color_tex);
   wnd->getRenderer().setAlphaTexture(alpha_tex);
   gl3::SceneInfo frame = locscene->GetSceneObjs();
   for (auto drawable_ptr : frame.needs_buffering)
   {
      wnd->getRenderer().buffer(drawable_ptr);
   }
   wnd->getRenderer().render(frame.queue);
}

void MyExpose()
{
   int w, h;
   wnd->getGLDrawSize(w, h);
   MyExpose(w, h);
   wnd->signalSwap();
}


thread_local mfem::Array<void (*)()> IdleFuncs;
thread_local int LastIdleFunc;
thread_local bool use_idle = false;

bool MainIdleFunc();

void InitIdleFuncs()
{
   IdleFuncs.SetSize(0);
   LastIdleFunc = 0;
   if (glvis_command)
   {
      wnd->setOnIdle(MainIdleFunc);
   }
}

bool CommunicationIdleFunc()
{
   int status = glvis_command->Execute();
   if (status < 0)
   {
      cout << "GLVisCommand signalled exit" << endl;
      wnd->signalQuit();
   }
   else if (status == 1)
   {
      // no commands right now - main loop should sleep
      return true;
   }
   return false;
}

bool MainIdleFunc()
{
   bool sleep = true;
#ifndef __EMSCRIPTEN__
   if (glvis_command && visualize == 1
       && !(IdleFuncs.Size() > 0 && use_idle))
   {
      // Execute the next event from the communication thread if:
      //  - a valid GLVisCommand has been set
      //  - the communication thread is not stopped
      //  - The idle function flag is not set, or no idle functions have been
      //    registered
      sleep = CommunicationIdleFunc();
      if (IdleFuncs.Size() > 0) { sleep = false; }
   }
   else if (IdleFuncs.Size() > 0)
   {
      LastIdleFunc = (LastIdleFunc + 1) % IdleFuncs.Size();
      if (IdleFuncs[LastIdleFunc])
      {
         (*IdleFuncs[LastIdleFunc])();
      }
      // Continue executing idle functions
      sleep = false;
   }
   use_idle = !use_idle;
#else
   if (IdleFuncs.Size() > 0)
   {
      LastIdleFunc = (LastIdleFunc + 1) % IdleFuncs.Size();
      if (IdleFuncs[LastIdleFunc])
      {
         (*IdleFuncs[LastIdleFunc])();
      }
      // Continue executing idle functions
      sleep = false;
   }
#endif
   return sleep;
   LastIdleFunc = (LastIdleFunc + 1) % IdleFuncs.Size();
   if (IdleFuncs[LastIdleFunc])
   {
      (*IdleFuncs[LastIdleFunc])();
   }
}

void AddIdleFunc(void (*Func)(void))
{
   IdleFuncs.Union(Func);
   wnd->setOnIdle(MainIdleFunc);
}

void RemoveIdleFunc(void (*Func)(void))
{
   IdleFuncs.DeleteFirst(Func);
   if (IdleFuncs.Size() == 0 && glvis_command == nullptr)
   {
      wnd->setOnIdle(NULL);
   }
}


thread_local double xang = 0., yang = 0.;
thread_local gl3::GlMatrix srot;
thread_local double sph_t, sph_u;
thread_local GLint oldx, oldy, startx, starty;

thread_local int constrained_spinning = 0;


void MainLoop()
{
   static int p = 1;
   if (locscene->spinning)
   {
      if (!constrained_spinning)
      {
         locscene->Rotate(xang, yang);
         SendExposeEvent();
      }
      else
      {
         locscene->PreRotate(xang, 0.0, 0.0, 1.0);
         SendExposeEvent();
      }
      std::this_thread::sleep_for(std::chrono::milliseconds{10}); // sleep for 0.01 seconds
   }
   if (locscene->movie)
   {
      char fname[20];
      snprintf(fname, 20, "GLVis_m%04d", p++);
      wnd->screenshot(fname);
   }
}

// Pressed mouse buttons events

inline double sqr(double t)
{
   return t*t;
}

inline void ComputeSphereAngles(int &newx, int &newy,
                                double &new_sph_u, double &new_sph_t)
{
   GLint viewport[4] = { 0, 0, 0, 0 };
   double r, x, y, rr;
   const double maxr = 0.996194698091745532295;

   wnd->getGLDrawSize(viewport[2], viewport[3]);
   r = sqrt(sqr(viewport[2])+sqr(viewport[3]))*M_SQRT1_2;

   x = double(newx-viewport[0]-viewport[2]/2) / r;
   y = double(-newy+viewport[1]+viewport[3]/2) / r;

   rr = sqrt(x*x+y*y);
   if (rr > maxr)
   {
      x *= maxr/rr, y *= maxr/rr, rr = maxr;
   }

   new_sph_u = 2.0 * acos(rr) - M_PI_2;
   new_sph_t = atan2(y, x);
}

void LeftButtonDown (EventInfo *event)
{
   locscene -> spinning = 0;
   RemoveIdleFunc(MainLoop);

   oldx = event->mouse_x;
   oldy = event->mouse_y;

   ComputeSphereAngles(oldx, oldy, sph_u, sph_t);

   srot.identity();
   srot.mult(locscene->cam.RotMatrix());
   srot.mult(locscene->rotmat);

   startx = oldx;
   starty = oldy;
}

void LeftButtonLoc (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;
   int sendexpose = 1;

   if (event->keymod & KMOD_CTRL)
   {
      if (event->keymod & KMOD_SHIFT)
      {
         locscene->PreRotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
      }
      else
      {
         double new_sph_u, new_sph_t;

         ComputeSphereAngles(newx, newy, new_sph_u, new_sph_t);

         gl3::GlMatrix newrot;
         newrot.identity();

         double scoord[3], ncoord[3], inprod, cross[3];
         scoord[0] = scoord[1] = cos(sph_u);     scoord[2] = sin(sph_u);
         ncoord[0] = ncoord[1] = cos(new_sph_u); ncoord[2] = sin(new_sph_u);
         scoord[0] *= cos(sph_t);     scoord[1] *= sin(sph_t);
         ncoord[0] *= cos(new_sph_t); ncoord[1] *= sin(new_sph_t);
         inprod = InnerProd(scoord, ncoord);
         CrossProd(scoord, ncoord, cross);

         newrot.mult(locscene->cam.TransposeRotMatrix());
         newrot.rotate(acos(inprod)*(180.0/M_PI), cross[0], cross[1], cross[2]);
         newrot.mult(srot.mtx);
         locscene->rotmat = newrot.mtx;
      }
   }
   else if (event->keymod & KMOD_ALT)
   {
      locscene->Rotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
   }
   else if (event->keymod & KMOD_SHIFT)
   {
      locscene->Rotate(double(newx-oldx)/2, double(newy-oldy)/2);
   }
   else
   {
      locscene->Rotate(double(newy-oldy)/2, 1.0, 0.0, 0.0);
      locscene->PreRotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
   }

   oldx = newx;
   oldy = newy;

   if (sendexpose)
   {
      SendExposeEvent();
   }
}

void LeftButtonUp (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;

   xang = (newx-startx)/5.0;
   yang = (newy-starty)/5.0;

   if ( (event->keymod & KMOD_SHIFT) && (xang != 0.0 || yang != 0.0) )
   {
      locscene -> spinning = 1;
      AddIdleFunc(MainLoop);
      if (xang > 20) { xang = 20; } if (xang < -20) { xang = -20; }
      if (yang > 20) { yang = 20; } if (yang < -20) { yang = -20; }

      if (event->keymod & KMOD_CTRL)
      {
         constrained_spinning = 1;
      }
      else
      {
         constrained_spinning = 0;
      }
   }
}

void MiddleButtonDown (EventInfo *event)
{
   startx = oldx = event->mouse_x;
   starty = oldy = event->mouse_y;
}

void MiddleButtonLoc (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;

   if ( !( event->keymod & KMOD_CTRL ) )
   {
      int w, h;
      double TrX, TrY, scale;

      if (locscene->OrthogonalProjection)
      {
         scale = locscene->ViewScale;
      }
      else
      {
         scale = 0.4142135623730950488/tan(locscene->ViewAngle*(M_PI/360));
      }
      wnd->getGLDrawSize(w, h);
      if (w < h)
      {
         scale *= w;
      }
      else
      {
         scale *= h;
      }
      TrX = 2.0*double(oldx-newx)/scale;
      TrY = 2.0*double(newy-oldy)/scale;
      locscene->ViewCenterX += TrX;
      locscene->ViewCenterY += TrY;
   }
   else
   {
      // locscene->Translate((double)(newx-oldx)/200,(double)(newy-oldy)/200);

      double dx = double(newx-oldx)/400;
      double dy = double(oldy-newy)/400;

      if (event->keymod & KMOD_SHIFT)  // ctrl + shift
      {
         double sx = double(newx-startx)/400;
         double sy = double(starty-newy)/400;

         locscene->cam.TurnLeftRight(dx-sx);
         locscene->cam.TurnUpDown(sy-dy);

         locscene->cam.TurnUpDown(-sy);
         locscene->cam.TurnLeftRight(sx);
      }
      else if (event->keymod & KMOD_ALT) // ctrl + alt
      {
         locscene->cam.MoveForwardBackward(dy);
         locscene->cam.TiltLeftRight(-dx);
      }
      else // ctrl
      {
         locscene->cam.MoveLeftRight(dx);
         locscene->cam.MoveUpDown(-dy);
      }
   }

   SendExposeEvent();

   oldx = newx;
   oldy = newy;
}

void MiddleButtonUp (EventInfo*)
{}

void RightButtonDown (EventInfo *event)
{
   startx = oldx = event->mouse_x;
   starty = oldy = event->mouse_y;
}

void RightButtonLoc (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;

   if (event->keymod & KMOD_SHIFT)
   {
      // glLoadIdentity();
      // GLfloat light[] = {newx,-newy, sqrt((float)(newx*newx+newy*newy)), 0.0 };
      newx -= startx;
      newy -= starty;
      double l, x, y, z;
      x =  (double)newx / 300;
      y = -(double)newy / 300;
      l = sqrt (x*x + y*y);
      if (l <= 1.)
      {
         z = sqrt (1. - l*l);
      }
      else if (l < 2.)
      {
         x *= (2./l-1);
         y *= (2./l-1);
         l = 2. - l;
         z = -sqrt (1. - l*l);
      }
      else
      {
         x = 0.; y = 0.; z = -1.;
      }
      cout << "(x,y,z) = (" << x << ',' << y << ',' << z << ')' << endl;
      locscene->SetLight0CustomPos({(float)x, (float)y, (float)z, 0.f});
   }
   else if ( !( event->keymod & KMOD_CTRL ) )
   {
      locscene -> Zoom (exp ( double (oldy-newy) / 100 ));
   }
   else
   {
      locscene -> Scale ( exp ( double (oldy-newy) / 50 ) );
   }

   SendExposeEvent();

   oldx = newx;
   oldy = newy;
}

void RightButtonUp (EventInfo*)
{}

void TouchPinch(SDL_MultiGestureEvent & e)
{
   // Scale or Zoom?
   locscene->Zoom(exp(e.dDist*10));
   SendExposeEvent();
}

#if defined(GLVIS_USE_LIBTIFF)
const char *glvis_screenshot_ext = ".tif";
#elif defined(GLVIS_USE_LIBPNG)
const char *glvis_screenshot_ext = ".png";
#else
const char *glvis_screenshot_ext = ".bmp";
#endif

// https://wiki.libsdl.org/SDL_CreateRGBSurfaceFrom
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
Uint32 rmask = 0xff000000;
Uint32 gmask = 0x00ff0000;
Uint32 bmask = 0x0000ff00;
Uint32 amask = 0x000000ff;
#else // little endian, like x86
Uint32 rmask = 0x000000ff;
Uint32 gmask = 0x0000ff00;
Uint32 bmask = 0x00ff0000;
Uint32 amask = 0xff000000;
#endif

// https://halfgeek.org/wiki/Vertically_invert_a_surface_in_SDL
#define SDL_LOCKIFMUST(s) (SDL_MUSTLOCK(s) ? SDL_LockSurface(s) : 0)
#define SDL_UNLOCKIFMUST(s) { if(SDL_MUSTLOCK(s)) SDL_UnlockSurface(s); }

int InvertSurfaceVertical(SDL_Surface *surface)
{
   Uint8 *t, *a, *b, *last;
   Uint16 pitch;

   if ( SDL_LOCKIFMUST(surface) < 0 )
   {
      return -2;
   }

   /* do nothing unless at least two lines */
   if (surface->h < 2)
   {
      SDL_UNLOCKIFMUST(surface);
      return 0;
   }

   /* get a place to store a line */
   pitch = surface->pitch;
   t = (Uint8*)malloc(pitch);

   if (t == NULL)
   {
      SDL_UNLOCKIFMUST(surface);
      return -2;
   }

   /* get first line; it's about to be trampled */
   memcpy(t,surface->pixels,pitch);

   /* now, shuffle the rest so it's almost correct */
   a = (Uint8*)surface->pixels;
   last = a + pitch * (surface->h - 1);
   b = last;

   while (a < b)
   {
      memcpy(a,b,pitch);
      a += pitch;
      memcpy(b,a,pitch);
      b -= pitch;
   }

   /* in this shuffled state, the bottom slice is too far down */
   memmove( b, b+pitch, last-b );

   /* now we can put back that first row--in the last place */
   memcpy(last,t,pitch);

   /* everything is in the right place; close up. */
   free(t);
   SDL_UNLOCKIFMUST(surface);

   return 0;
}

#ifdef GLVIS_USE_LIBPNG
int SaveAsPNG(const char *fname, int w, int h, bool is_hidpi, bool with_alpha,
              std::function<void(int,void*)> get_row)
{
   png_byte *pixels = new png_byte[(with_alpha ? 4 : 3)*w];
   if (!pixels)
   {
      return 1;
   }

   png_structp png_ptr =
      png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
   if (!png_ptr)
   {
      delete [] pixels;
      return 1;
   }
   png_infop info_ptr = png_create_info_struct(png_ptr);
   if (!info_ptr)
   {
      png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
      delete [] pixels;
      return 1;
   }

   FILE *fp = fopen(fname, "wb");
   if (!fp)
   {
      png_destroy_write_struct(&png_ptr, &info_ptr);
      delete [] pixels;
      return 2;
   }

   if (setjmp(png_jmpbuf(png_ptr)))
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      delete [] pixels;
      return 3;
   }

   png_uint_32 ppi = is_hidpi ? 144 : 72; // pixels/inch
   png_uint_32 ppm = ppi/0.0254 + 0.5;    // pixels/meter
   png_set_pHYs(png_ptr, info_ptr, ppm, ppm, PNG_RESOLUTION_METER);

   png_init_io(png_ptr, fp);
   png_set_IHDR(png_ptr, info_ptr, w, h, 8,
                with_alpha ? PNG_COLOR_TYPE_RGBA : PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT);


   png_write_info(png_ptr, info_ptr);
   for (int i = 0; i < h; i++)
   {
      if (!get_row)
      {
         glReadPixels(0, h-1-i, w, 1, with_alpha ? GL_RGBA : GL_RGB,
                      GL_UNSIGNED_BYTE, pixels);
      }
      else
      {
         get_row(i, pixels);
      }
      png_write_row(png_ptr, pixels);
   }
   png_write_end(png_ptr, info_ptr);

   fclose(fp);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   delete [] pixels;

   return 0;
}
#endif // GLVIS_USE_LIBPNG

int Screenshot(const char *fname, bool convert)
{
#ifdef GLVIS_DEBUG
   cout << "Screenshot: glFinish() ... " << flush;
#endif
   glFinish();
#ifdef GLVIS_DEBUG
   cout << "done." << endl;
#endif
#ifndef __EMSCRIPTEN__
   if (wnd->isExposePending())
   {
      MFEM_WARNING("Expose pending, some events may not have been handled." << endl);
   }
   string filename = fname;
   string convert_name = fname;
   bool call_convert = false;
   if (convert)
   {
      // check the extension of 'fname' to see if convert is needed
      size_t ext_size = strlen(glvis_screenshot_ext);
      if (filename.size() < ext_size ||
          filename.compare(filename.size() - ext_size,
                           ext_size, glvis_screenshot_ext) != 0)
      {
         call_convert = true;
         filename += glvis_screenshot_ext;
      }
   }
   else // do not call convert
   {
      filename += glvis_screenshot_ext;
   }

   int w, h;
   wnd->getGLDrawSize(w, h);
   if (wnd->isSwapPending())
   {
#ifdef GLVIS_DEBUG
      cerr << "Screenshot: reading image data from back buffer..." << endl;
#endif
      glReadBuffer(GL_BACK);
   }
   else
   {
#ifdef GLVIS_DEBUG
      cerr << "Screenshot: reading image data from front buffer..." << endl;
#endif
      MFEM_WARNING("Screenshot: Reading from the front buffer is unreliable. "
                   << " Resulting screenshots may be incorrect." << endl);
      glReadBuffer(GL_FRONT);
   }
#if defined(GLVIS_USE_LIBTIFF)
   // Save a TIFF image. This requires the libtiff library, see www.libtiff.org
   TIFF* image;

   // MyExpose(w,h);

   unsigned char *pixels = new unsigned char[3*w];
   if (!pixels)
   {
      return 1;
   }

   image = TIFFOpen(filename.c_str(), "w");
   if (!image)
   {
      delete [] pixels;
      return 2;
   }

   TIFFSetField(image, TIFFTAG_IMAGEWIDTH, w);
   TIFFSetField(image, TIFFTAG_IMAGELENGTH, h);
   TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
   TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
   TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
   TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 3);
   TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
   TIFFSetField(image, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB);
   TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
   for (int i = 0; i < h; i++)
   {
      glReadPixels(0, h-i-1, w, 1, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      if (TIFFWriteScanline(image, pixels, i, 0) < 0)
      {
         TIFFClose(image);
         delete [] pixels;
         return 3;
      }
   }

   TIFFFlushData(image);
   TIFFClose(image);
   delete [] pixels;

#elif defined(GLVIS_USE_LIBPNG)
   // Save as png image. Requires libpng.
   int status = SaveAsPNG(filename.c_str(), w, h, wnd->isHighDpi());
   if (status != 0) { return status; }

#else
   // use SDL for screenshots

   // https://stackoverflow.com/questions/20233469/how-do-i-take-and-save-a-bmp-screenshot-in-sdl-2
   unsigned char * pixels = new unsigned char[w*h*4]; // 4 bytes for RGBA
   glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

   SDL_Surface * surf = SDL_CreateRGBSurfaceFrom(pixels, w, h, 8*4, w*4, rmask,
                                                 gmask, bmask, amask);
   if (surf == nullptr)
   {
      std::cerr << "unable to take screenshot: " << SDL_GetError() << std::endl;
   }
   else
   {
      if (InvertSurfaceVertical(surf))
      {
         std::cerr << "failed to invert surface, your screenshot may be upside down" <<
                   std::endl;
      }
      SDL_SaveBMP(surf, filename.c_str());
      SDL_FreeSurface(surf);
      // automatically convert to png if not being used
      if (!call_convert)
      {
         call_convert = true;
         convert_name += ".png";
      }
   }
   delete [] pixels;
#endif

   if (call_convert)
   {
      ostringstream cmd;
      cmd << "convert " << filename << ' ' << convert_name;
      if (system(cmd.str().c_str()))
      {
         return 1;
      }
      remove(filename.c_str());
   }
   return 0;
#else
   cout << "Screenshots not yet implemented for JS" << endl;
   return 1;
#endif
}

void KeyS()
{
   static int p = 1;

   if (locscene -> spinning)
   {
      locscene -> movie = 1 - locscene -> movie;
      if (locscene -> movie)
      {
         cout << "Recording a movie (series of snapshots)..." << endl;
      }
      else
      {
         cout << endl;
      }
      // use (ImageMagik's) convert GLVis_m* GLVis.{gif,mpg}
   }
   else
   {
      cout << "Taking snapshot number " << p << "... ";
      char fname[20];
      snprintf(fname, 20, "GLVis_s%02d", p++);
      wnd->screenshot(fname);
      cout << "done" << endl;
   }
   SendExposeEvent();
}

inline GL2PSvertex CreatePrintVtx(gl3::FeedbackVertex v)
{
   return
   {
      {v.position.x, v.position.y, v.position.z},
      {v.color.r, v.color.g, v.color.b, v.color.a}
   };
}

void PrintCaptureBuffer(gl3::CaptureBuffer& cbuf)
{
   // print lines
   for (size_t i = 0; i < cbuf.lines.size(); i += 2)
   {
      GL2PSvertex lineOut[2] =
      {
         CreatePrintVtx(cbuf.lines[i]),
         CreatePrintVtx(cbuf.lines[i+1])
      };
      gl2psAddPolyPrimitive(GL2PS_LINE, 2, lineOut, 0, 0.f, 0.f,
                            0xFFFF, 1, 0.2, 0, 0, 0);
   }
   // print triangles
   for (size_t i = 0; i < cbuf.triangles.size(); i += 3)
   {
      GL2PSvertex triOut[3] =
      {
         CreatePrintVtx(cbuf.triangles[i]),
         CreatePrintVtx(cbuf.triangles[i+1]),
         CreatePrintVtx(cbuf.triangles[i+2])
      };
      gl2psAddPolyPrimitive(GL2PS_TRIANGLE, 3, triOut, 0, 0.f, 0.f,
                            0xFFFF, 1, 1, 0, 0, 0);
   }
   // print text
   for (const auto &entry : cbuf.text)
   {
      GL2PSvertex rpos = CreatePrintVtx({entry.offset, entry.color});
      gl2psForceRasterPos(&rpos);
      gl2psText(entry.text.c_str(), "Times", 12);
   }
}

void KeyCtrlP()
{
#ifdef __EMSCRIPTEN__
   cerr << "Printing in WebGL is not supported at this time." << endl;
#else
   cout << "Printing the figure to GLVis.pdf... " << flush;
   locscene -> print = 1;
   FILE * fp;
   fp = fopen("GLVis.pdf", "wb");
   GLint viewport[4] = { 0, 0, 0, 0 };
   wnd->getGLDrawSize(viewport[2], viewport[3]);
   {
      gl3::SceneInfo wnd_scn = locscene->GetSceneObjs();
      for (auto to_buf : wnd_scn.needs_buffering)
      {
         wnd->getRenderer().buffer(to_buf);
      }
      gl2psBeginPage ( "GLVis.pdf", "GLVis", viewport,
                       GL2PS_PDF, // or GL2PS_SVG, or GL2PS_EPS
                       GL2PS_BSP_SORT,
                       GL2PS_SIMPLE_LINE_OFFSET |
                       // GL2PS_NO_PS3_SHADING |
                       // GL2PS_NO_BLENDING |
                       // GL2PS_OCCLUSION_CULL |
                       // GL2PS_BEST_ROOT |
                       GL2PS_SILENT |
                       // GL2PS_DRAW_BACKGROUND |
                       GL2PS_NO_BLENDING |
                       GL2PS_NO_OPENGL_CONTEXT,
                       GL_RGBA, 0, NULL, 16, 16, 16, 0, fp, "a" );
      gl3::CaptureBuffer cbuf = wnd->getRenderer().capture(wnd_scn.queue);
      PrintCaptureBuffer(cbuf);
      gl2psEndPage();
   }
   locscene -> print = 0;
   fclose(fp);
   cout << "done" << endl;
   wnd->signalExpose();
#endif
}

void KeyQPressed()
{
   wnd->signalQuit();
   visualize = 0;
}

void ToggleThreads()
{
   static const char *state[] = { "running", "stopped" };
   if (visualize > 0 && visualize < 3)
   {
      visualize = 3 - visualize; //  1 <-> 2
      cout << "Communication thread(s): " << state[visualize-1] << endl;
   }
}

void ThreadsPauseFunc(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      glvis_command->ToggleAutopause();
   }
   else
   {
      ToggleThreads();
   }
}

void ThreadsStop()
{
   if (visualize == 1)
   {
      ToggleThreads();
   }
}

void ThreadsRun()
{
   if (visualize == 2)
   {
      ToggleThreads();
   }
}

void CheckSpin()
{
   if (fabs(xang) < 1.e-2)
   {
      xang = 0.;
   }
   if (xang != 0. || yang != 0.)
   {
      locscene->spinning = 1;
      AddIdleFunc(MainLoop);
   }
   else
   {
      locscene->spinning = 0;
      RemoveIdleFunc(MainLoop);
   }
   cout << "Spin angle: " << xang << " degrees / frame" << endl;
}

const double xang_step = 0.2; // angle in degrees

void Key0Pressed()
{
   if (!locscene -> spinning)
   {
      xang = 0;
   }
   xang -= xang_step;
   CheckSpin();
}

void KeyDeletePressed()
{
   if (locscene -> spinning)
   {
      xang = yang = 0.;
      locscene -> spinning = 0;
      RemoveIdleFunc(MainLoop);
      constrained_spinning = 1;
   }
   else
   {
      xang = xang_step;
      locscene -> spinning = 1;
      AddIdleFunc(MainLoop);
      constrained_spinning = 1;
   }
}

void KeyEnterPressed()
{
   if (!locscene -> spinning)
   {
      xang = 0;
   }
   xang += xang_step;
   CheckSpin();
}

void Key7Pressed()
{
   locscene->PreRotate(1.0, 0.0, -1.0, 0.0);
   SendExposeEvent();
}

void Key8Pressed()
{
   locscene->Rotate(0.0, -1.0);
   SendExposeEvent();
}

void Key9Pressed()
{
   locscene->PreRotate(-1.0, 1.0, 0.0, 0.0);
   SendExposeEvent();
}

void Key4Pressed()
{
   locscene->PreRotate(-1.0, 0.0, 0.0, 1.0);
   SendExposeEvent();
}

void Key5Pressed()
{
   if (locscene->view == 2)
   {
      locscene->CenterObject2D();
   }
   else
   {
      locscene->CenterObject();
   }
   SendExposeEvent();
}

void Key6Pressed()
{
   locscene->PreRotate(1.0, 0.0, 0.0, 1.0);
   SendExposeEvent();
}

void Key1Pressed()
{
   locscene->PreRotate(1.0, 1.0, 0.0, 0.0);
   SendExposeEvent();
}

void Key2Pressed()
{
   locscene->Rotate(1.0, 1.0, 0.0, 0.0);
   SendExposeEvent();
}

void Key3Pressed()
{
   locscene->PreRotate(1.0, 0.0, 1.0, 0.0);
   SendExposeEvent();
}

void ShiftView(double dx, double dy)
{
   double scale;
   if (locscene->OrthogonalProjection)
   {
      scale = locscene->ViewScale;
   }
   else
   {
      scale = 0.4142135623730950488/tan(locscene->ViewAngle*(M_PI/360));
   }
   locscene->ViewCenterX += dx/scale;
   locscene->ViewCenterY += dy/scale;
}

void KeyLeftPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(0.05, 0.);
   }
   else
   {
      locscene->Rotate(-5, 0.0, 1.0, 0.0);
   }
   SendExposeEvent();
}

void KeyRightPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(-0.05, 0.);
   }
   else
   {
      locscene->Rotate(5, 0.0, 1.0, 0.0);
   }
   SendExposeEvent();
}

void KeyUpPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(0., -0.05);
   }
   else
   {
      locscene->Rotate(-5, 1.0, 0.0, 0.0);
   }
   SendExposeEvent();
}

void KeyDownPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(0., 0.05);
   }
   else
   {
      locscene->Rotate(5, 1.0, 0.0, 0.0);
   }
   SendExposeEvent();
}

void KeyJPressed()
{
   locscene->OrthogonalProjection = !(locscene->OrthogonalProjection);
   SendExposeEvent();
}

void KeyMinusPressed()
{
   locscene -> Scale(1., 1., 1./1.1);
   SendExposeEvent();
}

void KeyPlusPressed()
{
   locscene -> Scale(1., 1., 1.1);
   SendExposeEvent();
}

void ZoomIn()
{
   locscene->Zoom(exp (0.05));
   SendExposeEvent();
}

void ZoomOut()
{
   locscene->Zoom(exp (-0.05));
   SendExposeEvent();
}

void ScaleUp()
{
   locscene->Scale(1.025);
   SendExposeEvent();
}

void ScaleDown()
{
   locscene->Scale(1.0/1.025);
   SendExposeEvent();
}

void LookAt()
{
   cout << "ViewCenter = (" << locscene->ViewCenterX << ','
        << locscene->ViewCenterY << ")\nNew x = " << flush;
   cin >> locscene->ViewCenterX;
   cout << "New y = " << flush;
   cin >> locscene->ViewCenterY;
   SendExposeEvent();
}

const double window_scale_factor = 1.1;

void ShrinkWindow()
{
   int w, h;
   wnd->getWindowSize(w, h);
   w = (int)ceil(w / window_scale_factor);
   h = (int)ceil(h / window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   ResizeWindow(w, h);
}

void EnlargeWindow()
{
   int w, h;
   wnd->getWindowSize(w, h);
   w = (int)ceil(w * window_scale_factor);
   h = (int)ceil(h * window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   ResizeWindow(w, h);
}

void MoveResizeWindow(int x, int y, int w, int h)
{
   wnd->setWindowSize(w, h);
   wnd->setWindowPos(x, y);
}

void ResizeWindow(int w, int h)
{
   wnd->setWindowSize(w, h);
}

void SetWindowTitle(const char *title)
{
   wnd->setWindowTitle(title);
}

int GetUseTexture()
{
   return locscene->palette.GetSmoothSetting();
}

void SetUseTexture(int ut)
{
   if (ut == 0)
   {
      locscene->palette.UseDiscrete();
   }
   else
   {
      locscene->palette.UseSmooth();
   }
}

int GetMultisample()
{
   return glvis_multisample;
}

void SetMultisample(int m)
{
   if (glvis_multisample > -2)
   {
      glvis_multisample = m;
   }
   else
   {
      cout << "Multisampling is disabled." << endl;
   }
}

void SetLineWidth(float width)
{
   line_w = width;
   if (wnd)
   {
      wnd->getRenderer().setLineWidth(line_w);
   }
}

void SetLineWidthMS(float width_ms)
{
   line_w_aa = width_ms;
   if (wnd)
   {
      wnd->getRenderer().setLineWidthMS(line_w_aa);
   }

}

float GetLineWidth()
{
   return line_w;
}

float GetLineWidthMS()
{
   return line_w_aa;
}


// Fontconfig patterns to use for finding a font file.
// Use the command:
//    fc-list "pattern" file style
// to find the font files that match the "pattern".
vector<string> fc_font_patterns =
{
   "Ubuntu Light:style=Regular",
   "Ubuntu:style=Regular:weight=80",
   "OpenSans:style=Regular",
   "DejaVu Sans:style=Book:width=Normal",
   "DejaVu LGC Sans:style=Book:width=Normal",
   "Bitstream Vera Sans:style=Roman",
   "FreeSans:style=Medium",
   "Ubuntu Mono:style=Regular",
   "DejaVu Sans Mono:style=Book",
   "DejaVu LGC Sans Mono:style=Book",
   "Helvetica:style=Regular",
   "Arial:style=Regular:weight=80"
};

#ifdef GLVIS_FONT_SIZE
constexpr int default_font_size = GLVIS_FONT_SIZE;
#else
constexpr int default_font_size = 12;
#endif
int font_size = default_font_size;

thread_local GlVisFont glvis_font;
std::string priority_font;

void InitFont()
{
   // This function is called after the window is created.
   GLenum alphaChannel =
      gl3::GLDevice::useLegacyTextureFmts() ? GL_ALPHA : GL_RED;
   glvis_font.setAlphaChannel(alphaChannel);
   bool try_fc_patterns = true;
   if (!priority_font.empty())
   {
      if (SetFont({priority_font}, font_size) ||
          glvis_font.LoadFont(priority_font, 0, font_size))
      {
         try_fc_patterns = false;
      }
      else
      {
         cerr << "InitFont(): Font not found: " << priority_font << endl;
      }
   }
   if (try_fc_patterns)
   {
      if (!SetFont(fc_font_patterns, font_size))
      {
         cerr <<
              "InitFont(): No fonts found matching the built-in patterns.\n"
              "Use the -fn option or edit 'fc_font_patterns' in lib/aux_vis.cpp"
              << endl;
      }
   }
   wnd->getRenderer().setFontTexture(glvis_font.getFontTex());
}

GlVisFont * GetFont()
{
   return &glvis_font;
}

bool SetFont(const vector<std::string>& font_patterns, int height)
{
#ifdef __EMSCRIPTEN__
   return glvis_font.LoadFont("OpenSans.ttf", 0, height);
#else
   if (!FcInit())
   {
      return false;
   }

   FcObjectSet * os = FcObjectSetBuild(FC_FAMILY, FC_STYLE, FC_FILE,
                                       FC_SCALABLE, FC_INDEX, FC_WEIGHT,
                                       nullptr);

   for (const string& pattern : font_patterns)
   {
      string patternScale = pattern + ":scalable=True";
      FcPattern * pat = FcNameParse((FcChar8*)patternScale.c_str());
      if (!pat)
      {
         continue;
      }

      FcFontSet * fs = FcFontList(0, pat, os);
      if (!fs)
      {
         FcPatternDestroy(pat);
         continue;
      }
#ifdef GLVIS_DEBUG
      if (fs->nfont >= 1)
      {
         cout << "Font pattern '" << pattern << "' matched fonts:\n";
      }
      else
      {
         cout << "Font pattern '" << pattern << "' matched no fonts.\n";
      }
#endif
      std::string font_file = "";
      std::string font_name = "";
      int font_index = 0;
      for (int match_idx = 0; match_idx < fs->nfont; match_idx++)
      {
         FcChar8 * s;
         FcBool scalable;
         int curr_font_idx;
         FcPatternGetBool(fs->fonts[match_idx], FC_SCALABLE, 0, &scalable);
         FcPatternGetInteger(fs->fonts[match_idx], FC_INDEX, 0, &curr_font_idx);
         FcResult res = FcPatternGetString(fs->fonts[match_idx], FC_FILE, 0, &s);
         FcChar8 * fnt = FcNameUnparse(fs->fonts[match_idx]);
#ifdef GLVIS_DEBUG
         cout << " - " << fnt << endl;
#endif
         if (res == FcResultMatch && s && font_file == std::string(""))
         {
            font_file = (char*) s;
            font_name = (char*) fnt;
            font_index = curr_font_idx;
         }
         free(fnt);
      }

      FcFontSetDestroy(fs);
      if (font_file != std::string(""))
      {
         if (glvis_font.LoadFont(font_file, font_index, height))
         {
            break;
         }
      }
   }

   if (os)
   {
      FcObjectSetDestroy(os);
   }

   FcFini();

   return glvis_font.isFontLoaded();
#endif
}

void SetFont(const std::string& fn)
{
   priority_font = fn;
   size_t pos = priority_font.rfind('-');
   if (pos != string::npos)
   {
      font_size = std::stoi(priority_font.substr(pos + 1));
      priority_font.erase(pos);
   }
#ifdef GLVIS_DEBUG
   cout << "SetFont: name = '"
        << (priority_font.empty() ? "(default)" : priority_font)
        << "', height = " << font_size << endl;
#endif
}

function<string(double)> NumberFormatter(int precision, char format,
                                         bool showsign)
{
   return [precision, format, showsign](double x) -> string
   {
      ostringstream oss;
      switch (format)
      {
         case 'f':
            oss << fixed;
            break;
         case 's':
            oss << scientific;
            break;
         case 'd':
            oss << defaultfloat;
            break;
         default:
            MFEM_WARNING("Unknown formatting type. Using default. "
                         << "Valid options include: ['f', 's', 'd']" << endl);
            oss << defaultfloat;
            break;
      };
      if (showsign)
      {
         oss << showpos;
      }
      oss << setprecision(precision) << x;
      return oss.str();
   };
}

function<string(double)> NumberFormatter(string formatting)
{
   if (!isValidNumberFormatting(formatting))
   {
      MFEM_WARNING("Invalid formatting string. Using default. " << endl);
      return NumberFormatter();
   }
   else
   {
      return [formatting](double x) -> string
      {
         char buf[64];
         snprintf(buf, sizeof(buf), formatting.c_str(), x);
         return string(buf);
      };
   }
}

bool isValidNumberFormatting(const string& formatting)
{
   regex rgx = regex(R"(%[\-+#0\s]?[0-9]{0,3}\.?[0-9]{0,3}[FfEeGg])");
   return regex_match(formatting, rgx);
}
