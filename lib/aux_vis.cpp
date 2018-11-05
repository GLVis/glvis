// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <ctime>

#include "mfem.hpp"
using namespace mfem;
#include "sdl.hpp"
#include "palettes.hpp"
#include "visual.hpp"
#include "gl2ps.h"
#include "gl3print.hpp"

#if defined(GLVIS_USE_LIBTIFF)
#include "tiffio.h"
#elif defined(GLVIS_USE_LIBPNG)
#include <png.h>
#endif

#include "font.hpp"
#ifndef __EMSCRIPTEN__
#include <fontconfig/fontconfig.h>
#endif

int visualize = 0;
VisualizationScene * locscene;

float MatAlpha = 1.0;
float MatAlphaCenter = 0.5;

#ifdef GLVIS_MULTISAMPLE
static int glvis_multisample = GLVIS_MULTISAMPLE;
#else
static int glvis_multisample = -1;
#endif

//TODO: anything but this
SdlWindow * wnd = nullptr;
GlState * state = nullptr;

SdlWindow * GetAppWindow()
{
    return wnd;
}

GlState * GetGlState()
{
    return state;
}

VisualizationScene * GetVisualizationScene()
{
    return locscene;
}

void MyExpose(GLsizei w, GLsizei h);
void MyExpose();

int InitVisualization (const char name[], int x, int y, int w, int h)
{

#ifdef GLVIS_DEBUG
   cout << "OpenGL Visualization" << endl;
#endif
   if (!wnd) {
       wnd = new SdlWindow();
       if (!wnd->createWindow(name, w, h)) {
           return 1;
       }

       state = new GlState();
       if (!state->compileShaders()) {
           return 1;
       }
       paletteInit();
   } else {
       wnd->clearEvents();
   }

   paletteInit();
   InitFont();

#ifdef GLVIS_DEBUG
   cout << "Window should be up" << endl;
#endif

   SetUseTexture(0);

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

   // auxKeyFunc (AUX_p, KeyCtrlP); // handled in vsdata.cpp
   wnd->setOnKeyDown (SDLK_s, KeyS);
   wnd->setOnKeyDown ('S', KeyS);

   wnd->setOnKeyDown (SDLK_q, KeyQPressed);
   //wnd->setOnKeyDown (SDLK_Q, KeyQPressed);

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
   //wnd->setOnKeyDown (AUX_J, KeyJPressed);

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
       delete locscene;
#endif
   locscene = nullptr;

   return 0;
}

void SendKeySequence(const char *seq) {
    for (const char* key = seq; *key != '\0'; key++) {
        if (*key == '~') {
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
               case '.': // Keypad ./Del
                  wnd->signalKeyDown(SDLK_PERIOD);
                  break;
               case 'E': // Keypad Enter
                  wnd->signalKeyDown(SDLK_RETURN);
                  break;
            }
            continue;
        } else if (*key == '*') {
            wnd->signalKeyDown(SDLK_KP_MULTIPLY);
        } else if (*key == '/') {
            wnd->signalKeyDown(SDLK_KP_DIVIDE);
        } else {
            if (*key == '('
                || *key == ')'
                || *key == '!'
                || isupper(*key)) {
                wnd->signalKeyDown(*key, KMOD_LSHIFT);
            }
            wnd->signalKeyDown(*key, KMOD_LSHIFT);
        }
    }
}


static bool disableSendExposeEvent = false;

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

   //auxMainLoop(NULL);
#ifndef __EMSCRIPTEN__
   wnd->mainLoop();
#endif
   InitIdleFuncs();
}

void KillVisualization()
{
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
   state->setViewport(w, h);

   state->projection.identity();

   double ViewCenterX = locscene->ViewCenterX;
   double ViewCenterY = locscene->ViewCenterY;

   if (locscene->OrthogonalProjection)
   {
      double scale = locscene->ViewScale;
      if (w <= h)
      {
         state->projection.ortho(-1.0, 1.0, -double(h)/w, double(h)/w, -10, 10);
      }
      else
      {
         state->projection.ortho(-double(w)/h, double(w)/h, -1, 1, -10, 10);
      }
      state->projection.scale(scale, scale, 1.0);
   }
   else
   {
      double ViewAngle = locscene->ViewAngle;

      if (w < h)
         ViewAngle = (360.0/M_PI)*atan( tan( ViewAngle*(M_PI/360.0) ) *
                                        double(h)/w );

      state->projection.perspective(ViewAngle, double(w)/h, 0.1, 5.0);
   }
   state->projection.translate(-ViewCenterX, -ViewCenterY, 0.0);
}

void MyExpose(GLsizei w, GLsizei h)
{
   MyReshape (w, h);
   locscene -> Draw();
}

void MyExpose()
{
   int w, h;
   wnd->getWindowSize(w, h);
   MyExpose(w, h);
}


Array<void (*)()> IdleFuncs;
int LastIdleFunc;

void InitIdleFuncs()
{
   IdleFuncs.SetSize(0);
   LastIdleFunc = 0;
   wnd->setOnIdle(NULL);
}

void MainIdleFunc()
{
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
   if (IdleFuncs.Size() == 0)
   {
      wnd->setOnIdle(NULL);
   }
}


double xang = 0., yang = 0.;
GlMatrix srot;
double sph_t, sph_u;
static GLint oldx, oldy, startx, starty;

int constrained_spinning = 0;


void MainLoop()
{
   static int p = 1;
   struct timespec req;
   if (locscene->spinning)
   {
      if (!constrained_spinning)
      {
         locscene->Rotate(xang, yang);
         locscene->Draw();
      }
      else
      {
         locscene->PreRotate(xang, 0.0, 0.0, 1.0);
         locscene->Draw();
      }
      req.tv_sec  = 0;
      req.tv_nsec = 10000000;
      nanosleep (&req, NULL);  // sleep for 0.01 seconds
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

   state->getViewport(viewport);
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

         state->modelView.identity();

         double scoord[3], ncoord[3], inprod, cross[3];
         scoord[0] = scoord[1] = cos(sph_u);     scoord[2] = sin(sph_u);
         ncoord[0] = ncoord[1] = cos(new_sph_u); ncoord[2] = sin(new_sph_u);
         scoord[0] *= cos(sph_t);     scoord[1] *= sin(sph_t);
         ncoord[0] *= cos(new_sph_t); ncoord[1] *= sin(new_sph_t);
         inprod = InnerProd(scoord, ncoord);
         CrossProd(scoord, ncoord, cross);

         state->modelView.mult(locscene->cam.TransposeRotMatrix());
         state->modelView.rotate(acos(inprod)*(180.0/M_PI), cross[0], cross[1], cross[2]);
         state->modelView.mult(srot.mtx);
         locscene->rotmat = state->modelView.mtx;
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
      GLint vp[4];
      double TrX, TrY, scale;

      if (locscene->OrthogonalProjection)
      {
         scale = locscene->ViewScale;
      }
      else
      {
         scale = 0.4142135623730950488/tan(locscene->ViewAngle*(M_PI/360));
      }
      state->getViewport(vp);
      if (vp[2] < vp[3])
      {
         scale *= vp[2];
      }
      else
      {
         scale *= vp[3];
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

void MiddleButtonUp (EventInfo *event)
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
      //glLoadIdentity();
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
      GLfloat light[] = { float(x), float(y), float(z), 0.0f };
      state->setLightPosition(0, light);
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

void RightButtonUp (EventInfo *event)
{}

#if defined(GLVIS_USE_LIBTIFF)
const char *glvis_screenshot_ext = ".tif";
#elif defined(GLVIS_USE_LIBPNG)
const char *glvis_screenshot_ext = ".png";
#else
const char *glvis_screenshot_ext = ".xwd";
#endif

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
   string filename = fname;
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
   wnd->getWindowSize(w, h);
   glReadBuffer(GL_FRONT);
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

   png_byte *pixels = new png_byte[3*w];
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

   FILE *fp = fopen(filename.c_str(), "wb");
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

   png_init_io(png_ptr, fp);
   png_set_IHDR(png_ptr, info_ptr, w, h, 8, PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT);

   png_write_info(png_ptr, info_ptr);
   for (int i = 0; i < h; i++)
   {
      glReadPixels(0, h-1-i, w, 1, GL_RGB, GL_UNSIGNED_BYTE, pixels);
      png_write_row(png_ptr, pixels);
   }
   png_write_end(png_ptr, info_ptr);

   fclose(fp);
   png_destroy_write_struct(&png_ptr, &info_ptr);
   delete [] pixels;

#elif defined(GLVIS_X11)
   // Use the external X Window Dump (xwd) tool.
   // Note that xwd does not work on OS X!
   ostringstream cmd;
   cmd << "xwd -silent -out " << filename << " -nobdrs -id " << wnd->getXWindow();
   if (system(cmd.str().c_str()))
   {
      return 1;
   }
   // View with xwud -in GLVis_s*.xwd, or use convert GLVis_s*.xwd
   // GLVis_s*.{jpg,gif}
#else
   cerr << "No method for taking screenshots detected!" << endl;
   return 0;
#endif

   if (call_convert)
   {
      ostringstream cmd;
      cmd << "convert " << filename << ' ' << fname;
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
}

void KeyCtrlP()
{
#ifndef __EMSCRIPTEN__
   if (!GetGlState()->renderToFeedback()) {
       cout << "Unable to initialize printing capture pipeline." << endl;
       return;
   }
   cout << "Printing the figure to GLVis.pdf... " << flush;
   locscene -> print = 1;
   FILE * fp;
   fp = fopen("GLVis.pdf", "wb");
   GLint viewport[4] = { 0, 0, 0, 0 };
   GetGlState()->getViewport(viewport);
   {
       gl3::GL2PSFeedbackHook fb_capture;
       gl3::GlDrawable::setDrawHook(&fb_capture);
       gl2psBeginPage ( "GLVis.pdf", "GLVis", viewport,
                        GL2PS_PDF, // or GL2PS_SVG, or GL2PS_EPS
                        GL2PS_BSP_SORT,
                        GL2PS_SIMPLE_LINE_OFFSET |
                        // GL2PS_NO_PS3_SHADING |
                        // GL2PS_NO_BLENDING |
                        // GL2PS_OCCLUSION_CULL |
                        // GL2PS_BEST_ROOT |
                        GL2PS_SILENT |
                        //GL2PS_DRAW_BACKGROUND |
                        GL2PS_NO_BLENDING |
                        GL2PS_NO_OPENGL_CONTEXT,
                        GL_RGBA, 0, NULL, 16, 16, 16, 0, fp, "a" );
       locscene -> Draw();
       int state = gl2psEndPage();
       gl3::GlDrawable::setDrawHook(nullptr);
       GetGlState()->renderToDefault();
   }
   locscene -> print = 0;
   fclose(fp);

   cout << "done" << endl;
   wnd->signalExpose();
#else
   cout << "Printing disabled" << endl;
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
}

void Key0Pressed()
{
   if (!locscene -> spinning)
   {
      xang = 0;
   }
   xang -= 0.2;
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
      xang = 0.2;
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
   xang += 0.2;
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
   GLint viewport[4];

   state->getViewport(viewport);
   int w = viewport[2], h = viewport[3];
   w = (int)ceil(w / window_scale_factor);
   h = (int)ceil(h / window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   ResizeWindow(w, h);
}

void EnlargeWindow()
{
   GLint viewport[4];

   state->getViewport(viewport);
   int w = viewport[2], h = viewport[3];
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


// Draw a cone of radius 1 with base in the x-y plane and center at (0,0,2)
void Cone()
{
   const int n = 8;
   const double step = 2*M_PI/n;
   const double nz = (1.0/4.0);
   double point = step;
   int i;

   glBegin(GL_TRIANGLE_FAN);
   glNormal3d(0.0, 0.0, 1.0);
   glVertex3d(0, 0, 0);
   glNormal3d(1.0, 0.0, nz);
   glVertex3d(1, 0, -4);
   for (i = 1; i < n; i++)
   {
      glNormal3d(cos(point), sin(point), nz);
      glVertex3d(cos(point), sin(point), -4);
      point += step;
   }
   glNormal3d(1.0, 0.0, nz);
   glVertex3d(1, 0, -4);
   glEnd();
}

int MySetColorLogscale = 0;
extern int RepeatPaletteTimes;
int UseTexture         = 0;

float MySetColor (double val, double min, double max, float (&rgba)[4])
{
   // static double eps = 1e-24;
   static const double eps = 0.0;
   if (MySetColorLogscale)
   {
      if (val < min)
      {
         val = min;
      }
      if (val > max)
      {
         val = max;
      }
      return MySetColor (log(fabs(val/(min+eps))) / (log(fabs(max/(min+eps)))+eps), rgba);
   }
   else
   {
      return MySetColor ((val-min)/(max-min), rgba);
   }
}

void MySetColor (gl3::GlBuilder& builder, double val, double min, double max) {
  static const double eps = 0.0;
  if (MySetColorLogscale)
  {
     if (val < min)
     {
        val = min;
     }
     if (val > max)
     {
        val = max;
     }
     return MySetColor (builder, log(fabs(val/(min+eps))) / (log(fabs(max/(min+eps)))+eps));
  }
  else
  {
     return MySetColor (builder, (val-min)/(max-min));
  }
}

float MySetColor (double val, float (&rgba)[4])
{
   int i;
   double t, *pal;

   int palSize = paletteGetSize();

   if (val < 0.0) { val = 0.0; }
   if (val > 1.0) { val = 1.0; }

   double malpha = MatAlpha;
   if (malpha < 1.0)
   {
      if (MatAlphaCenter > 1.0)
      {
         malpha *= exp(-(MatAlphaCenter)*fabs(val-1.0));
      }
      else if (MatAlphaCenter < 0.0)
      {
         malpha *= exp((MatAlphaCenter-1.0)*fabs(val-0.0));
      }
      else
      {
         malpha *= exp(-fabs(val-MatAlphaCenter));
      }
   }
   
   if (UseTexture)
   {
      //return 1-alpha since default attrib value is 0.0
      return MatAlpha < 1.0 ? 1.0 - malpha : 0.0;
   }

   val *= 0.999999999 * ( palSize - 1 ) * abs(RepeatPaletteTimes);
   i = (int) floor( val );
   t = val - i;

   if (((i / (palSize-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
       ((i / (palSize-1)) % 2 == 1 && RepeatPaletteTimes < 0))
   {
      pal = paletteGet() + 3 * ( i % (palSize-1) );
   }
   else
   {
      pal = paletteGet() + 3 * ( (palSize-2) -
                                i % (palSize-1) );
      t = 1.0 - t;
   }

   rgba[0] = (1.0 - t) * pal[0] + t * pal[3];
   rgba[1] = (1.0 - t) * pal[1] + t * pal[4];
   rgba[2] = (1.0 - t) * pal[2] + t * pal[5];
   rgba[3] = MatAlpha < 1.0 ? malpha : 1.0;
   return 0.0;
}

void MySetColor (gl3::GlBuilder& builder, double val)
{
   int i;
   double t, *pal;

   int palSize = paletteGetSize();

   if (val < 0.0) { val = 0.0; }
   if (val > 1.0) { val = 1.0; }

   double malpha = MatAlpha;
   if (malpha < 1.0)
   {
      if (MatAlphaCenter > 1.0)
      {
         malpha *= exp(-(MatAlphaCenter)*fabs(val-1.0));
      }
      else if (MatAlphaCenter < 0.0)
      {
         malpha *= exp((MatAlphaCenter-1.0)*fabs(val-0.0));
      }
      else
      {
         malpha *= exp(-fabs(val-MatAlphaCenter));
      }
   }

   if (UseTexture)
   {
      builder.glTexCoord2f(val, MatAlpha < 1.0 ? malpha : 1.0);
      return;
   }

   val *= 0.999999999 * ( palSize - 1 ) * abs(RepeatPaletteTimes);
   i = (int) floor( val );
   t = val - i;

   if (((i / (palSize-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
       ((i / (palSize-1)) % 2 == 1 && RepeatPaletteTimes < 0))
   {
      pal = paletteGet() + 3 * ( i % (palSize-1) );
   }
   else
   {
      pal = paletteGet() + 3 * ( (palSize-2) -
                                i % (palSize-1) );
      t = 1.0 - t;
   }
   float rgba[4] = {
       (float)((1.0 - t) * pal[0] + t * pal[3]),
       (float)((1.0 - t) * pal[1] + t * pal[4]),
       (float)((1.0 - t) * pal[2] + t * pal[5]),
       (float)(MatAlpha < 1.0 ? malpha : 1.0)
   };
   builder.glColor4fv(rgba);
}

int GetUseTexture()
{
   return UseTexture;
}

void SetUseTexture(int ut)
{
   if (UseTexture != ut)
   {
      UseTexture = ut;
      if (UseTexture == 1)
          paletteUseDiscrete();
      else
          paletteUseSmooth();
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


// Fontconfig patterns to use for finding a font file.
// Use the command:
//    fc-list "pattern" file style
// to find the font files that match the "pattern".
const char *fc_font_patterns[] =
{
   "Ubuntu Light:style=Regular",
   "Ubuntu:style=Regular:weight=80",
   "DejaVu Sans:style=Book:width=Normal",
   "DejaVu LGC Sans:style=Book:width=Normal",
   "Bitstream Vera Sans:style=Roman",
   "FreeSans:style=Medium",
   "Ubuntu Mono:style=Regular",
   "DejaVu Sans Mono:style=Book",
   "DejaVu LGC Sans Mono:style=Book"
};

const int num_font_patterns = sizeof(fc_font_patterns)/sizeof(char *);

#ifdef __EMSCRIPTEN__
int font_size = 12;
#else
int font_size = 10;
#endif

GlVisFont glvis_font;
std::string priority_font;

void InitFont()
{
    if (priority_font == "") {
        SetFont(fc_font_patterns, num_font_patterns, font_size);
    } else {
        if (!glvis_font.LoadFont(priority_font.c_str(), font_size))
           cout << "Font not found: " << priority_font << endl;
    }
}

GlVisFont * GetFont() {
   return &glvis_font;
}

bool SetFont(const char *font_patterns[], int num_patterns, int height) {
#ifdef __EMSCRIPTEN__
    return glvis_font.LoadFont("OpenSans.ttf", height);
#else
    if (!FcInit()) {
        return false;
    }

    FcObjectSet * os = FcObjectSetBuild(FC_FAMILY, FC_STYLE, FC_FILE, nullptr);

    for (int i = 0; i < num_patterns; i++) {
        FcPattern * pat = FcNameParse((FcChar8*)font_patterns[i]);
        if (!pat) {
            continue;
        }

        FcFontSet * fs = FcFontList(0, pat, os);
        if (!fs) {
            FcPatternDestroy(pat);
            continue;
        }
#ifdef GLVIS_DEBUG
        if (fs->nfont > 1) {
            cout << "Font pattern '" << font_patterns[i]
                 << "' matched multiple fonts:\n";
        }
#endif
        std::string font_file = "";
        std::string font_name = "";
        for (int fnt_idx = 0; fnt_idx < fs->nfont; fnt_idx++) {
            FcChar8 * s;
            FcResult res = FcPatternGetString(fs->fonts[fnt_idx], FC_FILE, 0, &s);
            FcChar8 * fnt = FcNameUnparse(fs->fonts[fnt_idx]);
            cout << fnt << endl;
            if (res == FcResultMatch && s && font_file == "") {
                font_file = (char*) s;
                font_name = (char*) fnt;
            }
            free(fnt);
        }
        
        FcFontSetDestroy(fs);
        if (font_file != "") {
            if (glvis_font.LoadFont(font_file.c_str(), height)) {
#ifdef GLVIS_DEBUG
                cout << "Using font: " << font_name << endl;
#endif
                break;
            }
        }
    }

    if (os)
        FcObjectSetDestroy(os);

    FcFini();

    return glvis_font.isFontLoaded();
#endif
}

void SetFont(const char *fn)
{
   priority_font = fn;
   size_t pos = priority_font.rfind('-');
   if (pos != string::npos)
   {
      font_size = std::stoi(priority_font.substr(pos + 1));
      priority_font.erase(pos);
   }
#ifdef GLVIS_DEBUG
   cout << "SetFont: name = '" << priority_font << "', height = " << font_size
        << endl;
#endif
}
