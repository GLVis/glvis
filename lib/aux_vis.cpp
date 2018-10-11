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
#include <X11/keysym.h>

#include "mfem.hpp"
using namespace mfem;

#include "palettes.hpp"
#include "gl2ps.h"
#include "visual.hpp"

#if defined(GLVIS_USE_LIBTIFF)
#include "tiffio.h"
#elif defined(GLVIS_USE_LIBPNG)
#include <png.h>
#endif

#ifdef GLVIS_USE_FREETYPE
#include <ft2build.h>
#include FT_FREETYPE_H
#include FT_GLYPH_H
#include <fontconfig/fontconfig.h>
#endif

string fontname;
GLuint fontbase;
int visualize = 0;
VisualizationScene * locscene;

float MatAlpha = 1.0;
float MatAlphaCenter = 0.5;

#ifdef GLVIS_MULTISAMPLE
static int glvis_multisample = GLVIS_MULTISAMPLE;
#else
static int glvis_multisample = -1;
#endif

void MyExpose(GLsizei w, GLsizei h);

int InitVisualization (const char name[], int x, int y, int w, int h)
{
   static int init = 0;

   if (!init)
   {
      Init_Palettes();
      init = 1;
   }

#ifdef GLVIS_DEBUG
   cout << "OpenGL Visualization" << endl;
#endif

   GLenum mode = AUX_DOUBLE | AUX_RGBA | AUX_DEPTH;
   // mode |= (AUX_ALPHA | AUX_ACCUM);
   auxInitDisplayMode(mode);
   auxInitPosition(x, y, w, h);
   if (auxInitWindow(name) == GL_FALSE)
   {
      return 1;
   }

   Set_Texture_Image();

#ifdef GLVIS_DEBUG
   cout << "Window should be up" << endl;
#endif

#ifndef GLVIS_USE_FREETYPE
   if (fontname.empty())
   {
      fontname = "lucidasanstypewriter-14";
      // "-adobe-times-medium-r-normal-*-*-*-*-p-54-*-*";
      // "-*-bitstream vera sans-medium-r-normal-*-30-*-*-*-*-*-*-*";
      // "8x13";
      // "-adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1";
   }
   fontbase = tkLoadBitmapFont(fontname.c_str());
   if (fontbase == 0)
   {
      cerr << "Error loading font '" << fontname << '\'' << endl;
   }
#endif

   SetUseTexture(0);

   // auxReshapeFunc (MyReshape); // not needed, MyExpose calls it
   auxReshapeFunc (NULL);
   auxExposeFunc (MyExpose);

   auxMouseFunc (AUX_LEFTBUTTON, AUX_MOUSEDOWN, LeftButtonDown);
   auxMouseFunc (AUX_LEFTBUTTON, AUX_MOUSEUP, LeftButtonUp);
   auxMouseFunc (AUX_LEFTBUTTON, AUX_MOUSELOC, LeftButtonLoc);
   auxMouseFunc (AUX_MIDDLEBUTTON, AUX_MOUSEDOWN, MiddleButtonDown);
   auxMouseFunc (AUX_MIDDLEBUTTON, AUX_MOUSEUP, MiddleButtonUp);
   auxMouseFunc (AUX_MIDDLEBUTTON, AUX_MOUSELOC, MiddleButtonLoc);
   auxMouseFunc (AUX_RIGHTBUTTON, AUX_MOUSEDOWN, RightButtonDown);
   auxMouseFunc (AUX_RIGHTBUTTON, AUX_MOUSEUP, RightButtonUp);
   auxMouseFunc (AUX_RIGHTBUTTON, AUX_MOUSELOC, RightButtonLoc);

   // auxKeyFunc (AUX_p, KeyCtrlP); // handled in vsdata.cpp
   auxKeyFunc (AUX_S, KeyS);

   auxKeyFunc (AUX_q, KeyQPressed);
   auxKeyFunc (AUX_Q, KeyQPressed);

   auxModKeyFunc (AUX_LEFT, KeyLeftPressed);
   auxModKeyFunc (AUX_RIGHT, KeyRightPressed);
   auxModKeyFunc (AUX_UP, KeyUpPressed);
   auxModKeyFunc (AUX_DOWN, KeyDownPressed);

   auxKeyFunc (XK_KP_0, Key0Pressed);
   auxKeyFunc (XK_KP_1, Key1Pressed);
   auxKeyFunc (XK_KP_2, Key2Pressed);
   auxKeyFunc (XK_KP_3, Key3Pressed);
   auxKeyFunc (XK_KP_4, Key4Pressed);
   auxKeyFunc (XK_KP_5, Key5Pressed);
   auxKeyFunc (XK_KP_6, Key6Pressed);
   auxKeyFunc (XK_KP_7, Key7Pressed);
   auxKeyFunc (XK_KP_8, Key8Pressed);
   auxKeyFunc (XK_KP_9, Key9Pressed);

   auxKeyFunc (XK_KP_Subtract, KeyMinusPressed);
   auxKeyFunc (XK_KP_Add, KeyPlusPressed);

   auxKeyFunc (XK_KP_Decimal, KeyDeletePressed);
   auxKeyFunc (XK_KP_Enter, KeyEnterPressed);

   auxKeyFunc (XK_period, KeyDeletePressed);
   auxKeyFunc (TK_RETURN, KeyEnterPressed);

   auxKeyFunc (XK_0, Key0Pressed);
   auxKeyFunc (XK_1, Key1Pressed);
   auxKeyFunc (XK_2, Key2Pressed);
   auxKeyFunc (XK_3, Key3Pressed);
   auxKeyFunc (XK_4, Key4Pressed);
   auxKeyFunc (XK_5, Key5Pressed);
   auxKeyFunc (XK_6, Key6Pressed);
   auxKeyFunc (XK_7, Key7Pressed);
   auxKeyFunc (XK_8, Key8Pressed);
   auxKeyFunc (XK_9, Key9Pressed);

   auxKeyFunc (XK_minus, KeyMinusPressed);
   auxKeyFunc (XK_plus, KeyPlusPressed);
   auxKeyFunc (XK_equal, KeyPlusPressed);

   auxKeyFunc (AUX_j, KeyJPressed);
   auxKeyFunc (AUX_J, KeyJPressed);

   auxKeyFunc (XK_KP_Multiply, ZoomIn);
   auxKeyFunc (XK_KP_Divide, ZoomOut);

   auxKeyFunc (XK_asterisk, ZoomIn);
   auxKeyFunc (XK_slash, ZoomOut);

   auxKeyFunc (XK_bracketleft, ScaleDown);
   auxKeyFunc (XK_bracketright, ScaleUp);
   auxKeyFunc (XK_at, LookAt);

   auxKeyFunc(XK_parenleft, ShrinkWindow);
   auxKeyFunc(XK_parenright, EnlargeWindow);

   locscene = NULL;

   return 0;
}

void SendKeyEvent (KeySym keysym, int Shift=0)
{
   XKeyEvent  xke;

   xke.display = auxXDisplay();
   xke.window  = auxXWindow();
   xke.type    = KeyPress;
   xke.state   = 0;

   if (!Shift)
   {
      xke.keycode = XKeysymToKeycode(xke.display, keysym);
      XSendEvent(auxXDisplay(), auxXWindow(), True, KeyPressMask, (XEvent*)&xke);
   }
   else
   {
      xke.keycode = XKeysymToKeycode(xke.display, XK_Shift_L);
      XSendEvent(auxXDisplay(), auxXWindow(), True, KeyPressMask, (XEvent*)&xke);
      xke.state |= ShiftMask;

      xke.keycode = XKeysymToKeycode(xke.display, keysym);
      XSendEvent(auxXDisplay(), auxXWindow(), True, KeyPressMask, (XEvent*)&xke);

      xke.type = KeyRelease;
      XSendEvent(auxXDisplay(), auxXWindow(), True, KeyPressMask, (XEvent*)&xke);

      xke.keycode = XKeysymToKeycode(xke.display, XK_Shift_L);
      XSendEvent(auxXDisplay(), auxXWindow(), True, KeyPressMask, (XEvent*)&xke);
   }
}

void SendKeySequence (const char *seq)
{
   const char *key = seq;

   for ( ; *key != '\0'; key++ ) // see /usr/include/X11/keysymdef.h
   {
      if ( ((*key - '0') < 10) && ((*key - '0') >= 0) ) // (keypad) number
      {
         SendKeyEvent(XK_0 + (*key) -'0');
         continue;
      }

      if ( ((*key - 'a') < 26) && ((*key - 'a') >= 0) ) // lowercase letter
      {
         SendKeyEvent(XK_a + (*key) -'a');
         continue;
      }

      if ( ((*key - 'A') < 26) && ((*key - 'A') >= 0) ) // uppercase letter
      {
         SendKeyEvent(XK_A + (*key) -'A',1);
         continue;
      }

      switch (*key)
      {
         case '+':
            SendKeyEvent(XK_plus);
            continue;
         case '-':
            SendKeyEvent(XK_minus);
            continue;
         case '*':
            SendKeyEvent(XK_KP_Multiply);
            continue;
         case '/':
            SendKeyEvent(XK_KP_Divide);
            continue;
         case '.':
            SendKeyEvent(XK_period);
            continue;
         case '[':
            SendKeyEvent(XK_bracketleft);
            continue;
         case ']':
            SendKeyEvent(XK_bracketright);
            continue;
         case '(':
            SendKeyEvent(XK_parenleft,1);
            continue;
         case ')':
            SendKeyEvent(XK_parenright,1);
            continue;
         case '!':
            SendKeyEvent(XK_exclam,1);
            continue;
         case '~': // special codes
            key++;
            switch (*key)
            {
               case 'e': // expose event
                  SendExposeEvent();
                  break;
               case 'l': // left arrow
                  SendKeyEvent(XK_Left);
                  break;
               case 'r': // right arrow
                  SendKeyEvent(XK_Right);
                  break;
               case 'u': // up arrow
                  SendKeyEvent(XK_Up);
                  break;
               case 'd': // down arrow
                  SendKeyEvent(XK_Down);
                  break;
               case '3': // F3
                  SendKeyEvent(XK_F3);
                  break;
               case '5': // F5
                  SendKeyEvent(XK_F5);
                  break;
               case '6': // F6
                  SendKeyEvent(XK_F6);
                  break;
               case '7': // F7
                  SendKeyEvent(XK_F7);
                  break;
               case '.': // Keypad ./Del
                  SendKeyEvent(XK_period);
                  break;
               case 'E': // Keypad Enter
                  SendKeyEvent(XK_Return);
                  break;
            }
            continue;
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
         auxCallKeyFunc(*key, 0);
      }
      else
      {
         switch (*key)
         {
            case 'l': // left arrow
               auxCallKeyFunc(XK_Left, 0);
               break;
            case 'r': // right arrow
               auxCallKeyFunc(XK_Right, 0);
               break;
            case 'u': // up arrow
               auxCallKeyFunc(XK_Up, 0);
               break;
            case 'd': // down arrow
               auxCallKeyFunc(XK_Down, 0);
               break;
            case '3': // F3
               auxCallKeyFunc(XK_F3, 0);
               break;
            case '5': // F5
               auxCallKeyFunc(XK_F5, 0);
               break;
            case '6': // F6
               auxCallKeyFunc(XK_F6, 0);
               break;
            case '7': // F7
               auxCallKeyFunc(XK_F7, 0);
               break;
            case '.': // Keypad ./Del
               auxCallKeyFunc(XK_period, 0);
               break;
            case 'E': // Keypad Enter
               auxCallKeyFunc(XK_Return, 0);
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

   auxMainLoop(NULL);

   InitIdleFuncs();
}

void KillVisualization()
{
   delete locscene;
#ifndef GLVIS_USE_FREETYPE
   if (fontbase)
   {
      tkUnloadBitmapFont(fontbase);
      fontbase = 0;
   }
#endif
   auxCloseWindow();
}

void SendExposeEvent()
{
   if (disableSendExposeEvent) { return; }
   XExposeEvent ev;
   ev.type = Expose;
   ev.count = 0;
   ev.x = 0;
   ev.y = 0;
   ev.width = 100;
   ev.height = 100;
   XSendEvent (auxXDisplay(), auxXWindow(), 1, Expose, (XEvent*)&ev);
}

void MyReshape(GLsizei w, GLsizei h)
{
   glViewport (0, 0, w, h);

   glMatrixMode (GL_PROJECTION);
   glLoadIdentity();

   double ViewCenterX = locscene->ViewCenterX;
   double ViewCenterY = locscene->ViewCenterY;

   if (locscene->OrthogonalProjection)
   {
      double scale = locscene->ViewScale;
      if (w <= h)
      {
         glOrtho (-1.0, 1.0, -double(h)/w, double(h)/w, -10, 10);
      }
      else
      {
         glOrtho (-double(w)/h, double(w)/h, -1, 1, -10, 10);
      }
      glScaled(scale, scale, 1.0);
   }
   else
   {
      double ViewAngle = locscene->ViewAngle;

      if (w < h)
         ViewAngle = (360.0/M_PI)*atan( tan( ViewAngle*(M_PI/360.0) ) *
                                        double(h)/w );

      gluPerspective(ViewAngle, double(w)/h, 0.1, 5.0);
   }
   glTranslated(-ViewCenterX, -ViewCenterY, 0.0);
}

void MyExpose(GLsizei w, GLsizei h)
{
   MyReshape (w, h);
   locscene -> Draw();
}

void MyExpose()
{
   XWindowAttributes wa;

   XGetWindowAttributes(auxXDisplay(), auxXWindow(), &wa);
   MyExpose(wa.width, wa.height);
}


Array<void (*)()> IdleFuncs;
int LastIdleFunc;

void InitIdleFuncs()
{
   IdleFuncs.SetSize(0);
   LastIdleFunc = 0;
   auxIdleFunc(NULL);
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
   auxIdleFunc(MainIdleFunc);
}

void RemoveIdleFunc(void (*Func)(void))
{
   IdleFuncs.DeleteFirst(Func);
   if (IdleFuncs.Size() == 0)
   {
      auxIdleFunc(NULL);
   }
}


double xang = 0., yang = 0.;
double srot[16], sph_t, sph_u;
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
      Screenshot(fname);
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
   GLint viewport[4];
   double r, x, y, rr;
   const double maxr = 0.996194698091745532295;

   glGetIntegerv(GL_VIEWPORT, viewport);
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

void LeftButtonDown (AUX_EVENTREC *event)
{
   locscene -> spinning = 0;
   RemoveIdleFunc(MainLoop);

   oldx = event->data[AUX_MOUSEX];
   oldy = event->data[AUX_MOUSEY];

   ComputeSphereAngles(oldx, oldy, sph_u, sph_t);
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();
   locscene->cam.GLMultRotMatrix();
   glMultMatrixd(locscene->rotmat);
   glGetDoublev(GL_MODELVIEW_MATRIX, srot);
   glPopMatrix();

   startx = oldx;
   starty = oldy;
}

void LeftButtonLoc (AUX_EVENTREC *event)
{
   GLint newx = event->data[AUX_MOUSEX];
   GLint newy = event->data[AUX_MOUSEY];
   int sendexpose = 1;

   if (event->data[2] & ControlMask)
   {
      if (event->data[2] & ShiftMask)
      {
         locscene->PreRotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
      }
      else
      {
         double new_sph_u, new_sph_t;

         ComputeSphereAngles(newx, newy, new_sph_u, new_sph_t);

         glMatrixMode (GL_MODELVIEW);
         glLoadIdentity();

         double scoord[3], ncoord[3], inprod, cross[3];
         scoord[0] = scoord[1] = cos(sph_u);     scoord[2] = sin(sph_u);
         ncoord[0] = ncoord[1] = cos(new_sph_u); ncoord[2] = sin(new_sph_u);
         scoord[0] *= cos(sph_t);     scoord[1] *= sin(sph_t);
         ncoord[0] *= cos(new_sph_t); ncoord[1] *= sin(new_sph_t);
         inprod = InnerProd(scoord, ncoord);
         CrossProd(scoord, ncoord, cross);

         locscene->cam.GLMultTransposeRotMatrix();
         glRotated(acos(inprod)*(180.0/M_PI), cross[0], cross[1], cross[2]);
         glMultMatrixd(srot);
         glGetDoublev(GL_MODELVIEW_MATRIX, locscene->rotmat);
      }
   }
   else if (event->data[2] & Mod1Mask)
   {
      locscene->Rotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
   }
   else if (event->data[2] & ShiftMask)
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

void LeftButtonUp (AUX_EVENTREC *event)
{
   GLint newx = event->data[AUX_MOUSEX];
   GLint newy = event->data[AUX_MOUSEY];

   xang = (newx-startx)/5.0;
   yang = (newy-starty)/5.0;

   if ( (event->data[2] & ShiftMask) && (xang != 0.0 || yang != 0.0) )
   {
      locscene -> spinning = 1;
      AddIdleFunc(MainLoop);
      if (xang > 20) { xang = 20; } if (xang < -20) { xang = -20; }
      if (yang > 20) { yang = 20; } if (yang < -20) { yang = -20; }

      if (event->data[2] & ControlMask)
      {
         constrained_spinning = 1;
      }
      else
      {
         constrained_spinning = 0;
      }
   }
}

void MiddleButtonDown (AUX_EVENTREC *event)
{
   startx = oldx = event->data[AUX_MOUSEX];
   starty = oldy = event->data[AUX_MOUSEY];
}

void MiddleButtonLoc (AUX_EVENTREC *event)
{
   GLint newx = event->data[AUX_MOUSEX];
   GLint newy = event->data[AUX_MOUSEY];

   if ( !( event->data[2] & ControlMask ) )
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
      glGetIntegerv(GL_VIEWPORT, vp);
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

      if (event->data[2] & ShiftMask)  // ctrl + shift
      {
         double sx = double(newx-startx)/400;
         double sy = double(starty-newy)/400;

         locscene->cam.TurnLeftRight(dx-sx);
         locscene->cam.TurnUpDown(sy-dy);

         locscene->cam.TurnUpDown(-sy);
         locscene->cam.TurnLeftRight(sx);
      }
      else if (event->data[2] & Mod1Mask) // ctrl + alt
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

void MiddleButtonUp (AUX_EVENTREC *event)
{}

void RightButtonDown (AUX_EVENTREC *event)
{
   startx = oldx = event->data[AUX_MOUSEX];
   starty = oldy = event->data[AUX_MOUSEY];
}

void RightButtonLoc (AUX_EVENTREC *event)
{
   GLint newx = event->data[AUX_MOUSEX];
   GLint newy = event->data[AUX_MOUSEY];

   if (event->data[2] & ShiftMask)
   {
      glLoadIdentity();
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
      glLightfv(GL_LIGHT0, GL_POSITION, light);
   }
   else if ( !( event->data[2] & ControlMask ) )
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

void RightButtonUp (AUX_EVENTREC *event)
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
   cout << "Screenshot: glXWaitX() ... " << flush;
#endif
   glXWaitX();
#ifdef GLVIS_DEBUG
   cout << "done." << endl;
#endif

#ifdef GLVIS_DEBUG
   cout << "Screenshot: glFinish() ... " << flush;
#endif
   glFinish();
   // events in the X event queue may not be complete
   // in particular ExposeEvents generated by SendExposeEvent()
   // or keys sent by SendKeySequence
#ifdef GLVIS_DEBUG
   cout << "done." << endl;
#endif

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

#if defined(GLVIS_USE_LIBTIFF)
   // Save a TIFF image. This requires the libtiff library, see www.libtiff.org
   TIFF* image;

   XWindowAttributes wa;
   XGetWindowAttributes(auxXDisplay(), auxXWindow(), &wa);
   int w = wa.width;
   int h = wa.height;
   // MyExpose(w,h);
   glReadBuffer(GL_FRONT);

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

   XWindowAttributes wa;
   XGetWindowAttributes(auxXDisplay(), auxXWindow(), &wa);
   int w = wa.width;
   int h = wa.height;
   glReadBuffer(GL_FRONT);

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

#else
   // Use the external X Window Dump (xwd) tool.
   // Note that xwd does not work on OS X!
   ostringstream cmd;
   cmd << "xwd -silent -out " << filename << " -nobdrs -id " << auxXWindow();
   if (system(cmd.str().c_str()))
   {
      return 1;
   }
   // View with xwud -in GLVis_s*.xwd, or use convert GLVis_s*.xwd
   // GLVis_s*.{jpg,gif}
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
      Screenshot(fname);
      cout << "done" << endl;
   }
}

void KeyCtrlP()
{
   int state, buffsize;
   FILE * fp;
   GLint viewport[4];

   cout << "Printing the figure to GLVis.pdf... " << flush;

   fp = fopen("GLVis.pdf", "wb");
   buffsize = 0;
   state = GL2PS_OVERFLOW;
   locscene -> print = 1;
   glGetIntegerv(GL_VIEWPORT, viewport);
   while (state == GL2PS_OVERFLOW)
   {
      buffsize += 1024*1024;
      gl2psBeginPage ( "GLVis.pdf", "GLVis", viewport,
                       GL2PS_PDF, // or GL2PS_SVG, or GL2PS_EPS
                       GL2PS_BSP_SORT,
                       GL2PS_SIMPLE_LINE_OFFSET |
                       // GL2PS_NO_PS3_SHADING |
                       // GL2PS_NO_BLENDING |
                       // GL2PS_OCCLUSION_CULL |
                       // GL2PS_BEST_ROOT |
                       GL2PS_SILENT |
                       GL2PS_DRAW_BACKGROUND,
                       GL_RGBA, 0, NULL, 16, 16, 16, buffsize, fp, "a" );
      gl2psPointSize(.4);
      gl2psLineWidth(.2);
      locscene -> Draw();
      state = gl2psEndPage();
   }
   locscene -> print = 0;
   fclose(fp);

   cout << "done" << endl;
   locscene -> Draw();
}

void KeyQPressed()
{
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
   if (state & ControlMask)
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
   if (state & ControlMask)
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
   if (state & ControlMask)
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
   if (state & ControlMask)
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
   if (state & ControlMask)
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

   glGetIntegerv(GL_VIEWPORT, viewport);
   int w = viewport[2], h = viewport[3];
   w = (int)ceil(w / window_scale_factor);
   h = (int)ceil(h / window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   XResizeWindow(auxXDisplay(), auxXWindow(), w, h);
}

void EnlargeWindow()
{
   GLint viewport[4];

   glGetIntegerv(GL_VIEWPORT, viewport);
   int w = viewport[2], h = viewport[3];
   w = (int)ceil(w * window_scale_factor);
   h = (int)ceil(h * window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   XResizeWindow(auxXDisplay(), auxXWindow(), w, h);
}

void MoveResizeWindow(int x, int y, int w, int h)
{
   XMoveResizeWindow(auxXDisplay(), auxXWindow(), x, y, w, h);
}

void ResizeWindow(int w, int h)
{
   XResizeWindow(auxXDisplay(), auxXWindow(), w, h);
}

void SetWindowTitle(const char *title)
{
   XSetStandardProperties(auxXDisplay(), auxXWindow(),
                          title, NULL, None, NULL, 0, NULL);
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
int RepeatPaletteTimes = 1;
int UseTexture         = 0;

void MySetColor (double val, double min, double max)
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
      MySetColor (log(fabs(val/(min+eps))) / (log(fabs(max/(min+eps)))+eps));
   }
   else
   {
      MySetColor ((val-min)/(max-min));
   }
}

void MySetColor (double val)
{
   int i;
   double t, r, g, b, *pal;

   if (UseTexture)
   {
      glTexCoord1d(val);
      return;
   }

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

   val *= 0.999999999 * ( RGB_Palette_Size - 1 ) * abs(RepeatPaletteTimes);
   i = (int) floor( val );
   t = val - i;

   if (((i / (RGB_Palette_Size-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
       ((i / (RGB_Palette_Size-1)) % 2 == 1 && RepeatPaletteTimes < 0))
   {
      pal = RGB_Palette + 3 * ( i % (RGB_Palette_Size-1) );
   }
   else
   {
      pal = RGB_Palette + 3 * ( (RGB_Palette_Size-2) -
                                i % (RGB_Palette_Size-1) );
      t = 1.0 - t;
   }

   r = (1.0 - t) * pal[0] + t * pal[3];
   g = (1.0 - t) * pal[1] + t * pal[4];
   b = (1.0 - t) * pal[2] + t * pal[5];

   if (MatAlpha < 1.0)
   {
      glColor4f ( r, g, b, malpha );
   }
   else
   {
      glColor3f ( r, g, b );
   }
}

// const int Max_Texture_Size = 512;
const int Max_Texture_Size = 4*1024;
int Texture_Size;
GLfloat Texture_Image[3*Max_Texture_Size];

void Make_Texture_From_Palette()
{
   int i, j;
   double t, *pal;

   Texture_Size = Max_Texture_Size;

   for (i = 0; i < Texture_Size; i++)
   {
      t = double(i) / (Texture_Size - 1);
      t *= 0.999999999 * ( RGB_Palette_Size - 1 ) * abs(RepeatPaletteTimes);
      j = (int) floor(t);
      t -= j;

      if (((j / (RGB_Palette_Size-1)) % 2 == 0 && RepeatPaletteTimes > 0) ||
          ((j / (RGB_Palette_Size-1)) % 2 == 1 && RepeatPaletteTimes < 0))
      {
         pal = RGB_Palette + 3 * ( j % (RGB_Palette_Size-1) );
      }
      else
      {
         pal = RGB_Palette + 3 * ( (RGB_Palette_Size-2) -
                                   j % (RGB_Palette_Size-1) );
         t = 1.0 - t;
      }

      Texture_Image[3*i+0] = (1.0 - t) * pal[0] + t * pal[3];
      Texture_Image[3*i+1] = (1.0 - t) * pal[1] + t * pal[4];
      Texture_Image[3*i+2] = (1.0 - t) * pal[2] + t * pal[5];
   }
}

void Make_Texture_From_Palette_2()
{
   if (RGB_Palette_Size > Max_Texture_Size)
   {
      Texture_Size = Max_Texture_Size;
   }
   else
   {
      Texture_Size = RGB_Palette_Size;
   }

   if (RepeatPaletteTimes > 0)
   {
      for (int i = 0; i < 3*Texture_Size; i++)
      {
         Texture_Image[i] = RGB_Palette[i];
      }
   }
   else
   {
      for (int i = 0; i < Texture_Size; i++)
      {
         Texture_Image[3*i+0] = RGB_Palette[3*(Texture_Size-1-i)+0];
         Texture_Image[3*i+1] = RGB_Palette[3*(Texture_Size-1-i)+1];
         Texture_Image[3*i+2] = RGB_Palette[3*(Texture_Size-1-i)+2];
      }
   }
}

void Write_Texture_To_File()
{
   const char ppm_fname[] = "GLVis_texture.ppm";
   cout << "Saving texture image --> " << flush;
   ofstream ppm_file(ppm_fname);
   ppm_file << "P3\n" << Texture_Size << " 1\n255\n";
   for (int i = 0; i < Texture_Size; i++)
      for (int j = 0; j < 3; j++)
      {
         ppm_file << ' ' << int(floor(255.*Texture_Image[3*i+j]+0.5));
      }
   ppm_file << endl;
   cout << ppm_fname << endl;
}

void Set_Texture_Image()
{
   if (UseTexture == 1)
   {
      Make_Texture_From_Palette_2();
   }
   else
   {
      Make_Texture_From_Palette();
   }

   glTexImage1D(GL_TEXTURE_1D, // GLenum target,
                0,             // GLint level,
                3,             // GLint internalFormat,
                Texture_Size,  // GLsizei width,
                0,             // GLint border,
                GL_RGB,        // GLenum format,
                GL_FLOAT,      // GLenum type,
                Texture_Image  // const GLvoid *pixels
               );
   glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
   // glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

   // glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   // glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

   // glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

   // glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

   // Write_Texture_To_File();
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
      Set_Texture_Image();
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


#ifdef GLVIS_USE_FREETYPE

class GLVisFont
{
private:
   int         init;

   FT_Library  library;
   FT_Face     face;

   FT_Bool     use_kerning;
   FT_Glyph   *glyph;
   FT_UInt     alloc_glyphs, num_glyphs;
   FT_BBox     bbox;

   int            image_width, image_height;
   unsigned char *image;

   void LoadSequence(const char *text)
   {
      int          err;
      FT_UInt      len, glyph_idx, prev_glyph, i;
      FT_Vector    pen, delta;

      FreeSeq();

      len      = strlen(text);
      glyph    = new FT_Glyph[len];
      alloc_glyphs = len;

      pen.x = pen.y = 0;
      prev_glyph = i = 0;
      for (FT_UInt ic = 0; ic < len; ic++)
      {
         glyph_idx = FT_Get_Char_Index(face, text[ic]);

         if (use_kerning && prev_glyph && glyph_idx)
         {
            FT_Get_Kerning(face, prev_glyph, glyph_idx,
                           FT_KERNING_DEFAULT, &delta);
            pen.x += delta.x;
         }

         err = FT_Load_Glyph(face, glyph_idx, FT_LOAD_DEFAULT);
         if (err) { continue; }

         err = FT_Get_Glyph(face->glyph, &glyph[i]);
         if (err) { continue; }

         FT_Glyph_Transform(glyph[i], 0, &pen);

         pen.x += face->glyph->advance.x;
         prev_glyph = glyph_idx;

         i++;
      }
      num_glyphs = i;
   }

   void GetBBox()
   {
      FT_BBox glyph_bbox;

      if (num_glyphs <= 0)
      {
         bbox.xMin = bbox.yMin = 0;
         bbox.xMax = bbox.yMax = 0;
         return;
      }

      FT_Glyph_Get_CBox(glyph[0], FT_GLYPH_BBOX_PIXELS, &bbox);
      for (FT_UInt i = 1; i < num_glyphs; i++)
      {
         FT_Glyph_Get_CBox(glyph[i], FT_GLYPH_BBOX_PIXELS, &glyph_bbox);

         if (glyph_bbox.xMin < bbox.xMin)
         {
            bbox.xMin = glyph_bbox.xMin;
         }
         if (glyph_bbox.xMax > bbox.xMax)
         {
            bbox.xMax = glyph_bbox.xMax;
         }
         if (glyph_bbox.yMin < bbox.yMin)
         {
            bbox.yMin = glyph_bbox.yMin;
         }
         if (glyph_bbox.yMax > bbox.yMax)
         {
            bbox.yMax = glyph_bbox.yMax;
         }
      }
   }

   void FreeSeq()
   {
      if (alloc_glyphs <= 0)
      {
         return;
      }

      for (FT_UInt i = 0; i < num_glyphs; i++)
      {
         FT_Done_Glyph(glyph[i]);
      }

      delete [] glyph;

      glyph = NULL;

      alloc_glyphs = num_glyphs = 0;
   }

   static int FindFontFile(const char *font_patterns[], int num_patterns,
                           string &font_file)
   {
      FcObjectSet *os  = NULL;
      FcFontSet   *fs  = NULL;
      FcPattern   *pat = NULL;
      FcChar8     *s   = NULL;
#ifdef GLVIS_DEBUG
      FcChar8     *fnt;
#endif
      FcResult     res;

      if (num_patterns <= 0)
      {
         return -2;
      }

      if (!FcInit())
      {
         return -1;
      }

      os = FcObjectSetBuild(FC_FAMILY, FC_STYLE, FC_FILE, (void *)(NULL));

      for (int i = 0; i < num_patterns && !s; i++)
      {
         pat = FcNameParse((FcChar8 *)font_patterns[i]);
         if (!pat)
         {
            continue;
         }

         fs = FcFontList(0, pat, os);
         FcPatternDestroy(pat);
         if (!fs)
         {
            continue;
         }

#ifdef GLVIS_DEBUG
         if (fs->nfont > 1)
         {
            cout <<
                 "Font pattern '" << font_patterns[i] << "' matched"
                 " multiple fonts:\n";
            for (int j = 0; j < fs->nfont; j++)
            {
               fnt = FcNameUnparse(fs->fonts[j]);
               cout << fnt << endl;
               free(fnt);
            }
            cout << "-----" << endl;
         }
#endif

         for (int j = 0; j < fs->nfont; j++)
         {
            res = FcPatternGetString(fs->fonts[j], FC_FILE, 0, &s);
            if (res == FcResultMatch && s)
            {
               font_file = (char *)s;
#ifdef GLVIS_DEBUG
               fnt = FcNameUnparse(fs->fonts[j]);
               cout << "Using font: " << fnt << endl;
               free(fnt);
#endif
               break;
            }
         }

         FcFontSetDestroy(fs);
      }

      if (os)
      {
         FcObjectSetDestroy(os);
      }

      FcFini();

      return s ? 0 : -2;
   }

public:
   GLVisFont() { init = 0; }

   int Initialized() const { return init; }

   void Init(const char *font_patterns[], int num_patterns, int height)
   {
      SetFont(font_patterns, num_patterns, height);
   }

   int SetFont(const char *font_patterns[], int num_patterns, int height)
   {
      string font_file;

      if (FindFontFile(font_patterns, num_patterns, font_file))
      {
         if (!init)
         {
            init = -4;
         }
         return -4;
      }

#ifdef GLVIS_DEBUG
      cout << "Using font file: " << font_file << endl;
#endif

      return SetFontFile(font_file.c_str(), height);
   }

   int SetFontFile(const char *font_file, int height)
   {
      int err;

      if (init == -1) // library init failed
      {
         return -1;
      }

      if (init <= 0)
      {
         err = FT_Init_FreeType(&library);
         if (err)
         {
            cout << "GLVis: Can not initialize FreeType library!" << endl;
            return (init = -1);
         }
      }
      else
      {
         FreeSeq();
         delete [] image;
         FT_Done_Face(face);
      }

      err = FT_New_Face(library, font_file, 0, &face);
      if (err)
      {
         cout << "GLVis: Can not open font file: " << font_file << endl;
         FT_Done_FreeType(library);
         return (init = -2);
      }

      if (1)
      {
         // set font height in points
         Screen *scr = tkGetXScreen();
         int ppi_w = (int)rint(25.4*WidthOfScreen(scr)/WidthMMOfScreen(scr));
         int ppi_h = (int)rint(25.4*HeightOfScreen(scr)/HeightMMOfScreen(scr));

         err = FT_Set_Char_Size(face, 0, height*64, ppi_w, ppi_h);
         if (err)
         {
            cout << "GLVis: Can not set font height: " << height << " pts"
                 << endl;
            FT_Done_Face(face);
            FT_Done_FreeType(library);
            return (init = -3);
         }
      }
      else
      {
         // set font height in pixels
         err = FT_Set_Pixel_Sizes(face, 0, height);
         if (err)
         {
            cout << "GLVis: Can not set font height: " << height << " pixels"
                 << endl;
            FT_Done_Face(face);
            FT_Done_FreeType(library);
            return (init = -3);
         }
      }

      glyph = NULL;
      use_kerning = FT_HAS_KERNING(face);

      image = NULL;

      init = 1;

      return 0;
   }

   int Render(const char *text)
   {
      int            err;
      FT_Glyph       a_glyph;
      FT_BitmapGlyph bit_glyph;
      FT_Bitmap     *bitmap;

      if (init <= 0)
      {
         return 1;
      }

      delete [] image;
      image = NULL;

      LoadSequence(text);

      if (num_glyphs <= 0)
      {
         return 2;
      }

      GetBBox();

      const int pad = 2;

      image_width  = (bbox.xMax - bbox.xMin) + 2*pad;
      image_height = (bbox.yMax - bbox.yMin) + 2*pad;

      if (image_width <= 0 || image_height <= 0)
      {
         return 3;
      }

      GLfloat col[4];
      glGetFloatv(GL_CURRENT_COLOR, col);

      image = new unsigned char[4*image_width*image_height];
      for (int i = 0; i < 4*image_width*image_height; i++)
      {
         image[i] = 0;
      }

      for (FT_UInt g = 0; g < num_glyphs; g++)
      {
         a_glyph = glyph[g];

         err = FT_Glyph_To_Bitmap(&a_glyph, FT_RENDER_MODE_NORMAL, 0, 0);
         if (err) { continue; }

         bit_glyph = (FT_BitmapGlyph)a_glyph;
         bitmap = &bit_glyph->bitmap;

         int off_i = bit_glyph->left - bbox.xMin + pad;
         int off_j = image_height - bit_glyph->top + bbox.yMin - pad;

         for (int j = 0; j < (int) bitmap->rows; j++)
         {
            int im_j = image_height - 1 - (j + off_j);
            if (im_j < 0 || im_j >= image_height)
            {
#ifdef GLVIS_DEBUG
               cout <<
                    "GLVisFont::Render : outside 'y' range!\n"
                    "   text            = " << text            << "\n"
                    "   j               = " << j               << "\n"
                    "   bitmap->width   = " << bitmap->width   << "\n"
                    "   bitmap->rows    = " << bitmap->rows    << "\n"
                    "   off_i           = " << off_i           << "\n"
                    "   off_j           = " << off_j           << "\n"
                    "   bit_glyph->left = " << bit_glyph->left << "\n"
                    "   bit_glyph->top  = " << bit_glyph->top  << endl;
#endif
               continue;
            }
            for (int i = 0; i < (int) bitmap->width; i++)
            {
               int im_i = i + off_i;
               if (im_i < 0 || im_i >= image_width)
               {
#ifdef GLVIS_DEBUG
                  cout <<
                       "GLVisFont::Render : outside 'x' range!\n"
                       "   text            = " << text            << "\n"
                       "   i               = " << i               << "\n"
                       "   j               = " << j               << "\n"
                       "   bitmap->width   = " << bitmap->width   << "\n"
                       "   bitmap->rows    = " << bitmap->rows    << "\n"
                       "   off_i           = " << off_i           << "\n"
                       "   off_j           = " << off_j           << "\n"
                       "   bit_glyph->left = " << bit_glyph->left << "\n"
                       "   bit_glyph->top  = " << bit_glyph->top  << endl;
#endif
                  continue;
               }
               int off = 4*(im_i + im_j*image_width);
               int val = bitmap->buffer[i + j*bitmap->width];
               image[off + 0] = (unsigned char)(col[0]*255);
               image[off + 1] = (unsigned char)(col[1]*255);
               image[off + 2] = (unsigned char)(col[2]*255);
               image[off + 3] = val;
            }
         }

         FT_Done_Glyph(a_glyph);
      }

      FreeSeq();

      return 0;
   }

   const unsigned char *GetImage() const { return image; }
   int GetImageWidth() const { return image_width; }
   int GetImageHeight() const { return image_height; }

   ~GLVisFont()
   {
      if (init <= 0)
      {
         return;
      }

      FreeSeq();

      delete [] image;

      FT_Done_Face(face);
      FT_Done_FreeType(library);
   }
};

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

int font_size = 12;

GLVisFont glvis_font;

int RenderBitmapText(const char *text, int &width, int &height)
{
   if (!glvis_font.Initialized())
   {
      if (fontname.empty())
      {
         glvis_font.Init(fc_font_patterns, num_font_patterns, font_size);
      }
      else
      {
         const char *fc_pat[1];
         fc_pat[0] = fontname.c_str();
         glvis_font.Init(fc_pat, 1, font_size);
         if (glvis_font.Initialized() == -4)
         {
            cout << "Font not found: " << fontname << endl;
            glvis_font.Init(fc_font_patterns, num_font_patterns, font_size);
         }
      }

      if (glvis_font.Initialized() == -4)
         cout <<
              "GLVis: No fonts found! Use the -fn option or"
              " edit 'fc_font_patterns' in lib/aux_vis.cpp" << endl;
   }

   int fail = glvis_font.Render(text);

   if (!fail)
   {
      width = glvis_font.GetImageWidth();
      height = glvis_font.GetImageHeight();
   }
   else
   {
      width = height = 0;
   }

   return !fail;
}

void DrawBitmapText()
{
   if (glvis_font.Initialized() > 0 && glvis_font.GetImage())
   {
      glPushAttrib(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glEnable(GL_BLEND);
      glDepthMask(GL_FALSE);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      glDrawPixels(
         glvis_font.GetImageWidth(),
         glvis_font.GetImageHeight(),
         GL_RGBA,
         GL_UNSIGNED_BYTE,
         glvis_font.GetImage());

      glPopAttrib();
   }
}

void DrawBitmapText(const char *text)
{
   int width, height;

   if (RenderBitmapText(text, width, height))
   {
      DrawBitmapText();
   }
}

int SetFontFile(const char *font_file, int height)
{
   return glvis_font.SetFontFile(font_file, height);
}

int SetFont(const char *font_patterns[], int num_patterns, int height)
{
   return glvis_font.SetFont(font_patterns, num_patterns, height);
}

#endif // GLVIS_USE_FREETYPE

void SetFont(const char *fn)
{
   fontname = fn;
#ifndef GLVIS_USE_FREETYPE
   if (visualize)
   {
      if (fontbase)
      {
         tkUnloadBitmapFont(fontbase);
      }
      fontbase = tkLoadBitmapFont(fontname.c_str());
      if (fontbase == 0)
      {
         cerr << "Error loading font '" << fontname << '\'' << endl;
      }
   }
#else
   size_t pos = fontname.rfind('-');
   if (pos != string::npos)
   {
      font_size = atoi(fontname.substr(pos + 1).c_str());
      fontname.erase(pos);
   }
#ifdef GLVIS_DEBUG
   cout << "SetFont: name = '" << fontname << "', height = " << font_size
        << endl;
#endif
   if (glvis_font.Initialized())
   {
      const char *fc_pat[1];
      fc_pat[0] = fontname.c_str();
      if (SetFont(fc_pat, 1, font_size) == -4)
      {
         cout << "Font not found: " << fontname << endl;
      }
   }
#endif
}
