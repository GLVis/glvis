// Copyright (c) 2010, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-443271. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see http://glvis.googlecode.com.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.

#include <iostream>
#include <sstream>
#include <math.h>

#include "mfem.hpp"

#include "palettes.hpp"
#include "aux_vis.hpp"
#include "gl2ps.h"
#include <X11/keysym.h>
#include <time.h>

#ifdef GLVIS_USE_LIBTIFF
#include "tiffio.h"
#endif

extern Window window;

int fontbase;
int visualize;
VisualizationScene * locscene;

float MatAlpha = 1.0;
float MatAlphaCenter = 0.5;

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
      return 1;

   Set_Texture_Image();

#ifdef GLVIS_DEBUG
   cout << "Window should be up" << endl;
#endif

   const char fontname[] = "lucidasanstypewriter-18";
   // "-adobe-times-medium-r-normal-*-*-*-*-p-54-*-*";
   // "-*-bitstream vera sans-medium-r-normal-*-30-*-*-*-*-*-*-*";
   // "8x13";
   // "-adobe-helvetica-medium-r-normal--12-120-75-75-p-67-iso8859-1";
   fontbase = tkLoadBitmapFont(fontname);
   if (fontbase == 0)
      cerr << "Error loading font '" << fontname << '\'' << endl;

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

   auxKeyFunc (AUX_P, KeyP);
   auxKeyFunc (AUX_S, KeyS);

   auxKeyFunc (AUX_q, KeyQPressed);
   auxKeyFunc (AUX_Q, KeyQPressed);

   auxKeyFunc (AUX_LEFT, KeyLeftPressed);
   auxKeyFunc (AUX_RIGHT, KeyRightPressed);
   auxKeyFunc (AUX_UP, KeyUpPressed);
   auxKeyFunc (AUX_DOWN, KeyDownPressed);

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

void SendKeySequence (char * seq)
{
   char * key = seq;

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

void InitIdleFuncs();

void SetVisualizationScene(VisualizationScene * scene, int view, char * keys)
{
   locscene = scene;

   locscene -> view = view;
   if (view == 2)
      scene -> CenterObject2D();
   else
      scene -> CenterObject();

   InitIdleFuncs();
   if (scene -> spinning)
      AddIdleFunc(MainLoop);

   if (keys)
      SendKeySequence(keys);

   auxMainLoop(NULL);

   InitIdleFuncs();
}

void KillVisualization()
{
   delete locscene;
   if (fontbase)
   {
      tkUnloadBitmapFont(fontbase);
      fontbase = 0;
   }
   auxCloseWindow();
}

void SendExposeEvent()
{
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
         glOrtho (-1.0, 1.0, -double(h)/w, double(h)/w, -10, 10);
      else
         glOrtho (-double(w)/h, double(w)/h, -1, 1, -10, 10);
      glScaled(scale, scale, 1.0);
      glTranslated(-ViewCenterX, -ViewCenterY, 0.0);
   }
   else
   {
      double ViewAngle = locscene->ViewAngle;
      double Distance  = 2.5;

      if (w < h)
         ViewAngle = (360.0/M_PI)*atan( tan( ViewAngle*(M_PI/360.0) ) *
                                        double(h)/w );

      gluPerspective(ViewAngle, double(w)/h, 0.4*Distance, 1.65*Distance);
      glTranslated(-ViewCenterX, -ViewCenterY, -Distance);
   }
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
      (*IdleFuncs[LastIdleFunc])();
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
      auxIdleFunc(NULL);
}


double xang, yang;
double srot[16], sph_t, sph_u;
static GLint oldx, oldy, startx, starty;

int constrained_spinning = 0;


void MainLoop()
{
   static int p = 1;
   struct timespec req;
   if (locscene -> spinning)
   {
      if (!constrained_spinning)
      {
         glMatrixMode (GL_MODELVIEW);
         glLoadIdentity();
         glRotatef(xang, 0.0f, 1.0f, 0.0f);
         glRotatef(yang, 1.0f, 0.0f, 0.0f);
         glMultMatrixd (locscene -> rotmat);
         glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
         locscene -> Draw();
      }
      else
      {
         glMatrixMode (GL_MODELVIEW);
         glLoadMatrixd (locscene -> rotmat);
         glRotated ( xang, 0.0, 0.0, 1.0 );
         glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
         locscene -> Draw();
      }
      req.tv_sec  = 0;
      req.tv_nsec = 10000000;
      nanosleep (&req, NULL);  // sleep for 0.01 seconds
   }
   if (locscene -> movie)
   {
      char fname[20];
      sprintf(fname, "GLVis_m%04d", p++);
      Screenshot(fname);
   }
}

//== PRESSED MOUSE BUTTONS EVENTS ============================================

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
      x *= maxr/rr, y *= maxr/rr, rr = maxr;

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
   double *rot = locscene->rotmat;
   for (int i = 0; i < 16; i++)
      srot[i] = rot[i];

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
         glMatrixMode (GL_MODELVIEW);
         glLoadMatrixd (locscene -> rotmat);
         glRotated ( double (newx-oldx) / 2, 0.0, 0.0, 1.0 );
         glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
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
         inprod = scoord[0]*ncoord[0]+scoord[1]*ncoord[1]+scoord[2]*ncoord[2];
         cross[0] =  scoord[1] * ncoord[2] - ncoord[1] * scoord[2];
         cross[1] =  ncoord[0] * scoord[2] - scoord[0] * ncoord[2];
         cross[2] =  scoord[0] * ncoord[1] - ncoord[0] * scoord[1];
         inprod = acos(inprod);

         glRotated(inprod*(180.0/M_PI), cross[0], cross[1], cross[2]);
         glMultMatrixd(srot);
         glGetDoublev (GL_MODELVIEW_MATRIX, locscene->rotmat);
      }
   }
   else if (event->data[2] & Mod1Mask)
   {
      glMatrixMode (GL_MODELVIEW);
      glLoadIdentity();

      glRotated(double (newx-oldx) / 2, 0.0, 0.0, 1.0);
      glMultMatrixd(locscene->rotmat);
      glGetDoublev (GL_MODELVIEW_MATRIX, locscene->rotmat);
   }
   else if (event->data[2] & ShiftMask)
   {
      locscene -> Rotate (double (newx-oldx) / 2, double (newy-oldy) / 2);
   }
   else
   {
      // locscene -> Rotate (double (newx-oldx) / 2, double (newy-oldy) / 2);
      glMatrixMode (GL_MODELVIEW);
      glLoadIdentity();
      glRotated ( double (newy-oldy) / 2, 1.0, 0.0, 0.0 );
      glMultMatrixd (locscene -> rotmat);
      glRotated ( double (newx-oldx) / 2, 0.0, 0.0, 1.0 );
      glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   }

   oldx = newx;
   oldy = newy;

   if (sendexpose)
      SendExposeEvent();
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
      if (xang > 20) xang = 20; if (xang < -20) xang = -20;
      if (yang > 20) yang = 20; if (yang < -20) yang = -20;

      if (event->data[2] & ControlMask)
         constrained_spinning = 1;
      else
         constrained_spinning = 0;
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
         scale = locscene->ViewScale;
      else
         scale = 0.4142135623730950488/tan(locscene->ViewAngle*(M_PI/360));
      glGetIntegerv(GL_VIEWPORT, vp);
      if (vp[2] < vp[3])
         scale *= vp[2];
      else
         scale *= vp[3];
      TrX = 2.0*double(oldx-newx)/scale;
      TrY = 2.0*double(newy-oldy)/scale;
      locscene->ViewCenterX += TrX;
      locscene->ViewCenterY += TrY;
   }
   else
      locscene -> Translate ((double)(newx-oldx)/200,(double)(newy-oldy)/200);

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
         z = sqrt (1. - l*l);
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
      GLfloat light[] = { x, y, z, 0.0 };
      glLightfv(GL_LIGHT0, GL_POSITION, light);
   }
   else if ( !( event->data[2] & ControlMask ) )
      locscene -> Zoom (exp ( double (oldy-newy) / 100 ));
   else
      locscene -> Scale ( exp ( double (oldy-newy) / 50 ) );

   SendExposeEvent();

   oldx = newx;
   oldy = newy;
}

void RightButtonUp (AUX_EVENTREC *event)
{}

int Screenshot(const char *fname)
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

#ifdef GLVIS_USE_LIBTIFF
   // Save a TIFF image. This requires the libtiff library, see www.libtiff.org
   TIFF* image;

   char filename[100];
   sprintf(filename,"%s.tif",fname);
   image = TIFFOpen(filename, "w");
   if (!image)
      return 1;

   XWindowAttributes wa;
   XGetWindowAttributes(auxXDisplay(), auxXWindow(), &wa);
   int w = wa.width;
   int h = wa.height;
   // MyExpose(w,h);
   glReadBuffer(GL_FRONT);

   unsigned char** pixels = new unsigned char*[ h ];
   int idx;
   for (idx=0; idx<h; idx++)
   {
      pixels[idx] = new unsigned char[w*3];
      glReadPixels(0, idx, w, 1, GL_RGB, GL_UNSIGNED_BYTE, pixels[idx]);
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

   for (idx=0; idx<h; idx++)
      if (TIFFWriteScanline(image, pixels[h-idx-1], idx, 0) < 0)
         return 2;

   for (idx=0; idx<h; idx++)
      delete[] pixels[idx];
   delete[] pixels;

   TIFFFlushData(image);
   TIFFClose(image);
   return 0;
#else
   // Use the external X Window Dump (xwd) tool.
   // Note that xwd does not work on OS X!
   ostringstream cmd;
   cmd << "xwd -silent -out " << fname << ".xwd -nobdrs -id " << window;
   return system(cmd.str().c_str());
   // View with xwud -in GLVis_s*.xwd, or use convert GLVis_s*.xwd
   // GLVis_s*.{jpg,gif}
#endif
}

void KeyS()
{
   static int p = 1;

   if (locscene -> spinning)
   {
      locscene -> movie = 1 - locscene -> movie;
      if (locscene -> movie)
         cout << "Recording a movie (series of snapshots)..." << endl;
      else
         cout << endl;
      // use (ImageMagik's) convert GLVis_m* GLVis.{gif,mpg}
   }
   else
   {
      cout << "Taking snapshot number " << p << "... ";
      char fname[20];
      sprintf(fname, "GLVis_s%02d", p++);
      Screenshot(fname);
      cout << "done" << endl;
   }
}

void KeyP()
{
   int state, buffsize;
   FILE * fp;
   GLint viewport[4];

   cout << "Printing the figure to GLVis.eps... " << flush;

   fp = fopen("GLVis.eps", "wb");
   buffsize = 0;
   state = GL2PS_OVERFLOW;
   locscene -> print = 1;
   glGetIntegerv(GL_VIEWPORT, viewport);
   while (state == GL2PS_OVERFLOW)
   {
      buffsize += 1024*1024;
      gl2psBeginPage ( "GLVis.eps", "GLVis", viewport,
                       /*
                         GL2PS_EPS, // or GL2PS_PDF,
                         GL2PS_BSP_SORT,
                         GL2PS_SIMPLE_LINE_OFFSET  |
                         GL2PS_NO_BLENDING |
                         GL2PS_SILENT |
                         GL2PS_OCCLUSION_CULL |
                         GL2PS_BEST_ROOT,
                         GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, "a" );
                       */
                       GL2PS_EPS,
                       GL2PS_BSP_SORT,
                       GL2PS_NO_PS3_SHADING | GL2PS_DRAW_BACKGROUND,
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

void CheckSpin()
{
   if (fabs(xang) < 1.e-2)
      xang = 0.;
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
   xang--;
   CheckSpin();
}

void KeyDeletePressed()
{
   if (locscene -> spinning)
   {
      xang = yang = 0.;
      locscene -> spinning = 0;
      RemoveIdleFunc(MainLoop);
   }
   else
   {
      xang = 1.;
      locscene -> spinning = 1;
      AddIdleFunc(MainLoop);
      constrained_spinning = 1;
   }
}

void KeyEnterPressed()
{
   xang++;
   CheckSpin();
}

void Key7Pressed()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadMatrixd (locscene -> rotmat);
   glRotated ( 1.0, 0.0, -1.0, 0.0 );
   glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   SendExposeEvent();
}

void Key8Pressed()
{
   locscene -> Rotate (0,-1);
   SendExposeEvent();
}

void Key9Pressed()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadMatrixd (locscene -> rotmat);
   glRotated ( -1.0, 1.0, 0.0, 0.0 );
   glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   SendExposeEvent();
}

void Key4Pressed()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadMatrixd (locscene -> rotmat);
   glRotated ( -1.0, 0.0, 0.0, 1.0 );
   glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   SendExposeEvent();
}

void Key5Pressed()
{
   if (locscene -> view == 2)
      locscene -> CenterObject2D();
   else
      locscene -> CenterObject();
   SendExposeEvent();
}

void Key6Pressed()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadMatrixd (locscene -> rotmat);
   glRotated ( 1.0, 0.0, 0.0, 1.0 );
   glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   SendExposeEvent();
}

void Key1Pressed()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadMatrixd (locscene -> rotmat);
   glRotated ( 1.0, 1.0, 0.0, 0.0 );
   glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   SendExposeEvent();
}

void Key2Pressed()
{
   locscene -> Rotate (0,1);
   SendExposeEvent();
}

void Key3Pressed()
{
   glMatrixMode (GL_MODELVIEW);
   glLoadMatrixd (locscene -> rotmat);
   glRotated ( 1.0, 0.0, 1.0, 0.0 );
   glGetDoublev (GL_MODELVIEW_MATRIX, locscene -> rotmat);
   SendExposeEvent();
}

void KeyLeftPressed()
{
   locscene -> Rotate (-5,0);
   SendExposeEvent();
}

void KeyRightPressed()
{
   locscene -> Rotate (5,0);
   SendExposeEvent();
}

void KeyUpPressed()
{
   locscene -> Rotate (0,-5);
   SendExposeEvent();
}

void KeyDownPressed()
{
   locscene -> Rotate (0,5);
   SendExposeEvent();
}

void KeyJPressed()
{
   locscene -> OrthogonalProjection = !(locscene -> OrthogonalProjection);
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
   w = (int)round(w / window_scale_factor);
   h = (int)round(h / window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   XResizeWindow(auxXDisplay(), auxXWindow(), w, h);
}

void EnlargeWindow()
{
   GLint viewport[4];

   glGetIntegerv(GL_VIEWPORT, viewport);
   int w = viewport[2], h = viewport[3];
   w = (int)round(w * window_scale_factor);
   h = (int)round(h * window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   XResizeWindow(auxXDisplay(), auxXWindow(), w, h);
}

void MoveResizeWindow(int x, int y, int w, int h)
{
   XMoveResizeWindow(auxXDisplay(), auxXWindow(), x, y, w, h);
}


// Draw a cone of radius 1 with base in the x-y plane and center at (0,0,2)
void Cone()
{
   const int n = 8;
   const double step = 2*M_PI/n;
   const double nz = (1.0/4.0);
   double point = step;
   int i;

   glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
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
   static double eps = 1e-24;
   if (MySetColorLogscale)
   {
      if (val < min)
         val = min;
      if (val > max)
         val = max;
      MySetColor (log(fabs(val/(min+eps))) / (log(fabs(max/(min+eps)))+eps));
   }
   else
      MySetColor ((val-min)/(max-min));
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

   if (val < 0.0) val = 0.0;
   if (val > 1.0) val = 1.0;

   double malpha = MatAlpha;
   if (MatAlphaCenter > 1.0)
      malpha *= exp(-(MatAlphaCenter)*fabs(val-1.0));
   else if (MatAlphaCenter < 0.0)
      malpha *= exp((MatAlphaCenter-1.0)*fabs(val-0.0));
   else
      malpha *= exp(-fabs(val-MatAlphaCenter));

   val *= 0.999999999 * ( RGB_Palette_Size - 1 ) * RepeatPaletteTimes;
   i = (int) floor( val );
   t = val - i;

   if ((i / (RGB_Palette_Size-1)) % 2 == 0)
      pal = RGB_Palette + 3 * ( i % (RGB_Palette_Size-1) );
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
      glColor4f ( r, g, b, malpha );
   else
      glColor3f ( r, g, b );
}

const int Max_Texture_Size = 512;
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
      t *= 0.999999999 * ( RGB_Palette_Size - 1 ) * RepeatPaletteTimes;
      j = (int) floor(t);
      t -= j;

      if ((j / (RGB_Palette_Size-1)) % 2 == 0)
         pal = RGB_Palette + 3 * ( j % (RGB_Palette_Size-1) );
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
      Texture_Size = Max_Texture_Size;
   else
      Texture_Size = RGB_Palette_Size;

   for (int i = 0; i < 3*Texture_Size; i++)
   {
      Texture_Image[i] = RGB_Palette[i];
   }
}

void Set_Texture_Image()
{
   if (UseTexture == 1)
      Make_Texture_From_Palette_2();
   else
      Make_Texture_From_Palette();

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
