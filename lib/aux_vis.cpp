// Copyright (c) 2010-2021, Lawrence Livermore National Security, LLC. Produced
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
#include <fstream>
#include <cmath>
#include <ctime>
#include <exception>

#include "mfem.hpp"
using namespace mfem;
#include "sdl.hpp"
#include "palettes.hpp"
#include "visual.hpp"
#include "gl2ps.h"

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

#ifdef GLVIS_MULTISAMPLE
static int glvis_multisample = GLVIS_MULTISAMPLE;
#else
static int glvis_multisample = -1;
#endif

float line_w = 1.f;
float line_w_aa = gl3::LINE_WIDTH_AA;

[[deprecated]] SdlWindow* wnd;
[[deprecated]] GLVisWindow * glvis_wnd = nullptr;

SdlWindow * GetAppWindow()
{
    return wnd;
}

GLVisWindow * GetGLVisWindow()
{
    return glvis_wnd;
}

VisualizationScene * GetVisualizationScene()
{
   return locscene;
}

struct GLVisWindow::RotationControl
{
    GLVisWindow* wnd;
    double xang = 0., yang = 0.;
    gl3::GlMatrix srot;
    double sph_t, sph_u;
    GLint oldx, oldy, startx, starty;

    bool constrained_spinning = 0;

    void LeftButtonDown  (EventInfo *event);
    void LeftButtonLoc   (EventInfo *event);
    void LeftButtonUp    (EventInfo *event);
    void MiddleButtonDown(EventInfo *event);
    void MiddleButtonLoc (EventInfo *event);
    void MiddleButtonUp  (EventInfo *event);
    void RightButtonDown (EventInfo *event);
    void RightButtonLoc  (EventInfo *event);
    void RightButtonUp   (EventInfo *event);

    void CheckSpin();
    void Key0Pressed();
    void KeyDeletePressed();
    void KeyEnterPressed();
    MouseDelegate CreateMouseEvent(void (GLVisWindow::RotationControl::*func)(EventInfo*))
    {
        return [this, func](EventInfo* ei) { (this->*func)(ei); };
    }

};

template<typename T>
KeyDelegate CreateKeyEvent(T* inst, void (T::*func)())
{
    return [inst, func](GLenum) { (inst->*func)(); };
}


GLVisWindow::GLVisWindow(std::string name, int x, int y, int w, int h, bool legacyGlOnly)
    : locscene(nullptr)
    , idle_funcs(0)
    , rot_data(new RotationControl)
{

#ifdef GLVIS_DEBUG
   cout << "OpenGL Visualization" << endl;
#endif
   rot_data->wnd = this;
   ::glvis_wnd = this;
   wnd.reset(new SdlWindow());
   if (!wnd->createWindow(name, x, y, w, h, legacyGlOnly))
   {
      throw std::runtime_error("Could not create an SDL window.");
   }
   ::wnd = wnd.get();
#ifdef GLVIS_DEBUG
   cout << "Window should be up" << endl;
#endif
   InitFont();
   wnd->getRenderer().setFont(&font);
   wnd->getRenderer().setLineWidth(line_w);
   wnd->getRenderer().setLineWidthMS(line_w_aa);

   // auxReshapeFunc (MyReshape); // not needed, MyExpose calls it
   // auxReshapeFunc (NULL);
   wnd->setOnExpose([this]() { MyExpose(); });
   auto LeftButtonDown   = rot_data->CreateMouseEvent(&RotationControl::LeftButtonDown);
   auto LeftButtonUp     = rot_data->CreateMouseEvent(&RotationControl::LeftButtonUp);
   auto LeftButtonLoc    = rot_data->CreateMouseEvent(&RotationControl::LeftButtonLoc);
   auto MiddleButtonDown = rot_data->CreateMouseEvent(&RotationControl::MiddleButtonDown);
   auto MiddleButtonUp   = rot_data->CreateMouseEvent(&RotationControl::MiddleButtonUp);
   auto MiddleButtonLoc  = rot_data->CreateMouseEvent(&RotationControl::MiddleButtonLoc);
   auto RightButtonDown  = rot_data->CreateMouseEvent(&RotationControl::RightButtonDown);
   auto RightButtonUp    = rot_data->CreateMouseEvent(&RotationControl::RightButtonUp);
   auto RightButtonLoc   = rot_data->CreateMouseEvent(&RotationControl::RightButtonLoc);

   auto Key0Pressed = CreateKeyEvent(rot_data.get(), &RotationControl::Key0Pressed);
   auto KeyEnterPressed = CreateKeyEvent(rot_data.get(), &RotationControl::KeyEnterPressed);
   auto KeyDeletePressed = CreateKeyEvent(rot_data.get(), &RotationControl::KeyDeletePressed);

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
   SetKeyEventHandler('A', &GLVisWindow::ToggleAntialiasing);

   // auxKeyFunc (AUX_p, KeyCtrlP); // handled in vsdata.cpp
   SetKeyEventHandler (SDLK_s, &GLVisWindow::Screenshot);
   SetKeyEventHandler ('S', &GLVisWindow::Screenshot);

   SetKeyEventHandler (SDLK_q, &GLVisWindow::Quit);
   // wnd->setOnKeyDown (SDLK_Q, KeyQPressed);

   SetKeyEventHandler (SDLK_LEFT, &GLVisWindow::KeyLeftPressed);
   SetKeyEventHandler (SDLK_RIGHT, &GLVisWindow::KeyRightPressed);
   SetKeyEventHandler (SDLK_UP, &GLVisWindow::KeyUpPressed);
   SetKeyEventHandler (SDLK_DOWN, &GLVisWindow::KeyDownPressed);

   wnd->setOnKeyDown (SDLK_KP_0, Key0Pressed);
   SetKeyEventHandler (SDLK_KP_1, &GLVisWindow::Key1Pressed);
   SetKeyEventHandler (SDLK_KP_2, &GLVisWindow::Key2Pressed);
   SetKeyEventHandler (SDLK_KP_3, &GLVisWindow::Key3Pressed);
   SetKeyEventHandler (SDLK_KP_4, &GLVisWindow::Key4Pressed);
   SetKeyEventHandler (SDLK_KP_5, &GLVisWindow::Key5Pressed);
   SetKeyEventHandler (SDLK_KP_6, &GLVisWindow::Key6Pressed);
   SetKeyEventHandler (SDLK_KP_7, &GLVisWindow::Key7Pressed);
   SetKeyEventHandler (SDLK_KP_8, &GLVisWindow::Key8Pressed);
   SetKeyEventHandler (SDLK_KP_9, &GLVisWindow::Key9Pressed);

   SetKeyEventHandler (SDLK_KP_MEMSUBTRACT, &GLVisWindow::KeyMinusPressed);
   SetKeyEventHandler (SDLK_KP_MEMADD, &GLVisWindow::KeyPlusPressed);

   wnd->setOnKeyDown (SDLK_KP_DECIMAL, KeyDeletePressed);
   wnd->setOnKeyDown (SDLK_KP_ENTER, KeyEnterPressed);

   wnd->setOnKeyDown (SDLK_PERIOD, KeyDeletePressed);
   wnd->setOnKeyDown (SDLK_RETURN, KeyEnterPressed);

   wnd->setOnKeyDown (SDLK_0, Key0Pressed);
   SetKeyEventHandler (SDLK_1, &GLVisWindow::Key1Pressed);
   SetKeyEventHandler (SDLK_2, &GLVisWindow::Key2Pressed);
   SetKeyEventHandler (SDLK_3, &GLVisWindow::Key3Pressed);
   SetKeyEventHandler (SDLK_4, &GLVisWindow::Key4Pressed);
   SetKeyEventHandler (SDLK_5, &GLVisWindow::Key5Pressed);
   SetKeyEventHandler (SDLK_6, &GLVisWindow::Key6Pressed);
   SetKeyEventHandler (SDLK_7, &GLVisWindow::Key7Pressed);
   SetKeyEventHandler (SDLK_8, &GLVisWindow::Key8Pressed);
   SetKeyEventHandler (SDLK_9, &GLVisWindow::Key9Pressed);

   SetKeyEventHandler (SDLK_MINUS, &GLVisWindow::KeyMinusPressed);
   SetKeyEventHandler (SDLK_PLUS, &GLVisWindow::KeyPlusPressed);
   SetKeyEventHandler (SDLK_EQUALS, &GLVisWindow::KeyPlusPressed);

   SetKeyEventHandler (SDLK_j, &GLVisWindow::KeyJPressed);
   // wnd->setOnKeyDown (AUX_J, KeyJPressed);

   SetKeyEventHandler (SDLK_KP_MULTIPLY, &GLVisWindow::ZoomIn);
   SetKeyEventHandler (SDLK_KP_DIVIDE, &GLVisWindow::ZoomOut);

   SetKeyEventHandler (SDLK_ASTERISK, &GLVisWindow::ZoomIn);
   SetKeyEventHandler (SDLK_SLASH, &GLVisWindow::ZoomOut);

   SetKeyEventHandler (SDLK_LEFTBRACKET, &GLVisWindow::ScaleDown);
   SetKeyEventHandler (SDLK_RIGHTBRACKET, &GLVisWindow::ScaleUp);
   SetKeyEventHandler (SDLK_AT, &GLVisWindow::LookAt);

   SetKeyEventHandler (SDLK_SPACE, &GLVisWindow::ThreadsPauseFunc);

#ifndef __EMSCRIPTEN__
   SetKeyEventHandler(SDLK_LEFTPAREN, &GLVisWindow::ShrinkWindow);
   SetKeyEventHandler(SDLK_RIGHTPAREN, &GLVisWindow::EnlargeWindow);
#endif
}

GLVisWindow::~GLVisWindow()
{
    if (glvis_command)
    {
        glvis_command->Terminate();
    }
}

void GLVisWindow::InitVisualization(int field_type, StreamState state,
                                    const mfem::Array<istream*>& input_streams,
                                    bool& keep_attr)
{
   prob_state = std::move(state);
   if (input_streams.Size() > 0)
   {
      glvis_command.reset(new GLVisCommand(this, prob_state, &keep_attr));
      comm_thread.reset(new communication_thread(glvis_command.get(), input_streams));
   }

   locscene = prob_state.CreateVisualizationScene(field_type);

   if (prob_state.mesh->SpaceDimension() == 2 && field_type == 2)
   {
      locscene->view = 2;
      locscene->CenterObject2D();
   }
   else
   {
      locscene->view = 3;
      locscene->CenterObject();
   }

   if (locscene->spinning)
   {
      AddIdleFunc(::MainLoop);
   }
}

void GLVisWindow::SetKeyEventHandler(int key, void (GLVisWindow::*handler)())
{
    auto handlerWrapper = [this, handler]() { (this->*handler)(); };
    wnd->setOnKeyDown(key, handlerWrapper);
}

void GLVisWindow::SetKeyEventHandler(int key, void (GLVisWindow::*handler)(GLenum))
{
    auto handlerWrapper = [this, handler](GLenum mod) { (this->*handler)(mod); };
    wnd->setOnKeyDown(key, handlerWrapper);
}

void SendKeySequence(const char *seq)
{
    glvis_wnd->SendKeySequence(seq);
}

void CallKeySequence(const char *seq)
{
    glvis_wnd->CallKeySequence(seq);
}

void GLVisWindow::SendKeySequence(const char *seq)
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


void GLVisWindow::CallKeySequence(const char *seq)
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

void GLVisWindow::RunVisualization()
{
   visualize = 1;
   wnd->setOnIdle([this](){return MainIdleFunc();});
#ifndef __EMSCRIPTEN__
   wnd->mainLoop();
#endif
}

void SendExposeEvent()
{
    glvis_wnd->SendExposeEvent();
}

void MyExpose()
{
    glvis_wnd->MyExpose();
}

void GLVisWindow::SendExposeEvent()
{
   if (disableSendExposeEvent) { return; }
   wnd->signalExpose();
}

void GLVisWindow::MyReshape(GLsizei w, GLsizei h)
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

void GLVisWindow::MyExpose(GLsizei w, GLsizei h)
{
   MyReshape (w, h);
   GLuint color_tex = locscene->palette.GetColorTexture();
   GLuint alpha_tex = locscene->palette.GetAlphaTexture();
   wnd->getRenderer().setColorTexture(color_tex);
   wnd->getRenderer().setAlphaTexture(alpha_tex);

   std::array<float, 4> bgcol = locscene->GetBackgroundColor();
   wnd->getRenderer().setClearColor(bgcol[0], bgcol[1], bgcol[2], bgcol[3]);

   gl3::SceneInfo frame = locscene->GetSceneObjs();
   for (auto drawable_ptr : frame.needs_buffering)
   {
      wnd->getRenderer().buffer(drawable_ptr);
   }
   wnd->getRenderer().render(frame.queue);
}

void GLVisWindow::MyExpose()
{
   int w, h;
   wnd->getGLDrawSize(w, h);
   MyExpose(w, h);
   wnd->signalSwap();
}

bool GLVisWindow::CommunicationIdleFunc()
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

bool GLVisWindow::MainIdleFunc()
{
   bool sleep = true;
   if (glvis_command && visualize == 1
       && !(idle_funcs.Size() > 0 && use_idle))
   {
       // Execute the next event from the communication thread if:
       //  - a valid GLVisCommand has been set
       //  - the communication thread is not stopped
       //  - The idle function flag is not set, or no idle functions have been
       //    registered
       sleep = CommunicationIdleFunc();
       if (idle_funcs.Size() > 0) { sleep = false; }
   }
   else if (idle_funcs.Size() > 0)
   {
       last_idle_func = (last_idle_func + 1) % idle_funcs.Size();
       if (idle_funcs[last_idle_func])
       {
           (*idle_funcs[last_idle_func])(this);
       }
       // Continue executing idle functions
       sleep = false;
   }
   use_idle = !use_idle;
   return sleep;
}

void GLVisWindow::AddIdleFunc(GLVisWindow::IdleFPtr Func)
{
   idle_funcs.Union(Func);
   use_idle = false;
   wnd->setOnIdle([this](){return MainIdleFunc();});
}

void GLVisWindow::RemoveIdleFunc(GLVisWindow::IdleFPtr Func)
{
   idle_funcs.DeleteFirst(Func);
   if (idle_funcs.Size() == 0)
   {
      use_idle = false;
      if (!glvis_command) { wnd->setOnIdle(nullptr); }
   }
}

void MainLoop(GLVisWindow* wnd)
{
    wnd->MainLoop();
}

void GLVisWindow::MainLoop()
{
   static int p = 1;
   struct timespec req;
   double xang = rot_data->xang;
   double yang = rot_data->yang;
   if (locscene->spinning)
   {
      if (!rot_data->constrained_spinning)
      {
         locscene->Rotate(xang, yang);
         SendExposeEvent();
      }
      else
      {
         locscene->PreRotate(xang, 0.0, 0.0, 1.0);
         SendExposeEvent();
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

inline void ComputeSphereAngles(int viewport_w, int viewport_h,
                                int &newx, int &newy,
                                double &new_sph_u, double &new_sph_t)
{
   GLint viewport[4] = { 0, 0, viewport_w, viewport_h };
   double r, x, y, rr;
   const double maxr = 0.996194698091745532295;

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

void GLVisWindow::RotationControl::LeftButtonDown (EventInfo *event)
{
   wnd->locscene -> spinning = 0;
   wnd->RemoveIdleFunc(::MainLoop);

   oldx = event->mouse_x;
   oldy = event->mouse_y;

   int vp_w, vp_h;
   wnd->wnd->getGLDrawSize(vp_w, vp_h);

   ComputeSphereAngles(vp_w, vp_h, oldx, oldy, sph_u, sph_t);

   srot.identity();
   srot.mult(wnd->locscene->cam.RotMatrix());
   srot.mult(wnd->locscene->rotmat);

   startx = oldx;
   starty = oldy;
}

void GLVisWindow::RotationControl::LeftButtonLoc (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;

   if (event->keymod & KMOD_CTRL)
   {
      if (event->keymod & KMOD_SHIFT)
      {
         wnd->locscene->PreRotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
      }
      else
      {
         double new_sph_u, new_sph_t;

         int vp_w, vp_h;
         wnd->wnd->getGLDrawSize(vp_w, vp_h);

         ComputeSphereAngles(vp_w, vp_h, newx, newy, new_sph_u, new_sph_t);

         gl3::GlMatrix newrot;
         newrot.identity();

         double scoord[3], ncoord[3], inprod, cross[3];
         scoord[0] = scoord[1] = cos(sph_u);     scoord[2] = sin(sph_u);
         ncoord[0] = ncoord[1] = cos(new_sph_u); ncoord[2] = sin(new_sph_u);
         scoord[0] *= cos(sph_t);     scoord[1] *= sin(sph_t);
         ncoord[0] *= cos(new_sph_t); ncoord[1] *= sin(new_sph_t);
         inprod = InnerProd(scoord, ncoord);
         CrossProd(scoord, ncoord, cross);

         newrot.mult(wnd->locscene->cam.TransposeRotMatrix());
         newrot.rotate(acos(inprod)*(180.0/M_PI), cross[0], cross[1], cross[2]);
         newrot.mult(srot.mtx);
         wnd->locscene->rotmat = newrot.mtx;
      }
   }
   else if (event->keymod & KMOD_ALT)
   {
      wnd->locscene->Rotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
   }
   else if (event->keymod & KMOD_SHIFT)
   {
      wnd->locscene->Rotate(double(newx-oldx)/2, double(newy-oldy)/2);
   }
   else
   {
      wnd->locscene->Rotate(double(newy-oldy)/2, 1.0, 0.0, 0.0);
      wnd->locscene->PreRotate(double(newx-oldx)/2, 0.0, 0.0, 1.0);
   }

   oldx = newx;
   oldy = newy;

   wnd->SendExposeEvent();
}

void GLVisWindow::RotationControl::LeftButtonUp (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;

   xang = (newx-startx)/5.0;
   yang = (newy-starty)/5.0;

   if ( (event->keymod & KMOD_SHIFT) && (xang != 0.0 || yang != 0.0) )
   {
      wnd->locscene -> spinning = 1;
      wnd->AddIdleFunc(::MainLoop);
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

void GLVisWindow::RotationControl::MiddleButtonDown (EventInfo *event)
{
   startx = oldx = event->mouse_x;
   starty = oldy = event->mouse_y;
}

void GLVisWindow::RotationControl::MiddleButtonLoc (EventInfo *event)
{
   GLint newx = event->mouse_x;
   GLint newy = event->mouse_y;

   if ( !( event->keymod & KMOD_CTRL ) )
   {
      int w, h;
      double TrX, TrY, scale;

      if (wnd->locscene->OrthogonalProjection)
      {
         scale = wnd->locscene->ViewScale;
      }
      else
      {
         scale = 0.4142135623730950488/tan(wnd->locscene->ViewAngle*(M_PI/360));
      }
      wnd->wnd->getGLDrawSize(w, h);
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
      wnd->locscene->ViewCenterX += TrX;
      wnd->locscene->ViewCenterY += TrY;
   }
   else
   {
      // wnd->locscene->Translate((double)(newx-oldx)/200,(double)(newy-oldy)/200);

      double dx = double(newx-oldx)/400;
      double dy = double(oldy-newy)/400;

      if (event->keymod & KMOD_SHIFT)  // ctrl + shift
      {
         double sx = double(newx-startx)/400;
         double sy = double(starty-newy)/400;

         wnd->locscene->cam.TurnLeftRight(dx-sx);
         wnd->locscene->cam.TurnUpDown(sy-dy);

         wnd->locscene->cam.TurnUpDown(-sy);
         wnd->locscene->cam.TurnLeftRight(sx);
      }
      else if (event->keymod & KMOD_ALT) // ctrl + alt
      {
         wnd->locscene->cam.MoveForwardBackward(dy);
         wnd->locscene->cam.TiltLeftRight(-dx);
      }
      else // ctrl
      {
         wnd->locscene->cam.MoveLeftRight(dx);
         wnd->locscene->cam.MoveUpDown(-dy);
      }
   }

   wnd->SendExposeEvent();

   oldx = newx;
   oldy = newy;
}

void GLVisWindow::RotationControl::MiddleButtonUp (EventInfo*)
{}

void GLVisWindow::RotationControl::RightButtonDown (EventInfo *event)
{
   startx = oldx = event->mouse_x;
   starty = oldy = event->mouse_y;
}

void GLVisWindow::RotationControl::RightButtonLoc (EventInfo *event)
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
      wnd->locscene->SetLight0CustomPos({(float)x, (float)y, (float)z, 0.f});
   }
   else if ( !( event->keymod & KMOD_CTRL ) )
   {
      wnd->locscene -> Zoom (exp ( double (oldy-newy) / 100 ));
   }
   else
   {
      wnd->locscene -> Scale ( exp ( double (oldy-newy) / 50 ) );
   }

   wnd->SendExposeEvent();

   oldx = newx;
   oldy = newy;
}

void GLVisWindow::RotationControl::RightButtonUp (EventInfo*)
{}

void GLVisWindow::ToggleAntialiasing()
{
   bool curr_aa = wnd->getRenderer().getAntialiasing();
   wnd->getRenderer().setAntialiasing(!curr_aa);
   const std::string strings_off_on[2] = { "off", "on" };

   cout << "Multisampling/Antialiasing: "
        << strings_off_on[!curr_aa ? 1 : 0] << endl;

   SendExposeEvent();
}

void TouchPinch(SDL_MultiGestureEvent & e)
{
   // Scale or Zoom?
   locscene->Zoom(exp(e.dDist*10));
   SendExposeEvent();
}

void GLVisWindow::Screenshot(std::string filename)
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
      if (filename == "")
      {
         char fname[20];
         snprintf(fname, 20, "GLVis_s%02d", p++);
         wnd->screenshot(fname);
      }
      else 
      {
          wnd->screenshot(filename);
      }
      cout << "done" << endl;
   }
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
   //print lines
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

void GLVisWindow::PrintToPDF()
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
                       //GL2PS_DRAW_BACKGROUND |
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

void GLVisWindow::Quit()
{
   wnd->signalQuit();
   visualize = 0;
}

void GLVisWindow::ToggleThreads()
{
   static const char *state[] = { "running", "stopped" };
   if (visualize > 0 && visualize < 3)
   {
      visualize = 3 - visualize; //  1 <-> 2
      cout << "Communication thread(s): " << state[visualize-1] << endl;
   }
}

void GLVisWindow::ThreadsPauseFunc(GLenum state)
{
#ifndef __EMSCRIPTEN__
   if (glvis_command)
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
#endif
}

void GLVisWindow::ThreadsStop()
{
   if (visualize == 1)
   {
      ToggleThreads();
   }
}

void GLVisWindow::ThreadsRun()
{
   if (visualize == 2)
   {
      ToggleThreads();
   }
}

void GLVisWindow::RotationControl::CheckSpin()
{
   if (fabs(xang) < 1.e-2)
   {
      xang = 0.;
   }
   if (xang != 0. || yang != 0.)
   {
      wnd->locscene->spinning = 1;
      wnd->AddIdleFunc(::MainLoop);
   }
   else
   {
      wnd->locscene->spinning = 0;
      wnd->RemoveIdleFunc(::MainLoop);
   }
   cout << "Spin angle: " << xang << " degrees / frame" << endl;
}

const double xang_step = 0.2; // angle in degrees

void GLVisWindow::RotationControl::Key0Pressed()
{
   if (!wnd->locscene -> spinning)
   {
      xang = 0;
   }
   xang -= xang_step;
   CheckSpin();
}

void GLVisWindow::RotationControl::KeyDeletePressed()
{
   if (wnd->locscene -> spinning)
   {
      xang = yang = 0.;
      wnd->locscene -> spinning = 0;
      wnd->RemoveIdleFunc(::MainLoop);
      constrained_spinning = 1;
   }
   else
   {
      xang = xang_step;
      wnd->locscene -> spinning = 1;
      wnd->AddIdleFunc(::MainLoop);
      constrained_spinning = 1;
   }
}

void GLVisWindow::RotationControl::KeyEnterPressed()
{
   if (!wnd->locscene -> spinning)
   {
      xang = 0;
   }
   xang += xang_step;
   CheckSpin();
}

void GLVisWindow::Key7Pressed()
{
   locscene->PreRotate(1.0, 0.0, -1.0, 0.0);
   SendExposeEvent();
}

void GLVisWindow::Key8Pressed()
{
   locscene->Rotate(0.0, -1.0);
   SendExposeEvent();
}

void GLVisWindow::Key9Pressed()
{
   locscene->PreRotate(-1.0, 1.0, 0.0, 0.0);
   SendExposeEvent();
}

void GLVisWindow::Key4Pressed()
{
   locscene->PreRotate(-1.0, 0.0, 0.0, 1.0);
   SendExposeEvent();
}

void GLVisWindow::Key5Pressed()
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

void GLVisWindow::Key6Pressed()
{
   locscene->PreRotate(1.0, 0.0, 0.0, 1.0);
   SendExposeEvent();
}

void GLVisWindow::Key1Pressed()
{
   locscene->PreRotate(1.0, 1.0, 0.0, 0.0);
   SendExposeEvent();
}

void GLVisWindow::Key2Pressed()
{
   locscene->Rotate(1.0, 1.0, 0.0, 0.0);
   SendExposeEvent();
}

void GLVisWindow::Key3Pressed()
{
   locscene->PreRotate(1.0, 0.0, 1.0, 0.0);
   SendExposeEvent();
}

void ShiftView(VisualizationScene* locscene, double dx, double dy)
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

void GLVisWindow::KeyLeftPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(locscene.get(), 0.05, 0.);
   }
   else
   {
      locscene->Rotate(-5, 0.0, 1.0, 0.0);
   }
   SendExposeEvent();
}

void GLVisWindow::KeyRightPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(locscene.get(), -0.05, 0.);
   }
   else
   {
      locscene->Rotate(5, 0.0, 1.0, 0.0);
   }
   SendExposeEvent();
}

void GLVisWindow::KeyUpPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(locscene.get(), 0., -0.05);
   }
   else
   {
      locscene->Rotate(-5, 1.0, 0.0, 0.0);
   }
   SendExposeEvent();
}

void GLVisWindow::KeyDownPressed(GLenum state)
{
   if (state & KMOD_CTRL)
   {
      ShiftView(locscene.get(), 0., 0.05);
   }
   else
   {
      locscene->Rotate(5, 1.0, 0.0, 0.0);
   }
   SendExposeEvent();
}

void GLVisWindow::KeyJPressed()
{
   locscene->ToggleProjectionMode();
   SendExposeEvent();
}

void GLVisWindow::KeyMinusPressed()
{
   locscene -> Scale(1., 1., 1./1.1);
   SendExposeEvent();
}

void GLVisWindow::KeyPlusPressed()
{
   locscene -> Scale(1., 1., 1.1);
   SendExposeEvent();
}

void GLVisWindow::ZoomIn()
{
   locscene->Zoom(exp (0.05));
   SendExposeEvent();
}

void GLVisWindow::ZoomOut()
{
   locscene->Zoom(exp (-0.05));
   SendExposeEvent();
}

void GLVisWindow::ScaleUp()
{
   locscene->Scale(1.025);
   SendExposeEvent();
}

void GLVisWindow::ScaleDown()
{
   locscene->Scale(1.0/1.025);
   SendExposeEvent();
}

void GLVisWindow::LookAt()
{
   cout << "ViewCenter = (" << locscene->ViewCenterX << ','
        << locscene->ViewCenterY << ")\nNew x = " << flush;
   cin >> locscene->ViewCenterX;
   cout << "New y = " << flush;
   cin >> locscene->ViewCenterY;
   SendExposeEvent();
}

const double window_scale_factor = 1.1;

void GLVisWindow::ShrinkWindow()
{
   int w, h;
   wnd->getWindowSize(w, h);
   w = (int)ceil(w / window_scale_factor);
   h = (int)ceil(h / window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   ResizeWindow(w, h);
}

void GLVisWindow::EnlargeWindow()
{
   int w, h;
   wnd->getWindowSize(w, h);
   w = (int)ceil(w * window_scale_factor);
   h = (int)ceil(h * window_scale_factor);

   cout << "New window size : " << w << " x " << h << endl;

   ResizeWindow(w, h);
}

void GLVisWindow::MoveResizeWindow(int x, int y, int w, int h)
{
   wnd->setWindowSize(w, h);
   wnd->setWindowPos(x, y);
}

void GLVisWindow::ResizeWindow(int w, int h)
{
   wnd->setWindowSize(w, h);
}

void GLVisWindow::SetWindowTitle(const char *title)
{
   wnd->setWindowTitle(title);
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
};

//constexpr int default_font_size = 12;
//int font_size = default_font_size;

//GlVisFont glvis_font;
//std::string priority_font;

void GLVisWindow::InitFont()
{
   // This function is called after the window is created.
   GLenum alphaChannel =
      gl3::GLDevice::useLegacyTextureFmts() ? GL_ALPHA : GL_RED;
   int ppi_w, ppi_h;
   wnd->getDpi(ppi_w, ppi_h);
   bool is_hidpi = wnd->isHighDpi();
   font.SetDPIParams(is_hidpi, ppi_w, ppi_h);
   font.setAlphaChannel(alphaChannel);
   bool try_fc_patterns = true;
   if (!priority_font.empty())
   {
      if (SetFont({priority_font}, font_size) ||
          font.LoadFont(priority_font, 0, font_size))
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
   wnd->getRenderer().setFontTexture(font.getFontTex());
}

bool GLVisWindow::SetFont(const vector<std::string>& font_patterns, int height)
{
#ifdef __EMSCRIPTEN__
   return font.LoadFont("OpenSans.ttf", 0, height);
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
         if (font.LoadFont(font_file, font_index, height))
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

   return font.isFontLoaded();
#endif
}

void GLVisWindow::SetFont(const std::string& fn)
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
