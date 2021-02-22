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

#ifndef GLVIS_AUX_VIS_HPP
#define GLVIS_AUX_VIS_HPP

#include "gl/platform_gl.hpp"
#include "gl/types.hpp"

#include "sdl.hpp"
#include "font.hpp"
#include "openglvis.hpp"

extern GLuint fontbase;
extern float MatAlpha;
extern float MatAlphaCenter;
extern int RepeatPaletteTimes;
extern int PaletteNumColors;

/// Initializes the visualization and some keys.
int InitVisualization(const char name[], int x, int y, int w, int h);

void SetVisualizationScene(VisualizationScene * scene,
                           int view = 3, const char *keys = NULL);

/// Start the infinite visualization loop.
void RunVisualization();

/// Send expose event. In our case MyReshape is executed and Draw after it.
void SendExposeEvent();

void MyExpose();

void MainLoop();

SdlWindow * GetAppWindow();
VisualizationScene * GetVisualizationScene();

void SetLegacyGLOnly(bool status);

void AddIdleFunc(void (*Func)(void));
void RemoveIdleFunc(void (*Func)(void));

void LeftButtonDown  (EventInfo *event);
void LeftButtonLoc   (EventInfo *event);
void LeftButtonUp    (EventInfo *event);
void MiddleButtonDown(EventInfo *event);
void MiddleButtonLoc (EventInfo *event);
void MiddleButtonUp  (EventInfo *event);
void RightButtonDown (EventInfo *event);
void RightButtonLoc  (EventInfo *event);
void RightButtonUp   (EventInfo *event);

void TouchPinch(SDL_MultiGestureEvent & e);

void KeyCtrlP();
void KeyS();
void KeyQPressed();
void ToggleThreads();
void ThreadsPauseFunc(GLenum);
void ThreadsStop();
void ThreadsRun();

void Key1Pressed();
void Key2Pressed();
void Key3Pressed();
void Key4Pressed();
void Key5Pressed();
void Key6Pressed();
void Key7Pressed();
void Key8Pressed();
void Key9Pressed();

void Key0Pressed();
void KeyDeletePressed();
void KeyEnterPressed();

void KeyLeftPressed(GLenum);
void KeyRightPressed(GLenum);
void KeyUpPressed(GLenum);
void KeyDownPressed(GLenum);
void KeyJPressed();
void KeyMinusPressed();
void KeyPlusPressed();

void ZoomIn();
void ZoomOut();
void ScaleUp();
void ScaleDown();
void LookAt();
void ShrinkWindow();
void EnlargeWindow();
void MoveResizeWindow(int x, int y, int w, int h);
void ResizeWindow(int w, int h);
void SetWindowTitle(const char *title);

/// Take a screenshot using libtiff, libpng or sdl2
int Screenshot(const char *fname, bool convert = false);

/// Send a sequence of keystrokes to the visualization window
void SendKeySequence(const char *seq);

// Directly call the functions assigned to the given keys. Unlike the above
// function, SendKeySequence(), this function does not send X events and
// actually disables the function SendExposeEvent() used by many of the
// functions assigned to keys. Call MyExpose() after calling this function to
// update the visualization window.
void CallKeySequence(const char *seq);

extern int MySetColorLogscale;
double GetColorCoord(double val, double min, double max);
void GetColorFromVal(double val, float * rgba);
void MySetColor(gl3::GlBuilder& builder, double val);
void MySetColor(gl3::GlBuilder& builder, double val, double min, double max);

void SetUseTexture(int ut);
int GetUseTexture();
void Set_Texture_Image();
int GetMultisample();
void SetMultisample(int m);

void SetLineWidth(float width);
float GetLineWidth();
void SetLineWidthMS(float width_ms);
float GetLineWidthMS();

void InitFont();
GlVisFont * GetFont();
bool SetFont(const vector<std::string>& patterns, int height);
void SetFont(const std::string& fn);

#endif
