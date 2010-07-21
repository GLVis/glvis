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

#ifndef GLVIS_AUX_VIS
#define GLVIS_AUX_VIS

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

#include "tk.h"
#include "aux_gl.hpp"

#include "openglvis.hpp"

extern float MatAlpha;
extern float MatAlphaCenter;

/// Initializes the visualization and some keys.
void InitVisualization (const char name[],
                        int x, int y, int w, int h);

/// Start the infinite visualization loop.
void SetVisualizationScene (VisualizationScene * scene,
                            int view = 3, char * keys = NULL);

void KillVisualization();

/// Send expose event. In our case MyReshape is executed and Draw after it.
void SendExposeEvent();

void MyExpose();

void MainLoop ();

void LeftButtonDown  (AUX_EVENTREC *event);
void LeftButtonLoc   (AUX_EVENTREC *event);
void LeftButtonUp    (AUX_EVENTREC *event);
void MiddleButtonDown(AUX_EVENTREC *event);
void MiddleButtonLoc (AUX_EVENTREC *event);
void MiddleButtonUp  (AUX_EVENTREC *event);
void RightButtonDown (AUX_EVENTREC *event);
void RightButtonLoc  (AUX_EVENTREC *event);
void RightButtonUp   (AUX_EVENTREC *event);

void KeyP ();
void KeyS ();
void KeyQPressed ();

void Key1Pressed ();
void Key2Pressed ();
void Key3Pressed ();
void Key4Pressed ();
void Key5Pressed ();
void Key6Pressed ();
void Key7Pressed ();
void Key8Pressed ();
void Key9Pressed ();

void Key0Pressed ();
void KeyDeletePressed ();
void KeyEnterPressed ();

void KeyLeftPressed ();
void KeyRightPressed();
void KeyUpPressed   ();
void KeyDownPressed ();
void KeyJPressed    ();
void KeyMinusPressed();
void KeyPlusPressed ();

void ZoomIn ();
void ZoomOut ();
void ScaleUp();
void ScaleDown();
void LookAt();
void ShrinkWindow();
void EnlargeWindow();
void MoveResizeWindow(int x, int y, int w, int h);

/// Take a screenshot using libtiff or xwd
int Screenshot(const char *fname);

/// Send a sequence of keystrokes to the visualization window
void SendKeySequence (char * seq);

void Cone();

extern int MySetColorLogscale;
void MySetColor (double val, double min, double max);
void MySetColor (double val);
void SetUseTexture(int ut);
int GetUseTexture();

#endif
