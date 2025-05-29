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

#ifndef GLVIS_AUX_VIS_HPP
#define GLVIS_AUX_VIS_HPP

#include "gl/platform_gl.hpp"
#include "gl/types.hpp"

#include "openglvis.hpp"
#include "glwindow.hpp"
#include "font.hpp"

#include <functional>

MainThread& GetMainThread(bool headless = false);
void MainThreadLoop(bool headless = false, bool server_mode = false);

class GLVisCommand;
void SetGLVisCommand(GLVisCommand *cmd);

/// Initializes the visualization and some keys.
GLWindow* InitVisualization(const char name[], int x, int y, int w, int h,
                            bool headless = false);

void SetVisualizationScene(VisualizationScene * scene,
                           int view = 3, const char *keys = NULL);

/// Start the infinite visualization loop.
void RunVisualization();

/// Send expose event. In our case MyReshape is executed and Draw after it.
void SendExposeEvent();

void MyExpose();

void MainLoop();

class SdlWindow;
SdlWindow* GetSdlWindow();
GLWindow* GetAppWindow();
VisualizationScene * GetVisualizationScene();

void SetLegacyGLOnly(bool status);

void AddIdleFunc(void (*Func)(void));
void RemoveIdleFunc(void (*Func)(void));

void LeftButtonDown  (GLWindow::MouseEventInfo *event);
void LeftButtonLoc   (GLWindow::MouseEventInfo *event);
void LeftButtonUp    (GLWindow::MouseEventInfo *event);
void MiddleButtonDown(GLWindow::MouseEventInfo *event);
void MiddleButtonLoc (GLWindow::MouseEventInfo *event);
void MiddleButtonUp  (GLWindow::MouseEventInfo *event);
void RightButtonDown (GLWindow::MouseEventInfo *event);
void RightButtonLoc  (GLWindow::MouseEventInfo *event);
void RightButtonUp   (GLWindow::MouseEventInfo *event);

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
#ifdef GLVIS_USE_LIBPNG
int SaveAsPNG(const char *fname, int w, int h, bool is_hidpi,
              bool with_alpha = false,
              std::function<void(int,void*)> get_row = nullptr);
#endif

/// Send a sequence of keystrokes to the visualization window
void SendKeySequence(const char *seq);

// Directly call the functions assigned to the given keys. Unlike the above
// function, SendKeySequence(), this function does not send X events and
// actually disables the function SendExposeEvent() used by many of the
// functions assigned to keys. Call MyExpose() after calling this function to
// update the visualization window.
void CallKeySequence(const char *seq);


void SetUseTexture(int ut);
int GetUseTexture();
int GetMultisample();
void SetMultisample(int m);

void SetLineWidth(float width);
float GetLineWidth();
void SetLineWidthMS(float width_ms);
float GetLineWidthMS();

void InitFont();
GlVisFont * GetFont();
bool SetFont(const std::vector<std::string>& patterns, int height);
void SetFont(const std::string& fn);

void SetUseHiDPI(bool status);
bool GetUseHiDPI();

std::function<std::string(double)> NumberFormatter(int precision=4,
                                                   char format='d', bool showsign=false);
std::function<std::string(double)> NumberFormatter(std::string formatting);
bool isValidNumberFormatting(const std::string& formatting);

// This is a helper function for prompting the user for inputs. The benefit
// over using just `cin >> input` is that you can specify a type and optionally
// a validator lambda. The a validator if not specified, it defaults to the
// True function. If the input cannot be type casted to the expected type, or
// if it fails the validation, the user is asked again for a new input.
template <typename T>
T prompt(const std::string question,
         const T* default_value = nullptr,
std::function<bool(T)> validator = [](T) { return true; })
{
   T input;
   std::string strInput;

   while (true)
   {
      std::cout << question << " ";
      std::getline(std::cin, strInput);
      std::stringstream buf(strInput);

      if (strInput.empty() && default_value != nullptr)
      {
         std::cout << "Input empty. Using default value: " << *default_value
                   << std::endl;
         return *default_value;
      }

      if (buf >> input)
      {
         if (validator(input))
         {
            return input;
         }
         else
         {
            std::cout << "Input is not valid. Please try again." << std::endl;
         }
      }
      else
      {
         std::cout << "Input can not be casted to expected type. Please try again."
                   << std::endl;
      }
   }
   return input;
}

#endif
