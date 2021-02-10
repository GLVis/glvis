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

#ifndef GLVIS_AUX_VIS_HPP
#define GLVIS_AUX_VIS_HPP

#include "mfem.hpp"

#include "gl/platform_gl.hpp"
#include "gl/types.hpp"

#include "sdl.hpp"
#include "font.hpp"
#include "openglvis.hpp"
#include "stream_reader.hpp"

class GLVisCommand;
class communication_thread;

class GLVisWindow
{
public:
    using IdleFPtr = void(*)(GLVisWindow* wnd);
    using KeyEvent = void(*)(GLVisWindow* wnd, int keystate);

    /// Initializes the visualization and some keys.
    GLVisWindow(std::string name, int x, int y, int w, int h, bool legacyGlOnly);

    ~GLVisWindow();

    void InitVisualization(int field_type, StreamState state,
                           bool& keep_attr,
                           const mfem::Array<istream*>& input_streams = {});

    void SetFont(const std::string& fn);

    VisualizationScene* getScene() { return locscene.get(); }

    GlVisFont* getFont() { return &font; }

    SdlWindow* getSdl() { return wnd.get(); }

    /// Start the infinite visualization loop.
    void RunVisualization();

    /// Send expose event. In our case MyReshape is executed and Draw after it.
    void SendExposeEvent();

    void MyExpose();

    /// Send a sequence of keystrokes to the visualization window
    void SendKeySequence(const char *seq);

    // Directly call the functions assigned to the given keys. Unlike the above
    // function, SendKeySequence(), this function does not send X events and
    // actually disables the function SendExposeEvent() used by many of the
    // functions assigned to keys. Call MyExpose() after calling this function to
    // update the visualization window.
    void CallKeySequence(const char *seq);

    void AddIdleFunc(IdleFPtr func);
    void RemoveIdleFunc(IdleFPtr func);

    void ThreadsStop();
    void ThreadsRun();
    void ToggleThreads();

    void MainLoop();

    void Quit();

    void ToggleAntialiasing();
    void Screenshot() { Screenshot(""); }
    void Screenshot(std::string filename);
    void PrintToPDF();

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

    template<typename... Args>
    void AddKeyEvent(int key, void(*eh)(Args...))
    {
        wnd->setOnKeyDown(key, eh);
    }

    template<typename TScene>
    void AddKeyEvent(int key, void (TScene::*eh)(), bool exposeAfter = true)
    {
       auto wrapped_eh = [this, eh, exposeAfter](GLenum e)
       {
           TScene* pScene = dynamic_cast<TScene*>(locscene.get());
           (pScene->*eh)();
           if (exposeAfter) { SendExposeEvent(); }
       };
       wnd->setOnKeyDown(key, wrapped_eh);
    }

    void AddKeyEvent(int key, void (*eh)(GLVisWindow*), bool exposeAfter = true)
    {
       auto wrapped_eh = [this, eh, exposeAfter](GLenum e)
       {
           (*eh)(this);
           if (exposeAfter) { SendExposeEvent(); }
       };
       wnd->setOnKeyDown(key, wrapped_eh);
    }

    void AddKeyEvent(int key, void (*eh)(GLVisWindow*, GLenum), bool exposeAfter = true)
    {
       auto wrapped_eh = [this, eh, exposeAfter](GLenum e)
       {
           (*eh)(this, e);
           if (exposeAfter) { SendExposeEvent(); }
       };
       wnd->setOnKeyDown(key, wrapped_eh);
    }
private:
    void InitFont();
    bool SetFont(const vector<std::string>& patterns, int height);

    void SetKeyEventHandler(int key, void (GLVisWindow::*handler)());
    void SetKeyEventHandler(int key, void (GLVisWindow::*handler)(GLenum));

    bool MainIdleFunc();
    bool CommunicationIdleFunc();

    void MyReshape(GLsizei w, GLsizei h);
    void MyExpose(GLsizei w, GLsizei h);

    // Internal event handlers for small scene rotations
    void Key1Pressed();
    void Key2Pressed();
    void Key3Pressed();
    void Key4Pressed();
    void Key5Pressed();
    void Key6Pressed();
    void Key7Pressed();
    void Key8Pressed();
    void Key9Pressed();

    // Internal event handlers for other scene rotations, transformations
    void KeyLeftPressed(GLenum);
    void KeyRightPressed(GLenum);
    void KeyUpPressed(GLenum);
    void KeyDownPressed(GLenum);
    void KeyJPressed();
    void KeyMinusPressed();
    void KeyPlusPressed();

    // Internal event handler for toggling state of threads
    void ThreadsPauseFunc(GLenum);

    std::unique_ptr<SdlWindow> wnd;
    std::unique_ptr<VisualizationScene> locscene;
    std::unique_ptr<GLVisCommand> glvis_command;
    std::unique_ptr<communication_thread> comm_thread;
    StreamState prob_state;
    bool use_idle = false;

    // Idle function callbacks
    mfem::Array<IdleFPtr> idle_funcs{};
    int last_idle_func = 0;

    int visualize = 0;

    bool disableSendExposeEvent = false;

    GlVisFont font;
    std::string priority_font;
    int font_size = 12;

    struct RotationControl;
    std::unique_ptr<RotationControl> rot_data;

};

/// Send expose event. In our case MyReshape is executed and Draw after it.
[[deprecated]] void SendExposeEvent();

[[deprecated]] void MyExpose();

void MainLoop(GLVisWindow* wnd);

[[deprecated]] SdlWindow * GetAppWindow();
[[deprecated]] GLVisWindow * GetGLVisWindow();
VisualizationScene * GetVisualizationScene();

void TouchPinch(SDL_MultiGestureEvent & e);



/// Take a screenshot using libtiff, libpng or sdl2
//int Screenshot(const char *fname, bool convert = false);

/// Send a sequence of keystrokes to the visualization window
void SendKeySequence(const char *seq);

// Directly call the functions assigned to the given keys. Unlike the above
// function, SendKeySequence(), this function does not send X events and
// actually disables the function SendExposeEvent() used by many of the
// functions assigned to keys. Call MyExpose() after calling this function to
// update the visualization window.
void CallKeySequence(const char *seq);


int GetMultisample();
void SetMultisample(int m);

void SetLineWidth(float width);
float GetLineWidth();
void SetLineWidthMS(float width_ms);
float GetLineWidthMS();

#endif
