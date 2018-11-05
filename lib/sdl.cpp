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

#include "platform_gl.hpp"
#include <iostream>
#include <chrono>
#include "sdl.hpp"
#include <SDL2/SDL_syswm.h>
#include "visual.hpp"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

using std::cerr;
using std::endl;

extern int GetMultisample();
extern int visualize;

struct SdlWindow::_SdlHandle {
    SDL_Window * hwnd;
    SDL_GLContext gl_ctx;
    _SdlHandle()
        : hwnd(nullptr)
        , gl_ctx(0) { }

    ~_SdlHandle() {
        if (gl_ctx)
            SDL_GL_DeleteContext(gl_ctx);
        SDL_DestroyWindow(hwnd);
    }
};

bool SdlWindow::isGlInitialized() {
    return (_handle->gl_ctx != 0);
}

SdlWindow::SdlWindow()
    : onIdle(nullptr)
    , onExpose(nullptr)
    , onReshape(nullptr)
    , ctrlDown(false)
    , requiresExpose(false)
    , takeScreenshot(false) {
}

bool SdlWindow::createWindow(const char * title, int w, int h) {
    if (!SDL_WasInit(SDL_INIT_VIDEO | SDL_INIT_EVENTS)) {
        if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0) {
            cerr << "Failed to initialize SDL: " << SDL_GetError() << endl;
            return false;
        }
    }

    //destroy any existing SDL window
    _handle.reset(new _SdlHandle);

    // technically, SDL already defaults to double buffering and a depth buffer
    // all we need is an alpha channel

    //SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 8);
    SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 24);
    if (GetMultisample() > 0) {
        SDL_GL_SetAttribute( SDL_GL_MULTISAMPLEBUFFERS, 1);
        SDL_GL_SetAttribute( SDL_GL_MULTISAMPLESAMPLES, GetMultisample());
    }

    _handle->hwnd = SDL_CreateWindow(title,
                                     SDL_WINDOWPOS_UNDEFINED,
                                     SDL_WINDOWPOS_UNDEFINED,
                                     w,
                                     h,
                                     SDL_WINDOW_OPENGL |
                                     SDL_WINDOW_ALLOW_HIGHDPI |
                                     SDL_WINDOW_RESIZABLE);

    SDL_GLContext context = SDL_GL_CreateContext(_handle->hwnd);
    if (!context) {
        cerr << "Failed to create an OpenGL context: " << SDL_GetError() << endl;
        return false;
    }
    _handle->gl_ctx = context;

#ifndef __EMSCRIPTEN__
    SDL_GL_SetSwapInterval(0);
    glEnable(GL_DEBUG_OUTPUT);
#endif

    GLenum err = glewInit();
    if (err != GLEW_OK) {
        cerr << "Failed to initialize GLEW: " << glewGetErrorString(err) << endl;
        return false;
    }
    cerr << glGetString(GL_VERSION) << "\n";

    if (!GLEW_ARB_vertex_shader ||
        !GLEW_ARB_fragment_shader ||
        !GLEW_ARB_shading_language_100) {
        cerr << "Shader support missing, failed to launch." << endl;
    }

    if (!GLEW_VERSION_3_0) {
        if (GLEW_EXT_transform_feedback) {
            glBindBufferBase            = glBindBufferBaseEXT;
            glTransformFeedbackVaryings = glTransformFeedbackVaryingsEXT;
            glBeginTransformFeedback    = glBeginTransformFeedbackEXT;
            glEndTransformFeedback      = glEndTransformFeedbackEXT;
        }
    }
    return true;
}

SdlWindow::~SdlWindow() {
}

void SdlWindow::windowEvent(SDL_WindowEvent& ew) {
    switch(ew.event) {
        case SDL_WINDOWEVENT_SIZE_CHANGED:
            if (onReshape)
                onReshape(ew.data1, ew.data2);
            break;
        case SDL_WINDOWEVENT_EXPOSED:
            if (onExpose)
                requiresExpose = true;
            break;
        default:
            break;
    }
}

void SdlWindow::motionEvent(SDL_MouseMotionEvent& em) {
    EventInfo info = {
        em.x, em.y,
        SDL_GetModState()
    };
    if (em.state & SDL_BUTTON_LMASK) {
        if (onMouseMove[SDL_BUTTON_LEFT]) {
            onMouseMove[SDL_BUTTON_LEFT](&info);
        }
    } else if (em.state & SDL_BUTTON_RMASK) {
        if (onMouseMove[SDL_BUTTON_RIGHT]) {
            onMouseMove[SDL_BUTTON_RIGHT](&info); 
        }
    } else if (em.state & SDL_BUTTON_MMASK) {
        if (onMouseMove[SDL_BUTTON_MIDDLE]) {
            onMouseMove[SDL_BUTTON_MIDDLE](&info); 
        }
    }
}

void SdlWindow::mouseEventDown(SDL_MouseButtonEvent& eb) {
    if (onMouseDown[eb.button]) {
        EventInfo info = {
            eb.x, eb.y,
            SDL_GetModState()
        };
        onMouseDown[eb.button](&info); 
    }
}

void SdlWindow::mouseEventUp(SDL_MouseButtonEvent& eb) {
    if (onMouseUp[eb.button]) {
        EventInfo info = {
            eb.x, eb.y,
            SDL_GetModState()
        };
        onMouseUp[eb.button](&info);
    }
}

bool SdlWindow::keyEvent(SDL_Keysym& ks) {
    bool handled = false;
    if (ks.sym > 128 || ks.sym < 32) {
        if (onKeyDown[ks.sym])
            onKeyDown[ks.sym](ks.mod);
        handled = true;
    } else if (ctrlDown == true) {
        onKeyDown[ks.sym](ks.mod);
        handled = true;
    }
    if (ks.sym == SDLK_RCTRL || ks.sym == SDLK_LCTRL) {
        ctrlDown = true;
    }
    return handled;
}

bool SdlWindow::keyEvent(char c) {
    if (onKeyDown[c]) {
        onKeyDown[c](SDL_GetModState());
        return true;
    }
    return false;
}

bool SdlWindow::mainIter() {
    SDL_Event e;
    static bool useIdle = false;
    bool needsSwap = false;
    while (SDL_PollEvent(&e)) {
        bool renderKeyEvent = false;
        switch(e.type) {
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
                    || e.key.keysym.sym == SDLK_RCTRL) {
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
            break;
    }
#ifndef __EMSCRIPTEN__
    if (onIdle) {
        if (glvis_command == NULL || visualize == 2 || useIdle) {
            onIdle();
            needsSwap = true;
        } else {
            if (glvis_command->Execute() < 0)
                running = false;
        }
        useIdle = !useIdle;
    }
#else
    if (onIdle) {
        onIdle();
        needsSwap = true;
    }
#endif
    if (requiresExpose) {
        onExpose();
        requiresExpose = false;
        needsSwap = true;
    }
    return needsSwap;
}

void SdlWindow::mainLoop() {
    running = true;
#ifdef __EMSCRIPTEN__
    emscripten_set_main_loop_arg([](void* arg) {
                                        ((SdlWindow*) arg)->mainIter();
                                    }, this, 0, 1);
#else
    visualize = 1;
    while (running) {
        bool glSwap = mainIter();
        if (glSwap) {
            SDL_GL_SwapWindow(_handle->hwnd);
        }
        if (takeScreenshot) {
            Screenshot(screenshot_file.c_str());
            takeScreenshot = false;
        }
        // sleep for 2 seconds to avoid pegging CPU at 100%
        SDL_Delay(2);
    }
#endif
}

void SdlWindow::getWindowSize(int& w, int& h) {
    if (_handle) {
#ifdef __EMSCRIPTEN__
        int is_fullscreen;
        emscripten_get_canvas_size(&w, &h, &is_fullscreen);
#else
        SDL_GL_GetDrawableSize(_handle->hwnd, &w, &h);
#endif
    }
}

void SdlWindow::getDpi(int& w, int& h) {
    if (_handle) {
        int disp = SDL_GetWindowDisplayIndex(_handle->hwnd);
        float f_w, f_h;
        SDL_GetDisplayDPI(disp, NULL, &f_w, &f_h);
        w = f_w;
        h = f_h;
    }
}

#ifdef GLVIS_X11
Window SdlWindow::getXWindow() {
    SDL_SysWMinfo info;
    SDL_VERSION(&info.version);

    if (SDL_GetWindowWMInfo(window, &info)) {
        if (info.subsystem == SDL_SYSWM_X11) {
            return info.x11.window;
        }
    }
}
#endif

void SdlWindow::setWindowTitle(std::string& title) {
    setWindowTitle(title.c_str());
}

void SdlWindow::setWindowTitle(const char * title) {
    if (_handle)
        SDL_SetWindowTitle(_handle->hwnd, title);
}

void SdlWindow::setWindowSize(int w, int h) {
    if (_handle)
        SDL_SetWindowSize(_handle->hwnd, w, h);
}

void SdlWindow::setWindowPos(int x, int y) {
    if (_handle)
        SDL_SetWindowPosition(_handle->hwnd, x, y);
}

void SdlWindow::signalKeyDown(SDL_Keycode k, SDL_Keymod m) {
    SDL_Event event;
    event.type = SDL_KEYDOWN;
    event.key.keysym.sym = k;
    event.key.keysym.mod = m;
    SDL_PushEvent(&event);
}

void SdlWindow::swapBuffer() {
    SDL_GL_SwapWindow(_handle->hwnd);
}

