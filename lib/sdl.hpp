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

#ifndef SDL_HPP
#define SDL_HPP
#include <string>
#include <memory>
#include <functional>
#include <map>
#include "platform_gl.hpp"


#define AUX_MOUSESTATUS		3

#define AUX_MOUSEDOWN	16
#define AUX_MOUSEUP	32
#define AUX_MOUSELOC	64

#define AUX_MOUSEX		0
#define AUX_MOUSEY		1

#define	AUX_LEFTBUTTON		1
#define	AUX_RIGHTBUTTON		2
#define	AUX_MIDDLEBUTTON	4

#define AUX_LEFT SDLK_LEFT
#define AUX_RIGHT SDLK_RIGHT

struct EventInfo {
    GLint event;
    GLint data[4];
};

typedef std::function<void(EventInfo*)> MouseDelegate;
typedef std::function<void(GLenum)> KeyDelegate;
typedef std::function<void(int, int)> WindowDelegate;
typedef std::function<void()> Delegate;


class SdlWindow
{
private:
    struct _SdlHandle;
    //Use shared_ptr to manage handle lifetimes
    std::shared_ptr<_SdlHandle> _handle;

    bool running;
    
    Delegate onIdle;
    Delegate onExpose;
    WindowDelegate onReshape;
    std::map<int, KeyDelegate> onKeyDown;
    std::map<int, MouseDelegate> onMouseDown;
    std::map<int, MouseDelegate> onMouseUp;
    std::map<int, MouseDelegate> onMouseMove;

    SDL_Keycode curr;
    bool keyDown;
   
    bool requiresExpose;
    bool takeScreenshot;
    std::string screenshot_file;
    //internal event handlers
    void windowEvent(SDL_WindowEvent& ew);
    void motionEvent(SDL_MouseMotionEvent& em);
    void mouseEventDown(SDL_MouseButtonEvent& eb);
    void mouseEventUp(SDL_MouseButtonEvent& eb);
    void keyEvent(SDL_Keysym& ks);
    void keyEvent(char c);
public:
    SdlWindow(const char * title, int w, int h);
    ~SdlWindow();

    /**
     * Creates a new OpenGL context with the window.
     * Returns false if OpenGL or GLEW intialization fails.
     */
    bool createGlContext();
    /**
     * Runs the window loop.
     */
    void mainLoop();
    bool mainIter();

    void setOnIdle(Delegate func) { onIdle = func; }
    void setOnExpose(Delegate func) { onExpose = func; }
    void setOnReshape(WindowDelegate func) { onReshape = func; }
    
    void setOnKeyDown(int key, Delegate func) {
        onKeyDown[key] = [func](GLenum e) { func(); };
    }
    void setOnKeyDown(int key, KeyDelegate func) { onKeyDown[key] = func; }
    
    void setOnMouseDown(int btn, MouseDelegate func) { onMouseDown[btn] = func; }
    void setOnMouseUp(int btn, MouseDelegate func) { onMouseUp[btn] = func; }
    void setOnMouseMove(int btn, MouseDelegate func) { onMouseMove[btn] = func; }
    
    void callKeyDown(SDL_Keycode k) { onKeyDown[k](0); } 

    void getWindowSize(int& w, int& h);
    void getDpi(int& wdpi, int& hdpi);
#ifdef GLVIS_X11
    int getXWindow();
#endif
    void setWindowTitle(std::string& title);
    void setWindowTitle(const char* title);
    void setWindowSize(int w, int h);
    void setWindowPos(int x, int y);

    void signalKeyDown(SDL_Keycode k, SDL_Keymod m = KMOD_NONE);
    void signalExpose() { requiresExpose = true; }
    void signalQuit() { running = false; }

    void screenshot(std::string filename) {
        takeScreenshot = true;
        screenshot_file = filename;
    }

    void swapBuffer();

    operator bool() { return (bool) _handle ; }
    bool isWindowInitialized() { return (bool) _handle; }
    /**
     * Returns true if the OpenGL context was successfully initialized.
     */
    bool isGlInitialized();
};

#endif
