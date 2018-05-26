#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>
#include <OpenGL/GL.h>
#include <iostream>
#include "sdl.hpp"

using std::cerr;
using std::endl;

struct SdlWindow::_SdlHandle {
    SDL_Window * hwnd;
    SDL_GLContext gl_ctx;
    _SdlHandle(SDL_Window * window, SDL_GLContext context)
        : hwnd(window)
        , gl_ctx(context) { }

    ~_SdlHandle() {
        if (context)
            SDL_GL_DeleteContext(gl_ctx);
        SDL_DestroyWindow(hwnd);
    }
};

bool SdlWindow::isGlInitialized() {
    return (_handle->gl_ctx != 0);
}

SdlWindow::SdlWindow(std::string& title, int w, int h) {

    if (!SDL_WasInit(SDL_INIT_VIDEO | SDL_INIT_EVENTS)) {
        if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS) != 0) {
            SDL_Log("Failed to initialize SDL: %s", SDL_GetError());
            return;
        }
    }

    _handle = std::make_shared(SDL_CreateWindow(title,
                                                SDL_WINDOWPOS_UNDEFINED,
                                                SDL_WINDOWPOS_UNDEFINED,
                                                w,
                                                h,
                                                SDL_WINDOW_OPENGL), 0);
}

bool SdlWindow::createGlContext() {
    if (!_handle) {
        cerr << "Can't initialize an OpenGL context without a valid window" << endl;
        return false;
    }
    if (_handle->gl_ctx) {
        // destroy existing opengl context
        SDL_GL_DeleteContext(_handle->gl_ctx);
        _handle->gl_ctx = 0;
    }
    // on OSX systems, only core profiles are available for OpenGL 3+, which
    // removes the fixed-function pipeline
    SDL_GL_SetAttribute( SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute( SDL_GL_CONTEXT_MINOR_VERSION, 1);
    
    // technically, SDL already defaults to double buffering and a depth buffer
    // all we need is an alpha channel

    //SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute( SDL_GL_ALPHA_SIZE, 1);
    //SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, 1);

    SDL_GLContext context = SDL_GL_CreateContext(hwnd);
    if (!context) {
        SDL_Log("Failed to create an OpenGL 2.1 context: %s", SDL_GetError());
        return false;
    }
    _handle->gl_ctx = context;
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        SDL_Log("Failed to initialize GLEW: %s", glewGetErrorString(err));
        return false;
    }
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
                onExpose();
            break;
    }
}

void SdlWindow::motionEvent(SDL_MouseMotionEvent& em) {
    EventInfo info;
    info.event = AUX_MOUSELOC;
    info.data[AUX_MOUSEX] = em.x;
    info.data[AUX_MOUSEY] = em.y;
    info.data[2] = SDL_GetModState();
    if (em.state | SDL_BUTTON_LMASK) {
        info.data[SDL_BUTTON_MOUSESTATUS] = AUX_LEFT;
        if (onMouseMove[SDL_BUTTON_LEFT])
            onMouseMove[SDL_BUTTON_LEFT](&info);
    } else if (em.state | SDL_BUTTON_RMASK) {
        info.data[SDL_BUTTON_MOUSESTATUS] = AUX_RIGHT
        if (onMouseMove[SDL_BUTTON_RIGHT])
            onMouseMove[SDL_BUTTON_RIGHT](&info);
    } else if (em.state | SDL_BUTTON_MMASK) {
        info.data[SDL_BUTTON_MOUSESTATUS] = AUX_MIDDLE
        if (onMouseMove[SDL_BUTTON_MIDDLE])
            onMouseMove[SDL_BUTTON_MIDDLE](&info);
    }
}

void SdlWindow::mouseEventDown(SDL_MouseButtonEvent& eb) {
    if (onMouseDown[eb.button]) {
        EventInfo info;
        info.event = AUX_MOUSEDOWN;
        info.data[AUX_MOUSEX] = x;
        info.data[AUX_MOUSEY] = y;
        info.data[2] = SDL_GetModState();
        info.data[AUX_MOUSESTATUS] = eb.button;
        result = onMouseDown[eb.button](&info);
    }
}

void SdlWindow::mouseEventUp(SDL_MouseButtonEvent& eb) {
    if (onMouseUp[eb.button]) {
        EventInfo info;
        info.event = AUX_MOUSEUP;
        info.data[AUX_MOUSEX] = x;
        info.data[AUX_MOUSEY] = y;
        info.data[2] = SDL_GetModState();
        info.data[AUX_MOUSESTATUS] = eb.button;
        result = onMouseUp[eb.button](&info);
    }
}

void SdlWindow::keyEvent(SDL_Keysym& ks) {
    if (e.key.keysym.mod & KMOD_SHIFT) {
        //check if separate caps handler exists
        if (onKeyDown[toupper(e.key.keysym.sym)]) {
            onKeyDown[toupper(e.key.keysym.sym)](e.key.keysym.mod);
            return;
        }
    }
    if (onKeyDown[e.key.keysym.sym])
        onKeyDown[e.key.keysym.sym](e.key.keysym.mod);
}

void SdlWindow::mainLoop() {
    bool running = true;
    SDL_Event e;
    while (running) {
        while (SDL_PollEvent(&e)) {
            switch(e.type) {
                case SDL_QUIT:
                    running = false;
                    break;
                case SDL_WINDOWEVENT:
                    windowEvent(e.window);
                    break;
                case SDL_KEYDOWN:
                    keyEvent(e.key.keysym);
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
        }
        if (onIdle)
            onIdle();
        SDL_GL_SwapWindow(_handle.hwnd);
    }
}

void getWindowSize(int& w, int& h) {
    if (_handle)
        SDL_GetWindowSize(_handle->hwnd, &w, &h);
}

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
        SDL_SetWindowPos(_handle->hwnd, x, y);
}

void SdlWindow::signalKeyDown(SDL_Keycode k, SDL_Keymod m) {
    SDL_Event event;
    event.type = SDL_KEYDOWN;
    event.key.keysym.sym = k;
    event.key.keysym.mod = m;
    SDL_PushEvent(&event);
}

void SdlWindow::signalExpose() {
   SDL_Event event;
   event.type = SDL_WINDOWEVENT;
   event.window.event = SDL_WINDOWEVENT_EXPOSED;
   SDL_PushEvent(&event);
}

