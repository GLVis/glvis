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

#ifndef GLVIS_SDL_HPP
#define GLVIS_SDL_HPP

#include <string>
#include <memory>
#include <functional>
#include <map>
#include <mutex>
#include <condition_variable>
#include <deque>
#include "gl/renderer.hpp"

struct EventInfo
{
   GLint mouse_x;
   GLint mouse_y;
   SDL_Keymod keymod;
};

using TouchDelegate = void (*)(SDL_MultiGestureEvent&);
using MouseDelegate = void (*)(EventInfo*);
using KeyDelegate = std::function<void(GLenum)>;
using WindowDelegate = void (*)(int, int);
using Delegate = void (*)();
using IdleDelegate = bool (*)();

class SdlMainThread;
SdlMainThread& GetMainThread();

class SdlWindow
{
private:
   friend class SdlMainThread;
   struct Handle
   {
      SDL_Window * hwnd{nullptr};
      SDL_GLContext gl_ctx{};

      Handle() = default;

      Handle(const std::string& title, int x, int y, int w, int h,
             Uint32 wndflags);

      ~Handle();

      Handle(Handle&& other)
         : hwnd(other.hwnd), gl_ctx(other.gl_ctx)
      {
         other.hwnd = nullptr;
         other.gl_ctx = 0;
      }

      Handle& operator= (Handle&& other)
      {
         std::swap(hwnd, other.hwnd);
         std::swap(gl_ctx, other.gl_ctx);
         return *this;
      }

      bool isInitialized() const
      {
         return (hwnd != nullptr && gl_ctx != 0);
      }
   };

   int window_id = -1;
   Handle handle;
   std::unique_ptr<gl3::MeshRenderer> renderer;
   static const int high_dpi_threshold = 144;
   // The display is high-dpi when:
   // - SDL's "screen coordinates" sizes are different from the pixel sizes, or
   // - either the horizontal or the vertical dpi, as returned by getDpi(),
   //   is >= high_dpi_threshold, defined above.
   bool high_dpi = false;
   // Ratio of SDL's "screen coordinates" to GLVis' "screen coordinates":
   // - set to 1 on non-high-dpi displays,
   // - set to 1 on high-dpi displays where SDL's "screen coordinates" sizes are
   //   different from the pixels sizes (e.g. Mac retina displays),
   // - set to 2 on other high-dpi displays, so that GLVis can always work with
   //   scaled "screen coordinates" on all high-dpi displays.
   float pixel_scale_x = 1.0f, pixel_scale_y = 1.0f;

   bool running;

   IdleDelegate onIdle{nullptr};
   Delegate onExpose{nullptr};
   WindowDelegate onReshape{nullptr};
   std::map<int, KeyDelegate> onKeyDown;
   std::map<int, MouseDelegate> onMouseDown;
   std::map<int, MouseDelegate> onMouseUp;
   std::map<int, MouseDelegate> onMouseMove;
   TouchDelegate onTouchPinch{nullptr};
   TouchDelegate onTouchRotate{nullptr};

#ifdef __EMSCRIPTEN__
   std::string canvas_id_;
#endif

   enum class RenderState
   {
      // window displayed is fully current (no events or backbuffer updates pending)
      Updated,
      // events issued which may require a call to MyExpose
      ExposePending,
      // back buffer updated by MyExpose, now awaiting swap to be displayed on window
      SwapPending
   };

   RenderState wnd_state{RenderState::Updated};

   bool update_before_expose{false};

   // bool requiresExpose;
   bool takeScreenshot{false};
   std::string screenshot_file;
   bool screenshot_convert;
   bool lastKeyDownProcessed;
   Uint16 lastKeyDownMods;
   char lastKeyDownChar;

   // internal event handlers
   void windowEvent(SDL_WindowEvent& ew);
   void motionEvent(SDL_MouseMotionEvent& em);
   void mouseEventDown(SDL_MouseButtonEvent& eb);
   void mouseEventUp(SDL_MouseButtonEvent& eb);
   void keyDownEvent(SDL_Keysym& ks);
   void textInputEvent(const SDL_TextInputEvent &e);
   void multiGestureEvent(SDL_MultiGestureEvent & e);

   bool is_multithreaded{true};

   bool call_idle_func{false};

   // Hand off events to the SdlWindow. Intended to be called by the main SDL
   // thread in MainThread::MainLoop().
   void queueEvents(std::vector<SDL_Event> events)
   {
      {
         std::lock_guard<std::mutex> evt_guard{event_mutex};
         waiting_events.insert(waiting_events.end(), events.begin(), events.end());
      }
      if (is_multithreaded)
      {
         events_available.notify_all();
      }
   }

   std::string saved_keys;

   std::condition_variable events_available;
   std::mutex event_mutex;
   // The window-specific events collected by the main event thread.
   std::deque<SDL_Event> waiting_events;
public:
   SdlWindow();
   ~SdlWindow();

   // Initializes SDL2 and starts the main loop. Should be called from the main
   // thread only.
   static void StartSDL(bool server_mode);

   /// Creates a new OpenGL window. Returns false if SDL or OpenGL initialization
   /// fails.
   bool createWindow(const char * title, int x, int y, int w, int h,
                     bool legacyGlOnly);

   /// Runs the window loop.
   void mainLoop();
   void mainIter();

   // Called by worker threads in GLVisCommand::signal()
   void signalLoop();

   void setOnIdle(IdleDelegate func) { onIdle = func; }
   void setOnExpose(Delegate func) { onExpose = func; }
   void setOnReshape(WindowDelegate func) { onReshape = func; }

   void setOnKeyDown(int key, Delegate func)
   {
      onKeyDown[key] = [func](GLenum) { func(); };
   }
   void setOnKeyDown(int key, KeyDelegate func) { onKeyDown[key] = func; }

   void setOnMouseDown(int btn, MouseDelegate func) { onMouseDown[btn] = func; }
   void setOnMouseUp(int btn, MouseDelegate func) { onMouseUp[btn] = func; }
   void setOnMouseMove(int btn, MouseDelegate func) { onMouseMove[btn] = func; }

   void setTouchPinchCallback(TouchDelegate cb) { onTouchPinch = cb; }
   void setTouchRotateCallback(TouchDelegate cb) { onTouchRotate = cb; }

   void clearEvents()
   {
      onIdle = nullptr;
      onExpose = nullptr;
      onReshape = nullptr;
      onKeyDown.clear();
      onMouseUp.clear();
      onMouseDown.clear();
      onMouseMove.clear();
   }

   void callKeyDown(SDL_Keycode k, Uint16 mod=0)
   {
      if (onKeyDown[k])
      {
         onKeyDown[k](mod);
      }
   }

   void getWindowSize(int& w, int& h);
   void getGLDrawSize(int& w, int& h);
   void getDpi(int& wdpi, int& hdpi);
   /// This property is set by createWindow().
   bool isHighDpi() const { return high_dpi; }

   gl3::MeshRenderer& getRenderer() { return *renderer.get(); }
   void setWindowTitle(std::string& title);
   void setWindowTitle(const char* title);
   void setWindowSize(int w, int h);
   void setWindowPos(int x, int y);

   void signalKeyDown(SDL_Keycode k, SDL_Keymod m = KMOD_NONE);
   void signalExpose() { wnd_state = RenderState::ExposePending; }
   void signalSwap() { wnd_state = RenderState::SwapPending; }
   void signalQuit() { running = false; }

   /// Returns the keyboard events that have been logged by the window.
   std::string getSavedKeys() const { return saved_keys; }

   /// Queues a screenshot to be taken.
   void screenshot(std::string filename, bool convert = false)
   {
      takeScreenshot = true;
      screenshot_file = filename;
      screenshot_convert = convert;
      // Queue up an expose, so Screenshot() can pull image from the back
      // buffer
      signalExpose();
   }

   void swapBuffer();

   operator bool() { return handle.isInitialized(); }
   bool isWindowInitialized() { return handle.isInitialized(); }
   /// Returns true if the OpenGL context was successfully initialized.
   bool isGlInitialized();

   bool isSwapPending() { return wnd_state == RenderState::SwapPending; }
   bool isExposePending() { return wnd_state == RenderState::ExposePending; }

#ifdef __EMSCRIPTEN__
   std::string getCanvasId() const { return canvas_id_; }
   void setCanvasId(std::string canvas_id) { canvas_id_ = canvas_id; }
#endif
};

#endif
