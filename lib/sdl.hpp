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
#include "gl/renderer.hpp"
#include "sdl_helper.hpp"

struct EventInfo
{
   GLint mouse_x;
   GLint mouse_y;
   SDL_Keymod keymod;
};

typedef void (*MouseDelegate)(EventInfo*);
typedef std::function<void(GLenum)> KeyDelegate;
typedef void (*WindowDelegate)(int, int);
typedef void (*Delegate)();

class SdlWindow
{
private:
   struct Handle;
   std::unique_ptr<Handle> handle;
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

   std::unique_ptr<SdlNativePlatform> platform;

   SDL_SysWMinfo sysinfo;
   static Uint32 glvis_event_type;

   bool running;

   Delegate onIdle;
   Delegate onExpose;
   WindowDelegate onReshape;
   std::map<int, KeyDelegate> onKeyDown;
   std::map<int, MouseDelegate> onMouseDown;
   std::map<int, MouseDelegate> onMouseUp;
   std::map<int, MouseDelegate> onMouseMove;

   bool ctrlDown;

   enum class RenderState
   {
      // window displayed is fully current (no events or backbuffer updates pending)
      Updated,
      // events issued which may require a call to MyExpose
      ExposePending,
      // back buffer updated by MyExpose, now awaiting swap to be displayed on window
      SwapPending
   };

   RenderState wnd_state;

   //bool requiresExpose;
   bool takeScreenshot;
   std::string screenshot_file;

   int probeGLContextSupport();
   // internal event handlers
   void windowEvent(SDL_WindowEvent& ew);
   void motionEvent(SDL_MouseMotionEvent& em);
   void mouseEventDown(SDL_MouseButtonEvent& eb);
   void mouseEventUp(SDL_MouseButtonEvent& eb);
   void keyEvent(SDL_Keysym& ks);
   void keyEvent(char c);

   std::string saved_keys;
public:
   SdlWindow();
   ~SdlWindow();

   /// Creates a new OpenGL window. Returns false if SDL or OpenGL intialization
   /// fails.
   bool createWindow(const char * title, int x, int y, int w, int h,
                     bool legacyGlOnly);
   /// Runs the window loop.
   void mainLoop();
   void mainIter();

   // Called by worker threads in GLVisCommand::signal()
   void signalLoop();

   void setOnIdle(Delegate func) { onIdle = func; }
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

   void callKeyDown(SDL_Keycode k)
   {
      if (onKeyDown[k])
      {
         onKeyDown[k](0);
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
   void screenshot(std::string filename)
   {
      takeScreenshot = true;
      screenshot_file = filename;
   }

   void swapBuffer();

   operator bool() { return (bool) handle ; }
   bool isWindowInitialized() { return (bool) handle; }
   /// Returns true if the OpenGL context was successfully initialized.
   bool isGlInitialized();

   bool isSwapPending() { return wnd_state == RenderState::SwapPending; }
   bool isExposePending() { return wnd_state == RenderState::ExposePending; }
};

#endif
