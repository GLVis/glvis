// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-443271.
//
// This file is part of the GLVis visualization tool and library. For more
// information and source code availability see https://glvis.org.
//
// GLVis is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef GLVIS_GLWINDOW_HPP
#define GLVIS_GLWINDOW_HPP

#include <memory>
#include <map>
#include <functional>
#include <string>

#include "gl/renderer.hpp"

class MainThread
{
public:
   virtual ~MainThread() { }

   // Handles all operations that are expected to be handled on the main
   // thread (i.e. events and window creation)
   virtual void MainLoop(bool server_mode) = 0;
};

class GLWindow
{
   friend class EglMainThread;
public:
   typedef Uint8 SDL_Mousebutton;
   struct MouseEventInfo
   {
      Sint32 mouse_x;
      Sint32 mouse_y;
      SDL_Keymod keymod;
   };

   using IdleDelegate = bool (*)();
   using Delegate = void (*)();
   using KeyDelegate = std::function<void(SDL_Keymod)>;
   using MouseDelegate = void (*)(MouseEventInfo*);

protected:
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

   std::unique_ptr<gl3::MeshRenderer> renderer;

   IdleDelegate onIdle{};
   Delegate onExpose{};
   std::map<SDL_Keycode, KeyDelegate> onKeyDown;
   std::map<SDL_Mousebutton, MouseDelegate> onMouseDown;
   std::map<SDL_Mousebutton, MouseDelegate> onMouseUp;
   std::map<SDL_Mousebutton, MouseDelegate> onMouseMove;

#ifdef __EMSCRIPTEN__
   std::string canvas_id;
#endif
   std::string saved_keys;

   bool initGLEW(bool legacyGlOnly);
   void recordKey(SDL_Keycode k, SDL_Keymod m);

public:
   virtual ~GLWindow() = default;

   /// Creates a new OpenGL window
   virtual bool createWindow(const char *title, int x, int y, int w, int h,
                             bool legacyGlOnly) = 0;

   /// Returns the renderer object
   inline gl3::MeshRenderer& getRenderer() { return *renderer.get(); }

   /// Runs the window loop.
   virtual void mainLoop() = 0;
   virtual void mainIter() = 0;

   /// Signals addition of a new event
   virtual void signalLoop() = 0;

   void setOnIdle(IdleDelegate func) { onIdle = func; }
   void setOnExpose(Delegate func) { onExpose = func; }

   void setOnKeyDown(SDL_Keycode key, Delegate func)
   {
      onKeyDown[key] = [func](SDL_Keymod) { func(); };
   }
   void setOnKeyDown(SDL_Keycode key, KeyDelegate func) { onKeyDown[key] = func; }

   void setOnMouseDown(SDL_Mousebutton btn, MouseDelegate func) { onMouseDown[btn] = func; }
   void setOnMouseUp(SDL_Mousebutton btn, MouseDelegate func) { onMouseUp[btn] = func; }
   void setOnMouseMove(SDL_Mousebutton btn, MouseDelegate func) { onMouseMove[btn] = func; }

   virtual void clearEvents()
   {
      onIdle = nullptr;
      onExpose = nullptr;
      onKeyDown.clear();
      onMouseUp.clear();
      onMouseDown.clear();
      onMouseMove.clear();
   }

   void callKeyDown(SDL_Keycode k, SDL_Keymod mod = KMOD_NONE)
   {
      if (onKeyDown[k])
      {
         onKeyDown[k](mod);
      }
   }

   /// Returns size of the window
   virtual void getWindowSize(int& w, int& h) const { w = h = 0; }

   /// Returns the drawable size
   virtual void getGLDrawSize(int& w, int& h) const { w = h = 0; }

   /// Returns the resolution (DPI) of the display
   virtual void getDpi(int& wdpi, int& hdpi) const { wdpi = hdpi = 0; }

   /// Checks if the display has high resolution (DPI)
   virtual bool isHighDpi() const { return false; }

   /// Set title of the window (string version)
   virtual void setWindowTitle(const std::string& title) { }

   /// Set title of the window (C-string version)
   virtual void setWindowTitle(const char* title) { }

   /// Set window size
   virtual void setWindowSize(int w, int h) { }

   /// Set window position
   virtual void setWindowPos(int x, int y) { }

   /// Returns true if the window has been succesfully initialized
   virtual bool isWindowInitialized() const { return false; }

   /// Returns true if the OpenGL context was successfully initialized
   virtual bool isGlInitialized() const { return false; }

   /// Signals key down event
   virtual void signalKeyDown(SDL_Keycode k, SDL_Keymod m = KMOD_NONE) { }

   /// Signals quit event
   virtual void signalQuit() { }

   /// Signals expose event when objects have been updated
   virtual void signalExpose() { wnd_state = RenderState::ExposePending; }

   /// Signals swap event when the back buffer is ready for swapping
   virtual void signalSwap() { wnd_state = RenderState::SwapPending; }

   /// Checks if the expose event is pending
   virtual bool isExposePending() const { return wnd_state == RenderState::ExposePending; }

   /// Checks if the swap event is pending
   virtual bool isSwapPending() const { return wnd_state == RenderState::SwapPending; }

   /// Returns the keyboard events that have been logged by the window.
   std::string getSavedKeys() const { return saved_keys; }

   /// Saves a screenshot ot the file, performing conversion optionally
   virtual void screenshot(std::string filename, bool convert = false) { }

#ifdef __EMSCRIPTEN__
   std::string getCanvasId() const { return canvas_id; }
   void setCanvasId(std::string canvas_id_) { canvas_id = canvas_id_; }
#endif
};

#endif //GLVIS_GLWINDOW_HPP
