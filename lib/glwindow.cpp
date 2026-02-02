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

#include "glwindow.hpp"
#include "aux_vis.hpp"
#include "gl/renderer_core.hpp"
#include "gl/renderer_ff.hpp"

#ifdef GLVIS_DEBUG
#define PRINT_DEBUG(s) std::cerr << s
#else
#define PRINT_DEBUG(s) {}
#endif

bool GLWindow::initGLEW(bool legacyGlOnly)
{
   GLenum err = glewInit();
#ifdef GLEW_ERROR_NO_GLX_DISPLAY
   // NOTE: Hacky workaround for Wayland initialization failure
   // See https://github.com/nigels-com/glew/issues/172
   if (err == GLEW_ERROR_NO_GLX_DISPLAY)
   {
      std::cerr << "GLEW: No GLX display found. If you are using Wayland this can "
                << "be ignored." << std::endl;
      err = GLEW_OK;
   }
#endif
   if (err != GLEW_OK)
   {
      std::cerr << "FATAL: Failed to initialize GLEW: "
                << glewGetErrorString(err) << std::endl;
      return false;
   }

   // print versions
   PRINT_DEBUG("Using GLEW " << glewGetString(GLEW_VERSION) << std::endl);
   PRINT_DEBUG("Using GL " << glGetString(GL_VERSION) << std::endl);

   renderer.reset(new gl3::MeshRenderer);
   renderer->setSamplesMSAA(GetMultisample());
#ifndef __EMSCRIPTEN__
   if (!GLEW_VERSION_1_1)
   {
      std::cerr << "FATAL: Minimum of OpenGL 1.1 is required." << std::endl;
      return false;
   }
   if (!GLEW_VERSION_1_3)
   {
      // Multitexturing was introduced into the core OpenGL specification in
      // version 1.3; for versions before, we need to load the functions from
      // the ARB_multitexture extension.
      if (GLEW_ARB_multitexture)
      {
         glActiveTexture = glActiveTextureARB;
         glClientActiveTexture = glClientActiveTextureARB;
         glMultiTexCoord2f = glMultiTexCoord2fARB;
      }
      else
      {
         std::cerr << "FATAL: Missing OpenGL multitexture support." << std::endl;
         return false;
      }
   }
   if (!GLEW_VERSION_3_0 && GLEW_EXT_transform_feedback)
   {
      glBindBufferBase            = glBindBufferBaseEXT;
      // Use an explicit typecast to suppress an error from inconsistent types
      // that are present in older versions of GLEW.
      glTransformFeedbackVaryings =
         (PFNGLTRANSFORMFEEDBACKVARYINGSPROC)glTransformFeedbackVaryingsEXT;
      glBeginTransformFeedback    = glBeginTransformFeedbackEXT;
      glEndTransformFeedback      = glEndTransformFeedbackEXT;
   }
   if (!legacyGlOnly && (GLEW_VERSION_3_0
                         || (GLEW_VERSION_2_0 && GLEW_EXT_transform_feedback)))
   {
      // We require both shaders and transform feedback EXT_transform_feedback
      // was made core in OpenGL 3.0
      PRINT_DEBUG("Loading CoreGLDevice..." << std::endl);
      renderer->setDevice<gl3::CoreGLDevice>();
   }
   else
   {
      PRINT_DEBUG("Shader support missing, loading FFGLDevice..." << std::endl);
      renderer->setDevice<gl3::FFGLDevice>();
   }

#else
   renderer->setDevice<gl3::CoreGLDevice>();
#endif

   if ((err = glGetError()) != GL_NO_ERROR)
   {
      std::cerr << "Renderer init OpenGL error: " << err << std::endl;
      return false;
   }

   return true;
}

void GLWindow::recordKey(SDL_Keycode sym, SDL_Keymod mod)
{
   // Record the key in 'saved_keys':
   bool isAlt = mod & (KMOD_ALT);
   bool isCtrl = mod & (KMOD_CTRL);
   if (isAlt || isCtrl)
   {
      saved_keys += "[";
   }
   if (isCtrl) { saved_keys += "C-"; }
   if (isAlt) { saved_keys += "Alt-"; }
   if (sym >= 32 && sym < 127)
   {
      saved_keys += (char)(sym);
   }
   else
   {
      saved_keys += SDL_GetKeyName(sym);
   }
   if (isAlt || isCtrl)
   {
      saved_keys += "]";
   }
}
