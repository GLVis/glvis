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

#ifndef GLVIS_DEPTH_PEEL_OIT_HPP
#define GLVIS_DEPTH_PEEL_OIT_HPP

#include "renderer.hpp"
#include "shader.hpp"
#include "framebuffer.hpp"

namespace gl3
{

class DepthPeeler : public IMainRenderPass
{
public:
   virtual void SetGLDevice(GLDevice* device);

   virtual bool Filter(const RenderParams& param)
   {
      return true;
   }

   virtual void PreRender();
   virtual void Render(const RenderQueue& queue);
   virtual void PostRender();
private:
   const double MAX_DEPTH = 1.0;
   const int NUM_PASSES = 4;

   void CreateScreenPeelObjs();

   void RenderOpaque(const RenderQueue& queue);

   void DoRenderPass(int i, const RenderQueue& queue);

   TextureHandle CreateScreenTexture(GLenum internalFmt,
                                     GLenum format,
                                     GLenum type);

   int screen_w, screen_h;

   ShaderProgram main_prgm;
   ShaderProgram blend_prgm;
   ShaderProgram finalize_prgm;

   TextureHandle depthTex[2];
   TextureHandle frontColorTex[2];
   TextureHandle backColorTex[2];
   TextureHandle opaqueColorTex, opaqueDepthTex;

   TextureHandle backBlendTex;

   Framebuffer main_peel_fbs[2];
   Framebuffer color_fbs[2];
   Framebuffer blend_back_fb;
   Framebuffer opaque_fb;

   // Drawing full-screen rectangles
   BufObjHandle rect_buf;
};

}

#endif
