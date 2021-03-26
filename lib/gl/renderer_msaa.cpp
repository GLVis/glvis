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

#include "renderer_msaa.hpp"

namespace gl3
{

void MultisamplePass::SetGLDevice(GLDevice* dev)
{
   IRenderPass::SetGLDevice(dev);
   int max_msaa_samples;
#ifdef __EMSCRIPTEN__
   const std::string versionString
      = reinterpret_cast<const char*>(glGetString(GL_VERSION));
   bool is_webgl2 = (versionString.find("OpenGL ES 3.0") != std::string::npos);
   feat_use_fbo_antialias = is_webgl2;
   if (feat_use_fbo_antialias)
   {
      glGetIntegerv(GL_MAX_SAMPLES, &max_msaa_samples);
   }
#else
   // TODO: we could also support ARB_framebuffer_object
   feat_use_fbo_antialias = GLEW_VERSION_3_0;
   glGetIntegerv(GL_MAX_SAMPLES, &max_msaa_samples);
#endif
   if (msaa_samples > max_msaa_samples)
   {
      std::cerr << "GL_MAX_SAMPLES = " << max_msaa_samples
                << " but requested " << msaa_samples << "x MSAA. ";
      std::cerr << "Setting antialiasing mode to "
                << max_msaa_samples << "x MSAA." << endl;
      msaa_samples = max_msaa_samples;
   }
   if (feat_use_fbo_antialias)
   {
      GLuint colorBuf, depthBuf;
      glGenRenderbuffers(1, &colorBuf);
      glGenRenderbuffers(1, &depthBuf);
      renderBufs[0] = RenderBufHandle(colorBuf);
      renderBufs[1] = RenderBufHandle(depthBuf);

      GLuint fbo;
      glGenFramebuffers(1, &fbo);

      int vp[4];
      device->getViewport(vp);
      int width = vp[2];
      int height = vp[3];
      glBindRenderbuffer(GL_RENDERBUFFER, colorBuf);
      glRenderbufferStorageMultisample(GL_RENDERBUFFER, msaa_samples,
                                       GL_RGBA8, width, height);
      glBindRenderbuffer(GL_RENDERBUFFER, depthBuf);
      glRenderbufferStorageMultisample(GL_RENDERBUFFER, msaa_samples,
                                       GL_DEPTH_COMPONENT24, width, height);
      glBindRenderbuffer(GL_RENDERBUFFER, 0);

      glBindFramebuffer(GL_FRAMEBUFFER, fbo);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                GL_RENDERBUFFER, colorBuf);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                                GL_RENDERBUFFER, depthBuf);

      if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      {
         cerr << "Unable to create multisampled renderbuffer." << flush;
         glDeleteFramebuffers(1, &fbo);
      }
      else
      {
         msaaFb = FBOHandle(fbo);
      }
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
   }
}

void MultisamplePass::PreRender()
{
   if (msaa_enable)
   {
      if (feat_use_fbo_antialias)
      {
         glBindFramebuffer(GL_FRAMEBUFFER, msaaFb);
      }
      else
      {
         glEnable(GL_LINE_SMOOTH);
         device->enableBlend();
      }
#ifndef __EMSCRIPTEN__
      glEnable(GL_MULTISAMPLE);
#endif
      device->setLineWidth(line_w_aa);
   }
   else
   {
      device->setLineWidth(line_w);
   }
}

void MultisamplePass::PostRender()
{
   if (msaa_enable && feat_use_fbo_antialias && msaaFb)
   {
      int vp[4];
      device->getViewport(vp);
      int width = vp[2];
      int height = vp[3];
      GLuint colorBufId;
      glGenRenderbuffers(1, &colorBufId);
      RenderBufHandle colorBuf(colorBufId);

      GLuint fboId;
      glGenFramebuffers(1, &fboId);
      FBOHandle resolveFb(fboId);

      glBindRenderbuffer(GL_RENDERBUFFER, colorBuf);
      glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);
      glBindRenderbuffer(GL_RENDERBUFFER, 0);

      glBindFramebuffer(GL_FRAMEBUFFER, resolveFb);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                                GL_RENDERBUFFER, colorBuf);

      if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      {
         cerr << "Unable to create resolve renderbuffer." << endl;
         glBindFramebuffer(GL_FRAMEBUFFER, 0);
      }

      // bind our draw framebuffer and blit the multisampled image
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, resolveFb);
      glBindFramebuffer(GL_READ_FRAMEBUFFER, msaaFb);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glBlitFramebuffer(0, 0, width, height,
                        0, 0, width, height,
                        GL_COLOR_BUFFER_BIT,
                        GL_NEAREST);
#ifndef __EMSCRIPTEN__
      glDisable(GL_MULTISAMPLE);
#endif
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
      glBindFramebuffer(GL_READ_FRAMEBUFFER, resolveFb);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      glBlitFramebuffer(0, 0, width, height,
                        0, 0, width, height,
                        GL_COLOR_BUFFER_BIT,
                        GL_LINEAR);
   }
   else if (msaa_enable && !feat_use_fbo_antialias)
   {
      glDisable(GL_MULTISAMPLE);
      glDisable(GL_LINE_SMOOTH);
      device->disableBlend();
   }
   device->setLineWidth(line_w);
}

void MultisamplePass::SetAntialiasing(bool aa_status)
{
   msaa_enable = aa_status;
}

void MultisamplePass::SetLineWidth(float w)
{
   line_w = w;
   if (device && !msaa_enable)
   {
      device->setLineWidth(line_w);
   }
}

void MultisamplePass::SetLineWidthMS(float w)
{
   line_w_aa = w;
   if (device && msaa_enable)
   {
      device->setLineWidth(line_w_aa);
   }
}

}
