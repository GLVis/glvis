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

#include "depth_peel_oit.hpp"
#include "renderer_core.hpp"
#include "shader.hpp"

const std::string DepthPeelVS[2] =
{
#include "shaders/default.vert"
   ,
#include "shaders/depth_peel_passthrough.vert"
};

const std::string DepthPeelFS[] =
{
#include "shaders/lighting.glsl"
#include "shaders/depth_peel.frag"
   ,
#include "shaders/depth_peel_blend_back.frag"
   ,
#include "shaders/depth_peel_finalize.frag"
};


namespace gl3
{

void DepthPeeler::SetGLDevice(GLDevice* dev)
{
   IMainRenderPass::SetGLDevice(dev);
   {
      GLuint vbo;
      glGenBuffers(1, &vbo);
      rect_buf = BufObjHandle{vbo};
      float quad_verts[] =
      {
         -1.f, 1.f, -1.f, -1.f, 1.f, -1.f,
            -1.f, 1.f, 1.f, -1.f, 1.f, 1.f,
         };
      glBindBuffer(GL_ARRAY_BUFFER, vbo);
      glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(float), quad_verts, GL_STATIC_DRAW);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
   }

   std::unordered_map<int, std::string> attribMap =
   {
      { CoreGLDevice::ATTR_VERTEX, "vertex"},
      { CoreGLDevice::ATTR_TEXT_VERTEX, "textVertex"},
      { CoreGLDevice::ATTR_NORMAL, "normal"},
      { CoreGLDevice::ATTR_COLOR, "color"},
      { CoreGLDevice::ATTR_TEXCOORD0, "texCoord0"}
   };
   if (!main_prgm.create(DepthPeelVS[0], DepthPeelFS[0], attribMap, 3))
   {
      std::cerr << "Unable to create depth peeling main program." << std::endl;
   }
   if (!blend_prgm.create(DepthPeelVS[1], DepthPeelFS[1], attribMap, 1))
   {
      std::cerr << "Unable to create depth peeling blending program." << std::endl;
   }
   if (!finalize_prgm.create(DepthPeelVS[1], DepthPeelFS[2], attribMap, 1))
   {
      std::cerr << "Unable to create depth peeling blending program." << std::endl;
   }

   GLint vp[4];
   device->getViewport(vp);
   int tex_w = vp[2];
   int tex_h = vp[3];

   {
      GLuint fb_id;
      GLuint tex_ids[2];
      glGenFramebuffers(1, &fb_id);
      glGenTextures(2, tex_ids);
      opaque_fb = fb_id;
      opaqueColorTex = tex_ids[0];
      opaqueDepthTex = tex_ids[1];

      glBindFramebuffer(GL_FRAMEBUFFER, opaque_fb);

      glBindTexture(GL_TEXTURE_2D, opaqueColorTex);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, tex_w, tex_h, 0, GL_RGBA,
                   GL_HALF_FLOAT, nullptr);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                             opaqueColorTex, 0);

      glBindTexture(GL_TEXTURE_2D, opaqueDepthTex);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, tex_w, tex_h, 0,
                   GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                             opaqueDepthTex, 0);
   }

   for (int i = 0; i < 2; i++)
   {
      {
         GLuint fb_ids[2], tex_ids[3];
         glGenFramebuffers(2, fb_ids);
         glGenTextures(3, tex_ids);
         main_peel_fbs[i] = fb_ids[0];
         color_fbs[i] = fb_ids[1];
         frontColorTex[i] = tex_ids[0];
         backColorTex[i] = tex_ids[1];
         depthTex[i] = tex_ids[2];
      }
      glBindFramebuffer(GL_FRAMEBUFFER, main_peel_fbs[i]);

      glBindTexture(GL_TEXTURE_2D, depthTex[i]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, tex_w, tex_h, 0, GL_RG, GL_FLOAT,
                   nullptr);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                             depthTex[i], 0);

      glBindTexture(GL_TEXTURE_2D, frontColorTex[i]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, tex_w, tex_h, 0, GL_RGBA,
                   GL_HALF_FLOAT, nullptr);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D,
                             frontColorTex[i], 0);

      glBindTexture(GL_TEXTURE_2D, backColorTex[i]);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, tex_w, tex_h, 0, GL_RGBA,
                   GL_HALF_FLOAT, nullptr);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D,
                             backColorTex[i], 0);

      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D,
                             opaqueDepthTex, 0);

      glBindFramebuffer(GL_FRAMEBUFFER, color_fbs[i]);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                             frontColorTex[i], 0);
      glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D,
                             backColorTex[i], 0);
   }

   {
      GLuint fb_id, tex_id;
      glGenFramebuffers(1, &fb_id);
      glGenTextures(1, &tex_id);
      blend_back_fb = fb_id;
      backBlendTex = tex_id;
   }
   glBindFramebuffer(GL_FRAMEBUFFER, blend_back_fb);
   glBindTexture(GL_TEXTURE_2D, backBlendTex);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, tex_w, tex_h, 0, GL_RGBA,
                GL_HALF_FLOAT, nullptr);
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                          backBlendTex, 0);

   glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


void DepthPeeler::PreRender()
{
   // Clear back blend buffer
   {
      GLenum blendBufs = GL_COLOR_ATTACHMENT0;
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, blend_back_fb);
      glDrawBuffers(1, &blendBufs);
      glClearColor(0, 0, 0, 0);
      glClear(GL_COLOR_BUFFER_BIT);
   }
   for (int i = 0; i < 2; i++)
   {
      // Initial depth texture should be set to [-MAX_DEPTH, MAX_DEPTH]
      GLenum depthBufs = GL_COLOR_ATTACHMENT0;
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, main_peel_fbs[i]);
      glDrawBuffers(1, &depthBufs);
      glClearColor(MAX_DEPTH, MAX_DEPTH, 0, 0);
      glClear(GL_COLOR_BUFFER_BIT);

      GLenum colorBufs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, color_fbs[i]);
      glDrawBuffers(2, colorBufs);
      glClearColor(0, 0, 0, 0);
      glClear(GL_COLOR_BUFFER_BIT);
   }
   device->enableBlend();
}

void DepthPeeler::RenderOpaque(const RenderQueue& queue)
{
   GLint vp[4];
   device->getViewport(vp);
   int tex_w = vp[2];
   int tex_h = vp[3];
   DefaultPass opaque_pass;
   opaque_pass.SetGLDevice(device);
   opaque_pass.SetTargetFramebuffer(opaque_fb);
   opaque_pass.setPalette(*palette);
   opaque_pass.setFontTexture(font_tex);
   // Render opaque pass objects - this fills in the correct depth attachment
   // texture for the depth testing of translucent peels.
   opaque_pass.PreRender();
   opaque_pass.Render(queue);
   opaque_pass.PostRender();
   // The opaque object depth data will be passed onto the peeling render
   // passes, so all we need to do now is to blit the color data to the output
   // framebuffer.
   glBindFramebuffer(GL_READ_FRAMEBUFFER, opaque_fb);
   glBindFramebuffer(GL_DRAW_FRAMEBUFFER, *target);
   glBlitFramebuffer(0, 0, tex_w, tex_h,
                     0, 0, tex_w, tex_h,
                     GL_COLOR_BUFFER_BIT,
                     GL_LINEAR);
}

void DepthPeeler::DoRenderPass(int i, const RenderQueue& queue)
{
   int src_i = i % 2;
   int dst_i = 1 - src_i;
   // Clear our target buffers for this pass
   {
      GLenum depthBufs = GL_COLOR_ATTACHMENT0;
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, main_peel_fbs[dst_i]);
      glDrawBuffers(1, &depthBufs);
      glClearColor(-10., -10., 0, 0);
      glClear(GL_COLOR_BUFFER_BIT);

      GLenum colorBufs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, color_fbs[dst_i]);
      glDrawBuffers(2, colorBufs);
      glClearColor(0, 0, 0, 0);
      glClear(GL_COLOR_BUFFER_BIT);
   }

   // Setup main peel program and framebuffer
   dynamic_cast<CoreGLDevice*>(device)->bindExternalProgram(main_prgm);
   main_prgm.setOutputFramebuffer(main_peel_fbs[dst_i]);

   // Bind source depth and front color texture
   glActiveTexture(GL_TEXTURE0 + 2);
   glBindTexture(GL_TEXTURE_2D, depthTex[src_i]);

   glActiveTexture(GL_TEXTURE0 + 3);
   glBindTexture(GL_TEXTURE_2D, frontColorTex[src_i]);

   glUniform1i(main_prgm.uniform("lastDepthTex"), 2);
   glUniform1i(main_prgm.uniform("lastFrontColorTex"), 3);

   int color_tex = palette->GetColorTexture();
   int alpha_tex = palette->GetAlphaTexture();
   glBlendEquation(GL_MAX);
   glDepthMask(GL_FALSE);
   // Render the geometry to peel
   for (auto& geom : queue)
   {
      const RenderParams& params = geom.first;
      device->setTransformMatrices(params.model_view.mtx, params.projection.mtx);
      device->setMaterial(params.mesh_material);
      device->setNumLights(params.num_pt_lights);
      for (int i = 0; i < params.num_pt_lights; i++)
      {
         device->setPointLight(i, params.lights[i]);
      }
      device->setAmbientLight(params.light_amb_scene);
      device->setStaticColor(params.static_color);
      device->setClipPlaneUse(params.use_clip_plane);
      device->setClipPlaneEqn(params.clip_plane_eqn);
      // aggregate buffers with common parameters
      std::vector<int> tex_bufs, no_tex_bufs;
      GlDrawable* curr_drawable = geom.second;
      auto buffers = curr_drawable->getArrayBuffers();
      for (const IVertexBuffer* buf : buffers)
      {
         if (buf->getVertexLayout() == LAYOUT_VTX_TEXTURE0
             || buf->getVertexLayout() == LAYOUT_VTX_NORMAL_TEXTURE0)
         {
            tex_bufs.emplace_back(buf->getHandle());
         }
         else
         {
            no_tex_bufs.emplace_back(buf->getHandle());
         }
      }
      device->attachTexture(GLDevice::SAMPLER_COLOR, color_tex);
      device->attachTexture(GLDevice::SAMPLER_ALPHA, alpha_tex);
      for (auto buf : tex_bufs)
      {
         device->drawDeviceBuffer(buf);
      }
      device->detachTexture(GLDevice::SAMPLER_COLOR);
      device->detachTexture(GLDevice::SAMPLER_ALPHA);
      for (auto buf : no_tex_bufs)
      {
         device->drawDeviceBuffer(buf);
      }
      device->attachTexture(1, font_tex);
      device->setNumLights(0);
      device->drawDeviceBuffer(curr_drawable->getTextBuffer());
   }
   glDepthMask(GL_TRUE);

   // Blend just-written back layer separately
   blend_prgm.bind();
   blend_prgm.setOutputFramebuffer(blend_back_fb);
   glBlendEquation(GL_FUNC_ADD);
   glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE,
                       GL_ONE_MINUS_SRC_ALPHA);

   glActiveTexture(GL_TEXTURE0 + 2);
   glBindTexture(GL_TEXTURE_2D, backColorTex[dst_i]);
   glUniform1i(blend_prgm.uniform("lastBackColor"), 2);

   glBindBuffer(GL_ARRAY_BUFFER, rect_buf);
   glEnableVertexAttribArray(CoreGLDevice::ATTR_VERTEX);
   glVertexAttribPointer(CoreGLDevice::ATTR_VERTEX,
                         2, GL_FLOAT, false, 0, 0);
   glDrawArrays(GL_TRIANGLES, 0, 6);
}

void DepthPeeler::PostRender()
{
   int src_i = (NUM_PASSES+1) % 2;

   finalize_prgm.bind();
   finalize_prgm.setOutputFramebuffer(*target);
   glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

   glActiveTexture(GL_TEXTURE0);
   glBindTexture(GL_TEXTURE_2D, frontColorTex[src_i]);

   glActiveTexture(GL_TEXTURE0 + 1);
   glBindTexture(GL_TEXTURE_2D, backBlendTex);

   glUniform1i(finalize_prgm.uniform("lastFrontColor"), 0);
   glUniform1i(finalize_prgm.uniform("lastBackColor"), 1);

   glBindBuffer(GL_ARRAY_BUFFER, rect_buf);
   glEnableVertexAttribArray(CoreGLDevice::ATTR_VERTEX);
   glVertexAttribPointer(CoreGLDevice::ATTR_VERTEX,
                         2, GL_FLOAT, false, 0, 0);
   glDrawArrays(GL_TRIANGLES, 0, 6);

   // Reset to the default program state
   device->initRenderMode();
   glBlendEquation(GL_FUNC_ADD);
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void DepthPeeler::Render(const RenderQueue& queue)
{
   // elements containing opaque objects should be rendered first
   RenderQueue sorted_queue = queue;
   auto begin_translucent =
      std::stable_partition(sorted_queue.begin(), sorted_queue.end(),
                            [](RenderQueue::value_type& renderPair)
   {
      return !renderPair.first.contains_translucent;
   });
   // Partition into two queues, one with opaque objects and one with
   // translucent objects
   RenderQueue opaque_queue(sorted_queue.begin(), begin_translucent);
   RenderQueue translucent_queue(begin_translucent, sorted_queue.end());

   RenderOpaque(opaque_queue);
   for (int i = 0; i < NUM_PASSES; i++)
   {
      DoRenderPass(i, translucent_queue);
   }
}

}


