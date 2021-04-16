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

#include "renderer.hpp"

namespace gl3
{

// Beginning in OpenGL 3.0, there were two changes in texture format support:
// - The older single-channel internal format GL_ALPHA was deprecated in favor
//   of GL_RED
// - New sized internal formats were introduced, e.g. GL_RGBA32F defines a 4-
//   channel texture with each channel holding a 32-bit floating point value
//
// An additional complication is introduced with OpenGL ES 3/WebGL 2 - the
// unsized formats like GL_RED and GL_RGBA no longer support floating-point
// data being passed in, so use of the sized internal formats is obligatory in
// WebGL 2.
bool GLDevice::useLegacyTextureFmts()
{
#ifdef __EMSCRIPTEN__
   const std::string versionString
      = reinterpret_cast<const char*>(glGetString(GL_VERSION));
   if (versionString.find("OpenGL ES 3.0") != std::string::npos)
   {
      return false;
   }
   else
   {
      return true;
   }
#else
   return !GLEW_VERSION_3_0;
#endif
}

void GLDevice::attachFramebuffer(const FBOHandle& fbo)
{
   glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
   GLenum drawOutput = GL_BACK;
   if (fbo != 0)
   {
      drawOutput = GL_COLOR_ATTACHMENT0;
   }
   glDrawBuffers(1, &drawOutput);
}

void DefaultPass::Render(const RenderQueue& queue)
{
   auto clear_color = device->getClearColor();
   glClearColor(clear_color[0], clear_color[1], clear_color[2], clear_color[3]);
   // elements containing opaque objects should be rendered first
   RenderQueue sorted_queue = queue;
   std::stable_partition(sorted_queue.begin(), sorted_queue.end(),
                         [](RenderQueue::value_type& renderPair)
   {
      return !renderPair.first.contains_translucent;
   });
   device->attachFramebuffer(*target);
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   int color_tex = palette->GetColorTexture();
   int alpha_tex = palette->GetAlphaTexture();
   bool always_blend = device->isBlendEnabled();
   for (auto& q_elem : sorted_queue)
   {
      const RenderParams& params = q_elem.first;
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
      GlDrawable* curr_drawable = q_elem.second;
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
      if (params.contains_translucent)
      {
         device->enableBlend();
      }
      else
      {
         device->enableDepthWrite();
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
      if (!params.contains_translucent)
      {
         device->enableBlend();
         device->disableDepthWrite();
      }
      device->attachTexture(1, font_tex);
      device->setNumLights(0);
      device->drawDeviceBuffer(curr_drawable->getTextBuffer());
      device->enableDepthWrite();
      if (!always_blend) { device->disableBlend(); }
   }
}

void MeshRenderer::render(const vector<IMainRenderPass*>& main_passes,
                          const vector<IRenderPass*>& extra_passes,
                          const RenderQueue& queued)
{
   for (IMainRenderPass* pass : main_passes)
   {
      pass->SetGLDevice(device.get());
   }
   for (IRenderPass* pass : extra_passes)
   {
      pass->SetGLDevice(device.get());
   }
   // Step 1: Match renderables in the queue with the *first* render pass that
   //         can handle them.
   std::vector<RenderQueue> matched_queues(main_passes.size());
   for (auto drawable : queued)
   {
      for (size_t ipass = 0; ipass < main_passes.size(); ipass++)
      {
         if (main_passes[ipass]->Filter(drawable.first))
         {
            matched_queues[ipass].emplace_back(drawable);
            break;
         }
      }
   }
   // Step 2: Setup the framebuffer with the first extra pass, and render the
   //         queue with the main passes.
   FBOHandle default_target(0);
   std::reference_wrapper<const FBOHandle> curr_out = default_target;
   if (extra_passes.size() > 0)
   {
      extra_passes[0]->PreRender();
      curr_out = extra_passes[0]->GetSourceFramebuffer();
   }
   for (size_t ipass = 0; ipass < main_passes.size(); ipass++)
   {
      main_passes[ipass]->SetTargetFramebuffer(curr_out);
      main_passes[ipass]->PreRender();
      main_passes[ipass]->Render(matched_queues[ipass]);
      main_passes[ipass]->PostRender();
   }

   if (extra_passes.size() > 0)
   {
      for (size_t ipass = 1; ipass < extra_passes.size(); ipass++)
      {
         // Finalize last stage's results onto next stage
         extra_passes[ipass]->PreRender();
         curr_out = extra_passes[ipass]->GetSourceFramebuffer();
         extra_passes[ipass-1]->PostRender();
      }
      extra_passes[extra_passes.size() - 1]->SetTargetFramebuffer(default_target);
      extra_passes[extra_passes.size() - 1]->PostRender();
   }
}

void MeshRenderer::buffer(GlDrawable* buf)
{
   auto buffers = buf->getArrayBuffers();
   for (IVertexBuffer* buf : buffers)
   {
      IIndexedBuffer* ibuf = dynamic_cast<IIndexedBuffer*>(buf);
      if (ibuf)
      {
         device->bufferToDevice(buf->getVertexLayout(), *ibuf);
      }
      else
      {
         device->bufferToDevice(buf->getVertexLayout(), *buf);
      }
   }
   device->bufferToDevice(buf->getTextBuffer());
}

void GLDevice::init()
{
   // enable depth testing
   glDepthFunc(GL_LEQUAL);
   glEnable(GL_DEPTH_TEST);
   // enable polygon offset to expose mesh lines
   glPolygonOffset(1,1);
   glEnable(GL_POLYGON_OFFSET_FILL);
   // use "over" blending equation
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
   // generate a white default texture modulation with default texture will just
   // pass through input color
   GLuint default_texture;
   glGenTextures(1, &default_texture);
   glBindTexture(GL_TEXTURE_2D, default_texture);
   int black_color = 0xFFFFFFFF;
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE,
                &black_color);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

   passthrough_texture = TextureHandle(default_texture);
}

void GLDevice::setViewport(GLsizei w, GLsizei h)
{
   vp_width = w;
   vp_height = h;
   glViewport(0, 0, w, h);
}

void GLDevice::getViewport(GLint (&vp)[4])
{
   vp[0] = vp[1] = 0;
   vp[2] = vp_width;
   vp[3] = vp_height;
}

void GLDevice::setTransformMatrices(glm::mat4 model_view, glm::mat4 projection)
{
   model_view_mtx = model_view;
   proj_mtx = projection;
}

}
