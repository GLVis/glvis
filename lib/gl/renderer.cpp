// Copyright (c) 2010-2020, Lawrence Livermore National Security, LLC. Produced
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

GLenum MeshRenderer::getDeviceAlphaChannel()
{
   if (!device) { return GL_NONE; }
#ifdef __EMSCRIPTEN__
   const std::string versionString
       = reinterpret_cast<const char*>(glGetString(GL_VERSION));
   if (versionString.find("OpenGL ES 3.0") != std::string::npos)
   {
       // WebGL 2 uses GL_RED as the single-channel texture
       return GL_RED;
   }
   else
   {
       // WebGL uses GL_ALPHA as the single-channel texture
       return GL_ALPHA;
   }
#else
   if (device->getType() == GLDevice::FF_DEVICE)
   {
      return GL_ALPHA;
   }
   else
   {
      return GL_RED;
   }
#endif
}

void MeshRenderer::setAntialiasing(bool aa_status)
{
   if (msaa_enable != aa_status)
   {
      msaa_enable = aa_status;
      if (msaa_enable)
      {
         device->enableMultisample();
         device->enableBlend();
         device->setLineWidth(line_w_aa);
      }
      else
      {
         device->disableMultisample();
         device->disableBlend();
         device->setLineWidth(line_w);
      }
   }
}

void MeshRenderer::setLineWidth(float w)
{
   line_w = w;
   if (device && !msaa_enable)
   {
      device->setLineWidth(line_w);
   }
}

void MeshRenderer::setLineWidthMS(float w)
{
   line_w_aa = w;
   if (device && msaa_enable)
   {
      device->setLineWidth(line_w_aa);
   }
}

void MeshRenderer::init()
{
}

void MeshRenderer::render(const RenderQueue& queue)
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   // elements containing opaque objects should be rendered first
   RenderQueue sorted_queue = queue;
   std::stable_partition(sorted_queue.begin(), sorted_queue.end(),
                         [](RenderQueue::value_type& renderPair)
   {
      return !renderPair.first.contains_translucent;
   });
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
      std::vector<TextBuffer*> text_bufs;
      GlDrawable* curr_drawable = q_elem.second;
      for (int i = 0; i < NUM_LAYOUTS; i++)
      {
         for (size_t j = 0; j < GlDrawable::NUM_SHAPES; j++)
         {
            if (curr_drawable->buffers[i][j])
            {
               if (i == LAYOUT_VTX_TEXTURE0 || i == LAYOUT_VTX_NORMAL_TEXTURE0)
               {
                  tex_bufs.emplace_back(curr_drawable->buffers[i][j].get()->getHandle());
               }
               else
               {
                  no_tex_bufs.emplace_back(curr_drawable->buffers[i][j].get()->getHandle());
               }
            }
            if (curr_drawable->indexed_buffers[i][j])
            {
               if (i == LAYOUT_VTX_TEXTURE0 || i == LAYOUT_VTX_NORMAL_TEXTURE0)
               {
                  tex_bufs.emplace_back(curr_drawable->indexed_buffers[i][j].get()->getHandle());
               }
               else
               {
                  no_tex_bufs.emplace_back(
                     curr_drawable->indexed_buffers[i][j].get()->getHandle());
               }
            }
         }
      }
      text_bufs.emplace_back(&curr_drawable->text_buffer);
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
      for (TextBuffer* buf : text_bufs)
      {
         device->drawDeviceBuffer(*buf);
      }
      device->enableDepthWrite();
      if (!msaa_enable) { device->disableBlend(); }
   }
}

CaptureBuffer MeshRenderer::capture(const RenderQueue& queue)
{
   CaptureBuffer cbuf;
   device->initXfbMode();
   for (auto& q_elem : queue)
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
      std::vector<TextBuffer*> text_bufs;
      GlDrawable* curr_drawable = q_elem.second;
      for (int i = 0; i < NUM_LAYOUTS; i++)
      {
         for (size_t j = 0; j < GlDrawable::NUM_SHAPES; j++)
         {
            if (curr_drawable->buffers[i][j])
            {
               if (i == LAYOUT_VTX_TEXTURE0 || i == LAYOUT_VTX_NORMAL_TEXTURE0)
               {
                  tex_bufs.emplace_back(curr_drawable->buffers[i][j].get()->getHandle());
               }
               else
               {
                  no_tex_bufs.emplace_back(curr_drawable->buffers[i][j].get()->getHandle());
               }
            }
            if (curr_drawable->indexed_buffers[i][j])
            {
               if (i == LAYOUT_VTX_TEXTURE0 || i == LAYOUT_VTX_NORMAL_TEXTURE0)
               {
                  tex_bufs.emplace_back(curr_drawable->indexed_buffers[i][j].get()->getHandle());
               }
               else
               {
                  no_tex_bufs.emplace_back(
                     curr_drawable->indexed_buffers[i][j].get()->getHandle());
               }
            }
         }
      }
      text_bufs.emplace_back(&curr_drawable->text_buffer);

      device->attachTexture(GLDevice::SAMPLER_COLOR, color_tex);
      device->attachTexture(GLDevice::SAMPLER_ALPHA, alpha_tex);
      for (auto buf : tex_bufs)
      {
         device->captureXfbBuffer(cbuf, buf);
      }
      device->detachTexture(GLDevice::SAMPLER_COLOR);
      device->detachTexture(GLDevice::SAMPLER_ALPHA);
      for (auto buf : no_tex_bufs)
      {
         device->captureXfbBuffer(cbuf, buf);
      }
      if (!params.contains_translucent)
      {
         device->enableBlend();
         device->disableDepthWrite();
      }
      device->attachTexture(1, font_tex);
      device->setNumLights(0);
      for (TextBuffer* buf : text_bufs)
      {
         device->captureXfbBuffer(cbuf, *buf);
      }
   }
   device->exitXfbMode();
   return cbuf;
}

void MeshRenderer::buffer(GlDrawable* buf)
{
   for (int i = 0; i < NUM_LAYOUTS; i++)
   {
      for (size_t j = 0; j < GlDrawable::NUM_SHAPES; j++)
      {
         if (buf->buffers[i][j])
         {
            device->bufferToDevice((array_layout) i, *(buf->buffers[i][j].get()));
         }
         if (buf->indexed_buffers[i][j])
         {
            device->bufferToDevice((array_layout) i, *(buf->indexed_buffers[i][j].get()));
         }
      }
   }
   device->bufferToDevice(buf->text_buffer);
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

void GLDevice::captureXfbBuffer(CaptureBuffer& capture, const TextBuffer& t_buf)
{
   for (const auto& entry : t_buf)
   {
      glm::vec3 raster = glm::project(
                            glm::vec3(entry.rx, entry.ry, entry.rz),
                            model_view_mtx,
                            proj_mtx,
                            glm::vec4(0, 0, vp_width, vp_height));
      capture.text.emplace_back(raster, glm::make_vec4(static_color.data()),
                                entry.text);
   }
}

}
