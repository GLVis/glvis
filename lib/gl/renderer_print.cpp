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

#include "renderer_print.hpp"
#include "renderer_ff.hpp"

namespace gl3
{

// Defined in renderer.hpp
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

// Defined in renderer_ff.hpp
void FFGLDevice::captureXfbBuffer(PaletteState& pal, CaptureBuffer& cbuf,
                                  int hnd)
{
   if (hnd == 0) { return; }
   if (disp_lists[hnd].count == 0) { return; }
   GLenum fbType;
   int fbStride;
   if (disp_lists[hnd].layout == VertexTex::layout
       || disp_lists[hnd].layout == VertexNormTex::layout)
   {
      //capture texture values too
      // [ X Y Z ] [ R G B A ] [ U V - - ]
      fbType = GL_3D_COLOR_TEXTURE;
      fbStride = 11;
   }
   else
   {
      // only capture pos and color
      // [ X Y Z ] [ R G B A ]
      fbType = GL_3D_COLOR;
      fbStride = 7;
   }
   // compute feedback buffer size
   int sizebuf = 0;
   if (disp_lists[hnd].shape == GL_LINES)
   {
      // for each line: LINE_TOKEN [Vtx] [Vtx]
      sizebuf = (disp_lists[hnd].count / 2) + disp_lists[hnd].count * fbStride;
   }
   else if (disp_lists[hnd].shape == GL_TRIANGLES)
   {
      // for each tri: POLY_TOKEN 3 [Vtx] [Vtx] [Vtx]
      // NOTE: when clip plane is enabled, we might get two triangles
      // or a quad for an input triangle. However, the other clipped
      // triangles get discarded, so this *should* be enough space.
      sizebuf = (disp_lists[hnd].count / 3) * (2 + fbStride * 4);
   }
   else
   {
      std::cerr << "Warning: unhandled primitive type in FFPrinter::preDraw()" <<
                std::endl;
      return;
   }
   // allocate feedback buffer
   vector<float> xfb_buf;
   xfb_buf.resize(sizebuf);
   glFeedbackBuffer(sizebuf, fbType, xfb_buf.data());
   // draw with feedback capture
   glRenderMode(GL_FEEDBACK);
   drawDeviceBuffer(hnd);
#ifndef GLVIS_DEBUG
   glRenderMode(GL_RENDER);
#else
   if (glRenderMode(GL_RENDER) < 0)
   {
      std::cerr << "Warning: feedback data exceeded available buffer size" <<
                std::endl;
   }
#endif
   size_t tok_idx = 0;
   // process feedback buffer
   while (tok_idx < xfb_buf.size())
   {
      switch ((GLuint)xfb_buf[tok_idx])
      {
         case GL_LINE_TOKEN:
         case GL_LINE_RESET_TOKEN:
         {
            tok_idx++;
            glm::vec3 coord0 = glm::make_vec3(&xfb_buf[tok_idx]),
                      coord1 = glm::make_vec3(&xfb_buf[tok_idx + fbStride]);
            glm::vec4 color0 = glm::make_vec4(&xfb_buf[tok_idx + 3]),
                      color1 = glm::make_vec4(&xfb_buf[tok_idx + 3 + fbStride]);
            if (fbStride == 11)
            {
               // get texture
               pal.GetColorFromVal(xfb_buf[tok_idx + 7], glm::value_ptr(color0));
               pal.GetColorFromVal(xfb_buf[tok_idx + 7 + fbStride], glm::value_ptr(color1));
            }
            cbuf.lines.emplace_back(coord0, color0);
            cbuf.lines.emplace_back(coord1, color1);
            tok_idx += fbStride * 2;
         }
         break;
         case GL_POLYGON_TOKEN:
         {
            int n = xfb_buf[tok_idx + 1];
            tok_idx += 2;
            // get vertex 0, 1
            glm::vec3 coord0 = glm::make_vec3(&xfb_buf[tok_idx]),
                      coord1 = glm::make_vec3(&xfb_buf[tok_idx + fbStride]);
            glm::vec4 color0 = glm::make_vec4(&xfb_buf[tok_idx + 3]),
                      color1 = glm::make_vec4(&xfb_buf[tok_idx + 3 + fbStride]);
            if (fbStride == 11)
            {
               // get texture
               pal.GetColorFromVal(xfb_buf[tok_idx + 7], glm::value_ptr(color0));
               pal.GetColorFromVal(xfb_buf[tok_idx + 7 + fbStride], glm::value_ptr(color1));
            }
            // decompose polygon into n-2 triangles [0 1 2] [0 2 3] ...
            for (int i = 0; i < n-2; i++)
            {
               // get last vertex of current triangle
               int vtxStart = fbStride * (2 + 3*i);
               glm::vec3 coord2 = glm::make_vec3(&xfb_buf[tok_idx + vtxStart]);
               glm::vec4 color2 = glm::make_vec4(&xfb_buf[tok_idx + 3 + vtxStart]);
               if (fbStride == 11)
               {
                  pal.GetColorFromVal(xfb_buf[tok_idx + 7 + vtxStart], glm::value_ptr(color2));
               }
               cbuf.triangles.emplace_back(coord0, color0);
               cbuf.triangles.emplace_back(coord1, color1);
               cbuf.triangles.emplace_back(coord2, color2);
               // last vertex becomes second vertex of next triangle
               coord1 = coord2;
               color1 = color2;
            }
            tok_idx += n * fbStride;
         }
         break;
         case GL_POINT_TOKEN:
         case GL_BITMAP_TOKEN:
         case GL_DRAW_PIXEL_TOKEN:
         case GL_COPY_PIXEL_TOKEN:
         default:
            // commands containing the token + a single vertex ignore for now
            tok_idx += 1 + fbStride;
            break;
      }
   }
}

void CapturePass::Render(const RenderQueue& queue)
{
   int color_tex = palette->GetColorTexture();
   int alpha_tex = palette->GetAlphaTexture();
   bool always_blend = device->isBlendEnabled();
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

      device->attachTexture(GLDevice::SAMPLER_COLOR, color_tex);
      device->attachTexture(GLDevice::SAMPLER_ALPHA, alpha_tex);
      for (auto buf : tex_bufs)
      {
         device->captureXfbBuffer(*palette, cbuf, buf);
      }
      device->detachTexture(GLDevice::SAMPLER_COLOR);
      device->detachTexture(GLDevice::SAMPLER_ALPHA);
      for (auto buf : no_tex_bufs)
      {
         device->captureXfbBuffer(*palette, cbuf, buf);
      }
      if (!params.contains_translucent)
      {
         device->enableBlend();
         device->disableDepthWrite();
      }
      device->attachTexture(1, font_tex);
      device->setNumLights(0);
      device->captureXfbBuffer(cbuf, curr_drawable->getTextBuffer());
      if (!always_blend) { device->disableBlend(); }
   }
}

}
