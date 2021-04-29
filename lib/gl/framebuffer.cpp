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

#include "framebuffer.hpp"
#include "renderer.hpp"
#include "shader.hpp"

static const std::string BlitVertexShader =
#include "shaders/depth_peel_passthrough.vert"
   ;

static const std::string BlitFragShader =
#include "shaders/fb_blit.frag"
   ;

namespace gl3
{

class ShaderBasedBlit
{
public:
   ShaderBasedBlit()
   {
      if (!passthrough_shader.create(BlitVertexShader, BlitFragShader,
   {{ 0, "vertex" }}, 1))
      {
         std::cerr << "Unable to create blit main program." << std::endl;
      }
      GLuint vbo;
      glGenBuffers(1, &vbo);
      quad_buffer = BufObjHandle{vbo};
      float quad_verts[] =
      {
         -1.f, 1.f, -1.f, -1.f, 1.f, -1.f,
            -1.f, 1.f, 1.f, -1.f, 1.f, 1.f,
         };
      glBindBuffer(GL_ARRAY_BUFFER, vbo);
      glBufferData(GL_ARRAY_BUFFER, 12 * sizeof(float), quad_verts, GL_STATIC_DRAW);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
   }
   ShaderProgram& GetPassthroughShader() { return passthrough_shader; }
   BufObjHandle& GetQuadBuffer() { return quad_buffer; }
private:
   ShaderProgram passthrough_shader;
   BufObjHandle quad_buffer;
};

static ShaderBasedBlit& GetShaderBasedBlit()
{
   static ShaderBasedBlit blit;
   return blit;
}

void Framebuffer::BlitFrom(const Framebuffer &fb_from, int w, int h,
                           GLenum filter) const
{
   if (GLDevice::isOpenGL3())
   {
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, handle);
      glBindFramebuffer(GL_READ_FRAMEBUFFER, fb_from);
      glBlitFramebuffer(0, 0, w, h,
                        0, 0, w, h,
                        GL_COLOR_BUFFER_BIT, filter);
   }
   else
   {
      if (!(filter == GL_NEAREST || filter == GL_LINEAR))
      {
         std::cerr << "Blit error: no texture bound to color attachment 0 "
                   << "of source framebuffer" << std::endl;
         return;
      }
      ShaderProgram& blit_shader = GetShaderBasedBlit().GetPassthroughShader();
      BufObjHandle& quad_buffer = GetShaderBasedBlit().GetQuadBuffer();

      blit_shader.bind();
      Bind();
      GLuint texture_src = fb_from.color_attached_textures[0];
      if (texture_src == 0)
      {
         std::cerr << "Blit error: no texture bound to color attachment 0 "
                   << "of source framebuffer" << std::endl;
         return;
      }
      glUniform1i(blit_shader.uniform("sourceColor"), 0);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, texture_src);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);

      glBindBuffer(GL_ARRAY_BUFFER, quad_buffer);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 2, GL_FLOAT, false, 0, 0);
      glDrawArrays(GL_TRIANGLES, 0, 6);
   }
}

}
