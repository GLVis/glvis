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

#ifndef GLVIS_FRAMEBUFFER_HPP
#define GLVIS_FRAMEBUFFER_HPP

#include "types.hpp"
#include <vector>

namespace gl3
{

class Framebuffer
{
public:
   Framebuffer() = default;

   explicit operator bool() const
   {
      return handle;
   }

   bool IsComplete() const
   {
      glBindFramebuffer(GL_FRAMEBUFFER, handle);
      return glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
   }

   void Init()
   {
      GLuint fb_id;
      glGenFramebuffers(1, &fb_id);
      handle = fb_id;
   }

   void Detach(GLuint attach_point)
   {
      if (handle == 0)
      {
         std::cerr << "Can't attach textures from a default framebuffer."
                   << std::endl;
         return;
      }
      if (attach_point >= GL_COLOR_ATTACHMENT0 &&
          attach_point < GL_COLOR_ATTACHMENT0 + NUM_ATTACHMENTS)
      {
         int attach_idx = attach_point - GL_COLOR_ATTACHMENT0;
         color_attach_active[attach_idx] = false;
         color_attached_textures[attach_idx] = 0;
      }
   }

   void Attach(GLuint attach_point,
               GLuint texture_binding,
               const resource::TextureHandle& tex_handle)
   {
      if (handle == 0)
      {
         std::cerr << "Can't attach textures to a default framebuffer."
                   << std::endl;
         return;
      }
      if (attach_point >= GL_COLOR_ATTACHMENT0 &&
          attach_point < GL_COLOR_ATTACHMENT0 + NUM_ATTACHMENTS)
      {
         int attach_idx = attach_point - GL_COLOR_ATTACHMENT0;
         color_attach_active[attach_idx] = true;
         color_attached_textures[attach_idx] = tex_handle;
      }
      glBindFramebuffer(GL_FRAMEBUFFER, handle);
      glFramebufferTexture2D(GL_FRAMEBUFFER, attach_point,
                             texture_binding, tex_handle, 0);
   }


   void Attach(GLuint attach_point,
               const resource::RenderBufHandle& renderbuf_handle)
   {
      if (handle == 0)
      {
         std::cerr << "Can't attach renderbuffers to a default framebuffer."
                   << std::endl;
         return;
      }
      if (attach_point >= GL_COLOR_ATTACHMENT0 &&
          attach_point < GL_COLOR_ATTACHMENT0 + NUM_ATTACHMENTS)
      {
         int attach_idx = attach_point - GL_COLOR_ATTACHMENT0;
         color_attach_active[attach_idx] = true;
         color_attached_textures[attach_idx] = 0;
      }
      glBindFramebuffer(GL_FRAMEBUFFER, handle);
      glFramebufferRenderbuffer(GL_FRAMEBUFFER, attach_point,
                                GL_RENDERBUFFER, renderbuf_handle);
   }

   void Bind(int nbufs = 0) const
   {
      glBindFramebuffer(GL_FRAMEBUFFER, handle);
      if (handle == 0)
      {
          if (nbufs > 1)
          {
              std::cerr << "Default framebuffer only has one valid output "
                        << "buffer." << std::endl;
          }
          glDrawBuffer(GL_BACK);
          return;
      }
      // Set active color attachments in order
      std::vector<GLenum> output_bufs;
      for (int ibuf = 0; ibuf < NUM_ATTACHMENTS; ibuf++)
      {
         if (color_attach_active[ibuf])
         {
            output_bufs.push_back(GL_COLOR_ATTACHMENT0 + ibuf);
         }
      }
      if (nbufs == 0)
      {
         nbufs = output_bufs.size();
      }
      while (output_bufs.size() < nbufs)
      {
         output_bufs.push_back(GL_NONE);
      }
      glDrawBuffers(nbufs, output_bufs.data());
   }

   void Bind(const std::vector<GLenum>& drawbufs) const
   {
      glBindFramebuffer(GL_FRAMEBUFFER, handle);
      glDrawBuffers(drawbufs.size(), drawbufs.data());
   }

   // Unbinds the current framebuffer and binds the default framebuffer.
   void Unbind() const
   {
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      GLenum drawOutput = GL_BACK;
      glDrawBuffers(1, &drawOutput);
   }

   // Blits an image located on GL_COLOR_ATTACHMENT0 of a source framebuffer
   // to all the color attachments of this framebuffer.
   void BlitFrom(const Framebuffer& fb_from,
                 int w, int h,
                 GLenum filter = GL_NEAREST) const;

private:
   resource::FBOHandle handle;
   static constexpr int NUM_ATTACHMENTS=8;
   bool color_attach_active[NUM_ATTACHMENTS] = {};
   GLuint color_attached_textures[NUM_ATTACHMENTS] = {};
};

}

#endif
