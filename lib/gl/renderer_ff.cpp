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

#include "renderer_ff.hpp"
#include "attr_traits.hpp"
#include "../aux_vis.hpp"

namespace gl3
{

template<typename TVtx>
void setupFFVertexArray(TVtx* buffer)
{
   static_assert(AttrCoord<TVtx>::exists,
                 "Invalid vertex type, requires at least TVtx::coord to be present.");
   AttrCoord<TVtx>::setupLegacy(buffer);
   AttrNormal<TVtx>::setupLegacy(buffer);
   AttrColor<TVtx>::setupLegacy(buffer);
   glClientActiveTexture(GL_TEXTURE0);
   AttrTexcoord<TVtx>::setupLegacy(buffer);
   glClientActiveTexture(GL_TEXTURE1);
   AttrTexcoord<TVtx>::setupLegacy(buffer);
}

template<typename TVtx>
void clearFFVertexArray()
{
   AttrCoord<TVtx>::clearLegacy();
   AttrNormal<TVtx>::clearLegacy();
   AttrColor<TVtx>::clearLegacy();
   glClientActiveTexture(GL_TEXTURE0);
   AttrTexcoord<TVtx>::clearLegacy();
   glClientActiveTexture(GL_TEXTURE1);
   AttrTexcoord<TVtx>::clearLegacy();
}

template<typename TVtx>
void FFGLDevice::bufferFFDeviceImpl(const VertexBuffer<TVtx>& buf)
{
   glNewList(disp_lists[buf.getHandle()].list, GL_COMPILE);
   setupFFVertexArray((TVtx*)buf.getData());
   glDrawArrays(buf.getShape(), 0, buf.count());
   glEndList();
   clearFFVertexArray<TVtx>();
}

template<typename TVtx>
void FFGLDevice::bufferFFDeviceImpl(const IndexedVertexBuffer<TVtx>& buf)
{
   glNewList(disp_lists[buf.getHandle()].list, GL_COMPILE);
   setupFFVertexArray((TVtx*)buf.getData());
   glDrawElements(buf.getShape(), buf.getIndices().size(), GL_UNSIGNED_INT,
                  buf.getIndices().data());
   glEndList();
   clearFFVertexArray<TVtx>();
}

void FFGLDevice::init()
{
   GLDevice::init();
   // Fixed-function pipeline parameters
   glEnable(GL_NORMALIZE);
   glShadeModel(GL_SMOOTH);
   glEnable(GL_COLOR_MATERIAL);
   glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
   // Texturing is set up such that the output color is computed as follows:
   // - color_out.rgb = color_in.rgb * tex0.rgb
   // - color_out.a = tex1.a
   // Texture unit 0 should contain the color palette, while texture unit 1
   // contains either the transparency alpha channel, or the font texture
   glActiveTexture(GL_TEXTURE0);
   glEnable(GL_TEXTURE_2D);
   glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
   glActiveTexture(GL_TEXTURE1);
   glEnable(GL_TEXTURE_2D);
   // With a GL_ALPHA format texture loaded, GL_MODULATE will pass through the
   // fragment rgb value.
   glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
}

void FFGLDevice::setTransformMatrices(glm::mat4 model_view,
                                      glm::mat4 projection)
{
   GLDevice::setTransformMatrices(model_view, projection);
   glMatrixMode(GL_MODELVIEW);
   glLoadMatrixf(glm::value_ptr(model_view));
   glMatrixMode(GL_PROJECTION);
   glLoadMatrixf(glm::value_ptr(projection));
}

void FFGLDevice::setNumLights(int i)
{
   if (i == 0)
   {
      glDisable(GL_LIGHTING);
      return;
   }
   glEnable(GL_LIGHTING);
   for (int light_id = 0; light_id < i; light_id++)
   {
      glEnable(GL_LIGHT0 + light_id);
   }
   for (int light_id = i; light_id < LIGHTS_MAX; light_id++)
   {
      glDisable(GL_LIGHT0 + light_id);
   }
}

void FFGLDevice::setMaterial(Material mat)
{
   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat.diffuse.data());
   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat.ambient.data());
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat.specular.data());
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, mat.shininess);

   GLfloat memis[] = { 0.0, 0.0, 0.0, 1.0 };
   glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, memis);
}

void FFGLDevice::setPointLight(int i, Light lt)
{
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();
   glLightfv(GL_LIGHT0 + i, GL_POSITION, lt.position.data());

   GLfloat lambi[] = { 0.0, 0.0, 0.0, 1.0 };
   glLightfv(GL_LIGHT0 + i, GL_AMBIENT, lambi);

   glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, lt.diffuse.data());
   glLightfv(GL_LIGHT0 + i, GL_SPECULAR, lt.specular.data());
   glPopMatrix();
}

void FFGLDevice::setAmbientLight(const std::array<float, 4>& amb)
{
   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, &amb[0]);
}

void FFGLDevice::setClipPlaneUse(bool enable)
{
   if (enable) { glEnable(GL_CLIP_PLANE0); }
   else { glDisable(GL_CLIP_PLANE0); }
}

void FFGLDevice::setClipPlaneEqn(const std::array<double, 4>& eqn)
{
   glClipPlane(GL_CLIP_PLANE0, eqn.data());
}

void FFGLDevice::bufferToDevice(ArrayLayout layout, IVertexBuffer& buf)
{
   if (buf.getHandle() == 0)
   {
      if (buf.count() == 0) { return; }
      GLuint new_hnd = glGenLists(1);
      buf.setHandle(disp_lists.size());
      disp_lists.emplace_back(DispListData_{new_hnd, buf.getShape(), buf.count(), layout});
   }
   else
   {
      disp_lists[buf.getHandle()].count = buf.count();
   }

   switch (layout)
   {
      case Vertex::layout:
         bufferFFDeviceImpl(static_cast<const VertexBuffer<Vertex>&>(buf));
         break;
      case VertexColor::layout:
         bufferFFDeviceImpl(static_cast<const VertexBuffer<VertexColor>&>(buf));
         break;
      case VertexTex::layout:
         bufferFFDeviceImpl(static_cast<const VertexBuffer<VertexTex>&>(buf));
         break;
      case VertexNorm::layout:
         bufferFFDeviceImpl(static_cast<const VertexBuffer<VertexNorm>&>(buf));
         break;
      case VertexNormColor::layout:
         bufferFFDeviceImpl(static_cast<const VertexBuffer<VertexNormColor>&>(buf));
         break;
      case VertexNormTex::layout:
         bufferFFDeviceImpl(static_cast<const VertexBuffer<VertexNormTex>&>(buf));
         break;
      default:
         cerr << "WARNING: Unhandled vertex layout " << layout << endl;
   }
}

void FFGLDevice::bufferToDevice(ArrayLayout layout, IIndexedBuffer& buf)
{
   if (buf.getHandle() == 0)
   {
      if (buf.count() == 0) { return; }
      GLuint new_hnd = glGenLists(1);
      buf.setHandle(disp_lists.size());
      disp_lists.emplace_back(DispListData_{new_hnd, buf.getShape(), buf.getIndices().size(), layout});
   }
   else
   {
      disp_lists[buf.getHandle()].count = buf.getIndices().size();
   }

   switch (layout)
   {
      case Vertex::layout:
         bufferFFDeviceImpl(static_cast<const IndexedVertexBuffer<Vertex>&>(buf));
         break;
      case VertexColor::layout:
         bufferFFDeviceImpl(static_cast<const IndexedVertexBuffer<VertexColor>&>(buf));
         break;
      case VertexTex::layout:
         bufferFFDeviceImpl(static_cast<const IndexedVertexBuffer<VertexTex>&>(buf));
         break;
      case VertexNorm::layout:
         bufferFFDeviceImpl(static_cast<const IndexedVertexBuffer<VertexNorm>&>(buf));
         break;
      case VertexNormColor::layout:
         bufferFFDeviceImpl(static_cast<const IndexedVertexBuffer<VertexNormColor>&>
                            (buf));
         break;
      case VertexNormTex::layout:
         bufferFFDeviceImpl(static_cast<const IndexedVertexBuffer<VertexNormTex>&>(buf));
         break;
      default:
         cerr << "WARNING: Unhandled vertex layout " << layout << endl;
   }
}

void FFGLDevice::bufferToDevice(TextBuffer&)
{
   // we can't really do anything here can only compute offset matrix at draw
}

void FFGLDevice::drawDeviceBuffer(int hnd)
{
   if (hnd == 0) { return; }
   if (disp_lists[hnd].count == 0) { return; }
   if (disp_lists[hnd].layout == Vertex::layout
       || disp_lists[hnd].layout == VertexNorm::layout)
   {
      glColor4fv(static_color.data());
   }
   else
   {
      glColor4f(1.f, 1.f, 1.f, 1.f);
   }
   if (!( disp_lists[hnd].layout == VertexNorm::layout
          || disp_lists[hnd].layout == VertexNormColor::layout
          || disp_lists[hnd].layout == VertexNormTex::layout))
   {
      glNormal3f(0.f, 0.f, 1.f);
   }
   glCallList(disp_lists[hnd].list);
   // reset texturing parameters
   // glMultiTexCoord2f(GL_TEXTURE0, 0.f, 0.f);
   // glMultiTexCoord2f(GL_TEXTURE1, 0.f, 0.f);
}

void FFGLDevice::drawDeviceBuffer(const TextBuffer& buf)
{
   glColor4fv(static_color.data());
   glNormal3f(0.f, 0.f, 1.f);
   glMultiTexCoord2f(GL_TEXTURE0, 0.f, 0.f);
   float tex_w = GetFont()->getAtlasWidth();
   float tex_h = GetFont()->getAtlasHeight();
   // Model-view transform:
   // - scale bounding boxes to relative clip-space/NDC coords
   // - add z-offset of -0.005 to reduce text hiding by mesh
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();
   glScalef(2.f / vp_width, 2.f / vp_height, 0.f);
   glTranslatef(0.f, 0.f, -0.005);
   glMatrixMode(GL_PROJECTION);
   for (const TextBuffer::Entry& e : buf)
   {
      glm::vec4 pos(e.rx, e.ry, e.rz, 1.0);
      // transform text starting position into NDC
      pos = model_view_mtx * pos;
      pos = proj_mtx * pos;
      pos = pos / pos.w;
      // Projection transform:
      // - add starting offset in NDC
      glPushMatrix();
      glLoadIdentity();
      glTranslatef(pos.x, pos.y, pos.z);
      float pen_x = e.ox, pen_y = e.oy;
      char prev_c = '\0';
      for (char c : e.text)
      {
         const GlVisFont::glyph &g = GetFont()->GetTexChar(c);
         pen_x += GetFont()->GetKerning(prev_c, c);
         // note: subtract 1 to account for the padding in the texture glyphs
         float cur_x = pen_x + g.bear_x - 1;
         float cur_y = -pen_y - g.bear_y - 1;
         pen_x += g.adv_x;
         pen_y += g.adv_y;
         if (!g.w || !g.h)
         {
            continue;
         }
         glBegin(GL_TRIANGLE_STRIP);
         glMultiTexCoord2f(GL_TEXTURE1, g.tex_x, 0);
         glVertex2f(cur_x, -cur_y);
         glMultiTexCoord2f(GL_TEXTURE1, g.tex_x + g.w / tex_w, 0);
         glVertex2f(cur_x + g.w, -cur_y);
         glMultiTexCoord2f(GL_TEXTURE1, g.tex_x, g.h / tex_h);
         glVertex2f(cur_x, -cur_y - g.h);
         glMultiTexCoord2f(GL_TEXTURE1, g.tex_x + g.w / tex_w, g.h / tex_h);
         glVertex2f(cur_x + g.w, -cur_y - g.h);
         glEnd();
      }
      glPopMatrix();
   }
   glMatrixMode(GL_MODELVIEW);
   glPopMatrix();
}

}
